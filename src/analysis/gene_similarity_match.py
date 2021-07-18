import heapq
import multiprocessing
import os
import re
from random import random

import pandas as pd
from typing import List, Mapping
from concurrent.futures import ThreadPoolExecutor

from analysis.models.match_candidate import MatchCandidate
from analysis.models.similarity_type import SimilarityType
from analysis.models.order_type import OrderType
from analysis.gene_location_analysis import GeneLocationAnalysis
from analysis.similarities.pattern_similarity import MatchPattern
from analysis.similarities.similarity_factory import SimilarityFactory
from utils.factories.logger_factory import LoggerFactory
from utils.ncbi_database import NCBIDatabase
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter
from collections import deque

CandidateClearSize = 10000


class GeneSimilarityMatch:
    def __init__(self,
                 gene_path: str,
                 data_path: str,
                 output_directory: str,
                 top_k: int = 20,
                 candidate_distance: int = 5,
                 patience: int = 0,
                 weighted: Mapping[SimilarityType, int] = None,
                 conditions: dict = None,
                 continuous_mismatch_limit: int = None,
                 order_type: OrderType = OrderType.Decrement,
                 gene_name_filter=None):
        self.gene_path = gene_path
        self.data_path = data_path
        self.output_directory = output_directory
        self.top_k = top_k
        self.candidate_distance = candidate_distance
        self.patience = patience
        self.weighted = {k: v for k, v in weighted.items() if v > 0}
        self.conditions = conditions
        self.continuous_mismatch_limit = continuous_mismatch_limit
        self.order_type = order_type
        self.gene_name_filter = gene_name_filter
        self.data_name = os.path.basename(self.data_path)
        self.dna_code = None
        self.rev_dna_code = None

        file_name = os.path.basename(self.gene_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory,
                                        '%s_match_result.txt' % file_prefix)
        self.gene_reader = NCBIDatabase(self.data_path)
        self.logger = LoggerFactory()
        self.weighted_sum = sum([v for k, v in self.weighted.items()])
        assert self.weighted_sum > 0
        self.initialize()

    def initialize(self):
        self.dna_code = self.gene_reader.dna_code
        self.rev_dna_code = get_opposite_dna(self.gene_reader.dna_code[::-1])

    def run(self, gene_name_filter: GeneLocationAnalysis = None):
        self.gene_name_filter = gene_name_filter
        with open(self.result_path, 'w', encoding='utf8') as fw:
            gene_datas = pd.read_csv(self.gene_path, sep='\t')
            records = [record for _, record in gene_datas.iterrows()]
            solved = 0
            total = len(gene_datas)
            self.logger.info_with_expire_time(
                'Doing Similarity Matching: %d/%d(%.2f%%)' % (
                    solved, total, solved * 100.0 / total), solved, total)
            with multiprocessing.Pool(min(multiprocessing.cpu_count(), 2)) as p:
                for ret in p.imap(self.find_candidate_for_gene, records):
                    fw.write(ret)
                    fw.flush()
                    solved += 1
                    self.logger.info_with_expire_time(
                        'Doing Similarity Matching: %d/%d(%.2f%%)' % (
                            solved, total, solved * 100.0 / total), solved, total)

    def find_candidate_for_gene(self, record: pd.Series):
        def next_interval(size, thread_cnt):
            batch_size = size // thread_cnt
            last_pos = 0
            while last_pos < size:
                next_pos = min(last_pos + batch_size, size)
                yield last_pos, next_pos
                last_pos = next_pos

        name, gene = record['name'], record['gene'].lower()
        candidates = []
        with ThreadPoolExecutor() as executor:
            tasks = []
            for start, end in next_interval(len(self.dna_code), 32):
                tasks.append(executor.submit(self.match_gene, name, gene, False, start, end))
                tasks.append(executor.submit(self.match_gene, name, gene, True, start, end))
            for task in tasks:
                candidates.extend(task.result())
        candidates = list(candidates)
        candidates.sort(key=lambda arg: -arg.weighted_similarity)
        candidates = candidates[:self.top_k]
        if self.order_type == OrderType.Increment:
            for candidate in candidates:
                candidate.weighted_similarity = -candidate.weighted_similarity
        results = self.render_similarity_for_candidates(gene, candidates[:self.top_k])

        headers = [
            'name',
            'direction',
            'weighted_similarity'
        ]
        for similarity_name, weight in self.weighted.items():
            headers.append(similarity_name.name.lower() + '_similarity')
        headers.append('original      :')
        sequence_headers = [
            'gene_format   :',
            'target_format :',
            'match_format  :']

        content = ''
        idx = 1
        for candidate_result in results:
            candidate = candidate_result[0]
            content += '(%d)\n' % idx
            attribute = {
                'name': name,
                'direction': '-' if candidate.is_reverse else '+',
                'weighted_similarity': '%.2f' % candidate.weighted_similarity,
                'original      :': gene
            }
            sequence_content = []
            offset = 1
            for similarity_name, weight in sorted(self.weighted.items(), key=lambda arg: arg[0]):
                attribute[similarity_name.name.lower() + '_similarity'] = '%.2f' % candidate.similarity_dict[
                    similarity_name]
                for sequence_header, value in zip(sequence_headers, candidate_result[offset:offset + 3]):
                    value = ''.join(value)
                    sequence_content.append(similarity_name.name.lower() + "_" + sequence_header + '=' + value)
                offset += 3

            content += '>%s/%s-%s\t%s,%s\n\n' % (
                self.data_name.replace(".txt", ''),
                candidate.start,
                candidate.end,
                ','.join(['%s=%s' % (key, attribute[key]) for key in headers if key in attribute]),
                ','.join(sequence_content)
            )
            idx += 1
        return content

    def match_gene(self, name, gene, is_reverse, start, end):
        database = self.rev_dna_code if is_reverse else self.dna_code
        candidates: List[MatchCandidate] = []
        gene_length = len(gene)
        min_weighted_similarity_in_candidates = 0.0
        database_length = len(database)
        tot = end - start
        similarity_heap = []
        buff = deque()
        match_pattern = MatchPattern(gene, self.conditions) if self.conditions else None
        current_logger = LoggerFactory(int(10 + random() * 120))
        solved = 0
        msg = 'Analysis [%s][%s] %d/%d(%.2f%%)' % (
            name,
            '-' if is_reverse else '+',
            solved,
            tot,
            solved * 100.0 / tot)
        current_logger.info_with_expire_time(msg,
                                             solved,
                                             tot)
        end = min(database_length - gene_length + 1, end)
        for offset in range(start, end):
            weighted_similarity, similarity_dict = count_similarity(
                weighted=self.weighted,
                gene=gene,
                database=database,
                offset=offset,
                max_patience=self.patience,
                match_pattern=match_pattern,
                continuous_mismatch_limit=self.continuous_mismatch_limit)
            if self.order_type == OrderType.Increment:
                weighted_similarity = -weighted_similarity
            new_candidate = MatchCandidate(
                left=offset,
                right=offset + gene_length - 1,
                is_reverse=is_reverse,
                database_length=database_length,
                weighted_similarity=weighted_similarity,
                similarity_dict=similarity_dict)

            added_flag = update_candidate_list(new_candidate,
                                               buff,
                                               candidates,
                                               self.candidate_distance)
            if added_flag:
                heapq.heappush(similarity_heap, candidates[-1])
                if len(similarity_heap) > self.top_k:
                    heapq.heappop(similarity_heap)
                    top = similarity_heap[0]
                    min_weighted_similarity_in_candidates = max(min_weighted_similarity_in_candidates,
                                                                top.weighted_similarity)

            solved += 1
            msg = 'Analysis for %s[%s](%d~%d): %d/%d(%.2f%%) ' \
                  '--top_k=%d --top_similarity_info=[%s] ' \
                  '--gene_length=%d --candidates_num=%d' % (
                      name,
                      '-' if is_reverse else '+',
                      start,
                      end,
                      solved,
                      tot,
                      solved * 100.0 / tot,
                      self.top_k,
                      similarity_heap[0].get_similarity_str() if len(similarity_heap) > 0 else 'None',
                      gene_length,
                      len(candidates)
                  )
            current_logger.info_with_expire_time(msg,
                                                 solved,
                                                 tot)

            if len(candidates) > CandidateClearSize:
                candidates.sort(key=lambda arg: -arg.weighted_similarity)
                candidates = candidates[:self.top_k]
        while len(buff) > 0:
            update_candidate_list(None, buff, candidates, 1)
        return candidates

    def render_similarity_for_candidates(self, gene, candidates):
        result = []
        for candidate in candidates:
            database = self.rev_dna_code if candidate.is_reverse else self.dna_code
            candidate_result = [candidate]
            for similarity_type, weight in sorted(self.weighted.items(), key=lambda arg: arg[0]):
                candidate_result.extend(
                    self.render_target_dna_sequence(similarity_type, gene, database, candidate.original_match_left))
            result.append(candidate_result)
        return result

    def render_target_dna_sequence(self, similarity_type: SimilarityType, gene, database, offset):
        similarity_render = SimilarityFactory.get_similarity(
            similarity_type=similarity_type,
            continuous_mismatch_limit=self.continuous_mismatch_limit,
            max_patience=self.patience,
            match_pattern=MatchPattern(gene, self.conditions) if self.conditions else None
        )
        sequence_gene, sequence_target, sequence = similarity_render.rendering_sequence(gene, database, offset)
        return sequence_gene, sequence_target, sequence


def update_candidate_list(new_candidate: MatchCandidate, buff: deque, candidate_result: list,
                          keep_size: int):
    added = False
    while (len(buff) >= keep_size) or (
            len(buff) > 0 and new_candidate is not None and abs(buff[0].start - new_candidate.start) >= keep_size):
        old_candidate = buff.popleft()
        if not old_candidate.should_ignore:
            candidate_result.append(old_candidate)
            added = True
    if new_candidate is not None:
        if new_candidate.weighted_similarity <= 0.0:
            new_candidate.should_ignore = True
        for candidate in buff:
            if candidate.weighted_similarity > new_candidate.weighted_similarity:
                new_candidate.should_ignore = True
            elif candidate.weighted_similarity < new_candidate.weighted_similarity:
                candidate.should_ignore = True
        buff.append(new_candidate)
    return added


def fast_skip(gene_dict, gene_length, database, offset, cut_same, pat):
    if pat is not None:
        if not re.match(pat, database[offset:offset + gene_length]):
            return True
    gene_dict_database = count_acgt(database[offset:offset + gene_length])
    same_count = 0
    for x, cnt in gene_dict.items():
        if x in gene_dict_database:
            same_count += min(cnt, gene_dict_database[x])
            if same_count >= cut_same:
                return False
    return True


def count_acgt(gene):
    gene_dict = {}
    for x in gene:
        if x not in gene_dict:
            gene_dict[x] = 1
        else:
            gene_dict[x] += 1
    return gene_dict


def count_similarity(weighted: Mapping[SimilarityType, int],
                     gene: str,
                     database: str,
                     offset: int,
                     max_patience=2,
                     match_pattern=None,
                     continuous_mismatch_limit=None):
    weighted_similarity = 0.0
    total_weight = 0.0
    similarity = {}
    for similarity_type, weight in weighted.items():
        similarity_compute = SimilarityFactory.get_similarity(
            similarity_type=similarity_type,
            continuous_mismatch_limit=continuous_mismatch_limit,
            max_patience=max_patience,
            match_pattern=match_pattern,
            mid_limit=10,
            end_limit=2
        )
        score, info = similarity_compute.get_similarity(gene, database, offset)
        similarity[similarity_type] = score
        weighted_similarity += score * weight
        total_weight += weight
    weighted_similarity = weighted_similarity / total_weight
    return weighted_similarity, similarity


def update_or_add_min_val(dp, i, j, update_val):
    if i not in dp:
        dp[i] = {}
    if j not in dp[i]:
        dp[i][j] = update_val
        return True
    if dp[i][j] > update_val:
        dp[i][j] = update_val
        return True
    return False
