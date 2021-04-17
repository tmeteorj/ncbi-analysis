import heapq
import os
import random
import re
import threading
import time
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Mapping

from analysis.enum_types import SimilarityType
from analysis.gene_location_analysis import GeneLocationAnalysis
from utils.factories.logger_factory import LoggerFactory
from utils.gene_file_util import GeneFileReader
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter
from collections import deque

CandidateClearSize = 10000



@dataclass
class MatchCandidate:
    left: int
    right: int
    is_reverse: bool
    database_length: int
    weighted_similarity: float
    similarity_dict: Mapping[SimilarityType, float] = field(default_factory=dict)

    def __post_init__(self):
        if self.is_reverse:
            self.start = self.database_length - self.left
            self.end = self.database_length - self.right
        else:
            self.start = self.left + 1
            self.end = self.right + 1
        self.should_ignore = False
        self.original_match_left = self.left
        self.original_match_right = self.right

    def get_similarity_str(self):
        result = 'weighted=%.2f' % self.weighted_similarity
        for key, value in self.similarity_dict.items():
            result += ', %s=%.2f' % (key.name, value)
        return result

    def __str__(self):
        return '[%d-%d], (%s), %s, ' % (self.start, self.right, self.get_similarity_str(), self.should_ignore)

    def __le__(self, other):
        return self.weighted_similarity <= other.weighted_similarity

    def __lt__(self, other):
        return self.weighted_similarity < other.weighted_similarity

    def __ge__(self, other):
        return self.weighted_similarity >= other.weighted_similarity

    def __gt__(self, other):
        return self.weighted_similarity > other.weighted_similarity


class HasReturnThread(threading.Thread):
    def __init__(self, func, args, name=''):
        threading.Thread.__init__(self)
        self.name = name
        self.func = func
        self.args = args
        self.result = None
        self.local_thread = None

    def start(self) -> None:
        new_args = [self.func]
        new_args.extend(self.args)
        new_args = tuple(new_args)
        self.local_thread = threading.Thread(target=self.trace_func, args=new_args)
        self.local_thread.start()

    def join(self) -> None:
        self.local_thread.join()

    def trace_func(self, func, *args, **kwargs):
        ret = func(*args, **kwargs)
        self.result = ret

    def get_result(self):
        try:
            return self.result
        except Exception:
            return None


@dataclass
class GeneSimilarityMatch:
    gene_path: str
    data_path: str
    output_directory: str
    top_k: int = 20
    candidate_distance: int = 5
    batch_size: int = 5
    patience: int = 0
    weighted: List[int] = field(default_factory=list)
    conditions: dict = None
    continuous_mismatch_limit: int = None
    order_type: OrderType = OrderType.Decrement
    dna_code = None
    rev_dna_code = None
    gene_name_filter = None

    def __post_init__(self):
        self.data_name = os.path.basename(self.data_path)
        file_name = os.path.basename(self.gene_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory,
                                        '%s_match_result.txt' % (file_prefix))
        self.gene_reader = GeneFileReader(self.data_path)
        self.dna_code = None
        self.rev_dna_code = None
        self.logger = LoggerFactory()

        self.lock = threading.Lock()
        self.solved = 0
        self.total = 0
        self.weighted_sum = sum(self.weighted)
        assert self.weighted_sum > 0 and len(self.weighted) == 5

    def run(self, gene_name_filter: GeneLocationAnalysis = None):
        self.gene_name_filter = gene_name_filter
        self.gene_reader.build_information()
        self.dna_code = self.gene_reader.dna_code
        self.rev_dna_code = get_opposite_dna(self.gene_reader.dna_code[::-1])
        with open(self.result_path, 'w', encoding='utf8') as fw:
            gene_sequences = open(self.gene_path, 'r', encoding='utf8').readlines()[1:]
            self.solved = 0
            self.total = len(self.gene_reader.dna_code) * len(gene_sequences) * 2
            self.logger.info_with_expire_time(
                'Doing Similarity Matching: %d/%d(%.2f%%)' % (
                    self.solved, self.total, self.solved * 100.0 / self.total), self.solved, self.total)
            pending_tasks = deque()
            running_tasks = []
            for gene_sequence in gene_sequences:
                items = gene_sequence.strip().split('\t')
                name, gene = items[0], items[1].lower()
                t = threading.Thread(target=self.find_candidate_for_gene, args=(name, gene, fw,))
                pending_tasks.append(t)
            while len(pending_tasks) > 0:
                running_tasks = [t for t in running_tasks if t.isAlive()]
                while len(running_tasks) < self.batch_size and len(pending_tasks) > 0:
                    t = pending_tasks.popleft()
                    t.start()
                    running_tasks.append(t)
                time.sleep(10)
            for t in running_tasks:
                t.join()

    def find_candidate_for_gene(self, name, gene, fw):

        t1 = HasReturnThread(func=self.match_gene,
                             args=(name, gene, self.dna_code, False,))
        t1.start()
        t2 = HasReturnThread(func=self.match_gene,
                             args=(name, gene, self.rev_dna_code, True,))
        t2.start()
        t1.join()
        t2.join()

        candidates = t1.get_result() + t2.get_result()
        candidates = list(candidates)
        candidates.sort(key=lambda arg: -arg.weighted_similarity)
        if self.order_type == OrderType.Increment:
            for candidate in candidates:
                candidate.weighted_similarity = -candidate.weighted_similarity
        results = self.render_similarity_for_candidates(gene, candidates[:self.top_k])
        self.lock.acquire()
        idx = 1
        headers = [
            'name',
            'direction',
            'weighted_similarity'
        ]
        for idx, similarity_name in enumerate(
                ['text_distance_similarity', 'direct_match_similarity', 'consistency_similarity',
                 'pattern_similarity', 'blat_similarity']):
            if self.weighted[idx] > 0:
                headers.append(similarity_name)
        headers.append('original      :')
        sequence_headers = [
            'gene_format   :',
            'target_format :',
            'match_format  :']
        for candidate_result in results:
            candidate = candidate_result[0]
            fw.write('(%d)\n' % idx)
            attribute = {
                'name': name,
                'direction': '-' if candidate.is_reverse else '+',
                'weighted_similarity': '%.2f' % candidate.weighted_similarity,
                'original      :': gene
            }
            for idx, similarity_name in enumerate(
                    ['text_distance_similarity', 'direct_match_similarity', 'consistency_similarity',
                     'pattern_similarity', 'blat_similarity']):
                if self.weighted[idx] > 0:
                    attribute[similarity_name] = '%.2f' % candidate.similarity_dict[
                        MatchAlgorithm.get_match_algorithm_by_name(similarity_name)]
            sequence_content = []
            offset = 1
            for idx, match_algorithm in enumerate(MatchAlgorithm.get_all_items()):
                if self.weighted[idx] > 0:
                    for sequence_header, value in zip(sequence_headers, candidate_result[offset:offset + 3]):
                        value = ''.join(value)
                        sequence_content.append(match_algorithm.name + "_" + sequence_header + '=' + value)
                    offset += 3

            fw.write('>%s/%s-%s\t%s,%s\n' % (
                self.data_name.replace(".txt", ''),
                candidate.start,
                candidate.end,
                ','.join(['%s=%s' % (key, attribute[key]) for key in headers if key in attribute]),
                ','.join(sequence_content)
            ))
            fw.write('\n')
            idx += 1
        self.lock.release()

    def match_gene(self, name, gene, database, is_reverse):
        candidates: List[MatchCandidate] = []
        gene_length = len(gene)
        min_weighted_similarity_in_candidates = 0.0
        database_length = len(database)
        limitation = database_length - gene_length + 1
        new_solved = 0
        similarity_heap = []
        buff = deque()
        match_pattern = MatchPattern(gene, self.conditions) if self.conditions else None
        for start in range(limitation):
            weighted_similarity, similarity_dict = count_similarity(weighted=self.weighted,
                                                                    gene=gene,
                                                                    database=database,
                                                                    offset=start,
                                                                    is_reverse=is_reverse,
                                                                    max_patience=self.patience,
                                                                    match_pattern=match_pattern,
                                                                    continuous_mismatch_limit=self.continuous_mismatch_limit,
                                                                    gene_name_filter=self.gene_name_filter)
            if self.order_type == OrderType.Increment:
                weighted_similarity = -weighted_similarity
            new_candidate = MatchCandidate(
                left=start,
                right=start + gene_length - 1,
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

            new_solved += 1
            if random.random() * 1000 < 1:
                self.lock.acquire()
                self.solved += new_solved
                self.logger.info_with_expire_time(
                    'Doing Similarity Matching for %s[%s]: %d/%d(%.2f%%) '
                    '--top_k=%d '
                    '--top_similarity_info=[%s] '
                    '--gene_length=%d '
                    '--candidates_num=%d' % (
                        name,
                        '-' if is_reverse else '+',
                        self.solved,
                        self.total,
                        self.solved * 100.0 / self.total,
                        self.top_k,
                        similarity_heap[0].get_similarity_str() if len(similarity_heap) > 0 else 'None',
                        gene_length,
                        len(candidates)
                    ),
                    self.solved,
                    self.total)
                self.lock.release()
                new_solved = 0

            if len(candidates) > CandidateClearSize:
                candidates.sort(key=lambda arg: -arg.weighted_similarity)
                candidates = candidates[:self.top_k]
        while len(buff) > 0:
            update_candidate_list(None, buff, candidates, 1)
        self.lock.acquire()
        self.solved += new_solved + gene_length - 1
        self.lock.release()
        return candidates

    def render_similarity_for_candidates(self, gene, candidates):
        result = []
        for candidate in candidates:
            database = self.rev_dna_code if candidate.is_reverse else self.dna_code
            candidate_result = [candidate]
            for idx, match_algorithm in enumerate(MatchAlgorithm.get_all_items()):
                if self.weighted[idx] > 0:
                    candidate_result.extend(
                        self.render_target_dna_sequence(match_algorithm, gene, database, candidate.original_match_left))
            result.append(candidate_result)
        return result

    def render_target_dna_sequence(self, match_algorithm: MatchAlgorithm, gene, database, offset):
        sequence_gene = []
        sequence_target = []
        sequence = []
        tot = len(gene)
        if match_algorithm == MatchAlgorithm.text_distance:
            score, dp = compute_text_distance_similarity(gene, database, offset)
            i, j = tot, tot
            while i > 0 or j > 0:
                gene_a, gene_b = gene[i - 1] if i > 0 else '.', database[j + offset - 1] if j > 0 else '.'
                if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + should_change(gene[i - 1],
                                                                                    database[j + offset - 1]):
                    sequence_gene.append(gene_a)
                    sequence_target.append(gene_b)
                    sequence.append('*' if should_change(gene[i - 1], database[j + offset - 1]) == 0 else '.')
                    i, j = i - 1, j - 1
                elif dp[i][j] == dp[i - 1][j] + 1:
                    sequence_gene.append(gene_a)
                    sequence_target.append('.')
                    sequence.append('.')
                    i -= 1
                elif dp[i][j] == dp[i][j - 1] + 1:
                    sequence_gene.append('.')
                    sequence_target.append(gene_b)
                    sequence.append('.')
                    j -= 1
                else:
                    raise ValueError('Should not go here!')
            sequence_gene.reverse()
            sequence_target.reverse()
            sequence.reverse()
        elif match_algorithm == MatchAlgorithm.direct_match:
            for i in range(tot):
                sequence_gene.append(gene[i])
                sequence_target.append(database[i + offset])
                if not should_change(gene[i], database[i + offset]):
                    sequence.append('*')
                else:
                    sequence.append('.')
        elif match_algorithm == MatchAlgorithm.consistency:
            score, score_queue, score_merge_idx = compute_consistency_similarity(gene, database, offset,
                                                                                 self.patience)
            sequence_gene.extend(gene[:])
            sequence_target.extend(database[offset:offset + tot])
            cur_pos = 0
            for idx, (same_cnt, same_end) in enumerate(score_queue):
                same_start = same_end - same_cnt
                while cur_pos < same_start:
                    if score_merge_idx[0] < idx <= score_merge_idx[1]:
                        sequence.append('-')
                    else:
                        sequence.append('.')
                    cur_pos += 1
                while cur_pos < same_end:
                    sequence.append('*')
                    cur_pos += 1
            while cur_pos < tot:
                sequence.append('.')
                cur_pos += 1
        elif match_algorithm == MatchAlgorithm.pattern:
            for i in range(tot):
                sequence_gene.append(gene[i])
                sequence_target.append(database[i + offset])
                if not should_change(gene[i], database[i + offset]):
                    sequence.append('*')
                else:
                    sequence.append('.')
        elif match_algorithm == MatchAlgorithm.blat:
            flag, pos_data_end = compute_blat_similarity(gene, database, offset)
            pos_data = offset
            pos_gene = 0
            while pos_gene < 4:
                if should_change(gene[pos_gene], database[pos_data]) > 0:
                    sequence_gene.append('-')
                    sequence_target.append(database[pos_data])
                    sequence.append('.')
                    pos_data += 1
                else:
                    sequence_gene.append(gene[pos_gene])
                    sequence_target.append(database[pos_data])
                    sequence.append('*')
                    pos_gene += 1
                    pos_data += 1
            rev_pos_gene = 7
            rev_pos_data = pos_data_end - 1
            rev_sequence_gene = []
            rev_sequence_target = []
            rev_sequence = []
            while rev_pos_gene > 3:
                if should_change(gene[rev_pos_gene], database[rev_pos_data]) > 0:
                    rev_sequence_gene.append('-')
                    rev_sequence_target.append(database[rev_pos_data])
                    rev_sequence.append('.')
                    rev_pos_data -= 1
                else:
                    rev_sequence_gene.append(gene[rev_pos_gene])
                    rev_sequence_target.append(database[rev_pos_data])
                    rev_sequence.append('*')
                    rev_pos_gene -= 1
                    rev_pos_data -= 1
            while pos_data <= rev_pos_data:
                sequence_gene.append('-')
                sequence_target.append(database[pos_data])
                sequence.append('.')
                pos_data += 1
            sequence_gene.extend(rev_sequence_gene[::-1])
            sequence_target.extend(rev_sequence_target[::-1])
            sequence.extend(rev_sequence[::-1])
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


def count_similarity(weighted, gene, database, offset, is_reverse, max_patience=2, match_pattern=None,
                     continuous_mismatch_limit=None,
                     gene_name_filter=None):
    tot = len(gene)
    weighted_similarity = 0.0
    similarity = {}
    for match_algorithm, weight in zip(
            MatchAlgorithm.get_all_items(),
            weighted):
        data_base_end = offset + len(gene)
        if weight == 0:
            score = 0.0
        elif match_algorithm == MatchAlgorithm.text_distance:
            score, dp = compute_text_distance_similarity(gene, database, offset, continuous_mismatch_limit)
        elif match_algorithm == MatchAlgorithm.direct_match:
            score = 0.0
            for i in range(tot):
                score += should_change(gene[i], database[i + offset])
            score = tot - score
        elif match_algorithm == MatchAlgorithm.consistency:
            score, score_queue, score_merge_idx = compute_consistency_similarity(gene, database, offset, max_patience)
        elif match_algorithm == MatchAlgorithm.pattern:
            score = compute_pattern_similarity(gene, database, offset, match_pattern)
        elif match_algorithm == MatchAlgorithm.blat:
            flag, data_base_end = compute_blat_similarity(gene, database, offset)
            if flag:
                score = 50 - (data_base_end - offset)
            else:
                score = 0
        if gene_name_filter and score > 0:
            if is_reverse:
                start = len(database) - offset
                end = len(database) - data_base_end + 1
            else:
                start = offset + 1
                end = data_base_end

            data = {
                'start': start,
                'end': end,
                'location_result': []
            }
            in_target = gene_name_filter.process_one_data(data)
            if not in_target:
                score = 0
        similarity[match_algorithm] = score
        weighted_similarity += score * weight
    weighted_similarity = weighted_similarity / sum(weighted)
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
