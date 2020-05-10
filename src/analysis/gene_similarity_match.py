import heapq
import os
import random
import re
import threading
import time

from utils.factories.logger_factory import LoggerFactory
from utils.gene_file_util import GeneFileReader
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter
from collections import deque


class GeneSimilarityMatch:
    def __init__(self, gene_path, data_path, output_directory, top_k=20, precision=1000,
                 match_algorithm='text_distance', candidate_distance=5, batch_size=5, min_similarity=0.3):
        self.data_path = data_path
        self.output_directory = output_directory
        self.data_name = os.path.basename(data_path)
        file_name = os.path.basename(gene_path)
        file_prefix = StrConverter.extract_file_name(file_name)
        self.result_path = os.path.join(self.output_directory,
                                        '%s_%s_match_result.txt' % (file_prefix, match_algorithm))
        self.gene_reader = GeneFileReader(self.data_path)
        self.gene_path = gene_path
        self.top_k = top_k
        self.precision = precision
        self.match_algorithm = match_algorithm
        self.candidate_distance = candidate_distance
        self.min_similarity = min_similarity
        self.batch_size = batch_size
        self.logger = LoggerFactory()

        self.lock = threading.Lock()
        self.solved = 0
        self.total = 0

    def run(self):
        self.gene_reader.build_information()
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
                pat = None if len(items) == 2 else items[2].lower()
                t = threading.Thread(target=self.find_candidate_for_gene, args=(name, gene, pat, fw,))
                pending_tasks.append(t)
            while len(pending_tasks) > 0:
                running_tasks = [t for t in running_tasks if t.isAlive()]
                while len(running_tasks) < self.batch_size and len(pending_tasks) > 0:
                    t = pending_tasks.popleft()
                    t.start()
                    running_tasks.append(t)
                time.sleep(60)
            for t in running_tasks:
                t.join()

    def find_candidate_for_gene(self, name, gene, pat, fw):
        candidates = [[], []]
        t1 = threading.Thread(target=self.match_gene,
                              args=(name, gene, self.gene_reader.dna_code, False, candidates[0], pat,))
        t1.start()
        rev_dna_code = get_opposite_dna(self.gene_reader.dna_code[::-1])
        t2 = threading.Thread(target=self.match_gene, args=(name, gene, rev_dna_code, True, candidates[1], pat,))
        t2.start()
        t1.join()
        t2.join()
        candidates = candidates[0] + candidates[1]
        candidates.sort(key=lambda arg: -arg.similarity)
        self.lock.acquire()
        for candidate in candidates[:self.top_k]:
            fw.write('>%s/%s-%s\tname=%s,similarity=%.2f%%,direction=%s\n' % (
                self.data_name.replace(".txt", ''),
                candidate.start,
                candidate.end,
                name,
                candidate.similarity,
                '-' if candidate.is_reverse else '+'
            ))
        self.lock.release()

    def match_gene(self, name, gene, database, is_reverse, candidates, pat):
        gene_length = len(gene)
        min_same = int(self.min_similarity * gene_length) - 1
        database_length = len(database)
        limitation = database_length - gene_length + 1
        new_solved = 0
        similarity_heap = []
        # gene_dict = count_acgt(gene)
        buff = deque()
        for start in range(limitation):
            # if fast_skip(gene_dict, gene_length, database, start, min_same, pat):
            #    continue
            similarity = count_similarity(self.match_algorithm, self.precision, gene, database, start, min_same)
            new_solved += 1
            if random.random() * 1000 < 1:
                self.lock.acquire()
                self.solved += new_solved
                self.logger.info_with_expire_time(
                    'Doing Similarity Matching for %s[%s]: %d/%d(%.2f%%) --top_k=%d --min_similarity= %.2f%% --min_same=%d --gene_length=%d --candidates_num=%d' % (
                        name,
                        '-' if is_reverse else '+',
                        self.solved,
                        self.total,
                        self.solved * 100.0 / self.total,
                        self.top_k,
                        similarity_heap[0] if len(similarity_heap) > 0 else 0.00,
                        min_same,
                        gene_length,
                        len(candidates)
                    ),
                    self.solved,
                    self.total)
                self.lock.release()
                new_solved = 0
            if similarity < 0:
                continue
            new_candidate = MatchCandidate(start, start + gene_length - 1, is_reverse, database,
                                           similarity * 100.0 / self.precision)
            added_flag = update_candidate_list(new_candidate,
                                               buff,
                                               candidates,
                                               self.candidate_distance)
            if added_flag:
                heapq.heappush(similarity_heap, candidates[-1].similarity)
                if len(similarity_heap) > self.top_k:
                    heapq.heappop(similarity_heap)
                    top = similarity_heap[0]
                    min_same = max(min_same, int(top / 100.0 * gene_length) - 1)
        while len(buff) > 0:
            update_candidate_list(None, buff, candidates, 0)
        self.lock.acquire()
        self.solved += new_solved + gene_length - 1
        self.lock.release()
        return min_same


def update_candidate_list(new_candidate, buff: deque, candidate_result: list,
                          keep_size: int):
    added = False
    while (len(buff) >= keep_size) or (
            len(buff) > 0 and new_candidate is not None and abs(buff[0].start - new_candidate.start) >= keep_size):
        old_candidate = buff.popleft()
        if not old_candidate.should_ignore:
            candidate_result.append(old_candidate)
            added = True
    if new_candidate is not None:
        for candidate in buff:
            if candidate.similarity > new_candidate.similarity:
                new_candidate.should_ignore = True
            elif candidate.similarity < new_candidate.similarity:
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


def count_similarity(match_algorithm, scalar, gene, database, offset, min_same):
    tot = len(gene)
    if match_algorithm == 'text_distance':
        max_step = tot - min_same
        dp = [[99999 for _ in range(tot + 1)] for _ in range(tot + 1)]
        dp[0][0] = 0
        for i in range(1, tot + 1):
            gene_a = gene[i - 1]
            min_step = 1000000
            for j in range(1, tot + 1):
                gene_b = database[j - 1 + offset]
                dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + (0 if gene_a == gene_b else 1))
                min_step = min(dp[i][j] + abs(i - j), min_step)
            if min_step >= max_step:
                return -1
        score = tot - dp[tot][tot]
    else:
        score = 0.0
        for i in range(tot):
            if gene[i] == database[i + offset]:
                score += 1.0
    return int(score * scalar / tot)


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


class MatchCandidate:
    def __init__(self, left, right, is_reverse, database, similarity):
        self.is_reverse = is_reverse
        if is_reverse:
            self.start = len(database) - left
            self.end = len(database) - right
        else:
            self.start = left + 1
            self.end = right + 1
        self.similarity = similarity
        self.should_ignore = False
