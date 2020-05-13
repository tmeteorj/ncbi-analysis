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
    def __init__(self, gene_path, data_path, output_directory, top_k=20, scalar=1000,
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
        self.scalar = scalar
        self.match_algorithm = match_algorithm
        self.candidate_distance = candidate_distance
        self.min_similarity = min_similarity
        self.batch_size = batch_size
        self.dna_code = None
        self.rev_dna_code = None
        self.logger = LoggerFactory()

        self.lock = threading.Lock()
        self.solved = 0
        self.total = 0

    def run(self):
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
                              args=(name, gene, self.dna_code, False, candidates[0], pat,))
        t1.start()
        t2 = threading.Thread(target=self.match_gene, args=(name, gene, self.rev_dna_code, True, candidates[1], pat,))
        t2.start()
        t1.join()
        t2.join()
        candidates = candidates[0] + candidates[1]
        candidates.sort(key=lambda arg: -arg.similarity)
        results = self.render_similarity_for_candidates(gene, candidates[:self.top_k])
        self.lock.acquire()
        idx = 1
        headers = ['name', 'direction', 'similarity', 'consistency', 'original', 'gene_format', 'target_format',
                   'match_format']
        for candidate, sequence_gene, sequence_target, sequence in results:
            fw.write('(%d)\n' % idx)
            attribute = {
                'name': name,
                'direction': '-' if candidate.is_reverse else '+'
            }
            if self.match_algorithm == 'consistency':
                attribute['similarity'] = '%.2f%%' % ((candidate.similarity % self.scalar) * 100.0 / len(gene))
                attribute['consistency'] = '%d' % (candidate.similarity // self.scalar)
            else:
                attribute['similarity'] = '%.2f%%' % candidate.similarity
            for key, value in zip(['origin', 'gene_format', 'target_format', 'match_format'],
                                  [gene, sequence_gene, sequence_target, sequence]):
                attribute[key] = value

            fw.write('>%s/%s-%s\t%s\n' % (
                self.data_name.replace(".txt", ''),
                candidate.start,
                candidate.end,
                ','.join(['%s=%s' % (key, attribute[key]) for key in headers if key in attribute])
            ))
            fw.write('\n')
            idx += 1
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
            similarity = count_similarity(self.match_algorithm, self.scalar, gene, database, start, min_same)
            new_solved += 1
            if random.random() * 1000 < 1:
                self.lock.acquire()
                self.solved += new_solved
                self.logger.info_with_expire_time(
                    'Doing Similarity Matching for %s[%s]: %d/%d(%.2f%%) --top_k=%d --min_score=%d --min_same=%d --gene_length=%d --candidates_num=%d' % (
                        name,
                        '-' if is_reverse else '+',
                        self.solved,
                        self.total,
                        self.solved * 100.0 / self.total,
                        self.top_k,
                        similarity_heap[0] if len(similarity_heap) > 0 else 0,
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
            if self.match_algorithm == 'text_distance':
                new_candidate = MatchCandidate(start, start + gene_length - 1, is_reverse, database,
                                               similarity * 100.0 / self.scalar)
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
            elif self.match_algorithm == 'consistency':
                if len(similarity_heap) == self.top_k and similarity < similarity_heap[0]:
                    continue
                heapq.heappush(similarity_heap, similarity)
                if len(similarity_heap) > self.top_k:
                    heapq.heappop(similarity_heap)
                new_candidate = MatchCandidate(start, start + gene_length - 1, is_reverse, database,
                                               similarity)
                candidates.append(new_candidate)
            else:
                raise ValueError('match algorithm not found')
        while len(buff) > 0:
            update_candidate_list(None, buff, candidates, 1)
        self.lock.acquire()
        self.solved += new_solved + gene_length - 1
        self.lock.release()
        return min_same

    def render_similarity_for_candidates(self, gene, candidates):
        result = []
        for candidate in candidates:
            database = self.rev_dna_code if candidate.is_reverse else self.dna_code
            if self.match_algorithm == 'text_distance':
                similarity, dp = count_similarity('text_distance',
                                                  self.scalar,
                                                  gene,
                                                  database,
                                                  candidate.original_match_left,
                                                  min_same=0,
                                                  return_dp=True)
                sequence_gene, sequence_target, sequence = render_dna_sequence(gene, database,
                                                                               candidate.original_match_left, dp)
                result.append([candidate, ''.join(sequence_gene), ''.join(sequence_target), ''.join(sequence)])
            elif self.match_algorithm == 'consistency':
                sequence_gene, sequence_target, sequence = render_dna_sequence(gene, database,
                                                                               candidate.original_match_left)
                result.append([candidate, ''.join(sequence_gene), ''.join(sequence_target), ''.join(sequence)])
        return result


def render_dna_sequence(gene, database, offset, dp=None):
    sequence_gene = []
    sequence_target = []
    sequence = []
    tot = len(gene)
    if dp is None:
        for i in range(tot):
            sequence_gene.append(gene[i])
            sequence_target.append(database[i + offset])
            if not should_change(gene[i], database[i + offset]):
                sequence.append('*')
            else:
                sequence.append('.')
    else:
        i, j = tot, tot
        while i > 0 or j > 0:
            gene_a, gene_b = gene[i - 1] if i > 0 else '.', database[j + offset - 1] if j > 0 else '.'
            if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + should_change(gene[i - 1], database[j + offset - 1]):
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
    return sequence_gene, sequence_target, sequence


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


def count_similarity(match_algorithm, scalar, gene, database, offset, min_same, return_dp=False):
    tot = len(gene)
    if match_algorithm == 'text_distance':
        max_step = tot - min_same
        dp = [[99999 for _ in range(tot + 1)] for _ in range(tot + 1)]
        pre = [[None for _ in range(tot + 1)] for _ in range(tot + 1)]
        dp[0][0] = 0
        for i in range(1, tot + 1):
            gene_a = gene[i - 1]
            min_step = 1000000
            for j in range(1, tot + 1):
                gene_b = database[j - 1 + offset]
                dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + should_change(gene_a, gene_b))
                min_step = min(dp[i][j] + abs(i - j), min_step)
            if min_step >= max_step:
                return -1
        score = tot - dp[tot][tot]
    elif match_algorithm == 'direct_match':
        score = 0.0
        for i in range(tot):
            score += should_change(gene[i], database[i + offset])
    elif match_algorithm == 'consistency':
        score = 0
        same = 0
        cur = 0
        for i in range(tot):
            if not should_change(gene[i], database[i + offset]):
                cur += 1
                same += 1
            else:
                cur = 0
            score = max(score, cur)
        return score * scalar + same
    if return_dp:
        return int(score * scalar / tot), dp
    return int(score * scalar / tot)


def should_change(a, b):
    if a == b:
        return 0
    elif a == 'c' and b == 't':
        return 0
    return 1


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
        self.original_match_left = left
        self.original_match_right = right
