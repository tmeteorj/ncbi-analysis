import os
import threading

from utils.factories.logger_factory import LoggerFactory
from utils.gene_file_util import GeneFileReader
from utils.gene_util import get_opposite_dna
from utils.str_util import StrConverter


class GeneSimilarityMatch:
    def __init__(self, gene_path, data_path, output_directory, top_k=20, precision=1000,
                 match_algorithm='text_distance'):
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
            ts = []
            for gene_sequence in gene_sequences:
                name, gene = gene_sequence.strip().split('\t')
                t = threading.Thread(target=self.find_candidate_for_gene, args=(name, gene, fw,))
                t.start()
                ts.append(t)
            for t in ts:
                t.join()

    def find_candidate_for_gene(self, name, gene, fw):
        gene = gene.strip().lower()
        candidates = []
        self.match_gene(gene, self.gene_reader.dna_code, False, candidates)
        rev_dna_code = get_opposite_dna(self.gene_reader.dna_code[::-1])
        self.match_gene(gene, rev_dna_code, True, candidates)
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

    def match_gene(self, gene, database, is_reverse, candidates):
        gene_length = len(gene)
        database_length = len(database)
        limitation = database_length - gene_length + 1
        upper = {i: 0 for i in range(self.precision + 1)}
        new_solved = 0
        for start in range(limitation):
            similarity = count_similarity(self.match_algorithm, self.precision, gene, database, start)
            new_solved += 1
            if (start + 1) % 1000 == 0:
                self.lock.acquire()
                self.solved += new_solved
                self.logger.info_with_expire_time(
                    'Doing Similarity Matching: %d/%d(%.2f%%)' % (
                        self.solved, self.total, self.solved * 100.0 / self.total), self.solved,
                    self.total)
                self.lock.release()
                new_solved = 0
            have = upper.get(similarity, 0)
            if have > self.top_k:
                continue
            for i in range(similarity):
                upper[i] += 1
            candidates.append(MatchCandidate(start, start + gene_length - 1, is_reverse, database,
                                             similarity * 100.0 / self.precision))
        self.lock.acquire()
        self.solved += new_solved + gene_length - 1
        self.lock.release()


def count_similarity(match_algorithm, scalar, gene, database, offset):
    tot = len(gene)
    if match_algorithm == 'text_distance':
        dp = [[999999 for i in range(tot + 1)] for j in range(tot + 1)]
        dp[0][0] = 0
        for i in range(1, tot + 1):
            gene_a = gene[i - 1]
            for j in range(1, tot + 1):
                gene_b = database[j - 1 + offset]
                match_cost = 1 if gene_a != gene_b else 0
                dp[i][j] = min(dp[i - 1][j - 1] + match_cost, dp[i - 1][j] + 1, dp[i][j - 1] + 1)
        score = tot - dp[tot][tot]
    else:
        score = 0.0
        for i in range(tot):
            if gene[i] == database[i + offset]:
                score += 1.0
    return int(score * scalar / tot)


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
