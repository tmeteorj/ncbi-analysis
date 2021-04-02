from analysis.similarities.base_similarity import BaseSimilarity

INF = 999999


class TextEditSimilarity(BaseSimilarity):
    def __init__(self, continuous_mismatch_limit: int = None):
        self.continuous_mismatch_limit = continuous_mismatch_limit

    def get_similarity(self, gene: str, database: str, offset: int):
        tot = len(gene)
        dp = [[INF for _ in range(tot + 1)] for _ in range(tot + 1)]
        dp[0][0] = 0
        for i in range(1, tot + 1):
            gene_a = gene[i - 1]
            min_step = INF
            for j in range(1, tot + 1):
                gene_b = database[offset + j - 1]
                dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1,
                               dp[i - 1][j - 1] + self.should_change(gene_a, gene_b))
                min_step = min(dp[i][j] + abs(i - j), min_step)
        score = float(tot - dp[tot][tot])
        if self.continuous_mismatch_limit is not None:
            i, j = tot, tot
            mismatch = 0
            while i > 0 or j > 0:
                gene_a, gene_b = gene[i - 1] if i > 0 else '.', database[j + offset - 1] if j > 0 else '.'
                if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + self.should_change(gene_a,
                                                                                         gene_b):
                    if self.should_change(gene[i - 1], database[j + offset - 1]) != 0:
                        mismatch += 1
                    else:
                        mismatch = 0
                    i, j = i - 1, j - 1
                elif dp[i][j] == dp[i - 1][j] + 1:
                    mismatch += 1
                    i -= 1
                elif dp[i][j] == dp[i][j - 1] + 1:
                    mismatch += 1
                    j -= 1
                else:
                    raise ValueError('Should not go here!')
                if mismatch >= self.continuous_mismatch_limit:
                    return 0, dp
        return score, dp
