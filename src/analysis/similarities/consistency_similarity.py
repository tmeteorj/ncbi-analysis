from analysis.similarities.base_similarity import BaseSimilarity


class ConsistencySimilarity(BaseSimilarity):
    def __init__(self, max_patience: int):
        self.max_patience = max_patience

    def get_similarity(self, gene: str, database: str, offset: int):
        tot = len(gene)
        score = 0.0
        same = 0
        cur_score = 0
        score_queue = []
        for i in range(tot):
            if not self.should_change(gene[i], database[i + offset]):
                cur_score += 1
                same += 1
                if i == tot - 1:
                    score_queue.append([cur_score, tot])
            else:
                score_queue.append([cur_score, i])
                cur_score = 0
            score = max(score, cur_score)
        score_merge_idx = [-1, -1]
        for idx in range(len(score_queue)):
            left = score_queue[idx][1] - score_queue[idx][0]
            total_score = 0
            for width in range(self.max_patience + 1):
                if width + idx < len(score_queue):
                    total_len = score_queue[idx + width][1] - left
                    total_score += score_queue[idx + width][0]
                    if total_len - total_score > self.max_patience:
                        break
                    if score < total_score:
                        score = total_score
                        score_merge_idx = [idx, idx + width]
        return score, (score_queue, score_merge_idx)
