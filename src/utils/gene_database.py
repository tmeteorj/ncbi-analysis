from typing import List


class GeneSegment:
    gene: str  # name
    left: int
    right: int


class GeneDatabase:
    gene_segments: List[GeneSegment]

    def find_first_greater_equal(self, pos):
        start, end = 0, len(self.gene_segments) - 1
        while start < end:
            mid = (start + end) // 2
            if self.gene_segments[mid].left < pos:
                start = mid + 1
            elif self.gene_segments[mid].left >= pos:
                end = mid
        if self.gene_segments[end].left >= pos:
            return end
        else:
            return end + 1
