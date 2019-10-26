# ncbi-analysis
This project is for analysis of DNA/RNA based on ncbi database

## Cluster Match
- Cluter match is to find the same dna sequence in the fna file and make them a cluster.
- When run a cluster match, you need to provide following parameters
  - rna_tag: the rna name in each element, like 7a_utr, 8a_utr
  - input_path: the fna file path
  - output_directory: the output directory path
- The output of cluster match contains 4 files as follows:
  - *_cluster_result.txt: the clusters of each element, the first column is the number of element, the second column is the group of elements separate by ','
  - *_sample_result.txt: the sample element in each cluster
  - *_all_result.txt: all elements will be list with format like following code
  ```
  >NameOfElement/StartPosition-EndPosition
  GeneStream
  ```
  - *_only_result.txt: similar with *_all_result.txt, but only keeps the one which only contains AUCG.

## Neighbor Analysis
- Neighbor analysis consumes result(*_all_result.txt or *_only_result.txt) from cluster match. For each element, it will find neighbor gene for the element. And also get the percentages of each next neighbor's gene or source.
- When run neighbor analysis, you need to provide following parameters:
  - input_path: result file from cluster match
  - download_directory: since this step will download the gene file from ncbi, you need to provide a directory for it
  - output_directory: the output directory path
  - keep_prefix_num: default 1, when name of source is too long, we can keep only few prefix words for clearly count
- The output of neighbor analysis contains 4 to 5 files as follows:
  - *_neighbor_result.txt : the neighbor gene for target input
  - *_next_neighbor_result.txt: the next neighbor gene for target input
  - *_source_count_result.txt: the name of source distribution for target input
  - *_gene_count_result.txt: the name of gene distribution for target input
  - *_error_result*: if there is any error, this file will log what happened

## Gene Extract
- For a given gene names and element path, it will extract gene's position and sequence from element path.
- You should provide following parameters:
  - data_path: the element path, which should be downloaded in neighbor analysis
  - rna_path: the file path which contains the names of gene
  - output_path: the result path for gene extract
