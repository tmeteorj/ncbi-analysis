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

## Gene Similarity Match
- Gene similarity analysis will get similarity gene sequence from a database like NC_000913.
- We provide three ways to compute the similarity as follows.
  - text_distance: means the minimal steps to convert one sequence to another sequence. The less step the better.
  - char_match: compare the sequence one by one, and count the number of same pairs. The more same pair the better.
  - consistency: compare the sequence one by one,  but allowed some of pairs not equal(called patience). The matched sequence who has longer length will get a higher score.
- When run the similarity analysis, you need to provide following parameters.
  - gene_path: the based gene you want to found, with tsv format and contains header sequence
  - data_path: the database you want to find sequence from. 
  - output_directory: where the result should be restored.
  - top_k: find best top_k number of sequence
  - scalar: not much important, to help the computing of similarity , default 1000
  - match_algorithm: match algorithm, which should comes from text_distance, char_match and consistency.
  - candidate_distance: the matched sequences' distance should at least be this, to avoid top similar sequence always come from a nearby of some position. 
  - batch_size: multi thread parameter, to find batch_size of sequence in the same time.
  - min_similarity: only store sequence with higher than this similarity.
  - patience: for consistency match, the maximal number of different steps allowed to ignore. 
- The output will include the top_k result of matched sequence in *_match_result.txt.

## Gene Location Analysis
- Gene location analysis based on the result of similarity match. It will generate more detail information like the position, promoter, product of matched sequence.
- The input include three parameters as follows.
  - input_file_path: means the match result of gene similarity match
  - ecocyc_file_path: this file contains gene information like promoter, product, start/end position, etc. And it is generated from ecocyc_analysis component.
  - output_directory: where the result should be restored.
- The output file include two parts.
  - *_location_result.txt: the detail information of whole matched sequence
  - *_sub_location_result.txt: the maximal consistency sub sequence's detail information
- The function format_data_to_tsv in the same file can convert location file into tsv format for better understanding.

## Location Reorder
- Location reorder is to reorder the location file from Gene Location Analysis. 
- You can give a list of orders, and this component will reorder the match information based on the order
- The input file contains the order list and location result file, while the output is what above described. 