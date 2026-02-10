# This directory contains scripts to analyze XBP1 splice site in the scRNA-seq data using Genotyping of the Transcriptomes (GoT)
## The core function employed here is IronThron-GoT: https://github.com/dan-landau/IronThrone-GoT
## IronThrone-GoT quantify wild-type (unspliced XBP1, uXBP1) and mutated transcripts (spliced XBP1, sXBP1) from fastq files
### Order of operations:
(1) Processing of sequencing reads with IronThrone-GoT. See example folder for file structure and config file layout. 
(2) Run the Combine_IronThrone_Parallel_Output.R. Modified from https://github.com/dan-landau/IronThrone-GoT/blob/master/Combine_IronThrone_Parallel_Output.R
(3) Run the analysis.R file to generate plots corresponding to the smoothened sXPB1 uXBP1 with pseudotime. 
