# A global survey of Salmonella plasmids and their associations with antimicrobial resistance
Plasmids are the primary vector for horizontal transfer of antimicrobial resistance (AMR) within bacterial populations, so it is critical to understand their population level distributions and associations with known AMR genes. Salmonella is a routinely sequenced bacterial pathogen in public health for surveillance and outbreak detection using Illumina-based whole genome sequencing (WGS). We applied the MOB-suite, a toolset for reconstructing and typing plasmids, to 150,767 publicly available Salmonella WGS samples covering 1,204 distinct serovars to produce a large scale population survey of plasmids based on the MOB-suite plasmid nomenclature.
<br>
<br>
This repository contains the analysis scripts used to analyse output files from MOB-suite and Abricate, along with metadata to reproduce 

## Requires
+ pandas~=1.3.4
+ plotly~=5.1.0
+ scipy~=1.7.1
+ scikit-bio~=0.5.7
+ scikit-learn~=1.0.1


##Prepare input files for use in visualization tools
Due to size limitations in github, the input files can be downloaded from zenodo here: https://doi.org/10.5281/zenodo.6617143

### Analysing tool outputs
For system-wide installation one can follow these commands on Ubuntu distro that includes Python
library dependencies and tools
```bash
python perform_analysis.py --abricate ./data/abricate-ncbi-amr-genes.txt.gz --contigs /data/salmonella-contig-reports.txt.gz --mobtyper ./data/2020-11-Salmonella-mobtyper_results.txt.gz --metadata ./data/metadata.txt.gz --outdir ./results
```

##Create visualizations
python create_plots.py
```bash
python perform_analysis.py --genes ./results/genes.summary.txt --serovar ./results/serovar.summary.txt --plasmid ./results/plasmid.summary.txt --outdir /Users/jrobertson/Desktop/results/plots --vignette ./results/specific.feature.data.txt"
```

## Output files
| file | Description |
| ------------ | ------------ |
| genes.summary.txt | Summarizes the counts of resistance genes and associations with serovars and MOB-clusters |
| plasmid.summary.txt | Summarizes different aspects of MOB-clusters including serovar, other metadata and resistance genes |
| serovar.summary.txt| Sumamrizes serovar information and associations with MOB-clusters and AMR |
| specific.feature.data.txt | Based on user input creates a hierarchical data structure of serovar and primary and secondary MOB-clusters  |


## genes.summary.txt field descriptions
| field  | Description |
| --------- |  --------- | 
| gene_id | Gene identifier |
| plasmid | Number of occurances on a plasmid contig |
| chromosome | Number of occurances on a chromosome contig |
| total | Total number of occurances |
| proportion | Number of plasmid occurances / Chromosome occurances |
| plasmid homogeneity | Homogenty score of MOB-clusters |
| serovar homogeneity | Homogenty score of Serovars |
| plasmid entropy | MOB-cluster entropy |
| serovar entropy | Serovar entropy |
| human proportion | Proportion of human associated samples |
| resistance | Resistance gene category |

## plasmid.summary.txt field descriptions
| field  | Description |
| --------- |  --------- | 
| plasmid_id | MOB-cluster primary ID |
| replicons | Counts of individual replicons associated with MOB-cluster |
| relaxases | Counts of individual relaxases associated with MOB-cluster |
| overall_mobility | Consensus mobility classification |
| mobility | Counts of the individual mobility classes |
| overall_serovar | Consensus serovar associated with MOB-cluster |
| serovar | Counts of individual serovars |
| continent | MOB-cluster entropy |
| country | Counts of individual continents |
| primary_sample_category | counts of primary sample categories |
| secondary_sample_category | counts of secondary sample categories |
| associated_taxa | Counts of taxa associated with the sample |
| earliest_year | earliest data for MOB-cluster |
| year | counts of individual years |
| samples | sample_ids with MOB-cluster |
| total_samples | number of samples positive for cluster |
| num_resistant | number of MOB-cluster members with >=1 resistance gene |
| proportion_resistant | num_resistant / total_samples |
| resistance_genes | counts of resistance genes found  |
| serovar_entropy | serovar entropy of MOB-cluster |
| serovar_shannon_index | serovar_shannon_index of MOB-cluster |
| serovar_simpson_index | serovar_simpson_index of MOB-cluster |
| serovar_simpson_index_e | serovar_simpson_index_e of MOB-cluster |
| serovar_chao1 | serovar_chao1 of MOB-cluster |
| num_serovars | total number of serovars found with MOB-cluster |
| poportion_human | Proportion of plasmid members associated with human isolates |
| taxa_entropy | taxa_entropy of MOB-cluster |
| taxa_shannon_index | taxa_shannon_index of MOB-cluster |
| taxa_simpson_index | taxa_simpson_index of MOB-cluster |
| taxa_simpson_index_e | taxa_simpson_index_e of MOB-cluster |
| taxa_chao1 | taxa_chao1 of MOB-cluster |

## serovar.summary.txt field descriptions
| field  | Description |
| --------- |  --------- | 
| serovar | Name of serovar |
| continent | Count of samples in each countinent |
| country | Count of samples in each country |
| primary_sample_category | Count of samples with primary sample category |
| secondary_sample_category | Count of samples with secondary sample category |
| associated_taxa | Count of samples with each taxa |
| year | Counts of years for each sample |
| samples | sample_ids with serovar |
| total_samples | total number of samples assigned to serovar |
| count_plasmid_negative_samples | Count of plasmid negative samples |
| count_plasmid_positive_samples | Count of plasmid positive samples |
| plasmids | Count of samples with specigic MOB-clusters |
| num_resistant | Count of samples with >= 1 resistance gene |
| proportion_resistant | num_resistant / total_samples |
| resistance_genes | Count of specific resistance genes |
| poportion_human | proportion of samples with a human label |
| taxa_shannon_index | shannon_index |
| taxa_simpson_index | simpson_index |
| taxa_simpson_index_e | simpson_index_e |
| taxa_chao1 | chao1 |
| taxa_entropy | entropy |
| plasmid_entropy | entropy of MOB-clusters |
| plasmid_shannon_index | shannon_index of MOB-clusters |
| plasmid_simpson_index | simpson_index of MOB-clusters |
| plasmid_simpson_index_e | simpson_index_e of MOB-clusters |
| plasmid_chao1 | chao1 of MOB-clusters |

## specific.feature.data.txt field descriptions
| field  | Description |
| --------- |  --------- | 
| serovar | Name of serovar |
| molecule_type | Chromosome or Plasmid |
| plasmid_id | Primary plasmid cluster |
| count | Count of members |