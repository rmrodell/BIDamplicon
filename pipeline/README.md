# File Prep

The point of file prep is to go from a compressed fastq file, named as a barcode from Nanopore demultiplexing, into a mapped bam file that is ready for downstream processing. All scripts are intended to work on Sherlock, Stanford's HPC.

For the sample map file, it is **strongly** recommended to set up you file names to follow the scheme of celltype_vector_rep_treat for integration with the rest of the pipeline.

Example sample map file:

``` bash
IV_noPUS_in_1:barcode66
IV_PUS7_BS_1:barcode67
IV_PUS7_BS_2:barcode68
IV_noPUS_BS_2:barcode70
```

## submit_file_prep.sh

This Bash script is a general-purpose SLURM job array submitter to launch a batch of jobs on an HPC cluster. This script does not perform the actual data processing itself.

Its steps are:
1. Parse Configuration: It reads all settings from command-line arguments. These include SLURM resource requests (CPU, memory, time) and the paths required by the pipeline script it will run (e.g., input/output directories, sample map).
2. Validate Inputs: It ensures all required arguments are provided and that the specified script, files, and directories exist and are accessible.
3. Submit Job Array: It counts the number of lines in the sample map file to determine the required number of jobs. It then uses the sbatch command to submit a job array to SLURM, where each task in the array is configured with the user-specified resources.
4. Launch Pipeline: Each job in the submitted array will execute the user-provided pipeline script. This submitter script passes the necessary file paths and settings on to that pipeline script, which then processes one sample per job.

One current limitation of this script is that it can only submit one primary executeable script at a time. To link together other scripts, they must be called in a single bash file.

Examples of how to use, referencing the scripts below:

```bash
bash $HOME/BIDamplicon/pipeline/submit_file_prep.sh \
    -u rodell@stanford.edu -t 4:00:00 \
    -s /home/users/rodell/BIDamplicon/pipeline/trim_map_dedup_mpra.sh \
    -a /scratch/users/rodell/20251106_pool1/invitro/barcodes_may.txt \
    -i /oak/stanford/groups/nicolemm/mmontes/053025_IVpsi_pool1_enzymeTestBS/RUN1/20250530_1516_P2S-01916-B_PBC44754_1d396595/fastq_pass \
    -o /scratch/users/rodell/20251106_pool1/invitro/maydata \
    -r /home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta

bash $HOME/BIDamplicon/pipeline/submit_file_prep.sh \
    -u rodell@stanford.edu -t 4:00:00 \
    -s /home/users/rodell/BIDamplicon/pipeline/trim_map_mpra.sh \
    -a /scratch/users/rodell/20251106_pool1/incell/barcodes.txt \
    -i /oak/stanford/groups/nicolemm/rodell/BIDamplicon/20241114_pool1/fastq \
    -o /scratch/users/rodell/20251106_pool1/incell/prep \
    -r /home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta
```


## trim_map_dedup.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for amplicons where there is only one set of adapters, the Nanopore adapters, to trim from the reads.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 and converts the resulting SAM output into a sorted BAM file with samtools.
    - TO DO: Update minimap2 command.
6. Deduplicate Reads: It uses umi_tools dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
7. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
8. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log. For BAM files, it reports the number of alignments. For fastq files, it reports the number of reads.

## trim_map_dedup_mpra.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for MPRA where there are two sets of adapters, the Nanopore adapters and the MPRA adapters, to trim from the reads. Nanopore adapters are trimmed prior to UMI extraction, MPRA adapters are trimmed following UMI extraction.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Second Trim: Performs a second cutadapt pass on the UMI-tagged reads to remove adapter sequences for the MPRA.
6. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 with default settings for short genomice reads and converts the resulting SAM output into a sorted BAM file with samtools.
7. Remove Multi-Aligners: Removes any reads that map to multiple locations within the reference, retaining only those that have a single primary alignment and a MAPQ score greater than or equal to 30.
8. Deduplicate Reads: It uses umi_tools dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
9. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
10. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log. For BAM files, it reports the number of alignments. For fastq files, it reports the number of reads.

## trim_bowite2_dedup_mpra.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for MPRA where there are two sets of adapters, the Nanopore adapters and the MPRA adapters, to trim from the reads. Nanopore adapters are trimmed prior to UMI extraction, MPRA adapters are trimmed following UMI extraction.

This script differs from the one above (trim_map_dedup_mpra.sh) through the use of bowtie2 as the aligner rather than minimap2. Preliminary insights indicate this is a better approach for highly degenerate sequences (i.e. mutagenesis pools) to prevent multi-mapping reads, which is a large problem with the minimap2 parameters in the above approach.

Due to the different arguments, this script needs to be submitted with submit_file_prep_bowtie.sh, which works similarly to submit_file_prep.sh.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Second Trim: Performs a second cutadapt pass on the UMI-tagged reads to remove adapter sequences for the MPRA.
6. Align Reads: It maps the UMI-tagged reads to a reference genome using bowtie2 with an end-to-end approach and converts the resulting SAM output into a sorted BAM file with samtools.
7. Deduplicate Reads: It uses umi_tools dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
8. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
9. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log. For BAM files, it reports the number of alignments. For fastq files, it reports the number of reads.

## trim_map_mpra.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for MPRA with the MPRA adapters, but where technical issues during library creation resulted in non-functional UMIs. Deduplication does not occur in this script.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the MPRA adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 and converts the resulting SAM output into a sorted BAM file with samtools.
    - TO DO: Update minimap2 command.
5. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
6. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log. For BAM files, it reports the number of alignments. For fastq files, it reports the number of reads.

KEY DIFFERNECES TO NOTE:
- Does not take a mapping file
- Input fastq files should already be named with the sample name.
- Input fastq files are not gzipped.
This truly is rather custom built for the mess I made of the Nov2024 Pool1 sequencing. Here's to learning and growing.

# Counting

Aggregates the nucleotide counts for all samples All scripts are intended to work on Sherlock, Stanford's HPC.

## BIDdetect.sh

This script automates the process of calculating and aggregating nucleotide counts from multiple bam files.

Required Inputs:
- -b, --bam_dir: A directory containing one or more sorted and indexed BAM files (.bam).
- -o, --output_dir: A path to a directory where all results and logs will be stored.
- -r, --ref_fasta: The path to the reference genome FASTA file.
- -e, --bed_file: The path to a BED file defining the genomic regions of interest.

Optional Input:
- -n, --col_names: A string defining column names for the final formatting step (Default: 'celltype_vector_rep_treat').

Key Outputs:

- Final Data File (${DEST_DIR}/BIDdetect_data.txt): The main result of the pipeline, containing the fully processed and formatted count data after all steps. This should be taken into further analyses.
- Aggregated Counts File (${DEST_DIR}/BIDdetect_counts.txt): A master table containing the raw, combined counts from all input BAM files before the final formatting step.
- Log File (${DEST_DIR}/logs/BIDdetect_...log): A timestamped log file that captures all screen output for debugging and record-keeping.
- Intermediate Files (${DEST_DIR}/intermediate_counts/): A directory holding the temporary count files generated for each individual sample before they are aggregated.

Workflow:

1. Argument Parsing & Validation: makes sure all required variables are present and valid
2. Environment Setup: loads appropriate modules, sets up directories, creates logging system
3. Main Processing Loop (Per-Sample Analysis): iterates through every (indexed) .bam file in the input directory
    1. Extracts a clean sample name from the filename. 
    2. Executes the bam_counts_fast.R script that performs a pileup on the BAM file at the positions defined in the BED file to generate nucleotide, deletion, and insertion counts. Calculates the deletion rate (deletions / total reads) for each site. 
    3. Saves the per-sample results to a temporary file in the intermediate_counts directory. 
    4. Appends the results from the temporary file to a single master file (BIDdetect_counts.txt), adding a new column containing the sample name.
4. Variable Extraction from File Name: calls sample_name.R to split the file name into 4 different columns, as defined in the --col_names argument.
    - default is celltype_vector_rep_treat, which is the recommended naming of bam files to ensure integration with analysis scripts
5. Final Processing & Cleanup: saves final data to BIDdetect_data.txt, outputs locations of final files

Use example:

```bash
bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
  --bam_dir /scratch/users/rodell/20251106_pool1/invitro/bam_links \
  --output_dir /scratch/users/rodell/20251106_pool1/invitro/counts_delpos \
  --ref_fasta /home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta \
  --bed_file /home/groups/nicolemm/rodell/pool1/pool1_cleaned_delpos_noadapters.bed \
  --col_names 'celltype_vector_treat_rep'

bash $HOME/BIDamplicon/pipeline/BIDdetect.sh \
  --bam_dir /scratch/users/rodell/20251106_pool1/invitro/bam_links \
  --output_dir /scratch/users/rodell/20251106_pool1/invitro/counts_full \
  --ref_fasta /home/groups/nicolemm/rodell/pool1/pool1_cleaned_noadapters.fasta \
  --bed_file /home/groups/nicolemm/rodell/pool1/pool1_cleaned_full_noadapters.bed \
  --col_names 'celltype_vector_treat_rep'
```

# Analysis

Perform analysis for statistical difference and equivalence of sites. At least **two replicates** are required for the statistical analysis to work properly, and more is of course better. 

## modification_analysis.R

This script performs statistical analysis on processed Nano-BID-Amp data to identify significantly modified sites. It takes a counts table (typically generated by BIDdetect.sh) with deletion rates as input and applies a mixed-effects model and equivalence testing to categorize each site.

Required Inputs:
- -i, --input: Path to the input data table containing per-replicate deletion rates.
- -o, --outdir: Parent directory where a new, prefixed output folder will be created.
- -p, --prefix: A unique name for the analysis run, used for the output folder and all file names.

Optional Inputs:
- --cores: Number of CPU cores to use for parallel processing. (Default: -1, all available cores).
- --sesoi: The Smallest Effect Size of Interest (SESOI) for equivalence testing. (Default: 0.05).
- --color: The hex code for coloring "modified" sites in plots.
- --plot_all_sites: A flag to generate a large PDF containing a barplot for every individual site. Include the flag to generate this plot.

Key Outputs:
- Final Results Table (data_summary/<prefix>_modification_significance.tsv): Contains the final category ('Modified', 'Unmodified', 'Inconclusive') and all statistical results for every site.
- Data Subsets (data_raw/ and data_summary/): Includes raw and summary tables for all sites in each category, as well as for quartiles of the modified sites.
- Plots (plots/): A folder containing summary plots, including a boxplot of average deletion rates, heatmaps of modified sites (with and without labels), and an optional PDF of all individual site plots.
- Log File (<prefix>_analysis_log.txt): A detailed log capturing all parameters, progress messages, and summary tables printed to the console.

Workflow:

1. Set Up: Before running, the script verifies that the R version is 4.3.0 or newer is available. If not, it stops with a helpful error message explaining how to load the correct modules on Sherlock. It parses all command-line arguments and creates the nested output directory structure.
2. Data Loading & Validation: It loads the input data table and verifies that all required columns are present. Required columns should be generated by the counting script if used properly. Required columns are:
    - chr, pos: extracted from BED file during counting, indicates the positon to analyze
    - treat: values should be in / input and BS / BID; column the statistics are performed between, **requires at least two data points in each condition** for the analysis to run, otherwise that site will be skipped
    - rep: indicates replicate number, need to have at least two for each treatment condition (BS, input)
    - delrate: deletion rate at a given sites
    - totalReads: total number of reads for a site
    - vector: (important) placeholder, indicates noPUS / PUS or WT / KD, more of a tool we'll use later but still required to have
3. Statistical Analysis (Parallelized): This is the core computational step.
    - It first pre-filters sites where all deletion rates are below the SESOI, categorizing them as 'Unmodified'.
    - It then sets up a parallel processing backend using the number of cores specified by the --cores argument. More cores means faster processing.
    - For each site, it fits a Bayesian mixed-effects model (bglmer) to test for significant differences and performs an equivalence test (TOST) to check for a lack of meaningful change.
4. Site Categorization: It uses the results from the statistical tests to assign a final category ('Modified', 'Unmodified', or 'Inconclusive') to every site.
    - Modified: adjusted p-value indicates a significant difference between BS-treated and input
    - Unmodified: equivalence between sites, or every data point for a given site is less than SESOI
    - Inconclusive: Did not meet statistical significance or equivalence, likely due to high variance in underlying data.
5. Generate Outputs: The script saves the final summary table, creates and saves the data subsets (including a quartile analysis for modified sites), and generates the summary plots (box plot, heat maps).

Use example:

```bash
Rscript $HOME/BIDamplicon/pipeline/modification_analysis.R \
  --input /scratch/users/rodell/20251106_pool1/invitro/counts_delpos/BIDdetect_data_invitro_delpos.txt \
  --outdir /scratch/users/rodell/20251106_pool1/invitro/analysis_delpos \
  --prefix invitro_delpos \
  --plot_all_sites
```