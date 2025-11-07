# File Prep

The point of file prep is to go from a compressed fastq file, named as a barcode from Nanopore demultiplexing, into a mapped bam file that is ready for downstream processing.

## submit_file_prep.sh

This Bash script is a general-purpose SLURM job array submitter to launch a batch of jobs on an HPC cluster. This script does not perform the actual data processing itself.

Its steps are:
1. Parse Configuration: It reads all settings from command-line arguments. These include SLURM resource requests (CPU, memory, time) and the paths required by the pipeline script it will run (e.g., input/output directories, sample map).
2. Validate Inputs: It ensures all required arguments are provided and that the specified script, files, and directories exist and are accessible.
3. Submit Job Array: It counts the number of lines in the sample map file to determine the required number of jobs. It then uses the sbatch command to submit a job array to SLURM, where each task in the array is configured with the user-specified resources.
4. Launch Pipeline: Each job in the submitted array will execute the user-provided pipeline script. This submitter script passes the necessary file paths and settings on to that pipeline script, which then processes one sample per job.

One current limitation of this script is that it can only submit one primary executeable script at a time. To link together other scripts, they must be called in a single bash file.


## trim_map_dedup.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for amplicons where there is only one set of adapters, the Nanopore adapters, to trim from the reads.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 and converts the resulting SAM output into a sorted BAM file with samtools.
6. Deduplicate Reads: It uses umi_tools dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
7. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
8. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log.

## trim_map_dedup_mpra.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for MPRA where there are two sets of adapters, the Nanopore adapters and the MPRA adapters, to trim from the reads. Nanopore adapters are trimmed prior to UMI extraction, MPRA adapters are trimmed following UMI extraction.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the Nanopore adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Extract UMIs: It uses umi_tools extract to pull 10-nucleotide Unique Molecular Identifiers (UMIs) from the 3' end of each read, embedding the UMI into the read's header for later use.
5. Second Trim: Performs a second cutadapt pass on the UMI-tagged reads to remove adapter sequences for the MPRA.
6. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 and converts the resulting SAM output into a sorted BAM file with samtools.
7. Deduplicate Reads: It uses umi_tools dedup to identify and remove PCR duplicates from the sorted BAM file, using the UMI and alignment position of each read.
8. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
9. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log.

## trim_map_mpra.sh

This script processes a single Nano-BID-Amp sample as part of a SLURM array job. This script is appropriate for MPRA with the MPRA adapters, but where technical issues during library creation resulted in non-functional UMIs. Deduplication does not occur in this script.

Its steps are:
1. Setup and Initialization: It parses command-line arguments, uses the SLURM array task ID to identify a specific sample from a map file, sets up output directories, and redirects all output to a log file.
2. Filter Reads: It uses cutadapt with linked adapters to process the compressed FASTQ file twice, creating two separate files: one containing "sense" reads and another with "antisense" reads, based on the MPRA adapter sequences.
3. Unify Read Orientation: It reverse-complements the "antisense" reads and then concatenates them with the "sense" reads. This creates a single FASTQ file where all reads are oriented in the sense direction.
4. Align Reads: It maps the UMI-tagged reads to a reference genome using minimap2 and converts the resulting SAM output into a sorted BAM file with samtools.
5. Finalize and Clean Up: It copies the final, deduplicated BAM file to a common output directory, verifies the copy's integrity with a checksum, indexes the file, and removes temporary intermediate files.
6. Reporting: Throughout the process, it tracks read counts and file sizes at each major step, printing a final summary table to the log.