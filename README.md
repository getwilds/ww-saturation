# ww-saturation: Saturation Mutagenesis Workflow

## Overview

This repository contains a WDL (Workflow Description Language) pipeline for performing saturation mutagenesis analysis on RNA-seq data. The workflow aligns paired-end sequencing reads to a reference genome using BWA, converts the resulting SAM files to BAM format, and performs saturation mutagenesis analysis using GATK's AnalyzeSaturationMutagenesis tool.

## Workflow Steps

1. **Alignment (BWA)**: Paired-end reads are aligned to a reference genome
2. **SAM to BAM Conversion**: The aligned reads are converted from SAM to BAM format and indexed
3. **Saturation Mutagenesis Analysis**: GATK's AnalyzeSaturationMutagenesis tool analyzes mutations across the specified ORF range

## Files in this Repository

- `ww-saturation.wdl`: The main workflow definition file
- `ww-saturation-inputs.json`: Example input parameters for the workflow
- `ww-saturation-options.json`: Runtime configuration options for Cromwell

## Requirements

- [Cromwell](https://github.com/broadinstitute/cromwell) or another WDL-compatible workflow execution engine
- Docker (all tools run in containers)
- Required Docker images:
  - getwilds/bwa:0.7.17
  - getwilds/samtools:1.11
  - getwilds/gatk:4.3.0.0

## Usage

1. **Prepare your input files**:
   - Paired-end FASTQ files
   - Reference genome (FASTA, index, and dictionary files)
   - Define your ORF range

2. **Modify the inputs JSON file**:
   ```json
   {
     "saturation_mutagenesis.fastq_1": "/path/to/your/reads_R1.fastq.gz",
     "saturation_mutagenesis.fastq_2": "/path/to/your/reads_R2.fastq.gz",
     ...
   }
   ```

3. **Run the workflow**:
   ```bash
   java -jar cromwell.jar run ww-saturation.wdl -i ww-saturation-inputs.json -o ww-saturation-options.json
   ```

## Input Parameters

| Parameter | Description |
|-----------|-------------|
| `fastq_1` | Path to first FASTQ file (R1) |
| `fastq_2` | Path to second FASTQ file (R2) |
| `reference_fasta` | Path to reference genome FASTA file |
| `reference_fasta_index` | Path to reference genome index (.fai) |
| `reference_dict` | Path to reference sequence dictionary (.dict) |
| `sample_name` | Sample identifier used for output file naming |
| `orf_range` | Open Reading Frame range (e.g., "122-781") |

## Resource Configuration

Task-specific resources can be configured in the inputs file:

- `AlignReads.threads`: Number of CPU threads for BWA alignment (default: 4)
- `AlignReads.memory_gb`: Memory allocation for alignment in GB (default: 8)
- `SamToBam.threads`: Number of CPU threads for SAM/BAM processing (default: 4)
- `SamToBam.memory_gb`: Memory allocation for SAM/BAM processing in GB (default: 8)
- `AnalyzeSaturationMutagenesis.memory_gb`: Memory allocation for GATK in GB (default: 16)

## Outputs

The workflow produces the following output files:

- `analysis_table`: Tab-delimited table of saturation mutagenesis results
- `analysis_plot`: PDF visualization of the saturation mutagenesis results
- `aligned_bam`: Sorted BAM file containing aligned reads
- `aligned_bai`: Index file for the aligned BAM

## Runtime Options

The `ww-saturation-options.json` file contains Cromwell runtime options:

- Continues execution of remaining tasks when possible if one task fails
- Enables workflow caching for faster reruns
- Configures task retry behavior
- Specifies output directory structure

## Contributing

Contributions to improve the workflow are welcome. Please submit pull requests or open issues to suggest enhancements or report bugs.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
