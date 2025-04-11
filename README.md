# Saturation Mutagenesis Workflow

## Overview

This repository contains a WDL (Workflow Description Language) pipeline for performing saturation mutagenesis analysis on RNA-seq data. The workflow processes multiple samples in parallel, aligning paired-end sequencing reads to a reference genome using BWA, converting the resulting SAM files to BAM format, and performing saturation mutagenesis analysis using GATK's AnalyzeSaturationMutagenesis tool.

## Workflow Steps

For each sample in parallel:

1. **Alignment (BWA)**: Paired-end reads are aligned to a reference genome
2. **SAM to BAM Conversion**: The aligned reads are converted from SAM to BAM format and indexed
3. **Saturation Mutagenesis Analysis**: GATK's AnalyzeSaturationMutagenesis tool analyzes mutations across the specified ORF range

## Files in this Repository

- `ww-saturation-parallel.wdl`: Multi-sample parallel workflow definition
- `ww-saturation-parallel-inputs.json`: Example input parameters for the workflow
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
   - Arrays of paired-end FASTQ files
   - Reference genome (FASTA, index, and dictionary files)
   - Define your ORF range
   - List of sample names

2. **Modify the inputs JSON file**:
   ```json
   {
     "saturation_mutagenesis.fastq_1_array": ["/path/to/sample1_R1.fastq.gz", "/path/to/sample2_R1.fastq.gz"],
     "saturation_mutagenesis.fastq_2_array": ["/path/to/sample1_R2.fastq.gz", "/path/to/sample2_R2.fastq.gz"],
     "saturation_mutagenesis.sample_names": ["sample1", "sample2"],
     ...
   }
   ```

3. **Run the workflow**:
   ```bash
   java -jar cromwell.jar run ww-saturation-parallel.wdl -i ww-saturation-parallel-inputs.json -o ww-saturation-options.json
   ```

## Input Parameters

| Parameter | Description |
|-----------|-------------|
| `fastq_1_array` | Array of paths to first FASTQ files (R1) for each sample |
| `fastq_2_array` | Array of paths to second FASTQ files (R2) for each sample |
| `sample_names` | Array of sample identifiers used for output file naming |
| `reference_fasta` | Path to reference genome FASTA file |
| `reference_fasta_index` | Path to reference genome index (.fai) |
| `reference_dict` | Path to reference sequence dictionary (.dict) |
| `orf_range` | Open Reading Frame range (e.g., "122-781") |

## Resource Configuration

Task-specific resources can be configured in the inputs file:

- `AlignReads.threads`: Number of CPU threads for BWA alignment (default: 4)
- `AlignReads.memory_gb`: Memory allocation for alignment in GB (default: 8)
- `SamToBam.threads`: Number of CPU threads for SAM/BAM processing (default: 4)
- `SamToBam.memory_gb`: Memory allocation for SAM/BAM processing in GB (default: 8)
- `AnalyzeSaturationMutagenesis.memory_gb`: Memory allocation for GATK in GB (default: 16)

## Outputs

The workflow produces arrays of output files:

- `analysis_tables`: Array of tab-delimited tables of saturation mutagenesis results
- `analysis_plots`: Array of PDF visualizations of the saturation mutagenesis results
- `aligned_bams`: Array of sorted BAM files containing aligned reads
- `aligned_bais`: Array of index files for the aligned BAMs

## Runtime Options

The `ww-saturation-options.json` file contains Cromwell runtime options:

- Continues execution of remaining tasks when possible if one task fails
- Enables workflow caching for faster reruns
- Configures task retry behavior
- Specifies output directory structure

## Parallelization

The workflow:

- Processes multiple samples simultaneously using WDL's scatter-gather mechanism
- Validates input arrays (checks that all arrays have the same length)
- Returns arrays of output files
- Can significantly reduce total execution time when processing multiple samples

By default, Cromwell will execute the scattered tasks in parallel up to the limit of your execution backend. You can configure the maximum concurrent tasks in your Cromwell configuration.

## Contributing

Contributions to improve the workflow are welcome. Please submit pull requests or open issues to suggest enhancements or report bugs.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
