version 1.0

workflow saturation_mutagenesis {
    input {
        Array[File] fastq_1_array
        Array[File] fastq_2_array
        Array[String] sample_names
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        String orf_range
    }

    # Scatter over the samples
    scatter (idx in range(length(sample_names))) {
        String sample_name = sample_names[idx]
        File fastq_1 = fastq_1_array[idx]
        File fastq_2 = fastq_2_array[idx]

        # Step 1: Align Nextera sequencing reads to reference using BWA
        call AlignReads {
            input:
                fastq_1 = fastq_1,
                fastq_2 = fastq_2,
                reference_fasta = reference_fasta,
                sample_name = sample_name
        }

        # Step 2: Convert SAM to BAM and sort
        call SamToBam {
            input:
                input_sam = AlignReads.aligned_sam,
                sample_name = sample_name
        }

        # Step 3: Analyze Saturation Mutagenesis
        call AnalyzeSaturationMutagenesis {
            input:
                input_bam = SamToBam.sorted_bam,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                orf_range = orf_range,
                sample_name = sample_name
        }
    }

    output {
        Array[File] variant_counts = AnalyzeSaturationMutagenesis.variant_counts
        Array[File] aa_counts = AnalyzeSaturationMutagenesis.aa_counts
        Array[File] aa_fractions = AnalyzeSaturationMutagenesis.aa_fractions
        Array[File] codon_counts = AnalyzeSaturationMutagenesis.codon_counts
        Array[File] codon_fractions = AnalyzeSaturationMutagenesis.codon_fractions
        Array[File] cov_length_counts = AnalyzeSaturationMutagenesis.cov_length_counts
        Array[File] read_counts = AnalyzeSaturationMutagenesis.read_counts
        Array[File] ref_coverage = AnalyzeSaturationMutagenesis.ref_coverage
    }
}

# Align Nextera sequencing reads to reference using BWA
task AlignReads {
    input {
        File fastq_1
        File fastq_2
        File reference_fasta
        String sample_name
        Int threads = 4
        Int memory_gb = 8
    }

    String output_sam = "~{sample_name}.aligned.sam"

    command <<<
        set -eo pipefail

        # Index the reference genome if not already indexed
        bwa index ~{reference_fasta}

        # Align paired-end reads using BWA MEM
        bwa mem \
            -t ~{threads} \
            -M \
            ~{reference_fasta} \
            ~{fastq_1} \
            ~{fastq_2} \
            > ~{output_sam}
    >>>

    runtime {
        docker: "getwilds/bwa:0.7.17"
        memory: "~{memory_gb} GB"
        cpu: threads
    }

    output {
        File aligned_sam = "~{output_sam}"
    }
}

# Convert SAM to BAM, sort, and index
task SamToBam {
    input {
        File input_sam
        String sample_name
        Int threads = 4
        Int memory_gb = 8
    }

    String output_bam = "~{sample_name}.sorted.bam"

    command <<<
        set -eo pipefail

        # Convert SAM to BAM
        samtools view -b -o ~{output_bam} ~{input_sam}
    >>>

    runtime {
        docker: "getwilds/samtools:1.11"
        memory: "~{memory_gb} GB"
        cpu: threads
    }

    output {
        File sorted_bam = "~{output_bam}"
    }
}

# Analyze Saturation Mutagenesis using GATK
task AnalyzeSaturationMutagenesis {
    input {
        File input_bam
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        String orf_range
        String sample_name
        Int memory_gb = 16
    }

    command <<<
        set -eo pipefail
        
        gatk AnalyzeSaturationMutagenesis \
            -R ~{reference_fasta} \
            -I ~{input_bam} \
            --orf ~{orf_range} \
            -O ~{sample_name}
    >>>

    runtime {
        docker: "getwilds/gatk:4.3.0.0"
        memory: "~{memory_gb} GB"
    }

    output {
        File aa_counts = "~{sample_name}.aaCounts"
        File aa_fractions = "~{sample_name}.aaFractions"
        File codon_counts = "~{sample_name}.codonCounts"
        File codon_fractions = "~{sample_name}.codonFractions"
        File cov_length_counts = "~{sample_name}.coverageLengthCounts"
        File read_counts = "~{sample_name}.readCounts"
        File ref_coverage = "~{sample_name}.refCoverage"
        File variant_counts = "~{sample_name}.variantCounts"
    }
}
