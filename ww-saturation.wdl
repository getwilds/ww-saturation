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

        # Step 1: Align RNA-seq reads to reference using BWA
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
                input_bai = SamToBam.sorted_bai,
                reference_fasta = reference_fasta,
                orf_range = orf_range,
                sample_name = sample_name
        }
    }

    output {
        Array[File] analysis_tables = AnalyzeSaturationMutagenesis.tablefile
        Array[File] analysis_plots = AnalyzeSaturationMutagenesis.plotfile
        Array[File] aligned_bams = SamToBam.sorted_bam
        Array[File] aligned_bais = SamToBam.sorted_bai
    }
}

# Align RNA-seq reads to reference using BWA
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
        # Convert SAM to BAM and sort
        samtools view -bS ~{input_sam} | samtools sort -@ ~{threads} -o ~{output_bam}
        
        # Index the BAM file
        samtools index ~{output_bam}
    >>>

    runtime {
        docker: "getwilds/samtools:1.11"
        memory: "~{memory_gb} GB"
        cpu: threads
    }

    output {
        File sorted_bam = "~{output_bam}"
        File sorted_bai = "~{output_bam}.bai"
    }
}

# Analyze Saturation Mutagenesis using GATK
task AnalyzeSaturationMutagenesis {
    input {
        File input_bam
        File input_bai
        File reference_fasta
        String orf_range
        String sample_name
        Int memory_gb = 16
    }

    String output_table = "~{sample_name}.mutagenesis_analysis.table"
    String output_plot = "~{sample_name}.mutagenesis_analysis.pdf"

    command <<<
        gatk AnalyzeSaturationMutagenesis \
            -R ~{reference_fasta} \
            -I ~{input_bam} \
            --orf ~{orf_range} \
            --output-table ~{output_table} \
            --output-plot ~{output_plot}
    >>>

    runtime {
        docker: "getwilds/gatk:4.3.0.0"
        memory: "~{memory_gb} GB"
    }

    output {
        File tablefile = "~{output_table}"
        File plotfile = "~{output_plot}"
    }
}
