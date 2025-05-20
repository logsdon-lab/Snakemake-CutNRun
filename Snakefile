from os.path import join, dirname


SAMPLES = config["samples"]
LOG_DIR = config.get("log_dir", "logs")
BMK_DIR = config.get("benchmark_dir", "benchmarks")
OUTPUT_DIR = config.get("output_dir", "results")
MIN_LEN, MAX_LEN = 100, 10_000

scattergather:
    split=4

rule convert_to_fq:
    input:
        ubam=lambda wc: SAMPLES[wc.sm]["path"]
    output:
        fastq=join(OUTPUT_DIR, "{sm}", "reads.fq"),
        fastq_fai=join(OUTPUT_DIR, "{sm}", "reads.fq.fai"),
    conda:
        "general"
    log:
        join(LOG_DIR, "{sm}_bam_to_fq.log")
    shell:
        """
        samtools bam2fq -T "*" {input.ubam} > {output.fastq} 2> {log}
        samtools faidx {output.fastq}
        """

# Cannot set buffersize in files with large number of reads for trim_galore so must partition.
rule partition_reads:
    input:
        rules.convert_to_fq.output.fastq
    output:
        scatter.split(
            join(OUTPUT_DIR, "{{sm}}", "reads_{scatteritem}.fq")
        )
    conda:
        "general"
    log:
        join(LOG_DIR, "{sm}_bam_to_fq.log")
    shell:
        """
        cat {input} | rustybam fastq-split {output}
        """
    
# Trim adapters
rule trim_adapter:
    input:
        fq=join(OUTPUT_DIR, "{sm}", "reads_{scatteritem}.fq")
    output:
        reads=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_{scatteritem}_trimmed.fq"),
    conda:
        "general"
    log:
        join(LOG_DIR, "{sm}_trim_adapters_{scatteritem}.log")
    threads:
        12
    params:
        min_len=MIN_LEN,
        stringency=5,
        outdir=lambda wc, output: dirname(output.reads),
    shell:
        """
        cutadapt \
        -j {threads} \
        -e 0.1 \
        -q 20 \
        -O 5 \
        -a AGATCGGAAGAGC \
        {input.fq} > {output} 2> {log}
        """

rule gather_trimmed_reads:
    input:
        gather.split(join(OUTPUT_DIR, "{{sm}}", "trimmed_reads", "reads_{scatteritem}_trimmed.fq"))
    output:
        reads=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_trimmed.fq"),
    shell:
        """
        cat {input} > {output}
        """

rule filter_reads:
    input:
        reads=rules.gather_trimmed_reads.output.reads,
    output:
        reads=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_trimmed_filtered.fq"),
        reads_fai=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_trimmed_filtered.fq.fai"),
    conda:
        "general"
    log:
        join(LOG_DIR, "{sm}_filter_reads.log")
    params:
        min_read_len=MIN_LEN,
        max_read_len=MAX_LEN,
    shell:
        """
        seqkit seq -m {params.min_read_len} -M {params.max_read_len} {input.reads} > {output.reads}
        samtools faidx {output.reads}
        """

rule get_number_of_reads:
    input:
        reads_fai=expand(rules.filter_reads.output.reads_fai, sm=SAMPLES),
    output:
        join(OUTPUT_DIR, "read_number.txt"),
    conda:
        "general"
    shell:
        """
        wc -l {input.reads_fai} | head -n-1 | awk '{{ print $1, $2 }}' | sort > {output}
        """

# Downsample to same number of reads across sample.
rule downsample_to_lowest_read_number:
    input:
        reads=rules.filter_reads.output.reads,
        read_nums=rules.get_number_of_reads.output
    output:
        downsampled_reads=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_trimmed_filtered_downsampled.fq"),
        downsampled_reads_fai=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_trimmed_filtered_downsampled.fq.fai"),
    params:
        seed=100
    conda:
        "general"
    log:
        join(LOG_DIR, "downsample_{sm}_reads_to_lowest_num.log")
    shell:
        """
        min_n_reads=$(head -n1 {input.read_nums} | cut -f 1)
        seqtk sample -s{params.seed} {input.reads} ${{min_n_reads}} > {output.downsampled_reads}
        samtools faidx {output.downsampled_reads}
        """

rule plot_reads:
    input:
        script="ont_stats",
        fai=rules.downsample_to_lowest_read_number.output.downsampled_reads_fai,
    output:
        join(OUTPUT_DIR, "{sm}", "cdf", "{sm}_read_length.pdf")
    params:
        outdir=lambda wc, output: os.path.dirname(str(output[0]))
    conda:
        "general"
    log:
        join(LOG_DIR, "plot_reads_{sm}.log")
    shell:
        """
        python {input.script} --fai {wildcards.sm}={input.fai} -p {params.outdir} --plot_ext pdf &> {log}
        """

rule index_asm:
    input:
        asm=config["ref"]
    output:
        multiext(
            config["ref"],
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        )
    conda:
        "general"
    log:
        join(LOG_DIR, "index_ref.log")
    shell:
        """
        bwa index {input.asm} &> {log}
        """

rule align_reads_bwa:
    input:
        reads=rules.downsample_to_lowest_read_number.output.downsampled_reads,
        reference=config["ref"],
        index=rules.index_asm.output,
    output:
        bam=join(OUTPUT_DIR, "{sm}", "reads_to_ref.bam"),
        bai=join(OUTPUT_DIR, "{sm}", "reads_to_ref.bam.bai"),
    conda:
        "general"
    params:
        k=50,
        c=1000000,
    threads:
        12
    log:
        join(LOG_DIR, "{sm}_align_reads_bwa.log")
    shell:
        """
        bwa mem -k {params.k} -c {params.c} -t {threads} {input.reference} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule all:
    input:
        expand(rules.align_reads_bwa.output, sm=SAMPLES),
        expand(rules.plot_reads.output, sm=SAMPLES)
    default_target:
        True
