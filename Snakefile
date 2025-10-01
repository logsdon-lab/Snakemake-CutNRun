from os.path import join, dirname


SAMPLES = config["samples"]
SAMPLES_NON_BASELINE = [sm for sm in SAMPLES if sm != config["sample_baseline"]]
LOG_DIR = config.get("log_dir", "logs")
BMK_DIR = config.get("benchmark_dir", "benchmarks")
OUTPUT_DIR = config.get("output_dir", "results")
MIN_LEN, MAX_LEN = config.get("min_length", 100), config.get("max_length", 10_000)

scattergather:
    split=4

rule convert_to_fq:
    input:
        ubam=lambda wc: SAMPLES[wc.sm]["path"]
    output:
        fastq=join(OUTPUT_DIR, "{sm}", "reads.fq"),
        fastq_fai=join(OUTPUT_DIR, "{sm}", "reads.fq.fai"),
    conda:
        "env.yaml"
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
        "env.yaml"
    log:
        join(LOG_DIR, "{sm}_bam_to_fq.log")
    shell:
        """
        cat {input} | rustybam fastq-split {output}
        """
    
# Trim adapters
# We don't need FastQC and frequently run into issues of trimgalore freezing
# See https://github.com/FelixKrueger/TrimGalore/issues/205
rule trim_adapter:
    input:
        fq=join(OUTPUT_DIR, "{sm}", "reads_{scatteritem}.fq")
    output:
        reads=join(OUTPUT_DIR, "{sm}", "trimmed_reads", "reads_{scatteritem}_trimmed.fq"),
    conda:
        "env.yaml"
    log:
        join(LOG_DIR, "{sm}_trim_adapters_{scatteritem}.log")
    threads:
        12
    params:
        error_rate=0.1,
        quality=20,
        stringency=5,
        motif="AGATCGGAAGAGC"
    shell:
        """
        cutadapt \
        -j {threads} \
        -e {params.error_rate} \
        -q {params.quality} \
        -O {params.stringency} \
        -a {params.motif} \
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
        "env.yaml"
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

rule plot_reads:
    input:
        script="ont_stats",
        fai=rules.filter_reads.output.reads_fai,
        # fai=rules.downsample_to_lowest_read_number.output.downsampled_reads_fai,
    output:
        join(OUTPUT_DIR, "{sm}", "cdf", "{sm}_read_length.pdf")
    params:
        outdir=lambda wc, output: os.path.dirname(str(output[0]))
    conda:
        "env.yaml"
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
        "env.yaml"
    log:
        join(LOG_DIR, "index_ref.log")
    shell:
        """
        bwa index {input.asm} &> {log}
        """

rule align_reads_bwa:
    input:
        reads=rules.filter_reads.output.reads,
        reference=config["ref"],
        index=rules.index_asm.output,
    output:
        bam=join(OUTPUT_DIR, "{sm}", "reads_to_ref.bam"),
        bai=join(OUTPUT_DIR, "{sm}", "reads_to_ref.bam.bai"),
    conda:
        "env.yaml"
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

rule convert_bam_to_bw:
    input:
        bam_1=rules.align_reads_bwa.output.bam,
        bam_2=expand(rules.align_reads_bwa.output.bam, sm=config["sample_baseline"])
    output:
        bw=join(OUTPUT_DIR, "{sm}", "reads_to_ref.bw"),
    params:
        # Noisy reads so MAPQ 1 filters everything.
        mapq=config.get("min_mapq", 0),
        binSize=config.get("binsize", 5000),
    conda:
        "env.yaml"
    threads:
        20
    log:
        join(LOG_DIR, "{sm}_convert_bam_to_bw.log")
    shell:
        """
        bamCompare \
        -b1 {input.bam_1} \
        -b2 {input.bam_2} \
        --operation ratio \
        --binSize {params.binSize} \
        --minMappingQuality {params.mapq} \
        -p {threads} \
        -o {output} 2> {log} 
        """

rule all:
    input:
        expand(rules.align_reads_bwa.output, sm=SAMPLES),
        expand(rules.convert_bam_to_bw.output, sm=SAMPLES_NON_BASELINE),
        expand(rules.plot_reads.output, sm=SAMPLES)
    default_target:
        True
