# CUT&RUN ONT workflow
Does the following:
1. Trims adapters
2. Filters reads between 100 - 10kbp
3. Aligns with bwa-mem
4. Normalizes to sample_baseline and gnera

Create conda environment.
```bash
conda env create --name cutnrun env.yaml
```

Modify to fit your use-case.
```yaml
ref: /project/logsdon_shared/projects/HGSVC3/new_65_asms_renamed/NA20355-asm-renamed-reort.fa
sample_baseline: IgG_2.5
samples:
  CENP-A_2.5:
    path: /project/logsdon_shared/projects/HGSVC3/CHM13_CUTNRUN/data/GM20355.bam
  IgG_2.5:
    path: /project/logsdon_shared/long_read_archive/unsorted/20250509_CNR_GM20355_SC_NBD114/BD1-3/20250509_1518_3C_PAU14051_56a49125/pod5/demux/56a49125-8307-45bf-94e1-c70aca65814b_SQK-NBD114-24_barcode03.bam
```

```bash
snakemake -np --configfile config.yaml -c 12
```
