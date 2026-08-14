# Check incorrect trimming effect on read length and CENP-A signal
New
```bash
git clone df14e05d229ad0d09c662bc4560c0092bfaaea62
snakemake -np --configfile config.yaml -c 12 --sdm conda
mv results_new
```

Old
```bash
git clone 81ab6cd650eadc4d1133ef3a7cee373409331504
snakemake -np --configfile config.yaml -c 12 --sdm conda
```

Plot difference in read lengths from new to old.
```bash
python experiments/incorrect_trimming/cmp.py
```

Also look at bigwigs in IGV
* New: `results/CENP-A/reads_to_ref.bw`
* Old: `results_new/NA20355/bam/reads_to_ref.bw`

See ppt:
https://docs.google.com/presentation/d/1iuILxFNxDyoZnGy3w7g2M_WqShA9ZQVsGW9YqkDJzkA/edit?usp=sharing
