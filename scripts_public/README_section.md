## Simulated Data Generation

We generate three simulated datasets using [SIMER](https://github.com/xiaolei-lab/SIMER) for large-scale benchmarking of genomic prediction methods.

### Scripts

| Script | Description | Scale |
|------|------|------|
| `sim_baseline_pig.R` | Pig growth traits, additive-only model, no fixed effects | 1.28M / 50k SNPs / 3 traits |
| `sim_cattle_dairy.R` | Dairy cattle pedigree with masked phenotypes (50% NA) | 5M / 50k SNPs / single trait |
| `sim_human_disease.R` | Human complex disease with fixed effects, three SNP panels | 1M / 500k–2M SNPs / 7 traits |

### Usage

```bash
Rscript sim_baseline_pig.R
Rscript sim_cattle_dairy.R
Rscript sim_human_disease.R
Rscript sim_human_disease.R --markers=2000000 --ncpus=32
```

### Dependencies

R >= 4.0, [SIMER](https://github.com/xiaolei-lab/SIMER), bigmemory, data.table, PLINK 1.9

All scripts create a `./output/` directory under the working directory to store the generated PLINK files.
