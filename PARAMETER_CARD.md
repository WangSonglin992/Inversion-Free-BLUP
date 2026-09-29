# Parameter Card Guide

IF-BLUP is configured with an INI-format **parameter card**. This document
describes every section and option. Ready-to-edit templates ship with the
repository: `config.ini` (minimal single-trait card wired to the example data
in `Data/`) and the documented samples in `examples/`.

General rules:

- The card is standard INI: `key = value` lines grouped under `[SECTION]`
  headers; end-of-line comments start with `;`.
- Column indices are **1-based** and refer to columns of the input files.
- All `[ADVANCED]` options are optional — leave the section out (or leave
  values empty) to use the built-in, benchmark-tuned defaults.

## Command line

```text
python ifblup.py --ini CARD [OPTIONS]

--ini CARD          Path to the parameter card (required)
--output DIR        Root output directory (overrides output_dir in the card)
--inbreeding PED    Standalone mode: compute inbreeding coefficients from a
                    pedigree file and exit (no MME solve)
--makeG             Standalone mode: build a G matrix from PLINK binary files
                    and export it as CSV
--bfile PREFIX      PLINK binary prefix for --makeG (without .bed/.bim/.fam)
--out CSV_FILE      Output CSV path for --makeG (default: G_matrix.csv)
--version, -v       Show version information
--about             Show about information
```

If `--output` is omitted and `output_dir` is empty, results are written to a
timestamped `<model>_output_<YYYYmmdd_HHMMSS>/` folder next to the card.

## [GENERAL] — run mode

| Option | Values | Description |
|---|---|---|
| `mode` | `single` \| `multi` | Single-trait or multi-trait evaluation. |
| `relation_model` | `auto` \| `ABLUP` \| `GBLUP` \| `ssGBLUP` | Relationship model. `auto` picks from the provided files: pedigree + genotype → `ssGBLUP`, genotype only → `GBLUP`, pedigree only → `ABLUP`. |
| `output_dir` | path | Optional default output root. |

## [PATHS] — input files

| Option | Required for | Description |
|---|---|---|
| `phenotype_file` | always | Phenotype file **with header**. |
| `pedigree_file` | ABLUP, ssGBLUP | Pedigree file, format `ID Sire Dam`; missing parents coded `0` / `NA`. |
| `bed_file_basename` | GBLUP, ssGBLUP | Prefix of PLINK binary files, **without** the `.bed`/`.bim`/`.fam` extensions. |

## [COLUMNS] — column indices (1-based)

| Option | Description |
|---|---|
| `ped_col_id` | Column of the individual ID in the pedigree file. |
| `ped_col_sire` | Column of the sire ID in the pedigree file. |
| `ped_col_dam` | Column of the dam ID in the pedigree file. |
| `pheno_col_id` | Column of the individual ID in the phenotype file. |

## [MODEL_PARAMS_SINGLE] — single-trait model (`mode = single`)

| Option | Description |
|---|---|
| `trait_col` | Column of the trait value in the phenotype file. |
| `fixed_effects` | Comma-separated fixed-effect column indices, e.g. `8,9`. |
| `covariates` | Comma-separated covariate column indices, e.g. `10`. |
| `h2` | Initial heritability, 0.0–0.99. For repeatability traits use `h2 = var_a / (var_a + var_e)` (PE variance excluded). |
| `random_effects` | Comma-separated columns of additional random-effect grouping factors (e.g. `1` = PE group column). |
| `random_ratios` | Variance ratios `var_re / var_e` of the corresponding random effects, same order as `random_effects`. For repeatability traits: `var_pe / var_e`. |

## [MODEL_PARAMS_MULTI] — multi-trait model (`mode = multi`)

| Option | Description |
|---|---|
| `G0` | Initial genetic (additive) covariance matrix, `n_traits × n_traits`, symmetric positive definite. Formats: `[[5.1, 1.6], [1.6, 3.2]]` or `5.1,1.6; 1.6,3.2`. |
| `R0` | Initial residual covariance matrix, same format. |
| `random_effect_names` | Names of additional random effects, comma-separated (e.g. `PE` or `Litter,PE`). |
| `random_covariances` | Covariance matrices of the additional random effects, e.g. `{'PE': [[1.9,0,0],[0,0.6,0],[0,0,1.6]]}`. |

## [T1], [T2], … — per-trait definitions (`mode = multi`)

One section per trait, in evaluation order.

| Option | Description |
|---|---|
| `trait_col` | Column of this trait in the phenotype file. |
| `fixed_effects` | Fixed-effect columns of this trait. |
| `covariates` | Covariate columns of this trait. |
| `random_effects` | Column of the additional random-effect grouping factor for this trait; must follow the order of `random_effect_names`. For a repeatability model set it to the animal-ID column (e.g. `1`). |

## [ADVANCED] — optional tuning

Omit the whole section to use the defaults below.

### Compute performance

| Option | Default | Description |
|---|---|---|
| `num_threads` | 4 | Parallel threads; recommended ≤ physical core count. |
| `block_size_G` | 110000 | Block size for G-matrix construction. |
| `snp_block_size` | 51000 | Block size for SNP reading. |

### Memory strategy

| Option | Default | Description |
|---|---|---|
| `g_mem_threshold_ratio` | 0.5 | RAM threshold for loading G (0.1–0.9). |
| `g_load_strategy` | auto | `auto` \| `in_memory` \| `memmap`. |
| `overwrite_g_memmap` | true | Overwrite an old G memmap cache. |

### Algorithm convergence

| Option | Default | Description |
|---|---|---|
| `max_solver_iterations` | 15000 | Iteration cap for the staged PCG solver. |
| `final_rtol` | 1e-6 | Final relative residual tolerance; a provided value takes effect as-is. |
| `final_atol` | 1e-6 | Final absolute residual tolerance; a provided value takes effect as-is. |

### Numerical stability

| Option | Default | Description |
|---|---|---|
| `apply_g_bending_build` | false | Apply diagonal bending while building G (benchmarks show bending reduces accuracy; keep disabled). |
| `g_bending_build` | 0.0 | Bending value used when building G. |
| `g_bending_ssgblup` | 0.05 | Bending value used during ssGBLUP integration. |
| `g_blending_alpha` | 0.05 | ssGBLUP blending ratio. |
| `scaling_mode_multi` | — | Multi-trait ssGBLUP scaling: `diag_inv_mean` \| `hutchinson`. |

### Solver strategy

| Option | Default | Description |
|---|---|---|
| `solver_strategy` | auto | `auto` \| `explicit` \| `implicit`. **explicit**: dense G with an APY-like preconditioner. **implicit**: M-matrix factorization + Woodbury (LMP) preconditioner, saves memory. |
| `matrix_storage_dtype` | float64 | `float64` \| `int2`. `int2` is a 2-bit genotype compression that saves ~50× memory with < 5% accuracy loss. |

### Coarse-space correction (augmentation)

| Option | Default | Description |
|---|---|---|
| `enable_augmentation` | false | Overlay a small Galerkin coarse-space correction before each PCG stage (main preconditioner unchanged); cheaply suppresses the slowest modes. |
| `augmentation_basis_mode` | ssgblup_blocks | `ssgblup_blocks` \| `blocks` \| `auto`. |
| `augmentation_max_vectors` | 6 | Maximum number of orthogonal coarse-space vectors (1–20). |
| `augmentation_stage_limit` | 3 | Apply the correction in stages 1–N only (0 disables). Recommended: first 2–3 stages, final stage clean. |

### Data quality / output

| Option | Default | Description |
|---|---|---|
| `write_inbreeding_file` | false | Also write `inbreeding_coefficients.csv` (ID, Inbreeding). |
| `standardize_covariates` | false | Standardize covariate columns before solving. |
| `missing_rate_threshold_pheno` | 0.5 | Warning threshold for the phenotype missing rate. |
| `calc_reliability` | false | Approximate breeding-value reliability. Exact reliability is incompatible with the inversion-free design and would add substantial overhead; only an approximation is computed. |

## Missing value recognition

Phenotype missing markers (case-insensitive): `NA`, `N/A`, `NaN`, `nan`,
`NAN`, `None`, `null`, `.`, `-999`, `-9999`, `-99999`, empty string.
Missing pedigree parents: `0`, `0.0` and `-999` all mean unknown parent.

## Recommended environment

| Item | Recommendation |
|---|---|
| OS | Linux x86_64 (CentOS 7.9+, Ubuntu 20.04+, Debian 10+) |
| CPU | 16+ cores (multithreaded PCG and G construction) |
| Memory | 20k individuals: ≥ 4 GB · 50k: ≥ 16 GB · 100k+: ≥ 64 GB |
| Storage | SSD recommended (streaming matrix I/O) |

## Worked example — multi-trait ssGBLUP

```ini
[GENERAL]
mode = multi
relation_model = ssGBLUP
output_dir =

[PATHS]
pedigree_file = ./data/pedigree.tsv
phenotype_file = ./data/phenotype.csv
bed_file_basename = ./data/genotype

[COLUMNS]
ped_col_id = 1
ped_col_sire = 2
ped_col_dam = 3
pheno_col_id = 1

[MODEL_PARAMS_MULTI]
G0 = [[0.50, 0.15, 0.10], [0.15, 0.45, 0.12], [0.10, 0.12, 0.40]]
R0 = [[1.00, 0.10, 0.05], [0.10, 0.90, 0.08], [0.05, 0.08, 0.80]]

[T1]
trait_col = 4
fixed_effects = 2,3

[T2]
trait_col = 5
fixed_effects = 2,3

[T3]
trait_col = 6
fixed_effects = 2,3

[ADVANCED]
num_threads = 8
```
