# IF-BLUP

IF-BLUP (Inversion-Free BLUP) is a Python solver for large-scale genetic
evaluation. It supports ABLUP, GBLUP and ssGBLUP with single-trait and
multi-trait models, repeatability random effects, CPU/GPU backends, and
memory-aware matrix storage with staged PCG solving.

This repository ships IF-BLUP as a **compiled distribution**: the command-line
entry points are plain Python, while all computational modules
(`core/`, `mme/`, `reporting/`, `utils/`) are precompiled binary extensions —
no build step is needed.

- Current version: **v1.3.6** (2026-08-25)
- Platform: **Linux x86_64**, GLIBC >= 2.28 (RHEL 8+, Ubuntu 20.04+)
- Python: **CPython 3.13.x** (the bundled extensions are CPython 3.13 ABI)

## Installation

Users install the runtime dependencies themselves; nothing needs to be
compiled.

```bash
git clone https://github.com/WangSonglin992/Inversion-Free-BLUP.git
cd Inversion-Free-BLUP
python3.13 -m venv .venv
source .venv/bin/activate
pip install -r requirements-cpu.txt

python ifblup.py --version
```

For GPU evaluation, additionally install the CuPy build matching your CUDA
runtime (see `requirements-gpu.txt`), for example:

```bash
pip install cupy-cuda12x
```

The CUDA backend is picked up automatically when CuPy is importable; CPU-only
environments need nothing extra.

## 中文快速开始

1. 克隆本仓库并创建 **Python 3.13** 虚拟环境；
2. `pip install -r requirements-cpu.txt`（GPU 环境另装与 CUDA 匹配的 CuPy，见 `requirements-gpu.txt`）；
3. 复制根目录 `config.ini` 或 `examples/` 下的参数卡模板到工作目录，按注释填写数据路径与列索引；
4. 运行：

```bash
python ifblup.py --ini /path/to/your-config.ini --output /path/to/output
```

5. 结果（EBV/可靠性/固定效应解）与诊断日志 `ifblup_diagnostic.log` 写入输出目录。

## Running

```bash
# full genetic evaluation
python ifblup.py --ini /path/to/config.ini --output /path/to/output

# standalone inbreeding coefficients from a pedigree file
python ifblup.py --inbreeding /path/to/pedigree.tsv --out /path/to/f.csv

# standalone G matrix from PLINK binary files
python ifblup.py --makeG --bfile /path/to/genotype-prefix --out /path/to/G.csv
```

`config.ini` at the repository root is a minimal single-trait template wired to
the example data in `Data/`; `examples/` documents the single-trait and
multi-trait parameter cards. **Every section and option of the parameter card
is described in [PARAMETER_CARD.md](PARAMETER_CARD.md).** Real evaluations must
use your own pedigree, phenotype and (optionally) PLINK binary files that
satisfy the column contract described in the cards.

## Repository layout

| Path | Content |
|---|---|
| `ifblup.py` | Python CLI entry point (implementation ships compiled) |
| `core/`, `mme/`, `reporting/`, `utils/` | Precompiled binary extensions (CPython 3.13 ABI) |
| `config.ini`, `examples/` | Parameter-card template and documented samples |
| `PARAMETER_CARD.md` | Full parameter-card reference (all sections and options) |
| `requirements-cpu.txt`, `requirements-gpu.txt` | Runtime dependencies |
| `Data/` | Minimal example data referenced by `config.ini` |
| `scripts_public/` | Public helper scripts |

## Notes

- The core algorithms are distributed as compiled extensions; the Python
  sources in this repository are limited to entry points and package markers.
- Version, release date and developer credits are embedded in the binary
  release info and shown by `python ifblup.py --version`.
- Do not commit run outputs, matrices caches or enterprise data back to this
  repository.

## License

Apache License 2.0 — see `LICENSE`.
