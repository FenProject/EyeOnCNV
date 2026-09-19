# Installation

| Component | Needed for | Notes |
|---|---|---|
| Python ≥ 3.10 with numpy, pandas, openpyxl | every Python step | installed by `pip install -e .` |
| bash, awk, samtools ≥ 1.10 | step 0 | macOS/Linux/WSL; Git Bash on Windows if samtools is on the PATH |
| [hgvs](https://github.com/biocommons/hgvs) + network access | step 3 only | macOS/Linux/WSL, not native Windows |
| xlrd | reading old `.xls` files | optional: `pip install -e ".[xls]"` |

## Option 1 - conda (recommended)

[Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Mamba provide samtools
and pre-built binaries of the libraries hgvs needs (`psycopg2`, `pysam`):

```bash
cd eyeoncnv
conda env create -f environment.yml
conda activate eyeoncnv
eyeoncnv --version
```

## Option 2 - pip

```bash
cd eyeoncnv
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e .              # steps 0-2 and 4
pip install -e ".[hgvs]"      # step 3 as well
```

samtools must then be installed separately:

- macOS: `brew install samtools`
- Debian/Ubuntu (and WSL): `sudo apt install samtools`
- any system with conda: `conda install -c conda-forge -c bioconda samtools`

If `pip install hgvs` fails while building `psycopg2`, install the binary wheel first
(`pip install psycopg2-binary`) or use conda.

## Windows

- **WSL (Ubuntu)** is the simplest option: every step works, follow the Linux instructions.
- **Native Windows**: steps 1, 2 and 4 work with pip. Step 0 needs
  bash (e.g. Git Bash) and a samtools executable on the PATH; Windows paths are written
  `/c/Users/...` in Git Bash. Step 3 does not work: `hgvs` depends on `pysam`, which cannot
  be built on Windows. Run step 3 on macOS, Linux or WSL - each step reads and writes plain
  Excel files, so the steps can be run on different computers.

## Step 3: UTA database and sequences

`eyeoncnv hgvs` maps positions with the [UTA](https://github.com/biocommons/uta) transcript
database. By default hgvs connects to the public server `uta.biocommons.org` on port
**5432** (PostgreSQL), and fetches reference bases from NCBI. Hospital networks often
block that port; in that case:

- run the step from another network, or
- run a local copy of UTA (a Docker image is described in the UTA README) and point to it
  with `--uta-url postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b`, or
  with the `UTA_DB_URL` environment variable;
- optionally install a local [SeqRepo](https://github.com/biocommons/biocommons.seqrepo)
  and set `HGVS_SEQREPO_DIR` to avoid calls to NCBI.

The UTA release in the URL above is an example: use the one given by the UTA documentation.

## Check the installation

```bash
pip install -e ".[dev]"
pytest
```

Tests needing samtools or several awk implementations are skipped when they are absent.
