AGENTS: BluePebble (BP1) – genetic epi / genomics

Cluster: BluePebble (BP1), University of Bristol
User work root: /user/work/fh6520
Default Slurm account: sscm013902

This file is generic – edit paths/modules per repo.

1. Core workflow rules

Never debug via Slurm if you can debug via bash.

When you add/change a script:

Run a quick bash test on the login node to catch missing modules, syntax errors, and path issues.

Only when the script works on a small test should you write/update a Slurm script and sbatch it.

Use timeouts here - first 30 seconds, then 1 minute, etc

Write lots of logging so it's easy to identify if a script is working but slow and will finish in 4m (keep in bash), or taking hours (we can identify at 1m should be a slurm job)

Everything must be scripted and reproducible.

Even one-off actions should be written as scripts and committed.

Use a numbered naming convention, e.g.:

01_munge_phenotypes.R

02_build_panel.R

03_run_gwas_chr22.R

Prefer R + tidyverse for data wrangling over ad-hoc shell or Python.

Use chr22 as the default test chromosome.

For new workflows, first implement and test on chr22 and/or a small N subset.

Only then extend to full chr1–22 or genome-wide runs.

2. Filesystems

Per-project root (example):

/user/work/fh6520/<project>/
  src/      # R / bash scripts
  slurm/    # sbatch scripts
  data/
  results/
  logs/
UK Biobank imputed arrays (read-only, do not copy):

/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/
  data/
    dosage_bgen/         # filtered dosage BGEN
    raw_downloaded/bgen/ # raw imputed BGEN (~50–190G/chr)
Use these shared paths as inputs only; don’t write into /bp1/mrcieu1/... and don’t copy full BGENs into your project.

3. Modules / tools

Use environment modules (versions can be adjusted per repo):

module purge

# R
module load languages/R/4.5.1

# Genomics tools
module load bcftools/1.19-openblas-5yp2
module load vcftools/0.1.16-ukwm
module load apps/plink2/2.00a68LM
# module load apps/regenie/<VERSION>  # if needed
Additional notes:

Use R (tidyverse) for wrangling by default.

ast-grep is available (installed via pip) and can be used for structured search/refactors in code.

4. Working with huge files

Principle: head/subset, don’t stream or copy everything.

Examples:

# gz text
zcat file.gz | head

# VCF/BCF header and small region
bcftools view -h file.vcf.gz
bcftools view -r 1:1-1e6 file.vcf.gz | head
UKB BGEN subset (small region, chr1 example):

plink2 \
  --bgen /bp1/mrcieu1/.../dosage_bgen/data.chr01.bgen ref-first \
  --sample /bp1/mrcieu1/.../dosage_bgen/data.chr1-22.sample \
  --chr 1 --from-bp 1 --to-bp 1e6 \
  --make-pgen \
  --out chr1_region_test
Do not:

zcat or otherwise stream entire 100+GB BGEN/VCF files.

Copy UKB BGENs from /bp1/mrcieu1/... into your own data/.

5. Slurm: minimal template

Partitions:

Use short for most jobs.

Use compute only when longer/larger jobs are actually needed.

Single-job template (adapt time/mem/cpus):

#!/bin/bash
#SBATCH --job-name=<name>
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=sscm013902
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load languages/R/4.5.1
module load bcftools/1.19-openblas-5yp2
module load vcftools/0.1.16-ukwm
module load apps/plink2/2.00a68LM

cd /user/work/fh6520/<project>

Rscript src/<script>.R
6. Testing pattern (chr22)

Default testing approach for new/changed workflows:

# On login node
cd /user/work/fh6520/<project>
module load languages/R/4.5.1

# Run a small test (e.g. chr22, reduced N)
Rscript src/run_gwas_chr22.R \
  --chr 22 \
  --test-mode TRUE
Once this works:

Generalise to 1–22 in the script.

Wire into a Slurm script (based on the template above) and submit via sbatch.

7. What agents may touch

OK to create/modify:

src/, R/, slurm/, configs/, results/, logs/, data/derived/.

Do not:

edit anything under /bp1/mrcieu1/....

create huge intermediates from UKB BGENs unless explicitly requested.

