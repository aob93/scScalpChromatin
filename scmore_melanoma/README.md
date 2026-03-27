# scMORE Melanoma Workflow

This is a standalone workflow scaffold for running melanoma GWAS integration with [`scMORE`](https://github.com/mayunlong89/scMORE) after the ArchR analysis completes.

It is designed to let you start immediately on:

- FUMA-ready summary-statistic normalization
- GWAS summary-statistic normalization
- MAGMA / FUMA gene-result normalization

and then later plug in:

- the final integrated ArchR project from the ATAC workflow

to build a `Seurat`/`Signac` object suitable for `scMORE`.

## Intended Analysis Scope

Recommended primary analysis:

- controls only
- disease-status filter: `C_SD,C_PB`
- cell type labels from RNA-guided integration, for example `NLabelClust_RNA`

Recommended sensitivity analyses:

- all samples
- ATAC-native labels instead of RNA-guided labels

## Expected Inputs

### Immediate Inputs

- raw melanoma GWAS summary statistics
- MAGMA gene results, or FUMA gene results

### Later Inputs

- final ArchR project directory
  - typically the `fine_clustered` ArchR project produced by the ATAC workflow

## Workflow Files

- [scmore.env.example](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/scmore.env.example)
- [run_scmore_workflow.sh](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/run_scmore_workflow.sh)
- [validate_inputs.R](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/validate_inputs.R)
- [00_prepare_fuma_sumstats.R](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/00_prepare_fuma_sumstats.R)
- [01_prepare_sumstats.R](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/01_prepare_sumstats.R)
- [02_prepare_gene_results.R](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/02_prepare_gene_results.R)
- [03_export_multimodal_object.R](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/03_export_multimodal_object.R)
- [04_run_scmore.R](/Users/obriena2/Documents/GitHub/scScalpChromatin/scmore_melanoma/04_run_scmore.R)

## Quick Start

1. Copy the env template:

```bash
cp scmore_melanoma/scmore.env.example scmore_melanoma/scmore.env
```

2. Fill in the raw melanoma GWAS path and MAGMA/FUMA path.

3. Create a FUMA-ready summary-statistics file now:

```bash
bash scmore_melanoma/run_scmore_workflow.sh prep-fuma
```

4. Normalize the scMORE inputs now:

```bash
bash scmore_melanoma/run_scmore_workflow.sh prep-gwas
bash scmore_melanoma/run_scmore_workflow.sh prep-gene
```

5. Once the ArchR workflow finishes, set `SCMORE_ARCHR_PROJECT_DIR` and build the multimodal object:

```bash
bash scmore_melanoma/run_scmore_workflow.sh build-object
```

6. Run scMORE:

```bash
bash scmore_melanoma/run_scmore_workflow.sh run
```

## scMORE Input Formats

The normalized summary-statistics output is written in scMORE's expected column layout:

- `CHR`
- `POS`
- `ES`
- `SE`
- `LP`
- `AF`
- `SZ`
- `SNP`

The normalized gene-result output is written with scMORE-compatible gene-level fields:

- `GENE`
- `CHR`
- `START`
- `STOP`
- `NSNPS`
- `NPARAM`
- `N`
- `ZSTAT`
- `P`
- optional `SYMBOL`

## FUMA Helper

The FUMA formatter writes:

- `SNP` if available
- `CHR`
- `BP`
- `A1`
- `A2`
- `P`
- optional `BETA`
- optional `OR`
- optional `SE`
- optional `N`

This is a safe generic format for FUMA's SNP2GENE upload page, where you can explicitly map the columns if needed. The FUMA tutorial states that the mandatory fields are a P-value plus either rsID or chromosome and position, while `OR`, `Beta`, and `SE` are optional.

Sources:

- [FUMA tutorial](https://fuma.ctglab.nl/tutorial)
- [FUMA user-group clarification on mandatory columns](https://groups.google.com/g/fuma-gwas-users/c/WqWzy6nUn18)

## Control-Only Recommendation

By default the multimodal object export is configured to keep control samples only:

- `SCMORE_KEEP_DISEASE_STATUS="C_SD,C_PB"`

That is deliberate, because melanoma GWAS should first be interpreted in the less perturbed control chromatin landscape. You can rerun later with:

- `SCMORE_KEEP_DISEASE_STATUS="AA,C_SD,C_PB"`

for a full-dataset sensitivity analysis.

## Install Notes

`scMORE` itself depends on several heavy R packages, including `ArchR`, `Pando`, `Seurat`, and `Signac`. The upstream package metadata currently lists:

- `ArchR`
- `BSgenome.Hsapiens.UCSC.hg38`
- `COSG`
- `GenomicRanges`
- `IRanges`
- `Pando`
- `Seurat`
- `Signac`
- `org.Hs.eg.db`

Source: [scMORE GitHub README](https://github.com/mayunlong89/scMORE) and [DESCRIPTION](https://raw.githubusercontent.com/mayunlong89/scMORE/main/DESCRIPTION).

I recommend using a dedicated R environment for this workflow rather than trying to reuse the RNA preprocessing environment.
