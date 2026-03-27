# Paper Reference Workflow

This reduced workflow starts from:

- RNA reference Seurat objects
- scATAC fragment files

It skips the raw scRNA preprocessing stack and runs the ArchR-based ATAC workflow plus RNA-guided integration used in the paper.

## What It Runs

The launcher wraps the existing scripts in this repo:

- `atac_preprocessing/01_prepare_archr_proj.R`
- `atac_preprocessing/02a_subgroup_archr.R`
- `atac_preprocessing/02b_recluster_*.R`
- `atac_preprocessing/03_call_peaks_archr.R`
- `atac_preprocessing/03b_p2gLinks_*.R`
- `atac_preprocessing/04_full_RNA_integration.R`

Optional downstream keratinocyte / hair-follicle analyses remain available as a separate stage:

- `atac_preprocessing/04b_Keratinocytes_analyses.R`
- `atac_preprocessing/04c_HF_Kc_analysis.R`

## Required Inputs

### ATAC

Set `SCSCALP_RAW_ATAC_FRAGMENT_DIR` to a directory containing one folder per sample. The workflow supports:

- `sample/fragments.tsv.gz`
- `sample/atac_fragments.tsv.gz`
- `sample/outs/fragments.tsv.gz`
- `sample/outs/atac_fragments.tsv.gz`

### RNA

Set `SCSCALP_RNA_REFERENCE_RDS` to the full RNA Seurat object used for cross-modality integration.

Set either:

- `SCSCALP_RNA_REFERENCE_SUBCLUSTER_DIR`

or subgroup-specific overrides:

- `SCSCALP_RNA_REFERENCE_SUBGROUP_ENDOTHELIAL_RDS`
- `SCSCALP_RNA_REFERENCE_SUBGROUP_FIBROBLASTS_RDS`
- `SCSCALP_RNA_REFERENCE_SUBGROUP_KERATINOCYTES_RDS`
- `SCSCALP_RNA_REFERENCE_SUBGROUP_LYMPHOID_RDS`
- `SCSCALP_RNA_REFERENCE_SUBGROUP_MYELOID_RDS`

If `SCSCALP_RNA_REFERENCE_SUBCLUSTER_DIR` is used, the workflow looks for either:

- `SUBCLUSTER_DIR/<Subgroup>/<Subgroup>.rds`
- `SUBCLUSTER_DIR/<Subgroup>.rds`

Expected subgroup names are `Endothelial`, `Fibroblasts`, `Keratinocytes`, `Lymphoid`, and `Myeloid`.

## Required RNA Metadata

The full RNA reference must contain:

- `NamedClust`
- `BroadClust` or a `NamedClust` column from which it can be derived
- `Sample` or `orig.ident`
- `diseaseStatus` or a `Sample` column from which it can be derived
- a `umap` reduction

Each subgroup RNA reference must contain:

- `FineClust`
- `Sample` or `orig.ident`
- `diseaseStatus` or a `Sample` column from which it can be derived
- a `umap` reduction

## Quick Start

1. Copy the env template:

```bash
cp paper_reference_workflow/reference.env.example paper_reference_workflow/reference.env
```

2. Edit the paths in `paper_reference_workflow/reference.env`.

3. Validate the inputs:

```bash
bash paper_reference_workflow/run_reference_workflow.sh validate
```

4. Run the core paper-style workflow:

```bash
bash paper_reference_workflow/run_reference_workflow.sh all
```

5. Optionally run the extra keratinocyte / hair-follicle analyses:

```bash
bash paper_reference_workflow/run_reference_workflow.sh advanced
```

## Stages

- `validate`: check fragment discovery and RNA reference objects
- `baseline`: create Arrow files, ArchR project, QC-filter, cluster the full ATAC dataset
- `subcluster`: create ATAC subgroup projects and RNA-guided subgroup re-clustering
- `peaks`: call peaks and subgroup-specific peak-to-gene analyses
- `integrate`: integrate the full ATAC project with the full RNA reference
- `advanced`: run the extra keratinocyte / hair-follicle downstream analyses
- `all`: `baseline` + `subcluster` + `peaks` + `integrate`

## Notes

- `03_call_peaks_archr.R` requires a working `macs2` installation visible to ArchR.
- The launcher writes per-script logs under `SCSCALP_RESULTS_ROOT/pipeline_logs/reference_workflow/`.
- This workflow reuses the existing project code; it does not duplicate the ArchR analysis logic.
