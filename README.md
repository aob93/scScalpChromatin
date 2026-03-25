# scScalpChromatin

Repository to host code for ["Integrated single-cell chromatin and transcriptomic analyses of human scalp identify gene-regulatory programs and critical cell types for hair and skin diseases" by Ober-Reynolds et al. Nature Genetics. 2023.](https://doi.org/10.1038/s41588-023-01445-4)

Scripts for performing initial scATAC data processing and subclustering are found in scATAC directory. Scripts for performing initial scRNA data processing and subclustering are found in the scRNA directory. Numbers preceeding these scripts indicate the order in which scripts should be used.

The preprocessing scripts support environment-based path configuration for cluster runs. Set `SCSCALP_PROJECT_ROOT` and optionally `SCSCALP_RESULTS_ROOT`, `SCSCALP_RAW_RNA_DIR`, `SCSCALP_RAW_ATAC_FRAGMENT_DIR`, and `SCSCALP_DEMUXLET_BEST` before running. Biowulf helpers are in `biowulf/run_pipeline.sh`, `biowulf/scscalp_pipeline.sbatch`, and `biowulf/biowulf.env.example`.

On Biowulf, copy `biowulf/biowulf.env.example` to `biowulf/scscalp.env`, edit the paths once, and then run either `bash biowulf/run_pipeline.sh all` interactively or `sbatch biowulf/scscalp_pipeline.sbatch` as a batch job. Both scripts will source `biowulf/scscalp.env` automatically if it exists.

If you want the local clone to remain the source of truth and only mirror files to Biowulf, use `biowulf/sync_to_biowulf.sh` with `BIOWULF_HOST` and `BIOWULF_REPO_DIR` set in your local shell.

If you want strict runtime package checks, export variables such as `SCSCALP_PKG_SEURAT='>=5.0.0'` or `SCSCALP_PKG_ARCHR='==1.0.3'` before launching. The scripts will stop early if the loaded package versions do not satisfy those requirements.

Fine-mapped SNP analyses (LDSR, intersection with peak-to-gene linkages, and GkmSVM models) are found in the GWAS directory.

Processed and annotated Seurat objects for scRNA-seq data can be found [here](https://drive.google.com/drive/folders/1klScH010LvYxdU-TZiSkGVip9IrilZN2?usp=sharing). This includes the objects for the subclustered keratinocytes, fibroblasts, endothelial, lymphoid, and myeloid groups. (We may change the location of these files in the future if we find a better place to host them, but will update this page with the location.)
