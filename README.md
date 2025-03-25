# 🧪 Minimal Example: `reduce` Operator error in Nextflow

This repository demonstrates a **minimal reproducible example** to highlight a potential issue with the `reduce` operator in [Nextflow](https://www.nextflow.io/).

## 🚀 Usage

To run the workflow:

```bash
nextflow run reduce.nf --samplesheet assets/Test1_Normal.csv --genome_version hg19
```

### ✅ First Run (Expected Behavior)

- The input sample sheet contains **only exome** data.
- The workflow runs as expected:
  - The `Exome_only` sub-workflow is triggered.
  - `work/` directory contains `Exome.csv`.
  - The `RNAseq_only` sub-workflow is **not triggered**.

## ⚠️ Demonstrating the Bug

1. Open `workflows/RNAseq_only.nf`.
2. **Uncomment lines 54–57** (the `reduce` function section).
3. Re-run the pipeline:

```bash
nextflow run reduce.nf --samplesheet assets/Test1_Normal.csv --genome_version hg19
```

### ❗ Unexpected Behavior

- Despite having **no RNA-seq input**, a part of `RNAseq_only` workflow still runs.
