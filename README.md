# 🧬 Create_Nextstrain_Tree_MOH

This pipeline generates the **Auspice dataset** for phylogenetic tree visualization on [auspice.us](https://auspice.us), starting from individual FASTA files or a multi-FASTA file.

---

## 📥 Required Inputs

You must provide either **individual FASTA files** *or* a **multi-FASTA file**, along with metadata and a GenBank reference.

### 🧾 Sequences (choose one of the following):
- **Individual FASTA files (as an array)**:
  - In **Terra.bio**:
    ```json
    ["path/to/seq1.fasta", "path/to/seq2.fasta", ...]
    ```
  - For **local use via `.json` file**:
    ```json
    {
      "NextstrainPhylogeny.individual_fasta_files": [
        "path/to/seq1.fasta",
        "path/to/seq2.fasta",
        ...
      ]
    }
    ```

- **Multi-FASTA file**:
  - Provide a single `.fasta` file containing all sequences.

### 🗂 Metadata file
- A `metadata.tsv` file with relevant fields (e.g., strain, date, location).

### 🧬 Reference file
- A GenBank file (`.gb`) containing the reference genome.

---

## 🔍 Optional Filtering Parameters

| Parameter              | Type     | Description                                      |
|------------------------|----------|--------------------------------------------------|
| `min_date`             | String?  | Filter sequences with collection date ≥ this.   |
| `max_date`             | String?  | Filter sequences with collection date ≤ this.   |
| `sequences_per_group` | Int?     | Subsample number of sequences per group.        |
| `group_by`             | String?  | Metadata field to group sequences (e.g., region).|
| `exclude_strains_file`| File?    | File listing strains to exclude.                |

---

## ⚙️ Runtime Parameters (Defaults)

| Parameter       | Default               |
|----------------|------------------------|
| `docker_image`  | `nextstrain/base:latest` |
| `cpu`           | `4`                    |
| `memory`        | `8 GB`                 |
| `disk_size`     | `50`                   |

---

## 📤 Output Files

| File                              | Description                                  |
|-----------------------------------|----------------------------------------------|
| `nextstrain_aligned_fasta`        | Aligned FASTA file                           |
| `nextstrain_phylogenetic_tree`    | Newick tree file                             |
| `nextstrain_ancestral_json`       | JSON with ancestral state reconstruction     |
| `nextstrain_traits_json`          | JSON with trait annotations                  |
| `nextstrain_branch_lengths`       | JSON with estimated branch lengths           |
| `nextstrain_auspice_file`         | Final `.json` file ready for Auspice         |

---
