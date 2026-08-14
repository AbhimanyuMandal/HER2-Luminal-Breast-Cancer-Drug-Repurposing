<p align="center">
  <img src="assets/banner.png" alt="HER2 Luminal Breast Cancer Drug Repurposing Banner" width="100%">
</p>

<h1 align="center">
🧬 HER2-Luminal Breast Cancer Drug Repurposing
</h1>

<p align="center">

Single-cell RNA-seq-driven disease characterization and LINCS-based computational drug repurposing

A reproducible computational biology workflow integrating single-cell RNA sequencing, subtype characterization, differential expression, disease-signature construction, LINCS perturbational screening, cross-cell-line validation, ChEMBL mechanism annotation, and external literature validation to prioritize candidate compounds for HER2/luminal breast cancer.

</p>

<p align="center">

<img src="https://img.shields.io/badge/R-4.x-blue?logo=r" />
<img src="https://img.shields.io/badge/Seurat-scRNA--seq-blueviolet" />
<img src="https://img.shields.io/badge/LINCS-Connectivity%20Analysis-orange" />
<img src="https://img.shields.io/badge/ChEMBL-Mechanism%20Annotation-green" />
<img src="https://img.shields.io/badge/Drug%20Repurposing-Computational-red" />
<img src="https://img.shields.io/badge/License-MIT-brightgreen" />

</p>

---

# Overview

Breast cancer is a molecularly heterogeneous disease in which distinct transcriptional states and molecular subtypes can exhibit different biological characteristics and therapeutic vulnerabilities.

This project develops an end-to-end computational framework to identify potential drug-repurposing candidates associated with the HER2/luminal breast cancer transcriptional phenotype.

The workflow begins with single-cell RNA-seq data and progressively narrows the analysis from cellular characterization to a disease-specific transcriptional signature and finally to candidate drug prioritization.

The pipeline integrates:

- Single-cell RNA-seq quality control
- Normalization and highly variable gene analysis
- PCA, RPCA integration, UMAP, and clustering
- Cell-type annotation and validation
- Breast cancer subtype characterization
- Normal breast tissue comparison
- Patient-level cell-type composition analysis
- Pseudobulk differential expression
- HER2/luminal disease-signature construction
- LINCS perturbational signature screening
- Connectivity-based compound prioritization
- Cross-cell-line robustness analysis
- ChEMBL compound and mechanism annotation
- External literature validation
- Final portfolio-oriented candidate ranking

The final analysis prioritized **66 robust candidate compounds**, including **3 Priority A candidates**, after integrating transcriptional reversal, cross-cell-line reproducibility, compound identity, mechanism evidence, and external validation.

> **Important:** This project represents computational drug-repurposing prioritization. The results do not establish clinical efficacy, safety, or therapeutic suitability.

---

# Project Objectives

The project aims to:

- Characterize breast cancer cellular populations using scRNA-seq.
- Identify and validate major epithelial and non-epithelial cell populations.
- Characterize breast cancer molecular subtypes.
- Compare breast cancer transcriptional states with normal breast tissue.
- Quantify patient-level cellular composition.
- Perform cell-type-aware pseudobulk differential expression analysis.
- Construct a disease-specific HER2/luminal transcriptional signature.
- Identify genes associated with the disease state.
- Map disease genes to LINCS perturbational data.
- Screen thousands of compound-induced transcriptional signatures.
- Quantify transcriptional connectivity between disease and drug perturbations.
- Identify compounds capable of reversing the disease-associated transcriptional state.
- Evaluate candidate reproducibility across multiple breast cancer cell lines.
- Annotate candidate compounds using ChEMBL.
- Integrate mechanism-of-action evidence.
- Perform external literature validation of high-ranked candidates.
- Produce a transparent final candidate ranking suitable for downstream investigation.

---

# Project Highlights

- ✅ End-to-End scRNA-seq Analysis
- ✅ Multiple Breast Cancer Cohorts
- ✅ HER2/Luminal Subtype Characterization
- ✅ Normal Breast Tissue Comparison
- ✅ Patient-Level Cell-Type Composition
- ✅ Pseudobulk Differential Expression
- ✅ Disease-Signature Construction
- ✅ LINCS Connectivity Screening
- ✅ 9,108 Perturbational Signatures Screened
- ✅ 6,720 LINCS Compounds Represented
- ✅ Five Breast Cancer Cell Lines
- ✅ Cross-Cell-Line Robustness Analysis
- ✅ 66 Robust Drug-Reversal Candidates
- ✅ ChEMBL Mechanism Annotation
- ✅ External Literature Validation
- ✅ Final Evidence-Based Candidate Prioritization
- ✅ Reproducible R-Based Workflow

---

# Repository Workflow

The complete computational workflow follows a progressive narrowing strategy:

<p align="center">
  <img src="assets/her2_workflow.png" width="90%" alt="HER2 Luminal Breast Cancer Drug Repurposing Workflow">
</p>

---

# Single-Cell RNA-seq analysis

The first stage of the project focuses on the characterization of breast cancer single-cell transcriptomic data. The workflow was designed to preserve biological structure while minimizing technical variation b/w samples.

The analysis includes:

- Raw data import
- Seurat object construction
- Quality-control assessment
- Filtering of low-quality cells
- Normalization
- Highly variable gene identification
- Scaling
- PCA
- Principal component diagnostics
- RPCA-based integration
- UMAP dimensionality reduction
- Graph-based clustering
- Cluster marker identification
- Cell-type annotation
- Annotation validation

---

# Breast Cancer Subtype Characterization

Following cellular annotation, epithelial populations were characterized to identify breast cancer molecular subtypes.

The analysis evaluates transcriptional patterns associated with:

- HER2-positive disease
- ER-positive/luminal disease
- Triple-negative breast cancer
- Other epithelial states

HER2/luminal characterization was subsequently used as the biological context for the downstream drug-repurposing analysis.

---

# Normal Breast Tissue Comparison

A separate normal breast single-cell dataset was processed using a compatible analytical framework.

The normal cohort was used to provide a non-malignant reference for:

- Cellular composition
- Cell-type annotation
- Marker expression
- UMAp/clustering structure
- Differential expression
- Disease-state comparison

This comparison provides additional biological context for identifying disease-associated transcriptional changes.

---

# Patient-Level Cell-Type Composition

Cell-type composition was evaluated at the patient level rather than relying exclusively on pooled cell counts.

This analysis helps distinguish:

- Differences in cellular abundance
- Patient-to-patient variability
- Subtype-associated composition changes
- Potential biological heterogeneity between samples

Patient-level summaries were retained as downstream results for visualization and interpretation.

---

# Pseudobulk Differential Expression

Cell-type-aware pseudobulk differential expression was performed to identify transcriptional differences between breast cancer groups and normal tissue.

Comparisons included:

- HER2-positive vs normal
- ER-positive vs normal
- TNBC vs normal
- HER2-positive vs ER-positive
- TNBC vs HER2-positive
- TNBC vs ER-positive

Pseudobulk analysis aggregates single-cell expression within biologically meaningful sample/cell-type units, providing a more appropriate framework for patient-level differential expression than treating individual cells as independent biological replicates.

---

# Disease Signature Construction

The HER2/luminal disease state was used to construct a disease-associated transcriptional signature.

The finalized extended disease signature contained:

### 454 genes

The signature was divided into:

- Disease-upregulated genes
- Disease-downregulated genes

A smaller high-confidence core signature was also retained for biological interpretation.

## Core Signature

The final core signature contained:

### 10 genes

The core genes included:

```text
DIO1
AC005013.5
SRPK3
NR4A1
ATP6V0A4
EGR4
TTC6
TCN1
LAMB3
FKBP1B
```

The contextual HER2/luminal markers included:

```text
ERBB2
GRB7
FOXA1
AR
```

These genes provide biological context for the transcriptional state used in the downstream connectivity analysis.

---

# LINCS Perturbational Screening

The disease signature was mapped to LINCS perturbational expression data to identify compounds that induce transcriptional states opposing the disease-associated transcriptional program.

The analysis used:

- **315 disease genes** successfully mapped to LINCS
- **9,108 perturbational signatures**
- **6,720 unique compounds**
- **5 breast cancer cell lines**

The five cell lines represented in the analysis were:

```text
BT20
HS578T
MCF7
MDAMB231
SKBR3
```

These LINCS signatures represented compound perturbations under the selected screening conditions.

---

# Connectivity Scoring

For each LINCS signature, a connectivity score was calculated to quantify the relationship between the disease transcriptional signature and the compound-induced transcriptional response.

Conceptually:
<p align="center">
  <img src="assets/connectivity_scoring.png" width="90%" alt="HER2 Luminal Breast Cancer Drug Repurposing Workflow">
</p>

Compounds producing reproducible positive reversal scores across relevant cell lines were prioritized for downstream analysis.

---

# Cross-Cell-Line Robustness

A major component of the prioritization strategy was cross-cell-line reproducibility.

Among the LINCS compounds:

**167 compounds** had perturbational signatures represented across all five breast cancer cell lines.

These compounds were evaluated based on:

- Number of positive cell lines
- Number of negative cell lines
- Median connectivity score
- Mean connectivity score
- Minimum connectivity score
- Reproducibility
- Positive connectivity fraction

**66 robust candidate compounds**

Candidates were divided into three tiers:

| Tier | Interpretation |
|-----------|------------|
| Tier 1 | Consistent reversal across 5/5 cell lines |
| Tier 2 | Strong reversal across >= 4/5 cell lines |
| Tier 3 | Moderate reversal across >= 3/5 cell lines | 

The final robust candidate set was then carried forward for compound annotation and evidence integration.

---

# ChEMBL Compound Annotation

The 66 robust candidates were mapped to ChEMBL to obtain additional chemical and pharmacological information.

The annotation workflow included:

- Compound identity matching
- ChEMBL identifier assignment
- Preferred compound name
- Clinical development phase where available
- Mechanism-of-action records
- Target identifiers
- Action types

The ChEMBL validation identified:

- **66 candidate compounds**
- **64 unique ChEMBL identifiers**
- **30 high-confidence compound identities**
- **36 identities requiring review**
- **23 candidates with mechanism evidence**
- **43 candidates without mechanism records**

ChEMBL annotation was used as supporting evidence rather than as proof of therapeutic activity.

---

# Mechansim Interpretation

Mechanism-supported candidates were grouped into broad mechanism categories.

The mechanism-supported candidate set contained:

| Mechanism Category | Candidates |
|-----------|------------|
| Kinase | 11 |
| Angiogenesis | 2 |
| Cytoskeleton | 1 | 
| DNA Damage Response | 1 |
| Other/Unclassified | 8 |

The dominance of kinase-associated candidates provides a biologically interpretable pattern within the computationally prioritized compound set.

However, mechanism annotation does not establish that a compound will therapeutically reverse the disease phenotype in patients.

---

# External Literature Validation

The top 20 computationally ranked candidates were subjected to an external literature screening step.

The screening evaluated evidence related to:

- Breast cancer
- HER2-associated research
- Oncology
- Published evidence for the candidate compound

Results:

- **20 candidates evaluated**
- **20 literature queries completed**
- **19 candidates with breast cancer/HER2-related evidence**
- **1 candidate with limited external evidence**

External evidence was incorporated into the final portfolio-oriented interpretation.

---




