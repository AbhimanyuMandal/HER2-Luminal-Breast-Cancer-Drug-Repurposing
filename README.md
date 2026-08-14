# HER2-Luminal Breast Cancer Drug Repurposing

### Single-cell RNA-seq driven disease characterization and LINCS-based drug repurposing

An end-to-end computational biology workflow integrating **single-cell RNA-seq analysis, subtype characterization, pseudobulk differential expression, disease-signature construction, pathway analysis, LINCS perturbational screening, cross-cell-line validation, ChEMBL mechanism annotation, and external literature validation** to prioritize potential drug-repurposing candidates for HER2/luminal breast cancer.

> **Important:** This project generates computational drug-repurposing hypotheses. The identified compounds are not experimentally or clinically validated therapeutic recommendations.

---

## Project Overview

The goal of this project was to develop a reproducible computational framework that moves from **single-cell disease characterization to drug-candidate prioritization**.

Rather than starting with a predefined list of breast cancer drugs, the workflow first identifies disease-associated transcriptional programs from breast cancer single-cell RNA-seq data and then uses those programs to search perturbational transcriptomic data for compounds capable of reversing the observed disease state.

The overall workflow is:

```text
Breast Cancer scRNA-seq
        │
        ▼
Quality Control
        │
        ▼
Normalization / HVG / PCA
        │
        ▼
RPCA Integration
        │
        ▼
UMAP / Clustering
        │
        ▼
Cell-type Annotation
        │
        ▼
HER2 / Luminal Characterization
        │
        ▼
Patient-level Cell Composition
        │
        ▼
Pseudobulk Differential Expression
        │
        ▼
Disease Signature Construction
        │
        ├───────────────┐
        ▼               ▼
Pathway Analysis    Signature QC
        │               │
        └───────┬───────┘
                ▼
       LINCS Perturbational Screening
                │
                ▼
       Connectivity Scoring
                │
                ▼
       Cross-cell-line Robustness
                │
                ▼
        66 Robust Candidates
                │
                ▼
       ChEMBL Mechanism Annotation
                │
                ▼
        Final Evidence Ranking
                │
                ▼
       External Literature Validation
                │
                ▼
        Final Candidate Prioritization
```
Key Results
Metric	Result
Extended disease-signature genes	454
Core disease-signature genes	10
Genes mapped to LINCS	315
LINCS signatures screened	9,108
LINCS compounds represented	6,720
Breast cancer cell lines	5
Compounds represented across all 5 cell lines	167
Robust LINCS candidates	66
Mechanism-supported candidates	23
Priority A candidates	3
Priority B candidates	8
Priority C candidates	29
Exploratory candidates	26
Candidates externally screened	20
Biological Question

The central question was:

Can disease-associated transcriptional programs identified from HER2/luminal breast cancer single-cell RNA-seq data be used to computationally prioritize compounds that may reverse those disease-associated expression patterns?

The project therefore combines two complementary perspectives:

Disease characterization using single-cell transcriptomics.
Drug perturbation analysis using LINCS transcriptional signatures.

This creates a data-driven bridge between disease biology and compound prioritization.

1. Single-cell RNA-seq Analysis

The project begins with breast cancer single-cell RNA-seq data and performs a complete single-cell analysis workflow.

Major steps
Quality control
Low-quality cell filtering
Normalization
Highly variable gene identification
PCA
RPCA integration
UMAP dimensionality reduction
Clustering
Cluster marker identification
Cell-type annotation
Annotation validation
Breast cancer subtype characterization

The project also processes a normal breast reference dataset to provide a comparison framework for downstream analyses.

2. Breast Cancer Subtype Characterization

Following cell-type annotation, epithelial populations were characterized across breast cancer subtypes.

The analysis included:

Cell-type composition
Patient-level cell-type proportions
Luminal epithelial populations
Basal epithelial populations
HER2-associated molecular features
Subtype-associated transcriptional differences

The resulting cell-level and patient-level analyses were used to identify transcriptional differences associated with the HER2/luminal disease state.

3. Pseudobulk Differential Expression

To reduce dependence on individual-cell observations, differential expression analysis was performed using patient-level pseudobulk profiles.

Comparisons included breast cancer subtypes and normal breast reference populations across relevant epithelial compartments.

This analysis provided the foundation for defining a disease-associated transcriptional signature.

4. Disease Signature Construction

A HER2/luminal disease-associated transcriptional signature was constructed from the differential expression results.

The final signature contained:

454 genes in the extended disease signature
10 genes in the high-confidence core signature
315 genes successfully mapped to the LINCS expression matrix

The core signature consisted of:

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

The disease signature was separated into directional components:

Disease UP genes
        +
Disease DOWN genes
        ↓
Disease transcriptional state

This directional signature was subsequently used for connectivity-based compound screening.

5. Pathway Analysis

The disease signature was interpreted using Gene Ontology pathway membership analysis.

Representative biological processes included:

Positive regulation of apoptotic process
Response to wounding
Wound healing
Chemotaxis
Epidermis development
Intrinsic apoptotic signaling
Leukocyte migration
Cell adhesion
Regulation of lymphocyte activation

The pathway analysis provides biological context for the disease-associated transcriptional program before drug-repurposing analysis.

6. LINCS Perturbational Screening

The disease signature was used to screen LINCS perturbational transcriptional profiles.

Screening scale

9,108 LINCS signatures

representing:

6,720 unique compounds

across five breast cancer cell lines:

MCF7
SKBR3
BT20
MDAMB231
HS578T

The LINCS expression matrix contained 315 disease-associated genes mapped to the available perturbational profiles.

For each perturbational signature, a connectivity score was calculated to estimate how strongly the compound-associated transcriptional response opposed the disease-associated expression state.

Conceptually:

Disease signature
      │
      ├── UP genes
      └── DOWN genes
             │
             ▼
    LINCS perturbational profile
             │
             ▼
      Connectivity score
             │
             ▼
Potential transcriptional reversal
7. Cross-Cell-Line Robustness

A major component of the workflow was testing whether candidate compounds showed reproducible reversal across multiple breast cancer cell lines.

Among the LINCS compounds:

167 compounds were represented across all five cell lines.

Cross-cell-line connectivity was then used to identify robust reversal candidates.

This resulted in:

66 robust LINCS candidates

Candidates were grouped into three robustness tiers based on the consistency of positive connectivity across cell lines.

This approach reduces dependence on a single cellular context and provides a more stringent candidate-selection strategy.

8. ChEMBL Mechanism Annotation

The robust candidate set was subsequently annotated using ChEMBL.

The annotation workflow evaluated:

Compound identity
ChEMBL identifiers
Mechanism of action
Target information
Action type
Mechanism evidence
Identity confidence

Among the 66 robust candidates:

23 candidates had mechanism evidence.

Candidates without available mechanism records were retained rather than automatically discarded.

Absence of a ChEMBL mechanism record does not mean that a compound has no biological mechanism. It indicates that mechanism evidence was not retrieved for that candidate through this annotation workflow.

9. Final Candidate Ranking

The final prioritization integrates multiple evidence layers:

LINCS connectivity
        +
Cross-cell-line reproducibility
        +
Positive cell-line consistency
        +
Compound identity confidence
        +
Mechanism evidence
        +
External evidence
        │
        ▼
Final candidate prioritization

The final candidate distribution was:

Priority	Number of candidates
Priority A	3
Priority B	8
Priority C	29
Exploratory	26
Priority A candidates
Rank	Compound
1	BMS-777607
2	Motesanib
3	OSI-930

These represent the highest-priority computational candidates according to the project's integrated scoring framework.

They should not be interpreted as clinically validated treatments.

10. External Literature Validation

The top 20 computational candidates were subjected to an external literature validation step.

The purpose was to determine whether the computationally prioritized compounds had previously reported evidence in breast cancer or HER2-related contexts.

The validation considered:

Literature presence
Breast cancer/HER2 evidence
Oncology evidence
Publication counts
External evidence classification
Validation summary
20 candidates evaluated
20 literature queries completed
19 candidates with breast cancer/HER2-related evidence
1 candidate with limited external evidence

External evidence was used as a supporting layer rather than as proof of therapeutic efficacy.

Top Portfolio Candidates

The top-ranked candidates from the final computational prioritization included:

Rank	Compound	Priority	Portfolio Evidence
1	BMS-777607	Priority A	High
2	Motesanib	Priority A	High
3	GSK-3-inhibitor-IX	Priority C	Secondary
4	Pazopanib	Priority B	Moderate
5	GW-5074	Priority C	Secondary
6	OSI-930	Priority A	High
7	PHA-665752	Priority B	Moderate
8	PLX-4720	Priority B	Moderate
9	GW-843682X	Priority C	Secondary
10	Sorafenib	Priority B	Moderate
Core Disease Biology

The high-confidence core signature contained ten genes:

Gene	Direction	logFC
DIO1	UP	9.49
AC005013.5	UP	5.92
SRPK3	UP	5.69
NR4A1	UP	5.46
ATP6V0A4	UP	4.97
EGR4	UP	4.82
TTC6	UP	3.60
TCN1	DOWN	-7.65
LAMB3	DOWN	-7.19
FKBP1B	DOWN	-7.00

The broader disease signature also contained contextual HER2/luminal markers including:

ERBB2
GRB7
FOXA1
AR
Mechanism Categories

Among mechanism-supported candidates, the retrieved mechanism annotations were broadly classified into:

Mechanism category	Candidates
Kinase	11
Other / Unclassified	8
Angiogenesis	2
Cytoskeleton	1
DNA damage response	1

This provides an additional biological perspective on the computationally prioritized candidate set.

Repository Structure
HER2-Luminal-Breast-Cancer-Drug-Repurposing/
│
├── README.md
│
├── assets/
│   └── project figures
│
├── scripts/
│   ├── 01–12  Single-cell breast cancer analysis
│   ├── 13–23  Normal breast reference analysis
│   ├── 24–26  Cohort harmonization and differential analysis
│   ├── 27–29  Disease signature and pathway analysis
│   ├── 30–31  LINCS extraction and connectivity scoring
│   ├── 31B–31C  LINCS QC and robust prioritization
│   ├── 32–32C  ChEMBL annotation and validation
│   ├── 33–33B  Final candidate ranking and QC
│   ├── 34  Biological interpretation
│   ├── 35  External validation
│   └── 36  Final portfolio summary
│
├── results/
│   ├── qc/
│   ├── pca/
│   ├── umap/
│   ├── subtype_characterization/
│   ├── patient_level_composition/
│   ├── pseudobulk_DE/
│   ├── pathway_analysis/
│   ├── disease_signature/
│   └── drug_repurposing/
│
├── docs/
│
└── assets/
Results Included

The repository contains selected analysis outputs from the major stages of the project.

Quality Control
results/qc/

Contains QC summaries and quality-control figures.

Dimensionality Reduction
results/pca/
results/umap/

Contains PCA and UMAP analysis outputs.

Subtype Characterization
results/subtype_characterization/

Contains breast cancer subtype-level characterization results.

Patient-level Composition
results/patient_level_composition/

Contains patient-level cell-type composition analyses.

Differential Expression
results/pseudobulk_DE/

Contains pseudobulk differential expression results and volcano plots.

Pathway Analysis
results/pathway_analysis/

Contains pathway-level interpretation of disease-associated genes.

Disease Signature
results/disease_signature/

Contains disease-signature construction and validation outputs.

Drug Repurposing
results/drug_repurposing/

Contains:

LINCS connectivity results
Robust compound rankings
Compound/cell-line connectivity
ChEMBL compound identities
Mechanism annotations
Final candidate rankings
Priority A/B/C candidate sets
External validation results
Portfolio summary outputs
Scripts

The analysis is organized as a sequential R-based workflow.

Scripts	Analysis stage
01–12	Breast cancer scRNA-seq processing and annotation
13–23	Normal breast reference analysis
24–26	Cohort harmonization, patient composition and pseudobulk DE
27–29	Disease signature and pathway analysis
30–31	LINCS extraction and connectivity scoring
31B–31C	LINCS QC and robust candidate prioritization
32–32C	ChEMBL compound and mechanism annotation
33–33B	Final candidate ranking and QC
34	Biological candidate interpretation
35	External candidate validation
36	Final portfolio summary

All scripts are included to document the analytical workflow and decision points used during the project.

Data Sources

The project uses publicly available transcriptomic datasets and pharmacogenomic resources.

Single-cell RNA-seq

GSE176078

Wu et al. breast cancer single-cell RNA-seq dataset.

GSE113196

Normal breast single-cell RNA-seq reference dataset.

Perturbational transcriptomics

LINCS perturbational expression profiles were used for connectivity-based compound screening.

Compound annotation

ChEMBL was used for compound identity and mechanism-of-action annotation.

External validation

Publicly available literature evidence was used to assess whether top computational candidates had previously reported breast cancer/HER2-related evidence.

Software and Tools

The analysis was primarily implemented in R.

Programming
R
Bash
Linux/macOS command line
Single-cell analysis
Seurat
dplyr
ggplot2
Matrix
Statistical analysis
Pseudobulk differential expression
PCA
RPCA integration
Gene-level differential expression
Multiple-testing correction
Functional analysis
Gene Ontology
Pathway analysis
Drug repurposing
LINCS
Connectivity scoring
ChEMBL
External literature validation
Reproducibility

The repository contains the R scripts used to perform the analysis.

The workflow is organized sequentially, with scripts numbered according to the analytical stage.

Example:

Rscript scripts/01_create_GSE176078_object.R
Rscript scripts/02_qc_GSE176078.R
Rscript scripts/03_apply_QC.R

The complete workflow continues through:

Rscript scripts/36_final_portfolio_summary.R

Large raw and intermediate single-cell objects are intentionally not included in the GitHub repository because of their size.

The datasets can be obtained from their corresponding public GEO accessions.

Why the Workflow Uses Multiple Validation Layers

A major design principle of this project was to avoid selecting compounds using a single metric.

A compound with a strong connectivity score in one cell line may not be robust across biological contexts.

Therefore, candidates were evaluated through multiple stages:

Single-cell disease biology
        ↓
Disease signature
        ↓
LINCS connectivity
        ↓
Cross-cell-line reproducibility
        ↓
Compound identity confidence
        ↓
Mechanism evidence
        ↓
External literature evidence
        ↓
Final prioritization

This does not eliminate biological uncertainty, but it provides a more transparent framework for ranking computational hypotheses.

Limitations

This project has several important limitations.

1. Computational evidence is not experimental validation

Connectivity-based reversal does not establish that a compound will kill tumor cells or produce a therapeutic response.

2. LINCS signatures are context-dependent

Perturbational responses can vary with:

Cell line
Dose
Treatment duration
Experimental conditions
3. Compound annotation is incomplete

Not every LINCS compound could be confidently mapped to a ChEMBL identity or mechanism.

4. Literature evidence can introduce confirmation bias

Previously studied compounds are more likely to have published evidence than poorly studied compounds.

5. Candidate prioritization is hypothesis-generating

The final candidates require experimental validation before any therapeutic interpretation.

Key Takeaway

This project demonstrates an end-to-end computational strategy for connecting single-cell disease biology with perturbational drug discovery.

The workflow progresses from:

Single-cell transcriptomics
        ↓
Cellular and subtype characterization
        ↓
Disease-associated transcriptional signature
        ↓
Perturbational transcriptomic screening
        ↓
Cross-cell-line robustness
        ↓
Mechanism annotation
        ↓
External validation
        ↓
Drug-repurposing candidate prioritization

The analysis ultimately identified:

66 robust computational drug-repurposing candidates, including 3 Priority A candidates, from 9,108 LINCS perturbational signatures representing 6,720 compounds.

The resulting candidate list provides a focused set of hypotheses for future experimental investigation.

Author

Abhimanyu Mandal

Integrated Dual Degree — Pharmaceutical Engineering & Technology
Indian Institute of Technology (BHU), Varanasi

Interests:

Computational Biology
Bioinformatics
Cancer Genomics
Single-cell RNA-seq
Drug Repurposing
Transcriptomics
Data-driven Drug Discovery
Disclaimer

This repository is intended for research, educational, and portfolio purposes.

The computational rankings presented here should not be interpreted as clinical recommendations or evidence of therapeutic efficacy.

No claim of treatment effectiveness is made for any compound identified by this workflow.

Citation

If you use or adapt this workflow, please cite the original data resources and databases used in the analysis, including the corresponding GEO datasets, LINCS/L1000 resources, and ChEMBL.


### One correction before you paste

There is one thing I deliberately **did not claim** in this README: that the three Priority A drugs are "the best drugs" or "effective against HER2 breast cancer." Your analysis doesn't establish that. Calling them **"highest-priority computational candidates"** is accurate and much safer for your CV/portfolio.

Also, I would **not upload the 62 GB `data/` directory or the multi-GB `.rds` objects**. Your current ~90 MB curated repository with scripts + selected results is much better for GitHub.
provide everthing under one. i just want to copy and paste.

Absolutely. Below is the entire README as one single copy-paste block. Paste everything from # HER2... through the end into your GitHub README.md.

# HER2-Luminal Breast Cancer Drug Repurposing

### Single-cell RNA-seq driven disease characterization and LINCS-based drug repurposing

An end-to-end computational biology workflow integrating **single-cell RNA-seq analysis, breast cancer subtype characterization, pseudobulk differential expression, disease-signature construction, pathway analysis, LINCS perturbational screening, cross-cell-line validation, ChEMBL mechanism annotation, and external literature validation** to prioritize potential drug-repurposing candidates for HER2/luminal breast cancer.

> **Important:** This project generates computational drug-repurposing hypotheses. The identified compounds are not experimentally or clinically validated therapeutic recommendations.

---

## Project Overview

The goal of this project was to develop a reproducible computational framework that moves from **single-cell disease characterization to drug-candidate prioritization**.

Rather than starting with a predefined list of breast cancer drugs, the workflow first identifies disease-associated transcriptional programs from breast cancer single-cell RNA-seq data and then uses those programs to search perturbational transcriptomic data for compounds capable of reversing the observed disease state.

The overall workflow is:

```text
Breast Cancer scRNA-seq
        │
        ▼
Quality Control
        │
        ▼
Normalization / HVG / PCA
        │
        ▼
RPCA Integration
        │
        ▼
UMAP / Clustering
        │
        ▼
Cell-type Annotation
        │
        ▼
HER2 / Luminal Characterization
        │
        ▼
Patient-level Cell Composition
        │
        ▼
Pseudobulk Differential Expression
        │
        ▼
Disease Signature Construction
        │
        ├───────────────┐
        ▼               ▼
Pathway Analysis    Signature QC
        │               │
        └───────┬───────┘
                ▼
       LINCS Perturbational Screening
                │
                ▼
       Connectivity Scoring
                │
                ▼
       Cross-cell-line Robustness
                │
                ▼
        66 Robust Candidates
                │
                ▼
       ChEMBL Mechanism Annotation
                │
                ▼
        Final Evidence Ranking
                │
                ▼
       External Literature Validation
                │
                ▼
        Final Candidate Prioritization
Key Results
Metric	Result
Extended disease-signature genes	454
Core disease-signature genes	10
Genes mapped to LINCS	315
LINCS signatures screened	9,108
LINCS compounds represented	6,720
Breast cancer cell lines	5
Compounds represented across all 5 cell lines	167
Robust LINCS candidates	66
Mechanism-supported candidates	23
Priority A candidates	3
Priority B candidates	8
Priority C candidates	29
Exploratory candidates	26
Candidates externally screened	20
Biological Question

The central question was:

Can disease-associated transcriptional programs identified from HER2/luminal breast cancer single-cell RNA-seq data be used to computationally prioritize compounds that may reverse those disease-associated expression patterns?

The project therefore combines two complementary perspectives:

Disease characterization using single-cell transcriptomics.
Drug perturbation analysis using LINCS transcriptional signatures.

This creates a data-driven bridge between disease biology and compound prioritization.

1. Single-cell RNA-seq Analysis

The project begins with breast cancer single-cell RNA-seq data and performs a complete single-cell analysis workflow.

Major steps
Quality control
Low-quality cell filtering
Normalization
Highly variable gene identification
PCA
RPCA integration
UMAP dimensionality reduction
Clustering
Cluster marker identification
Cell-type annotation
Annotation validation
Breast cancer subtype characterization

A normal breast reference dataset is also processed independently to provide a comparison framework for downstream analyses.

2. Breast Cancer Subtype Characterization

Following cell-type annotation, epithelial populations were characterized across breast cancer subtypes.

The analysis included:

Cell-type composition
Patient-level cell-type proportions
Luminal epithelial populations
Basal epithelial populations
HER2-associated molecular features
Subtype-associated transcriptional differences

The resulting cell-level and patient-level analyses were used to identify transcriptional differences associated with the HER2/luminal disease state.

3. Pseudobulk Differential Expression

To reduce dependence on individual-cell observations, differential expression analysis was performed using patient-level pseudobulk profiles.

Comparisons included breast cancer subtypes and normal breast reference populations across relevant epithelial compartments.

This analysis provided the foundation for defining a disease-associated transcriptional signature.

4. Disease Signature Construction

A HER2/luminal disease-associated transcriptional signature was constructed from the differential expression results.

The final signature contained:

454 genes in the extended disease signature
10 genes in the high-confidence core signature
315 genes successfully mapped to the LINCS expression matrix

The core signature consisted of:

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

The disease signature was separated into directional components:

Disease UP genes
        +
Disease DOWN genes
        ↓
Disease transcriptional state

This directional signature was subsequently used for connectivity-based compound screening.

5. Pathway Analysis

The disease signature was interpreted using Gene Ontology pathway membership analysis.

Representative biological processes included:

Positive regulation of apoptotic process
Response to wounding
Wound healing
Chemotaxis
Epidermis development
Intrinsic apoptotic signaling
Leukocyte migration
Cell adhesion
Regulation of lymphocyte activation

The pathway analysis provides biological context for the disease-associated transcriptional program before drug-repurposing analysis.

6. LINCS Perturbational Screening

The disease signature was used to screen LINCS perturbational transcriptional profiles.

Screening scale

9,108 LINCS signatures

representing:

6,720 unique compounds

across five breast cancer cell lines:

MCF7
SKBR3
BT20
MDAMB231
HS578T

The LINCS expression matrix contained 315 disease-associated genes mapped to the available perturbational profiles.

For each perturbational signature, a connectivity score was calculated to estimate how strongly the compound-associated transcriptional response opposed the disease-associated expression state.

Conceptually:

Disease signature
      │
      ├── UP genes
      └── DOWN genes
             │
             ▼
    LINCS perturbational profile
             │
             ▼
      Connectivity score
             │
             ▼
Potential transcriptional reversal
7. Cross-Cell-Line Robustness

A major component of the workflow was testing whether candidate compounds showed reproducible reversal across multiple breast cancer cell lines.

Among the LINCS compounds:

167 compounds were represented across all five cell lines.

Cross-cell-line connectivity was then used to identify robust reversal candidates.

This resulted in:

66 robust LINCS candidates

Candidates were grouped into three robustness tiers based on the consistency of positive connectivity across cell lines.

This approach reduces dependence on a single cellular context and provides a more stringent candidate-selection strategy.

8. ChEMBL Mechanism Annotation

The robust candidate set was subsequently annotated using ChEMBL.

The annotation workflow evaluated:

Compound identity
ChEMBL identifiers
Mechanism of action
Target information
Action type
Mechanism evidence
Identity confidence

Among the 66 robust candidates:

23 candidates had mechanism evidence.

Candidates without available mechanism records were retained rather than automatically discarded.

Absence of a ChEMBL mechanism record does not mean that a compound has no biological mechanism. It indicates that mechanism evidence was not retrieved for that candidate through this annotation workflow.

9. Final Candidate Ranking

The final prioritization integrates multiple evidence layers:

LINCS connectivity
        +
Cross-cell-line reproducibility
        +
Positive cell-line consistency
        +
Compound identity confidence
        +
Mechanism evidence
        +
External evidence
        │
        ▼
Final candidate prioritization

The final candidate distribution was:

Priority	Number of candidates
Priority A	3
Priority B	8
Priority C	29
Exploratory	26
Priority A candidates
Rank	Compound
1	BMS-777607
2	Motesanib
3	OSI-930

These represent the highest-priority computational candidates according to the project's integrated scoring framework.

They should not be interpreted as clinically validated treatments.

10. External Literature Validation

The top 20 computational candidates were subjected to an external literature validation step.

The purpose was to determine whether the computationally prioritized compounds had previously reported evidence in breast cancer or HER2-related contexts.

The validation considered:

Literature presence
Breast cancer/HER2 evidence
Oncology evidence
Publication counts
External evidence classification
Validation summary
20 candidates evaluated
20 literature queries completed
19 candidates with breast cancer/HER2-related evidence
1 candidate with limited external evidence

External evidence was used as a supporting layer rather than as proof of therapeutic efficacy.

Top Portfolio Candidates

The top-ranked candidates from the final computational prioritization included:

Rank	Compound	Priority	Portfolio Evidence
1	BMS-777607	Priority A	High
2	Motesanib	Priority A	High
3	GSK-3-inhibitor-IX	Priority C	Secondary
4	Pazopanib	Priority B	Moderate
5	GW-5074	Priority C	Secondary
6	OSI-930	Priority A	High
7	PHA-665752	Priority B	Moderate
8	PLX-4720	Priority B	Moderate
9	GW-843682X	Priority C	Secondary
10	Sorafenib	Priority B	Moderate
Core Disease Biology

The high-confidence core signature contained ten genes:

Gene	Direction	logFC
DIO1	UP	9.49
AC005013.5	UP	5.92
SRPK3	UP	5.69
NR4A1	UP	5.46
ATP6V0A4	UP	4.97
EGR4	UP	4.82
TTC6	UP	3.60
TCN1	DOWN	-7.65
LAMB3	DOWN	-7.19
FKBP1B	DOWN	-7.00

The broader disease signature also contained contextual HER2/luminal markers including:

ERBB2
GRB7
FOXA1
AR
Mechanism Categories

Among mechanism-supported candidates, the retrieved mechanism annotations were broadly classified into:

Mechanism category	Candidates
Kinase	11
Other / Unclassified	8
Angiogenesis	2
Cytoskeleton	1
DNA damage response	1

This provides an additional biological perspective on the computationally prioritized candidate set.

Repository Structure
HER2-Luminal-Breast-Cancer-Drug-Repurposing/
│
├── README.md
│
├── assets/
│   └── project figures
│
├── scripts/
│   ├── 01–12  Single-cell breast cancer analysis
│   ├── 13–23  Normal breast reference analysis
│   ├── 24–26  Cohort harmonization and differential analysis
│   ├── 27–29  Disease signature and pathway analysis
│   ├── 30–31  LINCS extraction and connectivity scoring
│   ├── 31B–31C  LINCS QC and robust prioritization
│   ├── 32–32C  ChEMBL annotation and validation
│   ├── 33–33B  Final candidate ranking and QC
│   ├── 34  Biological interpretation
│   ├── 35  External validation
│   └── 36  Final portfolio summary
│
├── results/
│   ├── qc/
│   ├── pca/
│   ├── umap/
│   ├── subtype_characterization/
│   ├── patient_level_composition/
│   ├── pseudobulk_DE/
│   ├── pathway_analysis/
│   ├── disease_signature/
│   └── drug_repurposing/
│
├── docs/
│
└── assets/
Results Included

The repository contains selected analysis outputs from the major stages of the project.

Quality Control
results/qc/

Contains QC summaries and quality-control figures.

Dimensionality Reduction
results/pca/
results/umap/

Contains PCA and UMAP analysis outputs.

Subtype Characterization
results/subtype_characterization/

Contains breast cancer subtype-level characterization results.

Patient-level Composition
results/patient_level_composition/

Contains patient-level cell-type composition analyses.

Differential Expression
results/pseudobulk_DE/

Contains pseudobulk differential expression results and volcano plots.

Pathway Analysis
results/pathway_analysis/

Contains pathway-level interpretation of disease-associated genes.

Disease Signature
results/disease_signature/

Contains disease-signature construction and validation outputs.

Drug Repurposing
results/drug_repurposing/

Contains:

LINCS connectivity results
Robust compound rankings
Compound/cell-line connectivity
ChEMBL compound identities
Mechanism annotations
Final candidate rankings
Priority A/B/C candidate sets
External validation results
Portfolio summary outputs
Scripts

The analysis is organized as a sequential R-based workflow.

Scripts	Analysis stage
01–12	Breast cancer scRNA-seq processing and annotation
13–23	Normal breast reference analysis
24–26	Cohort harmonization, patient composition and pseudobulk DE
27–29	Disease signature and pathway analysis
30–31	LINCS extraction and connectivity scoring
31B–31C	LINCS QC and robust candidate prioritization
32–32C	ChEMBL compound and mechanism annotation
33–33B	Final candidate ranking and QC
34	Biological candidate interpretation
35	External candidate validation
36	Final portfolio summary

All scripts are included to document the analytical workflow and decision points used during the project.

Data Sources

The project uses publicly available transcriptomic datasets and pharmacogenomic resources.

Single-cell RNA-seq

GSE176078

Wu et al. breast cancer single-cell RNA-seq dataset.

GSE113196

Normal breast single-cell RNA-seq reference dataset.

Perturbational transcriptomics

LINCS perturbational expression profiles were used for connectivity-based compound screening.

Compound annotation

ChEMBL was used for compound identity and mechanism-of-action annotation.

External validation

Publicly available literature evidence was used to assess whether top computational candidates had previously reported breast cancer/HER2-related evidence.

Software and Tools

The analysis was primarily implemented in R.

Programming
R
Bash
Linux/macOS command line
Single-cell analysis
Seurat
dplyr
ggplot2
Matrix
Statistical analysis
Pseudobulk differential expression
PCA
RPCA integration
Gene-level differential expression
Multiple-testing correction
Functional analysis
Gene Ontology
Pathway analysis
Drug repurposing
LINCS
Connectivity scoring
ChEMBL
External literature validation
Reproducibility

The repository contains the R scripts used to perform the analysis.

The workflow is organized sequentially, with scripts numbered according to the analytical stage.

Example:

Rscript scripts/01_create_GSE176078_object.R
Rscript scripts/02_qc_GSE176078.R
Rscript scripts/03_apply_QC.R

The workflow continues through:

Rscript scripts/36_final_portfolio_summary.R

Large raw and intermediate single-cell objects are intentionally not included in the GitHub repository because of their size.

The datasets can be obtained from their corresponding public GEO accessions.

Why the Workflow Uses Multiple Validation Layers

A major design principle of this project was to avoid selecting compounds using a single metric.

A compound with a strong connectivity score in one cell line may not be robust across biological contexts.

Therefore, candidates were evaluated through multiple stages:

Single-cell disease biology
        ↓
Disease signature
        ↓
LINCS connectivity
        ↓
Cross-cell-line reproducibility
        ↓
Compound identity confidence
        ↓
Mechanism evidence
        ↓
External literature evidence
        ↓
Final prioritization

This does not eliminate biological uncertainty, but it provides a more transparent framework for ranking computational hypotheses.

Limitations

This project has several important limitations.

1. Computational evidence is not experimental validation

Connectivity-based reversal does not establish that a compound will kill tumor cells or produce a therapeutic response.

2. LINCS signatures are context-dependent

Perturbational responses can vary with:

Cell line
Dose
Treatment duration
Experimental conditions
3. Compound annotation is incomplete

Not every LINCS compound could be confidently mapped to a ChEMBL identity or mechanism.

4. Literature evidence can introduce confirmation bias

Previously studied compounds are more likely to have published evidence than poorly studied compounds.

5. Candidate prioritization is hypothesis-generating

The final candidates require experimental validation before any therapeutic interpretation.

Key Takeaway

This project demonstrates an end-to-end computational strategy for connecting single-cell disease biology with perturbational drug discovery.

The workflow progresses from:

Single-cell transcriptomics
        ↓
Cellular and subtype characterization
        ↓
Disease-associated transcriptional signature
        ↓
Perturbational transcriptomic screening
        ↓
Cross-cell-line robustness
        ↓
Mechanism annotation
        ↓
External validation
        ↓
Drug-repurposing candidate prioritization

The analysis ultimately identified:

66 robust computational drug-repurposing candidates from 9,108 LINCS perturbational signatures representing 6,720 compounds, including 3 Priority A candidates and 20 candidates subjected to external validation.

The resulting candidate list provides a focused set of hypotheses for future experimental investigation.

Author

Abhimanyu Mandal

Integrated Dual Degree — Pharmaceutical Engineering & Technology
Indian Institute of Technology (BHU), Varanasi

Areas of Interest
Computational Biology
Bioinformatics
Cancer Genomics
Single-cell RNA-seq
Drug Repurposing
Transcriptomics
Data-driven Drug Discovery
Disclaimer

This repository is intended for research, educational, and portfolio purposes.

The computational rankings presented here should not be interpreted as clinical recommendations or evidence of therapeutic efficacy.

No claim of treatment effectiveness is made for any compound identified by this workflow.

Citation

If you use or adapt this workflow, please cite the original data resources and databases used in the analysis, including the corresponding GEO datasets, LINCS/L1000 resources, and ChEMBL.


**One recommendation:** after pasting this, don't spend more time adding more result files. Your repository is already sufficient for a portfolio. The next useful step is to add **4–5 strong figures to `assets/` and embed them near the relevant README sections** so recruiters can understand the project visually within 30 seconds.
