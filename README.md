# LoRNA-SH Manuscript

This repository contains pre-processing and downstream analysis scripts for the LoRNA-SH project. It includes scripts for processing raw RNA sequencing data and analyzing model predictions.

## Repository Structure

```
lorna-manuscript
├── pre-lornash/            # Pre-processing scripts
│   ├── mpaqt/             # MPAQT analysis (joint long+short reads)
│   ├── long-reads/         # Long-read RNA-seq processing
│   │   ├── bambu/         # Bambu pipeline
│   │   ├── espresso/      # ESPRESSO pipeline
│   │   ├── isoquant/      # IsoQuant pipeline
│   │   ├── pigeon/        # PIGEON pipeline
│   │   └── talon/         # TALON pipeline
│   └── short-reads/        # Short-read RNA-seq processing
│       ├── ggsashimi/     # Sashimi plots
│       ├── kallisto/      # Kallisto quantification
│       ├── miso/          # MISO analysis
│       ├── psi-sigma/     # PSI-Sigma analysis
│       ├── rmats/         # rMATS analysis
│       ├── star/          # STAR alignment
│       └── suppa/         # SUPPA analysis
└── post-lornash/          # Downstream analysis
    ├── probing/           # Model probing analyses
    │   ├── cell-line-specific-abundance/
    │   └── cell-line-specific-upregulation/
    └── zero-shot/         # Zero-shot prediction analyses
        ├── abundance/
        ├── clinvar/
        ├── depmap/
        ├── exon-trap/
        ├── motif-discovery/
        ├── rbp-binding-sites/
        └── splice-sites/
```

## Pre-LoRNA-SH Processing

### Long-read RNA-seq Processing
- **Bambu**: Transcript discovery and quantification
- **ESPRESSO**: Isoform detection and quantification
- **IsoQuant**: Long-read quantification
- **Pigeon**: Transcript annotation
- **TALON**: Transcript quantification and annotation

### Short-read RNA-seq Processing
- **Kallisto**: Transcript quantification
- **STAR**: Read alignment
- **rMATS**: Alternative splicing analysis
- **MISO**: Isoform analysis
- **SUPPA**: Alternative splicing analysis

## Post-LoRNA-SH Analysis

### Probing Analysis
- Cell-line-specific abundance prediction
- Cell-line-specific upregulation analysis

### Zero-shot Predictions
- Abundance prediction
- ClinVar non-coding variant analysis
- DepMap somatic mutations effect prediction
- Exon trap analysis
- Motif discovery analysis
- RBP binding site importance analysis
- Splice site analysis
