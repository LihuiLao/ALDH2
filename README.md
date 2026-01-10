# Python and R codes for: Aldehyde dehydrogenases are dual-function enzymes that oxidizes both alcohols and aldehydes

Mammalian ethanol metabolism is classically defined as a compartmentalized, two-step process: cytosolic alcohol dehydrogenases (ADH) oxidize ethanol to acetaldehyde, which is subsequently detoxified by mitochondrial aldehyde dehydrogenase 2 (ALDH2). Here, we challenge this canonical model by demonstrating that ALDH2 possesses robust, intrinsic alcohol dehydrogenase activity, mediating a previously unrecognized ethanol oxidation pathway within mammalian mitochondria. We show that this dual functionality is not unique to ALDH2 but is an evolutionarily conserved feature of the entire ALDH superfamily, effectively redefining these enzymes as broad-spectrum alcohol dehydrogenases rather than solely aldehyde processors. Consistent with this expanded role, we identified over 100 endogenous hydroxyl-containing metabolites that accumulate following the loss of ALDH2 activity. These findings provide a new mechanistic framework for understanding the pathology of the ALDH2*2 alleles that are linked to increased risks of cancer and cardiovascular disease.

This repository contains the Python and R scripts used to generate the figures in the manuscript.



## 1. System Requirements

### Operating System
The scripts have been tested on:
- Ubuntu 22.04.3

### Software Versions
To reproduce the figures exactly as shown in the paper, we recommend using the specific versions listed below, although newer versions may also work.

#### For R Scripts:
- **R Base**: 4.4.2
- **Key Libraries**:
  - ggplot2 3.5.2
  - dplyr 1.14
  - Biostrings 2.74.1
  - limma 3.62.2
  - ggseqlogo 0.2
  - pheatmap 1.0.12
  - psych 2.5.3
  - RColorBrewer 1.13
  - readxl 1.4.5
  - writexl 1.5.4
  - reshape2 1.4.4

#### For Python Scripts:
- **Python**: 3.12.11
- **Key Libraries**:
  - matplotlib 3.10.7
  - matplotlib-venn 1.1.2  
  - pandas 2.3.1
  - biopython 1.85
  - numpy 2.3.2
  - seaborn 0.13.2
  - rdkit 2025.9.1

*(No non-standard hardware is required.)*

## 2. Installation Guide

Since these are standalone scripts, no compilation or complex installation is required.

1. **Download**: Clone this repository.
   
   ```bash
   git clone https://github.com/LihuiLao/ALDH2

## 3. How to run

You can run the scripts directly in RStudio and via command line.

**Example:**

1. Open your terminal or command prompt.
2. Navigate to the directory containing the script and the dataset.
3. Run the following command

```bash
python OH_Peak_Pair_extract.py OH_Peak_Pair_Extract_Dataset.xlsx
```

### Expected Output

Upon successful execution, the script will generate:

1. **Processed Data**: Output files (CSV or Excel) containing the extracted peak pairs.