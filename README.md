# nSOS: Next-Generation Specific Oxidative Stress Index Framework

Built upon a non-linear, fold-change mathematical framework, **nSOS** resolves the long-standing issue of **directional cancellation** in multi-biomarker environmental indices. 

Unlike traditional linear indices (e.g., standard IBR) which mathematically mask high toxicity when opposing biological mechanisms (upregulation vs. downregulation) cancel each other out, nSOS integrates absolute directional deviations via a root-mean-square framework.

## Key Features
* **Directional Stability**: Preserves both enzyme induction and severe biochemical depletion/inhibition.
* **Noise Filtering**: Employs a quadratic integration method to allow major biological signal breaks to dominate over baseline environmental noise.
* **Biologically Validated**: Pre-loaded with synthetic tissue-specific validation parameters derived from empirical data for *Mytilus galloprovincialis*.

## Installation & Usage

To run the nSOS framework, source the core script into your R environment:

```R
library(dplyr)
library(tidyr)

# Source the core computational engine
source("nSOS_core.R")

# Generate the biologically validated mock validation dataset
raw_biomarker_data <- generate_mock_data()

# Compute the nSOS scores
results <- calculate_nsos(
  data = raw_biomarker_data, 
  metadata_cols = c("ID", "Site"), 
  biomarker_cols = c("SOD", "CAT", "GST", "LPO"), 
  control_site_name = "Reference"
)

print(head(results))
```

## References & Academic Anchorage
* **Core Preprint & Framework**: Yakimov, L. P. (2026). *The Next-Generation Specific Oxidative Stress Index (nSOS): A Non-Linear Fold-Change Framework to Resolve Directional Cancellation in Aquatic Biomonitoring*. Zenodo. Repository. DOI: [10.5281/zenodo.21535828](https://doi.org)
* **First-Generation Framework**: Yakimov, L., Tsvetanova, E., Georgieva, A., Petrov, L., & Alexandrova, A. (2018). Assessment of the oxidative status of Black Sea mussels (*Mytilus galloprovincialis* Lamark, 1819) from Bulgarian coastal areas with introduction of specific oxidative stress index. *Journal of Environmental Protection and Ecology*, 19(4), 1614–1622.
* **Physiological Baselines**: Regoli, F. (1995). Glutathione, glutathione-dependent and antioxidant enzymes in mussel, *Mytilus galloprovincialis*. *Aquatic Toxicology*, 31, 143–164.
* **Pathophysiological Rationale**: Vlahogianni, T., & Valavanidis, A. (2007). Heavy-metal effects on lipid peroxidation and antioxidant defence enzymes. *Chemical Ecology*, 23, 361–371.

## License
This project is open-source and licensed under the **MIT License**.