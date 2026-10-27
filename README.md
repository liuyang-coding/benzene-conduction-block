# Benzene exposure and cardiac conduction block

Code examples for the epidemiological and bioinformatics analyses.
Author: Yang Liu. Updated: 27 October 2026.

## Analysis scripts

| Script | Analysis |
|---|---|
| [epidemiological_analysis.R](epidemiological_analysis.R) | Baseline characteristics, Cox regression, RCS and subgroup analysis |
| [sensitivity_analysis.R](sensitivity_analysis.R) | Sensitivity analyses, MICE, Fine-Gray and SCAM |
| [genetic_susceptibility.R](genetic_susceptibility.R) | PRS associations, joint effects and interactions |
| [cmr_associations.R](cmr_associations.R) | Benzene-CMR and candidate omics-CMR associations |
| [nmr_analysis.R](nmr_analysis.R) | Metabolite screening, mediation and correlations |
| [olink_analysis.R](olink_analysis.R) | Protein screening and mediation |
| [functional_enrichment.R](functional_enrichment.R) | GO, KEGG and Reactome enrichment |
| [ppi_analysis.R](ppi_analysis.R) | STRING networks, modules and hub proteins |
| [single_cell_expression.py](single_cell_expression.py) | Expression in reference heart cell types |
| [spatial_localization.py](spatial_localization.py) | Spatial expression and paired tissue comparisons |
| [additional_statistical_methods.R](additional_statistical_methods.R) | CB-free life expectancy, PIF, weight-capped IPW and linear GPS adjustment |
| [e_value.do](e_value.do) | E-values |

## Use

The examples use prepared analysis data. `Benzene` is the exposure input and `PRS` is the genetic score. Follow the analysis sections and example calls in each script. Study data are not included.

Software: R 4.5.2, Python 3.13.11 and Stata 18. Required packages are listed in the scripts.
