# Species Sampling Problems with Ordering

This repository contains the code for the paper "A Bayesian Nonparametric Approach to Species Sampling Problems with Ordering" (2022) by Balocchi, Camerlenghi and Favaro [arxiv](https://arxiv.org/abs/2203.07342).

## Script to replicate numerical evaluation of the ordered PYP using synthetic and real data

All functions used for the analysis are contained in `scripts/funs.R` and `scripts/funs.cpp`. 
Additionally, `scripts/fun_lgfact.cpp` contains a function to compute log generalized factorial numbers, which are needed for the distribution of the number of distinct species; this function was taken from [here](https://github.com/danieledurante/ESBM/blob/master/Source/stirling.cpp) - many thanks to Tommaso Rigon! 

## Simulations

Simulations were ran and analyzed using different scripts:
- `scripts/simulation_from_model.R` calls `scripts/simulate_model.R` to generate sevaral synthetic datasets from the model; then analyzes and saves the results.
- `scripts/simulation_from_dir.R` calls `scripts/simulate_dir.R` to generate sevaral synthetic datasets with a DP clustering distribution; then analyzes and saves the results.
- `scripts/simulation_from_PYP.R` calls `scripts/simulate_pyp.R` to generate sevaral synthetic datasets with a PYP clustering distribution; then analyzes and saves the results.
- `scripts/simulation_from_zipf.R` calls `scripts/simulate_zipf.R` to generate sevaral synthetic datasets with a Zipf clustering distribution; then analyzes and saves the results.
	
The script `scripts/simulations_plots.R` is used to combine the results and produce plots and tables.

Note: the datasets with the compiled results are included in `results/`: `results_model.RData`, `results_DP.RData`, `results_PYP.RData`, `results_zipf.RData`. Thus `simulations_plots.R` can be ran without re-running all the simulations.

## Real data: 1000 Genomes and Human Genome Dating

We analyzed the variation for genes BRCA2 and EDAR combining datasets from the 1000 Genomes Project and the Human Genome Dating Projects.
To download the data, run (and follow instructions) in `scripts/get_data_clean_BRCA.R` and `scripts/get_data_clean_EDAR.R`.
The obtained cleaned data files (not the raw files) are also saved into `data/` (`data/output_cleandata_BRCA.Rdata`, `data/output_cleandata_EDAR.Rdata`, `data/Gen1000_cleaned_BRCA.Rdata` and `data/Gen1000_cleaned_EDAR.Rdata`). This allows to run the following scripts without going through the whole data cleaning process.

To analyze this data and produce the plots included in the manuscript, run the scripts `scripts/analysis_BRCA.R` and `scripts/analysis_EDAR.R`; the cleaned output files from these two scripts are saved in `results/` (`results/crossval_BRCA.Rdata` and `results/crossval_EDAR.Rdata`).
The script `scripts/plot_BRCA_EDAR.R` produces the plots reported in the paper.

## Real data: citation data

The citation data used are the ones described and analyzed in 

> Ji, Pengsheng, and Jiashun Jin. "Coauthorship and citation networks for statisticians." The Annals of Applied Statistics 10.4 (2016): 1779-1812.

Part of the dataset was available [here](https://www.stat.uga.edu/directory/people/pengsheng-ji), while the bibtex file was provided by the authors.

Script `scripts/get_data_clean_bibtex.R` imports the bib file, extracts information and creates the order between the paper. Then estracts data related to citations, and saves the cleaned datasets in `data/` (included in `data/bibtex_clean.rdata` and `data/ordered_cit_data.rdata`).

To analyze this data and produce the plots included in the manuscript, run the scripts `scripts/analysis_bibtex.R`.