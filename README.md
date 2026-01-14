# Sparse DIB Simulations

This repository includes the code related to "An information-theoretic algorithm for sparse clustering" by Costa, Papatsouma, and Markos (2026).

- `data`: Includes the cleaned version of bladder cancer data set.
- `data_res`: Includes results on the bladder cancer data set using Sparse DIB with varying sparsity constraint values, COSA/PAM, Sparse K-Means, VarSelLCM, PCA & Sparse PCA followed by K-Means and RPEClust.
- `res`: Results of simulations on synthetic data. The `.zip` files include the results for the 7 sparse clustering algorithms implemented. The `.RDS` files include summaries of the results; `results_df` for overall ARI/AMI performance, `results_best_vars` for the variable selection process of Sparse DIB & Sparse K-Means, `results_entropies` and `sparsity_param_df` include results relevant to the sparsity constraint value of Sparse DIB & Sparse K-Means.
- `src`: Directory of source files including R functions for constructing weighted perturbed similarity matrices (`coord_to_pxy_weights.R`), the Dykstra projection algorith (`dykstra.R`), estimating the mutual information between each feature and a clustering solution (`mi_comps.R`) and the main Sparse DIB implementation function (`sparse_DIB.R`).
- `gen_sparse_data.R` generates the synthetic data sets for the simulation study; these have not been added to the repository as they take up a lot of space. Please run this function first if you want to reproduce the synthetic simulation results.
- `_script.R` files (e.g. `RPEClust_script.R`) implement the clustering algorithms on the synthetic data sets.
- `_bladder.R` files implement the clustering algorithms on the cleaned bladder cancer data.
- `omics_preprocess.R` includes the pipeline for accessing the TCGA-BLCA data and the data cleaning process that has been used.
- `res_realdata.R` performs a results analysis in terms of number of features selected and ARI scores on the bladder cancer data.
- `sdib_bladder_meta.R` includes code for extracting the genes selected by Sparse DIB on the bladder cancer data and for reproducing Fig. 3 in the paper.
