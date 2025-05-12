# Running a spatialproteomics pipeline with snakemake
This repository shows how to go from raw CODEX TMA data to analysis-ready data using snakemake.
The `notebooks` folder contains explanations of how to perform the individual steps.

# Step 1: Cropping
Cropping of the TMAs requires some manual oversight. The notebook `notebooks/01_cropping_tmas_in_qupath.ipynb` details how this can be done using QuPath.

# Step 2: Running the pipeline
The entire snakemake pipeline can be found in the `pipeline` folder. Please refer to the notebook `notebooks/02_snakemake.ipynb` for details on how to run this pipeline.
