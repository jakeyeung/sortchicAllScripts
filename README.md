* Code for analysis of scChiC

** Structure of directory

- data/: Contains small files that might be relevant for the project like barcode IDs and other metadata

- notebooks/: Where I put .Rmd files to generate slides or documents.

- scripts/: Where I put all my analysis scripts
    - scripts/processing: All scripts that are run on the cluster (big jobs that don't require interactive or plotting)
    - scripts/Rfunctions: Rfunctions often used in scripts_analysis, but also can be used in processing
    - scripts/transfer_scripts: scripts to remember which files I download on my computer to explore
    - scripts/scripts_analysis: All analysis run on Rstudio, interactive plots. Paths here often are local on computer. Maybe we can fix this by getting an Rstudio server, but need support from IT department.
    - scripts/scripts_analysis/primetime: After exploring some datasets, I put analyses that generates "pretty figures" here.
    - scripts/scripts_analysis/primetime_from_server: Primetime scripts that can run from the server 

** How to run LDA

Pipelines are in processing. Usually a pipeline is put inside a subdirectory that can be run "as is". Example, the pipeline that runs LDA pipeline on 100 kb binned matrix is in "run_LDA_pipeline_on_binned_matrix"

Inside there are two numbered files. Sequence tells you the order you run scripts (more or less): 1-filter_count_mat.R, that creates a count matrix. Then 2-run.run_LDA_model.hiddenDomains.sh runs the LDA model on the count matrix output. 

The script "2-run.run_LDA_model.hiddenDomains.sh" is named to imply it's running run_LDA_model.R with the proper input parameters.

lib/run_LDA_model.R is the script that takes input arguments and runs the LDA model. 

To use run_LDA_model.R all you need is a count matrix (in normal or sparse format) put into an object that can be accessed by "count.dat$counts" (I use the $counts here because Rsubread package creates this object, but in theory it can be made without Rsubread). 
