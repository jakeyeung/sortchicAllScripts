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
