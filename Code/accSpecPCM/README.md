This folder contains code for the generation and analysis of accession-specific PCM versions.

- The pipeline starts in `araAccs.ipynb`. The fasta files and the localization prediction output data are used to extract fasta files of chloroplast proteins.
- DIAMOND is run against the genes in the union model using `blast.sh`
- This input is used by the `araAccsModelCreation.m` script to generate conventional accession-specific models
- Analysis of these is resumed in `araAccs.ipynb`
- The models are made enzyme-constrained in `Code/gecko/panAraEcPcm/`:
    - The adapter files for the different accessions are created by `generate_adapters.py` (run from within the PCM project)
    - The missing uniprot data was replaced by data calculated from the protein sequences, with `generate_uniprot.py`
    - The main function is in `code/makeAccSpecEcPcm.m`. Kcat tuning seems to vary between runs, and the resulting models might deviate slightly
- The kcats of these mdoels are recorded by the script `araAccsRecordKcats.m`
- The analysis, including pFBA predictions, is carried out in `araAccsGEM.ipynb`