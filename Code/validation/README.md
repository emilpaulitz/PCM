In this folder, find the code used to perform the KO study and comparison with experimental data, as well as the comparison of enzyme abundance prediction with measured data.

### KO study
- Formatting of GPR rules along with removal of non-chloroplastic genes is performed in orgSpecificChanges.m, writing v2 of the respective models
- The respective chloroplast models were plugged into conventional Aracore using plugIntoAracore.m
- The analysis and plotting of results is done with cobraPy in KoCompare.ipynb

### Enzyme abundance
- For this, we used the enzyme-constrained PCM version of Arabidopsis (´Data/pcm/ecPCM/Ath.tuned.xml´) that was generated in ´Code/gecko/ecAraPcm/makePcGemFinal.m´. The first time you run this, it will take a while to download KEGG data, which was too large to put on GitHub
- The enzyme-constrained pcm was plugged into an enzyme-constrained AraCore (From https://github.com/pwendering/AraTModel/tree/master/metabolic-models and updated with KEGG IDs) using the code in ´Code/validation/ecPcmIntoAraCore.m´
- The analysis is done in protDistr.ipynb