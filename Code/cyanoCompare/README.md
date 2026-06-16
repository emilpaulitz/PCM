This folder contains the Code to generate the result regarding the comparison to cyanobacterial metabolism. 

- `modelCuration.m` was kindly provided by Mauricio Alexander de Moura Ferreira to showcase how the model `iSynCJ816_curated.xml` was generated 
- `cyanoGenerate.m` curates the cyanobacterial model and brings it to a format so it can be plugged into a parent model. It also annotates it with KEGG IDs, which were extracted from KEGG by searching for the metabolite names. Then, the model is plugged into AraCore and post-processed
- The models are further processed in in `cyanoCompare.ipynb`, including the removal of transporters between the chloroplast and the cytosol that is not present in both models. All plots are generated in this notebook
- Flux coupling analysis is carried out in `cyanoCouplingAnalysis.m`, using F2C2 (Larhlimi et al 2012)