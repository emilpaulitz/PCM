Code for calculating the values presented in the first Results section, and generating the plots therein. 
- `circularplots.py` contains code for generating the circular pan.reactome plot
- `panReactome.m` reads all species-specific models and extracts a table of reaction presence/absence
- `plotsAndAnalysis.ipynb` is a Python notebook that contains code for figure generation and permutation analysis
- `speciesSpecModels.m` contains the code that was used for generating species-specific versions of the PCM

Because several incompatibility issues plague current BUSCO versions, it was run using a patched version available under https://github.com/emilpaulitz/busco-apptainer (commit ccd3dfff2ba81a8b610baffbed7812ae76cfa03e).