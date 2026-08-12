# PCM
Code and data relating to the pan-chloroplast model (PCM). 

## Getting started
- Clone the repository with `git clone git@github.com:emilpaulitz/PCM.git`
- Open the startup.m script in matlab
- Fill in the path to your cobratoolbox installation (see [Cobra toolbox Documentation](https://opencobra.github.io/cobratoolbox/stable/installation.html))
- If you have access to gurobi, indicate the path to the matlab folder of gurobi
- Save and run the startup.m script
- For running flux coupling analysis, please download [F2C2](https://doi.org/10.1186/1471-2105-13-57) into the Resources folder

## Generating a genotype-specific chloroplast model
- To generate a new genotype-specific chloroplast model, you need a set of chloroplast proteins as a `.fasta` file.
- Run BLAST against the genes in the union model. Example commands and the required output specifications can be found in `Code/accSpecPCM/blast.sh`. The resulting files should be named according to the pattern `${ID}.out.tsv`.
- Then, open the script `Code/carveNewPCM.m`, adapt the parameters in the ADAPT THESE section, and run. 

## PlugAndPlay
To plug the pcm version of your choice into a whole-cell model, apply the interactive function `plugAndPlay.m`. Examples for its usage can be found in `Code/validation/plugIntoAraCore.m`.

The script first deletes the whole-cell model's chloroplast and then renames metabolites that the models share across all compartments of the PCM. Keep that in mind when it comes to difficult decisions. So, say you integrate the PCM into a model from the MNXM namespace and the script asks you:
```
Rxn Tr_OPDA_c1 (974 / 1612) already exists. Please choose which version to keep (0 is from cyt, 1 from chl model):
(0) Tr_OPDA_c1 : 12-OPDA[hm] <==> MNXM729361[cy]
(1) Tr_OPDA_c1 : MNXM729361[hm] <==> MNXM729361[cy]
(2) Keep both, appending _chl to chl rxn ID
Index of reaction to keep (0|1|2):
```
The reaction (1) from the PCM uses the metabolite that was re-named using the name of the cytosol metabolite. However, given the old metabolite in the whole-cell model's chloroplast inner membrane (`12-OPDA[hm]`) had a different name, the script has to ask the user. Here, you should probably type 1, because `12-OPDA[hm]` is likely not connected anymore because the whole-cell model's chloroplast was removed. However, it can never hurt to choose 2 and check this manually after the procedure ran through.

## Structure of this repository
The Repository is stuctured in folders Code, Data, and Figures
- Code contains folders for the different analyses performed for the manuscript. Code and data for generating enzyme-constrained PCM models is found in the gecko folder. Further, Code contains functions used by multiple scripts, as well as: 
    - `plugAndPlay.m` for plugging a PCM model into a whole cell model
    - `carveNewPCM.m` for generating species-specific versions of the PCM
- Data contains several subfolders:
    - pcm contains all models generated in this work
    - comparison_models contain the corrected and annotated plant models used for comparison
    - sequences contains genomic data used for construction of the pcm
    - EC_predictions contains the EC annotations as presented in the manuscript
    - supplementary contains all supplementary data also available with the manuscript
    - analysis contains sub-folders for each of the analyses performed in the manuscript (see `Code` sub-folders)

## Feedback
If you spot a bug, problem, or possible improvement in the code or model, do not hesitate to get in contact, open an issue, or a pull request. Any feebdack is appreciated!

For improvements to the model, either suggest the changes, or: add a script that shows the changes made, and upload the new model with the version number updated. 

## How to cite
TBA

[![DOI](https://zenodo.org/badge/1246445440.svg)](https://doi.org/10.5281/zenodo.21031415)
