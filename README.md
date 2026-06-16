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
carvePCM: TODO

## PlugAndPlay
To plug the pcm version of your choice into a whole-cell model, apply the interactive function `plugAndPlay.m`. Examples for its usage can be found in `Code/validation/plugIntoAraCore.m`.

## Structure of this repository
The Repository is stuctured in folders Code, Data, and Figures
- Code contains functions used by multiple scripts, and is structured by the different analyses performed for the manuscript. Code and data for generating enzyme-constrained PCM models is found in the gecko folder
- Data contains several subfolders:
    - pcm contains all models generated in this work
    - comparison_models contain the corrected and annotated plant models used for comparison
    - sequences contains genomic data used for construction of the pcm
    - EC_predictions contains the EC annotations as presented in the manuscript
    - supplementary contains all supplementary data also available with the manuscript
    - analysis contains sub-folders for each of the analyses performed in the manuscript (see `Code` sub-folders)

## How to cite
TBA