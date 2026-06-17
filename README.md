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
To generate a new genotype-specific chloroplast model, you need a set of chloroplast proteins as a `.fasta` file. Run BLAST against the genes in the union model. Example commands and the required output specifications can be found in `Code/accSpecPCM/blast.sh`. The resulting files should be named as the pattern `${ID}.out.tsv`. Then, open the script `Code/carveNewPCM.m`, adapt the parameters in the ADAPT THESE section, and run. 

## PlugAndPlay
To plug the pcm version of your choice into a whole-cell model, apply the interactive function `plugAndPlay.m`. Examples for its usage can be found in `Code/validation/plugIntoAraCore.m`.

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
If you spot a bug, problem, or possible improvement in the code or model, do not hesitate to get in , open an issue, or a pull request. Any feebdack is appreciated!

For improvements to the model, either suggest the changes, or: add a script that shows the changes made, and upload the new model with the version number updated. 

## How to cite
TBA
