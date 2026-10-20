Code in this folder was used to curate the PCM from version 1 to version 2, unblocking all reactions.
- The main script is `updateToV2.m`
- `addKEGGReactionNew.m` is a function to automatically retrieve a reaction from KEGG with all its relevant data and metabolites' data, and integrate it into the PCM. 
- `addCytosolMet.m` adds the cytosol-version of a metabolite that is already present in the chloroplast of the PCM
- `automaticBalancing` attempts to automatically balance reactions by adding a small amount of water and proton metabolites, if that fully balances the reaction
- `get_unbalanced.m` gathers unbalanced reactions from the model and shows in which regard these are unbalanced
- `format_mnx.py` formats data tables downloaded from MetaNetX and subsets them to the contents of the PCM. The tables were written into `Data/*_translation.csv` and used to update Miriams in the PCM.

This folder is not automatically added to the path in `startup.m`. If you want to run `updateToV2.m` or use the scripts in this folder, add it to the path manually. 