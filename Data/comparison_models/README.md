This folder contains updated sbml files of plant and chloroplast models against which we conducted our comparison. The models were downloaded from their original paper, except for the following changes:

- AraGEM was downloaded from the following link, because the model in the supplementary does not conform to SBML: https://www.ebi.ac.uk/biomodels/MODEL1507180028#Files
- Mintz-Oron Arabidopsis leaf model: 10 GPR rules that had an incorrect format were corected manually to be able to read the model.
- AraCore was downloaded from https://github.com/pwendering/ArabidopsisCoreModel and anotated with KEGG identifiers by combining multiple sources in the Nikoloski lab, from: Ashwin Ananthanarayanan, Sebastian Huß, Philipp Wendering (Thank you!)
- panAlgae by Røkke et al 2020: KEGG IDs were extracted from the rxnNotes field to make them available in the .xml 
- CAM by Shameer et al 2018: The file could be read by the matlab COBRA toolbox, and was exported again to make it also readable by CobraPy. 
- Poplar: was read once with CobraPy and then exported, because the exported version could be read much faster by cobraPy.
- Soybean: The compartment information was not read properly. The model was read with cobraPy and a list of metabolites was generated to extract the compartment abbreviations from metabolite IDs into the comps field. Names for the compartments were extracted from the paper and stored in the compNames field.
- Maize: The Full-Maize-Model.sbml contained duplicate metabolite IDs, but could be read by the matlab COBRA toolbox. After exporting again, metabolite IDs have duplicate compartment tags but the model can be read by cobraPy.