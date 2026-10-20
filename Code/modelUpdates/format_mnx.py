'''Format chem_xref.tsv and reac_xref.tsv from a user-supplied folder.'''

import argparse
import os
import pandas as pd
import cobra

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Format the MetaNetX cross-reference tables.'
    )
    parser.add_argument(
        'folder',
        help='Folder containing chem_xref.tsv and reac_xref.tsv',
    )
    parser.add_argument(
        'outfolder',
        help='Folder to put the formatted tables',
    )
    parser.add_argument(
        'modelpath',
        help='Path to the cobra model file',
    )
    args = parser.parse_args()

    # check that the required files exist
    out_folder = args.outfolder + ('/' if args.outfolder[-1] != '/' else '')
    input_folder = args.folder + ('/' if args.folder[-1] != '/' else '')
    chem_xref = input_folder + 'chem_xref.tsv'
    chem_prop = input_folder + 'chem_prop.tsv'
    reac_xref = input_folder + 'reac_xref.tsv'
    reac_prop = input_folder + 'reac_prop.tsv'
    for table in (chem_xref, chem_prop, reac_xref, reac_prop):
        if not os.path.isfile(table):
            parser.error(f'Required table not found: {table}')
    if not os.path.isfile(args.modelpath):
        parser.error(f'Required model file not found: {args.modelpath}')

    print('Read the model...')
    model = cobra.io.read_sbml_model(args.modelpath)

    print('Formatting metabolite information...')
    model_mets = gather_model_mets(model)
    format_chem_xref(chem_xref, chem_prop, model_mets).to_csv(out_folder + 'met_translation.csv')

    print('Formatting reaction information...')
    model_rxns = gather_model_rxns(model)
    format_reac_xref(reac_xref, reac_prop, model_rxns).to_csv(out_folder + 'rxn_translation.csv')


# convert chem_xref.tsv to nice csv
def format_chem_xref(chem_xref_path: str, chem_prop_path: str, model_mets: dict) -> pd.DataFrame:

    with open(chem_xref_path, 'r') as file:
        for i, line in enumerate(file):
            if not line.startswith('#'):
                break
    chem = pd.read_csv(chem_xref_path, sep = '\t', skiprows = i, header = None)

    # format and filter table
    chem.columns = ['source', 'ID', 'description']
    mnx_IDs_in_model = chem.loc[chem['source'].isin(model_mets.values()), 'ID'].unique()
    assert len(mnx_IDs_in_model) == len(set(model_mets.values())), 'Not all model metabolites are in chem_xref.tsv'
    chem = chem.loc[(chem['description'] != 'secondary/obsolete/fantasy identifier') &
                    (chem['ID'].isin(mnx_IDs_in_model))].copy()

    # infer the type of identifier and extract only the raw identifier (seed.reaction:R12345 -> seed, R12345)
    tmp_d = {'MNXM1':'metanetx', 'MNXM01': 'metanetx', 'MNXM02': 'metanetx', 'MNXM03':'metanetx', 'BIOMASS': 'metanetx', 'WATER': 'metanetx', 
            'bigg.metabolite':'bigg', 'biggM':'bigg', 'kegg.compound':'kegg', 'keggC':'kegg','kegg.drug':'keggd', 'keggD':'keggd',
            'kegg.glycan':'keggg', 'keggG':'keggg', 'metacyc.compound':'biocyc', 'metacycM':'biocyc',
            'mnx':'metanetx', 'sabiork.compound':'sabiork', 'envipath': 'envipath', 'envipathM': 'envipath',
            'sabiorkM': 'sabiork', 'seed.compound':'seed', 'seedM':'seed', 'CHEBI': 'chebi', 'chebi': 'chebi',
            'reactome': 'reactome', 'reactomeM': 'reactome', 'lipidmaps': 'lipidmaps', 'lipidmapsM': 'lipidmaps',
            'rheaP': 'rhea', 'rheaG': 'rhea', 'slm': 'slm', 'SLM': 'slm', 'hmdb': 'hmdb', 'vmhM': 'vmh', 
            'vmhmetabolite': 'vmh', 'glycosphingo': 'glycosphingo'}
    chem['db'] = chem['source'].apply(lambda x: tmp_d[x.split(':')[0]])
    chem['xrefid'] = chem['source'].apply(lambda x: x.split(':')[1] if ':' in x else x)

    # drop unnecessary data
    chem = chem.loc[chem['db'] != 'metanetx', :]
    chem = chem.drop('source', axis=1).drop_duplicates()

    # namespaces that we can transport into xml: 
    # seed, inchi, chebi, hmdb, reactome, metanetx, bigg, biocyc
    chem = chem.loc[chem['db'].isin(['seed', 'inchi', 'chebi', 'hmdb',
                                    'reactome', 'metanetx', 'bigg', 'biocyc', 'kegg'])]

    # For inchis, we also need chem_prop.tsv
    with open(chem_prop_path, 'r') as infile:
        for i, line in enumerate(infile):
            if not line.startswith('#'):
                break
    chem_prop = pd.read_csv(chem_prop_path, sep='\t', skiprows=i, header = None)
    chem_prop.columns = 'ID     name    reference       formula charge  mass    InChI   InChIKey        SMILES'.split()
    chem_prop.set_index('ID', inplace=True)

    # Filter chem_prop to reduce computation time
    chem_prop = chem_prop.loc[chem_prop.index.intersection(chem['ID'])]

    # Format chem_prop
    chem_prop.loc[chem_prop['formula'] == '*', 'formula'] = float('nan')
    chem_prop = chem_prop[['InChI', 'SMILES']].dropna(axis=0, how='all')
    chem_prop = chem_prop.unstack().reset_index(0).sort_index()
    chem_prop.columns = ['db', 'xrefid']

    # concatenate all information collected and pivot wide
    chem_final = pd.concat([chem.drop('description', axis = 1), chem_prop.reset_index()], join='outer').sort_values('ID')
    return (chem_final.fillna('')
                      .groupby(['ID', 'db'])
                      .agg(lambda x: ';'.join(x))
                      .reset_index()
                      .pivot(index='ID', columns='db', values='xrefid'))

def gather_model_mets(model: cobra.Model) -> dict:
    # These have been mapped to multiple KEGG compounds, which leads to none being
    #   retained in the xml. Selected one for this purpose
    lost_in_translation = {'C00025[h]': 'C00025', 'C00085[h]': 'C00085', 'C00354[h]': 'C00354',
                        'C02051[h]': 'C02051', 'C04272[h]': 'C04272', 'C06010[h]': 'C06010',
                        'C10385[h]': 'C10385', 'C16695[h]': 'C16695', 'C00095[h]': 'C00095',
                        'DHL[h]': 'C15973'}

    # Retrieve KEGG compound IDs if only none is present (other cases are handled above) and add
    model_keggs = {m.id: ('keggC:' + m.annotation['kegg.compound']) for m in model.metabolites if m.annotation and 'kegg.compound' in m.annotation and not m.id in lost_in_translation}
    model_mets = model_keggs | {k: 'keggC:' + v for k, v in lost_in_translation.items()}

    # For the remaining metabolites, we have CHEBI IDs
    for m in model.metabolites:
        if m.annotation and not 'kegg.compound' in m.annotation and not m.id in lost_in_translation:
            if 'chebi' in m.annotation:
                model_mets[m.id] = m.annotation['chebi']

    # For (11R)-carlactonoate, KEGG is not mapped in MNX
    model_mets['C23076[h]'] = 'CHEBI:194504'

    return model_mets

def format_reac_xref(reac_xref_path: str, reac_prop_path: str, model_rxns: dict) -> pd.DataFrame:
    # reac_xref to my_rxn_translation
    # convert rxn lookup-table from MNXref to csv with EC numbers

    # read table and drop unnecessary first row
    with open(reac_xref_path, 'r') as file:
            for i, line in enumerate(file):
                if not line.startswith('#'):
                    break
    rxn_xref = pd.read_csv(reac_xref_path, sep = '\t', skiprows = i, header = None)
    rxn_xref.columns = ['source', 'ID', 'description']
    rxn_xref = rxn_xref.loc[(rxn_xref['ID'] != 'EMPTY') &
                            (rxn_xref['description'] != 'secondary/obsolete/fantasy identifier'), :]

    # infer the type of identifier and extract only the raw identifier (seed.reaction:R12345 -> seed, R12345)
    tmp_d = {'MNXR01':'metanetx', 'MNXR02': 'metanetx', 'MNXR03':'metanetx', 'bigg.reaction':'bigg', 'biggR':'bigg',
                'kegg.reaction':'kegg', 'keggR':'kegg', 'metacyc.reaction':'biocyc', 'metacycR':'biocyc',
                'mnx':'metanetx', 'rhea':'rhea', 'rheaR':'rhea', 'rh':'rhea', 'sabiork.reaction':'sabiork',
                'sabiorkR':'sabiork', 'seed.reaction':'seed', 'seedR':'seed', 'vmhR': 'vmhR', 'vmhreaction': 'vmhR'}
    rxn_xref['db'] = rxn_xref['source'].apply(lambda x: tmp_d[x.split(':')[0]])
    rxn_xref['xrefid'] = rxn_xref['source'].apply(lambda x: x.split(':')[1] if ':' in x else x)

    # drop unnecessary data
    rxn_xref = rxn_xref.loc[rxn_xref['db'] != 'metanetx', :]
    rxn_xref = rxn_xref.drop('source', axis=1).drop_duplicates()

    # there are 10 rxns that could not be found in the data table, and these are indeed not in MNX
    # print(len(set(model_rxns.values()) - set(rxn_xref['xrefid'])))

    # filter for KEGG reactions that are present in the model
    mnx_IDs_in_model = rxn_xref.loc[(rxn_xref['db'] == 'kegg') &
                                    (rxn_xref['xrefid'].isin(model_rxns.values())), 'ID']
    rxn_xref = rxn_xref.loc[rxn_xref['ID'].isin(mnx_IDs_in_model), :]

    # namespaces that we can transport into xml: 
    # reactions: rhea, kegg, seed, metanetx, bigg, reactome, ec-code, brenda, biocyc
    # reactome is not present in the databases mapped in metanetx
    # for ec-code, we need reac_prop
    rxn_xref = rxn_xref.loc[rxn_xref['db'].isin(['rhea', 'kegg', 'seed', 'metanetx', 'bigg', 'biocyc'])]

    # read EC number information
    reac_prop = pd.read_csv(reac_prop_path, sep = '\t', comment='#', header=None)
    reac_prop.columns='ID     mnx_equation    reference       ec        is_balanced     is_transport'.split()
    reac_prop = reac_prop.loc[~reac_prop['ec'].isna(), ['ID', 'ec']]
    reac_prop['db'] = 'EC'
    reac_prop.rename({'ec': 'xrefid'}, axis=1, inplace=True)
    reac_prop = reac_prop.loc[reac_prop['ID'].isin(rxn_xref['ID']), :]

    # merge with reac_xref
    reac_final = pd.concat([rxn_xref.drop('description', axis=1), reac_prop])

    # Pivot wide
    return (reac_final.groupby(['ID', 'db'])
                      .agg(lambda x: ';'.join(x))
                      .reset_index()
                      .pivot(index='ID', columns='db', values='xrefid'))

def gather_model_rxns(model: cobra.Model) -> dict:
    # For reactions, we only have KEGG IDs and EC numbers. From EC-numbers, we cannot infer the reaction
    # manual associations
    man_assoc = {'StP1_h': 'R01821', 'R01821_2': 'R01821', 'R03052': 'R03052', 'R00772': 'R00772'}
    return {r.id: r.annotation['kegg.reaction'] for r in model.reactions 
                if r.annotation and 'kegg.reaction' in r.annotation} | man_assoc
    
if __name__ == '__main__':
	main()

