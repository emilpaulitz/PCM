import os
import pandas as pd

projPath = '/home/emil/Desktop/PhD-Synch/pan_chloroplast/pan_chl_model/'
geckoProjPath = projPath + 'Code/gecko/panAraEcPcm/'

# get accession from the files in projPath + 'Data/analysis/acc-spec-analysis/04_models/'
accs = [fname.removesuffix('.pcm.v23.mat') 
        for fname in os.listdir(projPath + 'Data/analysis/acc-spec-analysis/04_models/')
        if fname.endswith('.pcm.v23.mat')]

# for calculating protein masses
aa_weights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
       'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
       'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
       'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06 }
def calc_weight(seq):
    total_weight = 0
    printed = False
    for aa in seq:
        if aa in aa_weights:
            total_weight += aa_weights[aa]
        else:
            if not printed:
                print(f'Warning: unknown amino acid {aa} in sequence {seq}')
                printed = True
            total_weight += sum(aa_weights.values()) / len(aa_weights)
     # subtract weight of water for each peptide bond
    return round(total_weight - (18.015 * (len(seq) - 1)))

# read the real arabidopsis uniprot to append to the accession-specific results
ara_uni = pd.read_csv(projPath + 'Code/gecko/ecAraPcm/data/uniprot.tsv', sep='\t', header=0)

for acc in accs:

    # get the corresponding fasta file from projPath + 'Data/analysis/acc-spec-analysis/02_chloroplast_fastas/'
    fasta_file = f'{projPath}Data/analysis/acc-spec-analysis/02_chloroplast_fastas/{acc}.faa'

    # extract gene names and sequences
    with open(fasta_file, 'r') as infile:
        gene_names = []
        sequences = []
        masses = []
        curr_seq = ''
        for line in infile:
            if line.startswith('>'):
                if curr_seq:
                    sequences.append(curr_seq)
                    curr_seq = ''
                gene_names.append(line[1:].strip().split()[0])
            else:
                curr_seq += line.strip()
        if curr_seq:
            sequences.append(curr_seq)

    for name, seq in zip(gene_names, sequences):
        if any(aa not in aa_weights for aa in seq):
            print(f'Warning: sequence for gene {name} contains unknown amino acids')

    # calculate protein weights
    masses = [calc_weight(sequence) for sequence in sequences]

    # write to a pandas df
    
    df = pd.DataFrame({
        'Entry': [f'G{str(i).zfill(5)}' for i in range(len(gene_names))],
        'Gene Names (ordered locus)': gene_names,
        'EC number': ['' for _ in gene_names],
        'Mass': masses,
        'Sequence': sequences
    })

    # append the data from general arabidopsis uniprot
    (pd.concat([df, ara_uni], ignore_index=True)
        .to_csv(f'{geckoProjPath}data/uniprot_{acc}.tsv',
                sep='\t', index=False, header=True))