# remember to activate the appropriate conda env with DIAMOND installed
proj_dir=${PWD%Code/accSpecPCM}
seqs_dir="${proj_dir}Data/analysis/accSpecPCM/02_chloroplast_fastas/"

diamond makedb --in "${proj_dir}Data/GenesInModel.fasta" --db "${proj_dir}Data/analysis/accSpecPCM/03_BLAST_files/union"
for fasta in $(ls $seqs_dir)
do
    # with sequences from multiple organisms in the model, we needed to increase the max-target-seqs parameter
    diamond blastp --query $seqs_dir$fasta \
        --db "${proj_dir}Data/analysis/accSpecPCM/03_BLAST_files/union.dmnd" \
        --threads 8 --evalue 1e-3 --id 80 --max-target-seqs 10 \
        --out "${proj_dir}Data/analysis/accSpecPCM/03_BLAST_files/${fasta%.faa}.out.tsv" \
        --outfmt 6 qseqid sseqid pident evalue bitscore score qlen slen length
done
