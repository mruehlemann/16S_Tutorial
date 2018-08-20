wosho=/home/qiime/Desktop/workshop
otu_folder=$wosho/output/otu_picking/VSEARCH_dnVS

# Assign taxonomy command
assign_taxonomy.py -o ${otu_folder}/uclust_assigned_taxonomy -i ${otu_folder}/vsearch_otus/All.VSEARCH_rep_set.fasta

# Make OTU table command
make_otu_table.py -i ${otu_folder}/vsearch_otus/All.VSEARCH_otus.txt -t ${otu_folder}/uclust_assigned_taxonomy/All.VSEARCH_rep_set_tax_assignments.txt -o ${otu_folder}/otu_table.biom

# Align sequences command
align_seqs.py -i ${otu_folder}/vsearch_otus/All.VSEARCH_rep_set.fasta -o ${otu_folder}/pynast_aligned_seqs

# Filter alignment command
filter_alignment.py -o ${otu_folder}/pynast_aligned_seqs -i ${otu_folder}/pynast_aligned_seqs/All.VSEARCH_rep_set_aligned.fasta

# Build phylogenetic tree command
make_phylogeny.py -i ${otu_folder}/pynast_aligned_seqs/All.VSEARCH_rep_set_aligned_pfiltered.fasta -o ${otu_folder}/rep_set.tre
