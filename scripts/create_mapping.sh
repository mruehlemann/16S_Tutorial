wosho=/home/qiime/Desktop/workshop
otu_folder=$wosho/output/otu_picking/QIIME_dnVS

for entry in $(grep '>' ${otu_folder}/vsearch_otus/dnVS.fa  | tr -d '>' | tr ' ' ';'); do
otu=$(echo $entry | cut -d ';' -f 1)
rep=$(echo $entry | cut -d ';' -f 2)
cc=$(grep -w $rep ${otu_folder}/vsearch_otus/dnVS.uc | grep "^H" | cut -f 9 | tr '\n' '\t')
echo "$otu\t$rep\t$cc" >> ${otu_folder}/vsearch_otus/All.QIIME_otus.txt
done
cut -d ' ' -f 1 ${otu_folder}/vsearch_otus/dnVS.fa > ${otu_folder}/vsearch_otus/All.QIIME_rep_set.fasta

# otu_folder2=$wosho/output/otu_picking/VSEARCH_dnVS
# for entry in $(grep '>' ${otu_folder2}/vsearch_otus/dnVS.fa  | tr -d '>' | tr ' ' ';'); do
# otu=$(echo $entry | cut -d ';' -f 1)
# rep=$(echo $entry | cut -d ';' -f 2)
# cc=$(grep -w $rep ${otu_folder2}/vsearch_otus/dnVS.uc | grep "^H" | cut -f 9 | tr '\n' '\t')
# echo "$otu\t$rep\t$cc" >> ${otu_folder2}/vsearch_otus/All.VSEARCH_otus.txt
# done
# cut -d ' ' -f 1  ${otu_folder2}/vsearch_otus/dnVS.fa > ${otu_folder2}/vsearch_otus/All.VSEARCH_rep_set.fasta
