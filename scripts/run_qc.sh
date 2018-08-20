wosho=/home/qiime/Desktop/workshop
out=$wosho/output/QIIME
raw=$wosho/raw_data

for sample in $(ls $raw/*gz | awk -F '/' '{print $NF}' | cut -d '_' -f 1 | uniq | grep "^B"); do
echo "Running sample $sample"
samplenew=$(echo $sample | tr '-' '.')
$wosho/software/vsearch --fastq_filter $raw/${sample}_*_R1_001.fastq.gz --fastq_trunclen 220 --fastqout $out/trimmed/${samplenew}_R1.fq
$wosho/software/vsearch --fastq_filter $raw/${sample}_*_R2_001.fastq.gz --fastq_trunclen 220 --fastqout $out/trimmed/${samplenew}_R2.fq
join_paired_ends.py --min_overlap=150 -f $out/trimmed/${samplenew}_R1.fq -r $out/trimmed/${samplenew}_R2.fq -o $out/join_pairs_output/${samplenew}
split_libraries_fastq.py -i $out/join_pairs_output/${samplenew}/fastqjoin.join.fastq -o $out/slout_q20/${samplenew} -q 19 --barcode_type 'not-barcoded' --sample_ids ${samplenew}
identify_chimeric_seqs.py -m usearch61 -i $out/slout_q20/${samplenew}/seqs.fna -o $out/usearch61_chimera/${samplenew} -r $wosho/resources/gold.fa
filter_fasta.py -f $out/slout_q20/${samplenew}/seqs.fna -o $out/usearch61_chimera/${samplenew}/seqs_clean.fna -s $out/usearch61_chimera/${samplenew}/chimeras.txt -n
done
echo "Done"
