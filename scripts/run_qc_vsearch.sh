wosho=/home/qiime/Desktop/workshop
out=$wosho/output/QIIME
out2=$wosho/output/VSEARCH
raw=$wosho/raw_data

for sample in $(ls $raw/*gz | awk -F '/' '{print $NF}' | cut -d '_' -f 1 | uniq | grep "^B"); do
echo "Running sample $sample"
samplenew=$(echo $sample | tr '-' '.')
mkdir -p $out2/merged/${samplenew}
$wosho/software/vsearch --fastq_mergepairs $out/trimmed/${samplenew}_R1.fq --reverse $out/trimmed/${samplenew}_R2.fq --fastq_maxee 2 --fastq_minmergelen 220 --fastq_maxmergelen 280 --fastaout $out2/merged/${samplenew}/${samplenew}.fna
mkdir -p $out2/nonchimeras/${samplenew}
$wosho/software/vsearch --derep_fulllength $out2/merged/${samplenew}/${samplenew}.fna --sizeout --output /dev/stdout | vsearch --uchime_denovo /dev/stdin --nonchimeras /dev/stdout | vsearch --rereplicate /dev/stdin --output $out2/nonchimeras/${samplenew}/${samplenew}_clean.fna --relabel ${samplenew}_
done
echo "Done"
