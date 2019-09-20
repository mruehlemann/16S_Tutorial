print 'Please wait while the files are being downloaded. This will take a few minutes'

cd raw_data
for pair in $(cut -f 10,12 PRJEB12201.txt | tr '\t' ';' | grep 'fastq/B' ); do
files=($(echo $pair | tr ';' ' '))
f3name=$(echo ${files[2]} | awk -F'/' '{print $NF}')
f4name=$(echo ${files[3]} | awk -F'/' '{print $NF}')
wget -O $f3name ${files[0]}
wget -O $f4name ${files[1]}
done

cd ..
