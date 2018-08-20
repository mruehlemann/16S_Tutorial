print 'Please wait while the files are being downloaded. This will take a few minutes'

cd raw_data
for file in $(cut -f 12 | grep 'fastq/B' | tr ';' '\n'); do wget $file; done
cd ..
