# 16S Amplicon Sequencing Workshop Tutorial #

Written by Malte Rühlemann

Hosted by: Malte Rühlemann & Felix Sommer

## 1. Introduction ##

This workshop is intended to give a broad overview of the processing and analysis of sequencing data from 16S rRNA gene amplicon-based experiments aiming to survey bacterial communities. We will start with command line basics and handling of sequencing data in FastQ format, to introduce what to expect when performing such experiments. After this we will move through the different data processing steps that are necessary when dealing with amplicon data. Finally, we will have a look at a final dataset and compare the outcomes of different processing possibilities, using QIIME output and also the statistical software R. The latter we will also use to discuss properties of microbiome data that need to be taken into account when applying statistical test. 

Please feel free to explore any output file in between as much as you like. We will also put this tutorial online so you can redo the analysis again and as a general resource for your own analyses. What we provide here is one or a few of the possible ways to process and analyze microbiome data, however there is not a single “right” way to do it and looking through the literature, you will easily find people doing something different from what you have been introduced to today. That is why we are being extensive in the explanations in this tutorial, so you can really understand what you are doing and why you are doing it. 

As we will be using the bash command line a lot, you can go here for a general resource about the commands used: https://ss64.com/bash/

The VirtualBox on the provided computers should already be configured and working. Please start it up by clicking the funny looking blue icon in the taskbar on the left.

The dataset we will use in the tutorial is a subset of the data from the article ‚The Gut Microbiota Modulates Energy Metabolism in the Hibernating Brown Bear Ursus arctos’ by Sommer et al. (Cell Reports, 2016, doi.org/10.1016/j.celrep.2016.01.026).

About the dataset: this dataset consists of a total of 22 samples generated from brown bear fecal samples. The region used for amplification and sequencing is the V4 region of the 16S rRNA gene. The primers used were 515F and 806R, creating an amplicon of  around 250bp. Sequencing was performed in 2x250 bp paired-end mode on an Illumina MiSeq machine. This means, the amplicon is basically completely covered twice by each read-pair, which is important to assure high quality 16S amplicon sequences. 

## 2. The FastQ format – working with the command line and Illumina sequencing data ##

After sending your samples out for sequencing, you will very likely not receive the very raw data from the sequencer, but a bunch of files for each sample you handed over for sequencing. These files carry the raw, already base-called and demultiplexed, sequences in the so-called FastQ format. The following steps will introduce you to this format and explain how to proceed with this kind of data. But first, let’s get a little familiar with the command line.

Open a Terminal by clicking on the icon in the sidebar on the left. The terminal should open. You will see something like:
```bash
qiime@qiime-190-virtual-box:~$ 
```
and a ‚flashing’ box. This is where you can enter your commands. Every line with a command will now be written starting simply with a dollar sign $. For example:
```bash
$ pwd
```

Please do not include the $ when copying the command. When you type this command, the terminal will put out the folder you are currently working in (pwd is short for ‚print work directory’). You can learn more about commands if you enter man followed by the command:
```bash
$ man pwd
```
This shows you information about the pwd command. You can close the manual by pressing q. You can do this for many of the default commands if you want to have further information.

We first need to download the data needed for the tutorial. Type

```bash
$ cd Desktop
$ git clone http://github.com/mruehlemann/16S_Tutorial
$ sh scripts/get_data.sh
```

Let’s switch to the folder where we can access our sequencing data (cd; change directory) and show the contents (ls; list directory content):
```bash
$ cd /home/qiime/Desktop
$ ls
```
We can see that in this folder there is another folder called raw_data. We now change into this directory, and list the first 10 files (head –n 10):
```bash
$ cd raw_data
$ ls | head -n 10
```
Between the commands ls and head, we see a vertical bar. This bar is called the pipe and is used to pass the output of one command (also called stdout) to another command, in this case the head command. The head command shows only a specified number of lines of the input. The number of lines to be displayed can be adjusted using the -n parameter.

What you see as output is the typical data you will get from your sequencing facility. Let us dissect the filename of the first file to make sense of it:

<p align="center"><tt>B-1104-S_S70_L001_R1_001.fastq.gz</tt></p>
<table>
<tr>
<td><tt>B-1104-S<td>This is the sample name and can be anything you want it to be. In this case this name is composed of ‘B’ for ‘Bear’, ‘1104’ is the ID of the bear, and ‘S’ which means that the sample was taken in summer. You will also see “B-1104-W” a few lines further down, which is the data of the same individual but taken in winter.
<tr>
<td><tt>S70<td>This part of the filename simply tells you that it is sample 70 of this specific sequencing run. Among other things to expect at this position of the filename are for example the index-sequences used in multiplexing.

<tr>
<td><tt>L001<td>This stands for ‘Lane 1’. When working with MiSeq data, this will always be Lane 1, as a MiSeq sequencer only has a single lane. When working with data from e.g. a HiSeq sequencer, this might be something between L001 and L008.

<tr>
<td><tt>R1<td>This tells you that this file contains the forward read (Read 1) of the sequencing data. When you look at the second file in the list, you will see that the filename is exactly the same apart from this position, where it says ‘R2’. This means that these two files belong together and contain the read pairs of a single sample.

<tr>
<td><tt>001<td>This is always 001. Don’t ask why.

<tr>
<td><tt>fastq.gz<td>This is the filename ending, indicating the file format. fastq indicates that the file is in FastQ format, the gz at the end means, that the file is compressed using the Gzip compression. More about this later.
</table>

Okay, now that we know what all those filename mean, let’s have a look at what is in the file. For this we need to first decompress the file (gunzip), pass the output to stdout (-c), pipe it to the next command and display the first four lines in the file:

```bash
$ gunzip -c B-1104-S_S70_L001_R1_001.fastq.gz | head -n 4
```

What you get will look something like this:
```
@M02765:9:000000000-A8H6D:1:1101:15136:1407 1:N:0:70
TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGG
+
ABBBBFFBFFFFGGGGFGGEGFHGGEFFHGGHHGHHGGGGHHHGGGEGHHHH
```
This output is the standard format of the FastQ sequencing data format. One sequence in FastQ format always consists of four lines:
<table>
<tr>
<td>Line 1:<td>This line is the start of the entry and always starts with the @ symbol, followed by a unique sequence ID. This ID is the same in both read files (R1+R2) and defines which sequences belong together. When looking at Illumina sequencing data, this ID consists of the unique instrument ID (M02765), the run ID (9), the flowcell ID (A8H6D) and the lane (1), followed by coordinates of the sequence cluster on the flowcell. After the space there are more informations on the sequence, however these do not belong to the unique sequence ID.

<tr>
<td>Line 2:<td>This is the actual sequence obtained from basecalling of the raw data. In our case this is the 16S amplicon sequence.

<tr>
<td>Line 3:<td>this line can contain additional information about the sequence, however usually only contains a mandatory ‚+’ symbol. 

<tr>
<td widht=3>Line 4:<td>In this line the quality of every single base in the sequence is encoded. The first position of this line corresponds to the first position of the sequence, the second to the second, and so on. Each ASCII-character endcodes a special value between 0 and 41 (in the most widely used Illumina 1.8 Phred+33 format; +33 stands for an offset of 33 in the ASCII characters).  The highest Phred score a base can get is in this case 41, encoded by the letter J (ASCII: 74), the lowest is 0, encoded by ! (ASCII: 33). The Phred score (Q) is defined as the probability (P), that the base at this position is incorrect: P(Q)=10-Q10. This corresponds to 0.01% error probability at Q=40 and 10% at Q=10.
</table>

Phred+33 encoding of quality scores:
```
Symbol	 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
        |                                        |          
ASCII	  33                                      73
Q Score	0.2......................26...31........41                      
```


Of course it is not possible to check each single sequence by ‘hand’, however there are a lot of very handy tools to check if the sequencing in general yielded satisfactory data. One tool we want to introduce is the very easy to use and freely available FastQC. Please go the website (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) using the browser in the VirtualBox, download the Linux version of the tool and unpack it. For this, please create a folder in the "16S_Tutorial" folder called "software", double-click on the downloaded archive and move the ‘FastQC’ folder into the newly created ‘software’ folder. Go to the terminal again and change to the FastQC software folder:
```
$ cd /home/qiime/Desktop/16S_Tutorial/software/FastQC
```

To use fastqc, we first have to make the binary file executable. For this type this into the terminal:

```
$ chmod +x fastqc
```

The command chmod can be used to change the mode or access permission of a file. Here we add permission to execute (+x) to the fastqc binary. No we can start FastQC:

```
$ ./fastqc
```

This should open the FastQC software. In the Menu bar go to: “File > Open …” and navigate to the raw_data folder in Shared_Folder. By holding the Shift key you can select multiple files at once. Please select the R1 and R2 read of the Sample “B-1105-W” and click the “OK” button. You will see the files being processed and quickly the reports should open, each as a tab in the FastQC window. On the left side you will see several green, orange and red symbols. These are the different quality criteria assessed by FastQC. However, let’s start with the first one called “Basic statistic”. Here you can already see some basic properties of the opened files, for example that each of the files holds 77868 sequences (as they are the paired reads of the same sample this should be the same; if it is not, there is something wrong), the sequence lengths are between 248 and 251 nt and the GC content differs sligthly between R1 and R2, however both being a little above 50%.

Next please have a look at the “Per base sequence quality” tab for both files. What you see here is a pretty typical (though not perfect) example of how amplicon sequencing data normally looks like:

the R1 is relatively good quality all along, only dropping a little at the end of the read, however mostly staying in the green area until ~180-200 nt.
the R2 is usually of worse quality than the R1. This is  is basically what you can expect and mainly due to technical reasons, so don’t worry.  Only if you see really strange patterns of very drastic drops in quality (also maybe only at a single position across all samples) things might have to be investigated more deeply.

Move on to the tab called “Per base sequence content”. This tab will likely be marked with a red cross, indicating a problem. When you look at the image in this tab, you can see the per-position composition of the four nucleotides and notice, that there is only very little variation in the first 8-9 nt. This is just as it is supposed to be, as the amplification primers have to bind in a non- or at least less-variable region of the 16S gene to work. As FastQC is intended to be used for any kind of sequencing data, it assumes that this is an error, however in our case, again, this is what we are expecting.

The same rules apply for GC content (of course there is a bias in GC content when only looking at a specific sequence) and duplicated/overrepresented sequences: we expect to have many duplicated or similar sequences, as this is exactly the purpose of using an amplicon.

In general we can say: our data look just like we expect it to look!

One more very convenient thing about FastQC is, that you can also use it for batch processing of FastQ file. Close the FastQC window and type:

```bash
$ wosho=/home/qiime/Desktop/16S_Tutorial
$ echo $wosho
$ raw=$wosho/raw_data
$ echo $raw
$ ./fastqc $raw/B-1105-W_* --outdir /home/qiime/Desktop
```

In the first line, we set a variable named “wosho”. This is very easy in the bash environment and might look a little different in other programming languages. In general, many things can be stored in variables, however, here we use it for a path to a folder. In the third line, we set another variable “raw”, which has the same content as “wosho”, plus “/raw_data” at the end. As you can see, the content of the variables can be show by using “echo” followed by the variable with a dollar sign prepended. The last line starts the fastqc binary, applied to all files in the folder that is stored in the “raw” variable that start with “B-1105-W_”. The asterisk “*” is a so-called “wildcard” and matches every character when searching for filenames.

Now look at your desktop. There should be two html files you can open with your browser. These file hold all the information we saw before in the FastQC window and this way of processing can be easily applied to a whole batch of sequencing files for fast processing to get an idea if your data looks usable.

Another software that is very nice for getting an overview of sequencing data in general is MultiQC (http://multiqc.info), which is able to aggregate the output of different QC tools to create reports. However, this is currently more aiming at RNAseq or RRBS experiments, but worth keeping an eye on. 

## 3. Step-by-step processing of amplicon sequencing data ##

Now that we have an idea how our sequencing data looks like and that it seems to be alright to work with, we can start to process the data. There are several popular tools for processing amplicon sequencing data, among them mothur (https://mothur.org), QIIME (http://qiime.org) and USEARCH (https://drive5.com/usearch/), to only name a few of them. Most pipelines in use are based on one of these three frameworks, however partly they are also overlapping (e.g. USEARCH optionally being used in some parts of QIIME). While all of them have their advantages and drawbacks, none of them is particularly better or worse than the other, it more depends on what one is looking for (though it is likely that once you are accustomed to one, you will develop strong feelings against any other of them). 

As we are using the QIIME VirtualBox there is little to be installed additionally. One software we will need is the usearch binaries, which are external dependencies not included in the QIIME VirtualBox due to licensing reasons (more info on that here). We need to add a folder where we can put binaries that should be accessible by the system:

```bash
$ echo "export PATH=$PATH:$wosho/software" >> $HOME/.bashrc
$ source $HOME/.bashrc
```
These commands add the folder ``$wosho/software`` to the PATH variable. The PATH variable tells the system, where it should look for executable binaries that can be called directly in the terminal. 

Please now register for your personal copies of USEARCH v5.2.32 and USEARCH 6.1.544 on https://www.drive5.com/usearch/download.html and save downloaded binaries in the software folder as ``usearch`` and ``usearch61``, respectively.

In addition to the usearch and usearch61 binaries already provided in the  folder, we need VSEARCH (https://github.com/torognes/vsearch). VSEARCH is a free, open-source, 64-bit alternative to USEARCH. Download the latest version from here (https://github.com/torognes/vsearch/releases) and unpack it to the “Downloads” folder. Then run the following command to get a shortcut and can then run the VSEARCH software:

```bash
$ vsearch
$ cd $HOME/Downloads/vsearch-2.8.1-linux-x86_64/ 
$ cp bin/vsearch $wosho/software
$ vsearch
```

The first command calling “vsearch” will fail. Only after copying (this is the cp command) the binary to the folder we just added to the PATH variable it is possible to run “vsearch” directly from the command line.

### 3.1 QIIME - the ‘out-of-the-box’ solution ###

QIIME (short for Quantitative Insights Into Microbial Ecology; pronounced ‘chime’) is a more or less complete suite for microbiome analysis. Earlier this year QIIME2 was released, however we will only look at the ‘final’ Version of QIIME1 (v.1.9) for now. In QIIME you can find lots of handy scripts and an active and helpful community (https://groups.google.com/forum/#!forum/qiime-forum). Downsides are however that it is a bit like a ‘blackbox’, where finding out what exactly it does might sometimes be a little complicated and a bit harder to tweak if you have special requirements. Also installation can be quite painful, however using the VirtualBox this has become a lot easier.

Let’s get started by creating an output folder for our QIIME data processing (mkdir) and defining our raw data folder:
```bash
$ cd /home/qiime/Desktop
$ out=$wosho/output/QIIME
$ mkdir -p $out
$ raw=$wosho/raw_data
```
The first step we need to do is trimming bad quality ends from forward and reverse read. As mentioned in the introduction, the amplicon is almost fully covered twice. This can cause trouble in the read-joining, as the amplicon might be shorter than the sequence reads, leading sequences running over the “edges” of the amplicon. Thus, we truncate both reads to 220bp using VSEARCH (unfortunately, QIIME has no option to do this). We can use a backslash ‘\’ to enter commands that are too long for a single line:
	
```bash
$ mkdir -p $out/trimmed
$ vsearch --fastq_filter $raw/B-1104-S_S70_L001_R1_001.fastq.gz --fastq_trunclen 220 --fastqout $out/trimmed/B.1104.S_R1.fq
$ vsearch --fastq_filter $raw/B-1104-S_S70_L001_R2_001.fastq.gz --fastq_trunclen 220 --fastqout $out/trimmed/B.1104.S_R2.fq
```

Now we can assemble the reads trimmed into a single sequence per read-pair (as mentioned in the introduction, the amplicon is almost fully covered twice). Typing the command in the terminal without any further option will give you information about basic usage.:
```bash
$ join_paired_ends.py
```
Adding “--help” or “-h” to the command will give detailed information about the options that can be specified for data processing (By the way: this is applicable for many other programs as well, e.g. also for VSEARCH):

```bash
$ join_paired_ends.py -h
```
To now join the sequences, we expect a sequence of ~250bp, so we set the minimum overlap between R1 and R2 to 150bp. we run:

```bash
$ mkdir -p $out/join_pairs_output
$ join_paired_ends.py --min_overlap=150 -f $out/trimmed/B.1104.S_R1.fq -r $out/trimmed/B.1104.S_R2.fq -o $out/join_pairs_output/B.1104.S
```
In the output folder we now have one file with the joined forward and reverse sequences (fastqjoin.join.fastq), and two (smaller) files, where merging was not successful, likely due to low quality and/or mismatches in the sequencing. This output can now be used for quality control:

```bash
$ split_libraries_fastq.py -i $out/join_pairs_output/B.1104.S/fastqjoin.join.fastq -o $out/slout_q20/B.1104.S -q 19 --barcode_type 'not-barcoded' --sample_ids B.1104.S
```

The -q 19 option defines, that sequences should be truncated at the first base with a quality score of <= 19. You can have a look at the output files in the final folder to get an Idea what your data looks like now. Especially of interest is the file seqs.fna, as this one now holds the 16S sequences in FastA format. This format is much simpler than the FastQ format. Each sequence has a header containing the sequence ID (the lines that start with the “>” symbol), followed by the sequence itself. Different from the FastQ format, in FastA format the sequence can be running over several lines, ending when the next header starts with “>” or at file end.

One last thing that is critical before we can proceed is the identification and removal of chimeric sequences. Chimeric sequences arise as a PCR artifact due to the high similarity of the amplicons which sometimes causes unspecific pairings of incomplete sequences in PCR, which are then written into a single molecule by the polymerase. The removal is a two step process in QIIME using a reference database to identify chimeras:

```bash
$ identify_chimeric_seqs.py -m usearch61 -i $out/slout_q20/B.1104.S/seqs.fna -o $out/usearch61_chimera/B.1104.S -r $wosho/resources/gold.fa
$ filter_fasta.py -f $out/slout_q20/B.1104.S/seqs.fna -o $out/usearch61_chimera/B.1104.S/seqs_clean.fna -s $out/usearch61_chimera/B.1104.S/chimeras.txt -o $out/usearch61_chimera/B.1104.S -s $out/usearch61_chimera/B.1104.S/chimeras.txt
```

These sequences can now be used for OTU picking and/or sequence annotation. However, since we have more than this sample, we need to do the same for all other samples to have a consistent annotation of OTU sequences. For this we provided a small script to run the previously performed trimming, merging, quality control and chimera removal  for all samples. Please run:

```bash
$ cd /home/qiime/Desktop
$ sh $wosho/scripts/run_qc.sh
```

This will take some minutes, while this is running you may want to proceed reading the next section about OTU picking approaches (Section 4). 

### 3.2 USEARCH/VSEARCH pipeline ###

The USEARCH pipeline also provides options for all steps you can perform in QIIME for preprocessing, however with even more filtering and quality control options. As mentioned before, VSEARCH is a free alternative to USEARCH that is open source and has no memory restrictions, so we will use VSEARCH instead. Let’s first create an output folder:

```bash
$ cd /home/qiime/Desktop
$ out2=$wosho/output/VSEARCH
$ mkdir -p $out2
```

We already used VSEARCH for the read trimming of the QIIME input and we will proceed with that data. We merge the reads and directly apply estimated-error filtering on the merged reads (--fastq_maxee). The estimated error is a composite score for a read calculated from the FastQ quality scores. Additionally we define our merged reads to be between 220 and 280 nt in length:

```bash
$ mkdir -p $out2/merged/B.1104.S
$ vsearch --fastq_mergepairs $out/trimmed/B.1104.S_R1.fq --reverse $out/trimmed/B.1104.S_R2.fq --fastq_maxee 2 --fastq_minmergelen 220 --fastq_maxmergelen 280 --fastaout $out2/merged/B.1104.S/B.1104.S.fna
```

As you see, VSEARCH also directly provides statistics about your reads, making it a little more convenient to keep track of your data. Like in QIIME we will now also filter out chimeric reads. Instead of a reference based filtering we will now use a de-novo approach. This approach uses the provided sequences themselves and differential abundance to identify chimeras, as these are likely to be the product of two other sequences seen in the dataset. ‘On-the-fly’ with this analysis we first dereplicate our sequences to only unique sequences to reduce the computational burden, then pass the sequences through the filtering, and ‘re-replicate’ them, to return to the initial number of sequences (minus the chimeras). This is all done using pipes and the input (/dev/stdin) and output (/dev/stdout) streams :

```bash
$ mkdir -p $out2/nonchimeras/B.1104.S
$ vsearch --derep_fulllength $out2/merged/B.1104.S/B.1104.S.fna --sizeout --output /dev/stdout | vsearch --uchime_denovo /dev/stdin --nonchimeras /dev/stdout | vsearch --rereplicate /dev/stdin --output $out2/nonchimeras/B.1104.S/B.1104.S_clean.fna --relabel B.1104.S_
```
As you see, this is not more complicated than using QIIME, however gives you more flexibility in controlling different parameters. We also run this for all the other samples:

```bash
$ cd /home/qiime/Desktop
$ sh $wosho/scripts/run_qc_vsearch.sh
```
We will later use the different approaches to compare if and how much they differ.

## 4. OTU picking approaches ##

Generally underlying the analysis of marker genes is the assumption, that all living organisms share a single common ancestor - this is called being monophyletic - and that (dis-)similarities in the marker gene are due to shared/different evolutionary trajectories, meaning more similar sequences are a sign if more closely related individuals. The term OTU is an acronym for Operational Taxonomic Unit, which in the end is simply a term for grouping a number of related individuals without further knowledge. In the analysis of marker genes, this grouping can happen in different ways: either based on a collection of annotated reference sequences or without any prior assumptions, solely based on the similarity of the sequences to one another. QIIME includes several different approaches.

### 4.1 Closed-reference OTU picking

Closed-reference OTU picking is entirely based on the similarity of a sequence under investigation to a set of reference sequences of know origin. This approach has the advantage that it is computationally relatively cheap, as the computational cost scales linearly with the number of sequences under investigation (N) and in the reference database (R), meaning that for each sequence a maximum of R comparisons have to be made, adding up to a total of N x R comparisons. Memory consumptions is more or less only dependent on database size. Generally complexity can be reduced by reducing database size based on sequence similarity (see de-novo clustering). Another advantage is, that datasets classified using the same database are comparable to each other. The major disadvantage of closed-reference OTU picking, is that query sequences that are not present in the database will not be identified and discarded from the analysis. Also incorrectly annotated databases (which is certainly not uncommon) and/or missing data can lead to wrong results. Especially species level assignments of sequences should be seen with extreme caution.

### 4.2 De-novo-clustering

De-novo-clustering is the absolute counterpart of reference-based approaches. Instead of comparing sequences to set of known, annotated sequences, this approach uses only the pairwise similarities of the query sequences to each other. Sequences exceeding a pre-defined similarity threshold are collapsed into a single cluster. This threshold is traditionally defined and widely used as 97% similarity (3% divergence) for the same bacterial species, the so-called “species-level OTU” (however, this threshold is permanently debated, very recently for example here: https://doi.org/10.1093/bioinformatics/bty113). Irrespective of the validity of a 3% threshold to define a species, this approach groups similar sequences which are very likely coming from closely related bacterial taxa and thus also likely sharing functional properties, making it a very valid approach for unbiased analysis of marker gene data. Additionally, de-novo approaches are very valuable in situations, where no reference database is available. A downsides of the de-novo clustering is the higher computational demand, as here for a number of N query sequences a total of N(N-1)2comparisons are necessary. Also in terms of memory this is more demanding, as basically all pairwise differences need to be stored while calculating. However, efficient algorithms have been implemented, most notably the UCLUST and UPARSE algorithms, both part of the USEARCH software (details: https://drive5.com). 

### 4.3 Open-reference OTU picking

The open-reference OTU picking is a hybrid of the two approaches that were just introduced. It takes those sequences that could not be matched to a reference database and subjects them to de-novo clustering. While this can eliminate one of the downsides of the closed-reference approach, the quality of the reference database remains a crucial factor influencing the final outcome. 

### 4.4 Amplicons sequence variants (ASVs) / zero OTUs (zOTUs)

Recently, along the discussion of what similarity thresholds to use for OTU picking, efforts have been undertaken to move away from hard cutoffs for clustering to the identification of exact sequences. The most notable contributions to these approaches are the UNOISE pipeline by the USEARCH developer (again: https://drive5.com) and DADA2 (http://benjjneb.github.io/dada2/index.html). Breifly, DADA2 tries to estimate error frequencies to identify sequencing error and correct them in the sequencing data to get really exact amplicon sequences for a finer scale resolution of the data. DADA2 is a core piece of the new QIIME2 toolbox. 

### 4.5 Direct classification of the unclustered sequences

One option is also to entirely skip OTU clustering and directly assign taxonomic labels to the sequences. This of course can take a longer than the other approaches. Also often assignment at low taxonomic levels is not often possible, so this option might not always give the finest resolution of the data. In the end, the classification of the sequences does not differ from e.g. taxonomic annotation of de-novo OTUs and is usually performed using so-called (Naïve) Bayesian Classifiers. One of the most used classifiers is the RDP (Ribosomal Database Project; http://rdp.cme.msu.edu) Classifier, however basically the same concept is also used by the SINTAX classifier (again part of the USEARCH pipeline: https://drive5.com) or the built-in classifier of the DADA2 package in R, to name just some. Briefly, such classifier break down the query sequences into shorter fragments (words) of a specific size (default in RDP:  w=8) and compare them to a database that underwent the same processing. The final annotation is then the assignment with the highest probability (irrespective of how high this probability is or how close other assignment probabilities are). Afterwards, a so-called ‘bootstrapping’ is used to calculate the confidence of this assignment on all levels of the taxonomic annotation. Details can be found in the original publication introducing the RDP classifier: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950982/

### 4.6 Excursion: 16S Sequence Databases

Next to the algorithm used for classification, the reference database is one of the most important things to consider when working with 16S data (or also in general). In 16S analysis, there is the luxury of three major databases to be available to choose from. While none of them is particularly bad, especially the GreenGenes database will likely soon become obsolete, unless anybody revives it with a major update:

GreenGenes: the GreenGenes database (http://greengenes.secondgenome.com/) is very widely used and is the default database for QIIME. Unfortunately, some of the data are wrongly annotated and the last update was in 2013.
SILVA: the SILVA database (https://www.arb-silva.de/) is a large database that is supported by the German network for Bioinformatics (de.NBI) and elixir. It is updated regularly and provides several tools for the analysis and annotation of 16S sequence data.
Ribosomal Database Project: the RDP differs a little from the other two databases, as it follows a slightly different taxonomy (Bergey’s manual) and it is also entirely hand curated. This means that the sequence quality is likely the best, however with the trade-off that updates (new/updated taxa, etc.) might come a little later than with the others.

### 4.7 Processing the data

We will now process the data using the different algorithm. The de-novo approach will be performed twice, once with the QIIME script (using the UCLUST algorithm: http://drive5.com/usearch/manual/uclust_algo.html) and once with VSEARCH for comparison. But first let us get a quick overview of the data again using QIIMEs “count_seqs.py” script:

```bash
$ count_seqs.py -i "$out/usearch61_chimera/*/*clean.fna"
$ count_seqs.py -i "$out2/nonchimeras/*/*clean.fna"
```

What we see is, that we have on average slightly more sequences in the QIIME based quality control, which might come from the stricter estimated-error filtering using VSEARCH. Also we see a little larger mean size in the VSEARCH output. This is likely due to the fact, that in the QIIME analysis, bad quality reads are truncated and in VSEARCH they are discarded. However, both approaches show that more than enough sequences remain after Quality Control, irrespective of the used processing pipeline. For now, we will only work with the QIIME data and come back to the VSEARCH based data once we start working with R. In fact, we will now subsample these to 10,000 sequences per sample, in order to speed-up the OTU picking steps. 

We will program a small loop for that:

```bash
$ for sample in $(ls $out/usearch61_chimera/); 
do
vsearch --fastx_subsample  $out/usearch61_chimera/${sample}/seqs_clean.fna --sample_size 10000 --fastaout $out/usearch61_chimera/${sample}/${sample}.clean.10k.fna
done
```

Once we have done this we collect all sequences from all samples into a single FastA file per approach:

```bash
$ cat $out/usearch61_chimera/*/*.clean.10k.fna > $out/All.QIIME.fna
```

The cat command stands for concatenate. This basically just writes all the files that match the given pattern to the standard output (again, the * are wildcards, matching any name). The “>” symbol redirects everything that is written to the standard output into a file, in this case into “$out/All.QIIME.fna”. So in short, we simply create this file, which holds the data of all cleaned sequences.

First start the closed-ref picking:

```bash
$ pick_closed_reference_otus.py -i $out/All.QIIME.fna -o $wosho/output/otu_picking/QIIME_closed_ref
```
Now start the open-ref picking:
```bash
$ pick_open_reference_otus.py -i $out/All.QIIME.fna -o $wosho/output/otu_picking/QIIME_open_ref -m usearch61
```
 ...and the de-novo picking:

```bash
$ pick_de_novo_otus.py -i $out/All.QIIME.fna -o $wosho/output/otu_picking/QIIME_denovo
```

All OTU picking scripts of QIIME are so-called ‘wrappers’, meaning scripts, which pass the options on to other scripts.  Looking into the log-files in the output folders can give you more information about what the scripts actually do. Use your file browser to go to the output folder of the open-reference picking of the QIIME output data (output/otu_picking/QIIME_open_ref) and open the log file (a text file with the name log_2018[a lot of numbers here].txt). 

The file starts with information about the QIIME version, followed by parameters used in the running of the script, e.g. the paths to the default database used by QIIME (GreenGenes 13.5 release). A little further down you see the commands that were used to generate the data:

pick_otus.py: reference-based OTU picking agains the GreenGenes Database
filter_fasta.py: get sequences that did not match the database (“failures”)
pick_rep_set.py: get reference sequences of those database sequences that  were found in the dataset
pick_otus.py: de-novo OTU picking on “failure” sequences
combine reference-based and de-novo OTUs
this is then followed by a few processing steps and taxonomic annotations.

Looking at the log of the the de-novo picking, you will see, that it is roughly the same starting at step 4. This is actually very nice, as this allows us to use the processing scripts also for the VSEARCH-based OTU-picking. However, let’s do the picking first:

```bash
$ mkdir -p $wosho/output/otu_picking/QIIME_dnVS/vsearch_otus
$ vsearch --cluster_size $wosho/output/QIIME/All.QIIME.fna --relabel OTU_ --uc $wosho/output/otu_picking/QIIME_dnVS/vsearch_otus/dnVS.uc --centroids $wosho/output/otu_picking/QIIME_dnVS/vsearch_otus/dnVS.fa --id 0.97 --relabel_keep
```

Now start a small script which will create the sequence-to-OTU mapping file needed to use the QIIME scripts:

```bash
$ sh $wosho/scripts/create_mapping.sh
```

And now we can follow the same procedure as given in the log-file:

```bash
$ otu_folder=$wosho/output/otu_picking/QIIME_dnVS
$ assign_taxonomy.py -o ${otu_folder}/uclust_assigned_taxonomy -i ${otu_folder}/vsearch_otus/All.QIIME_rep_set.fasta 
$ make_otu_table.py -i ${otu_folder}/vsearch_otus/All.QIIME_otus.txt -t ${otu_folder}/uclust_assigned_taxonomy/All.QIIME_rep_set_tax_assignments.txt -o ${otu_folder}/otu_table.biom 
$ align_seqs.py -i ${otu_folder}/vsearch_otus/All.QIIME_rep_set.fasta -o ${otu_folder}/pynast_aligned_seqs
$ filter_alignment.py -o ${otu_folder}/pynast_aligned_seqs -i ${otu_folder}/pynast_aligned_seqs/All.QIIME_rep_set_aligned.fasta 
$ make_phylogeny.py -i ${otu_folder}/pynast_aligned_seqs/All.QIIME_rep_set_aligned_pfiltered.fasta -o ${otu_folder}/rep_set.tre 
```

While this is running, there is one thing we have not yet talked about, which is … 

## 4.8 Excursion: Alignment and phylogenies

As mentioned before, one of the basic concepts in the analysis of microbial amplicon data, is that all the analyzed organisms share a common ancestor, making them monophyletic. This concept is the basis of the creation of so-called phylogenetic trees, which are created to visualize relations between e.g. organisms. One famous example is Charles Darwin’s Tree of Life:

![alt text](/static/tol.png "Logo Title Text 1")

This tree depicts each organism under study as a tip or ‘leaf’ in the tree, and all tips are connected by the branches in the tree, the closer the organisms are to each other, the ‘earlier’ they are connected at an internal node. While in general this is also what we talk about when looking at bacteria, one property that he did not cover, is that we assume, that all splits in the branches of the tree are bifurcating, meaning that there are exactly two branches going out of each internal node. 

It is possible to reconstruct such trees from amplicon sequencing data, which can then be used in the analysis, as relations between the bacteria can be an important property, assuming that closer related organisms are also likely to share e.g. functions or niche-specificity.

Again, there a multiple ways to go from individual OTU sequences to a tree, each with advantages and disadvantages, and again, these approaches can be placed somewhere between entirely reference-based and completely de-novo. Usually the whole process can be split into two parts: the alignment and the tree construction. 

Running closed-reference OTU picking also includes the most likely easiest way to generate such a tree. The GreenGenes database which is used by default by QIIME provides both, a reference alignment and a reference tree. After having identified the sequences that occur in the sample, these sequences simple need to be ‘extracted’ from the reference. 

However, when running de-novo OTU picking, also these steps can be a little more computationally expensive. For the alignment there are again two choices: either trying to align to a reference alignment or creating this alignment from scratch. For an amplicon based on the well-studied 16S gene, an alignment against a reference is usually the better choice, as this is faster, computationally not very expensive and even if a sequence has no perfect match in the reference alignment, it usually still can be placed close to the closest relative in the database. The algorithm used for this task is usually NAST (Nearest Alignment Space Termination). In QIIME, this alignment is done by the PyNAST (https://doi.org/10.1093/bioinformatics/btp636) tool, and mothur also has its own implementation of the algorithm. One tool that is mostly used for de-novo alignment of sequences is MUSCLE (http://www.drive5.com/muscle/), but there are literally hundreds of different multiple-sequence aligners available. An advantage of de-novo alignment is that it might resolve relations between closely related sequences a bit better, however at the cost of much longer computational time. Also in large datasets, large alignments can always also introduce errors, thus a (possibly even curated) template alignment is always a good place to start with.

Once the alignment is completed, this has to be translated into a tree. This step is really a science of its own with many different model assumptions being possible for the construction (for further reading and to get an overview you can start here: http://www.iqtree.org/doc/Command-Reference#partition-model-options). The default software used in QIIME is FastTree (http://www.microbesonline.org/fasttree/), others are IQ-TREE (http://www.iqtree.org), which provides automatic model selection, or RAxML (https://github.com/amkozlov/raxml-ng/releases). Most important to know about this is the Newick output format, which is mostly used for text-based tree representation (and can be read by most visualization softwares). Here is an example tree in Newick format:

((A:0.1,B:0.2):0.2,(F:0.2,(C:0.1,D:0.1):0.1):0.2);

Which translates to this tree:


Explained: One split in the tree is always separated by a comma, the numbers given after the colon are the branch lengths.

How these trees are part of the analysis will be explained later when talking about alpha- and beta-diversities. 
5. Analysis of amplicon-based microbiome survey data
5.1 QIIME’s ‘core_diversity_analysis.py’ wrapper

QIIME provides a great wrapper script for getting a broad initial overview of the data. All this script needs is a BIOM-file (which is generated by any of the previously presented picking approaches), a phylogenetic tree (same) and a mapping file, which includes ‘meta’-data of all the samples in the analysis. The mapping file (map.tsv) can be found in $wosho/resources. This mapping file can carry information about the samples. There is a single line per sample with a few fields that are mandatory (see here) and additional fields that can be specified by the user. In this example we have “Bear” and “Season” as variables in the mapping file.

We will run this analysis for the four QIIME-based parts of the analysis. Starting with the closed-reference picking (ignore the warnings), specifying to use ‘Season’ as an explanatory variable and subsampling to 5000 sequences per sample:

$ qcr=$wosho/output/otu_picking/QIIME_closed_ref
$ core_diversity_analyses.py -i $qcr/otu_table.biom \ 
-m $wosho/resources/map.tsv -c Season -t $qcr/97_otus.tree \
-o $wosho/output/analysis/QIIME_closed_ref --sampling_depth 5000

Followed by the open-ref picking and the QIIME de-novo approach (this one also for the VSEARCH pipeline):

$ qor=$wosho/output/otu_picking/QIIME_open_ref
$ core_diversity_analyses.py -i \
$qor/otu_table_mc2_w_tax_no_pynast_failures.biom \
-m $wosho/resources/map.tsv -c Season -t $qor/rep_set.tre \
-o $wosho/output/analysis/QIIME_open_ref --sampling_depth 5000
$ qdn=$wosho/output/otu_picking/QIIME_denovo
$ core_diversity_analyses.py -i $qdn/otu_table.biom \
-m $wosho/resources/map.tsv -c Season -t $qdn/rep_set.tre \
-o $wosho/output/analysis/QIIME_denovo --sampling_depth 5000
$ qdnvs=$wosho/output/otu_picking/QIIME_dnVS
$ core_diversity_analyses.py -i $qdnvs/otu_table.biom \
-m $wosho/resources/map.tsv -c Season -t $qdnvs/rep_set.tre \
-o $wosho/output/analysis/QIIME_dnVS --sampling_depth 5000

While these are running, we can have a look at the output of the closed-reference core diversity analyses, which gives an opportunity to explain fundamental terms and concepts in the analysis of amplicon data. Please navigate to the output folder (workshop/output/analysis/QIIME_closed_ref).
5.2 Taxa summaries

Open the “index.html” file in the Chrome browser to see an overview of all available subcategories that were created by the script. Click on “bar_charts.html” in the “Taxonomic summary results” (NOT by Season!). What you see on top is a bar-chart summarizing the microbial communities in the samples on the phylum level. The phylum is the broadest taxonomic level below the “kingdom” or “domain” level. If you compare your results with those of your neighbors, you will see that they might look slightly different from yours. That is, because we randomly subsampled twice in the process of getting to these results (once before OTU picking and then again now before running the diversity analyses). Below the plot you can see a table with all the phyla present in the dataset and their respective relative abundance in the samples. You will notice that there are six phyla which make up the majority of the community in all samples: Firmicutes, Bacteroidetes, Proteobacteria, Actinobacteria, Fusobacteria and Cyanobacteria. And you will also notice, that even at this level, the distributions of these six major groups already differ quite a lot between samples.

One Sample that looks different from the others is “B.1104.S”. In this sample you see that the bacterial community consists almost entirely of bacteria of the phylum “Cyanobacteria” and you might think “what the hell are blue-green algae doing in brow-bear poop?”...and by this you are totally right. If you scroll down a little further you will see similar plots for more fine-scale taxonomic levels from class to genus. Having a look at the class level plot for the same sample will show you, that these sequences are all assigned to the class “Chloroplast” in the phylum “Cyanobacteria”. Maybe you recall from your microbiology lecture that Chloroplasts are organelles that originate from the endosymbiosis of cyanobacteria into other cells. Through millions of years of evolution, these endosymbionts still kept some of their genomic DNA and thus still have a 16S rRNA gene, so these sequences are very likely coming from food remnants and undigested plant material in the bears diet (and not from a secret ocean in the bears gut...). The same can actually be seen for mitochondria, which very likely derive from a bacterium of the phylum Proteobacteria and can be found in the family level table in the order Rickettsiales.

You can go through these plots a little more if you like, but in the end you will see, that the finer the classification gets, the harder it gets to grasp what is going on in the samples. As mentioned before, we defined “Season” as an explanatory variable in the processing script Go back to index page and again open the “bar_charts.html” file, now by Season. Now you see the same plot succession as before, however as summaries for the seasons summer (S) and winter (W). This can give you an overview about differences between the groups. Going back to the index page, you will also find the section “Group significance results” which includes a very large text-file of all results from the differential abundance testing for all OTUs (more info here: http://qiime.org/scripts/group_significance.html). By default QIIME uses a G-test (aka goodness of fit, log-likelihood ratio test) to test for deviances in the frequencies between the groups under investigation. 

5.3 Alpha diversity

Alpha diversity is a name for the within-sample diversity. Many different indices exist to describe specific properties of alpha-diversity, most of them are designed to have a higher value for samples with higher diversity. Automatically provided by QIIME are Faith’s Phylogenetic Diversity (PD_whole_tree), Chao1 estimator (chao1) and the number of observed OTUs (observed_otus). 



5.3.1 Number of observed OTUs

The number of observed OTUs is probably the easiest to calculate measure of alpha diversity. For this, you simply count the number of OTU per sample, that have an abundance greater than 0. QIIME automatically provides a boxplot of the values and a test statistic which can be accessed over the index.html file (this is also true for the other indices provided by QIIME).

5.3.2 Chao1 estimator

The Chao1 estimator is a descriptor of alpha-diversity that tries to estimate (well...obviously) the number of OTUs including those, that were not actually observed. The calculation is based on the number of low-abundant OTUs observed in the sample and follows the rationale: if you observe many low-abundant OTUs, it is likely that there are even more that were not observed due to too low sequencing depth or too shallow subsampling. If none or only a few are observed, it is more likely that the whole community was sampled at sufficient depth.

5.3.3 Faith’s Phylogenetic Diversity (PD)
 
This measure of alpha diversity not only takes into account the OTU abundances, but also the relationship of the OTUs to each other by using the previously mentioned phylogenetic tree constructed from the OTU sequences. Briefly, the PD calculates the branchlengths of the tree spanning all OTUs observed in a sample.


5.3.4 Shannon Diversity

An widely used index for alpha diversity not covered in QIIME by default is the Shannon index (Shannon entropy; H’), which comes from the field of information theory and is basically a measure of uncertainty. It is calculated as:
 H'= -i=1Rpi ln pi 
With pi being the relative abundance of organism i. The Shannon entropy is highest with a high number of entries that are equally distributed. In some articles the Shannon entropy is not used directly, but as the ‘species equivalent’, which is calculated as eH'. This can be translated to the number of equally abundant species needed to have this entropy.

5.3.5 Rarefaction curves

In the QIIME output index you will see the link to “rarefaction_plots.html” in the alpha diversity subsection. Here you can choose one of the presented indices and a category to use for comparison. Try it with “PD_whole_tree” and “Season” as category. You are then shown a plot with a red and a blue line, representing the Seasons, and the x-Axis called “Sequences Per Sample”. From the rarefaction curves you can estimate, if the sampling depth (5000 in our case) is enough to cover the full diversity estimated to be present in the samples. The rationale behind this analysis is, that if the diversity is fully covered (saturated), this curves go into a plateau, which means, that the inclusion of additional sequencing data will not change the diversity anymore. You will see, that 5000 sequences are actually not enough to analyze the bear data to its full extend, however you also see, that the “Winter” curve is very likely to stay below the “Summer” curve.
5.4 Beta diversity

The beta diversity is a measure for diversity not within samples, but between samples (inter-sample diversity). This already implies, that there is not a single value of beta diversity for each sample, but rather for each pair of samples, resulting in a NxN distance matrix when analyzing N samples. As for alpha diversity, there are a lot of different measures for beta diversity, each with its own assumptions and possibilities of application.
5.4.1 Jaccard index

The Jaccard index is one of the most basic measures of beta diversity. It does not take abundances into account, but only the presence (and absence) of the taxa under investigation. As mentioned before, beta diversity indices are always calculated pairwise for each samples, and the Jaccard index can be summarized as the proportion of shared taxa between two samples. Let’s assume an abundance table of three samples (S1, S2 and S3; each with 100 sequence counts) that looks like this:



OTU_1
OTU_2
OTU_3
OTU_4
OTU_5
OTU_6
S1
20
10
0
0
50
20
S2
30
10
10
0
0
50
S3
0
0
30
70
0
0




Now we want to calculate the Jaccard index for Samples S1 and S2:

Present OTUs in S1: OTU_1, OTU_2, OTU_5 and OTU_6 → OS1	 = 4
Present OTUs in S2: OTU_1, OTU_2, OTU_3 and OTU_6 → OS2 	= 4
Present OTUs in S1 and S2: OTU_1, OTU_2 and OTU_6 → OS1,S2	 = 3 

J(S1,S2) = OS1,S2 / (OS1 + OS2  - OS1,S2 ) = 3 / (4 + 4 - 3) = 3 / 5 = 0.6

This simply means, that 0.6 = 60% of the community seen in S1 and S2 together is shared by both samples. As you see, this is actually a measure of similarity! For the Jaccard dissimilarity, we simply need to subtract this value from 1:

dJ(S1,S2) = 1 - J(S1,S2) = 1 - 0.6 = 0.4

We can calculate this for the pairs (S1,S3) and (S2,S3) as well:

J(S1,S3) = 0 / (4 + 2 - 0) = 0 / 6 = 0
⇒ dJ(S1,S3) = 1 - J(S1,S3) = 1 - 0 = 1

J(S2,S3) = 1 / (4 + 2 - 0) = 1 / 6 = 0.1666
⇒ dJ(S2,S3) = 1 - J(S2,S3) = 1 - 0.1666 = 0.8333

You see, Sample S1 and S3 do not share any OTUs, thus their dissimilarity is 1. Samples sharing every OTU would have the value 0 . This is a basic property of most dissimilarity measures, taking values between 0 (= being the same) and 1 (=being absolutely different from another).

5.4.2 Bray-Curtis dissimilarity

The Bray-Curtis dissimilarity (not distance!) is another widely used measure for describing beta diversity. In general it follows a similar pattern as the Jaccard dissimilarity, however not only looking at presence/absence, but also incorporating the abundances of the OTUs / taxonomic groups.  It is defined for samples A and B as:

BC(A,B) = 1 - 2 i=1nmin(CA,i,CB,i)CA + CB

with CA and CB being the total counts of samples A and B, respectively, and CA,i and CB,i being the counts of OTU i for samples A and B, respectively. 

For samples S1, S2 and S3 from the table in the previous section this is then:

BC(S1,S2) = 1 - 2  (20 + 10 + 0 + 0 + 0 + 20)100 + 100= 1 - 2  50200= 1 - 0.5 = 0.5

BC(S1,S3) = 1 - 2  (0 + 0 + 0 + 0 + 0 + 0)100 + 100= 1 - 2  0200= 1 - 0 = 1

BC(S2,S3) = 1 - 2  (0 + 10 + 0 + 0 + 0 + 0)100 + 100= 1 - 2  10200= 1 - 0.1 = 0.9

You will mostly see the pairwise differences as distance matrices, e.g.:



S1
S2
S3
S1
0




S2
0.5
0


S3
1
0.9
0

5.4.3 UniFrac distances

The weighted and unweighted UniFrac (https://en.wikipedia.org/wiki/UniFrac) distances are the default distances used by QIIME when following the default analysis path implemented in the “core_diversity_analysis.py” wrapper. UniFrac uses the phylogenetic tree and thus the relationship of OTUs for beta diversity calculation similar as the Phylogenetic Diversity does it for alpha diversity. The unweighted UniFrac uses presence/absence values only for the distance calculation, like in the Jaccard distance. The weighted UniFrac weighs branches according to the abundance of the tips, so this is more comparable to the Bray-Curtis dissimilarity. Introduced in 2005, the UniFrac distances are still debated, with critics stating they put too much weight on certain branches (rare and abundant, respectively, depending on which measure is used). In 2012 a generalized UniFrac (gUniFrac) was introduced, which is supposed to create a more balance weighing by combining the weighted and unweighted UniFrac into a single distance measure.
5.4.4 Beta diversity in QIIME

QIIME calculates UniFrac distances (all accessible through the index page) and subjects these to PCoA (Principal Coordinate Analysis). PCoA is a dimension-reduction technique that can be applied to distance matrices, trying to project as much variation possible onto a set of orthogonal axes, meaning that in the resulting plot (3D-space in QIIME by default) a point is placed into space for each sample. Similar samples (so samples with less dissimlarity/distance) are closer to one another in space, than samples that are more dissimilar. In QIIME these plots can be accessed through the “Emperor” tool, which is an interactive viewer of such plots. Please open the PCoA Plot “index.html” for weighted UniFrac from the index website in the browser. You will see a plot with red spheres on black background that you can move by using the mouse. On the right hand side, you can choose the tab “Colors” and choose “Season” in the second DropDown menu. Now you can see, that in the plot your Winter and Summer samples are colored differently, and it should be relatively obvious, that they make up two separate clusters, meaning there is a difference in the microbial communities of brown bears in Summer and Winter while hibernation. 
5.5 QIIME Summary

As you see, QIIME provides a good framework to get a broad analysis of microbiome data very easily. A lot more analyses are certainly possible and the (very nice) documentation of the available scripts can be accessed via the QIIME Scripts website (http://qiime.org/scripts/index.html) and, as mentioned earlier, by the --help argument, which usually provides a good overview of the functions directly within the commline window.

To conclude the QIIME based analysis, we will use the “compare_taxa_summaries.py” to get an idea, how and if the approaches differ in their outcomes (L5 are family level abundance tables; L6 genus level):

$ cd $wosho
$ compare_taxa_summaries.py -i \
output/analysis/QIIME_closed_ref/taxa_plots/table_mc5000_sorted_L5.txt,\
output/analysis/QIIME_open_ref/taxa_plots/table_mc5000_sorted_L5.txt \
-m paired -o output/comp/fam_qcr_qor --perform_detailed_comparison
$ compare_taxa_summaries.py -i \
output/analysis/QIIME_closed_ref/taxa_plots/table_mc5000_sorted_L5.txt,\
output/analysis/QIIME_denovo/taxa_plots/table_mc5000_sorted_L5.txt \
-m paired -o output/comp/fam_qcr_qdn --perform_detailed_comparison
$ compare_taxa_summaries.py -i \
output/analysis/QIIME_closed_ref/taxa_plots/table_mc5000_sorted_L6.txt,\
output/analysis/QIIME_denovo/taxa_plots/table_mc5000_sorted_L6.txt \
-m paired -o output/comp/gen_qcr_qdn --perform_detailed_comparison
$ compare_taxa_summaries.py -i \
output/analysis/QIIME_dnVS/taxa_plots/table_mc5000_sorted_L5.txt,\
output/analysis/QIIME_denovo/taxa_plots/table_mc5000_sorted_L5.txt \
-m paired -o output/comp/fam_vdn_qdn --perform_detailed_comparison
$ compare_taxa_summaries.py -i \
output/analysis/QIIME_dnVS/taxa_plots/table_mc5000_sorted_L6.txt,\
output/analysis/QIIME_denovo/taxa_plots/table_mc5000_sorted_L6.txt \
-m paired -o output/comp/gen_vdn_qdn --perform_detailed_comparison

Looking through the output files created by these comparisons (in the folders in output/comp), we can conclude, that irrespective of the used quality control and OTU picking approaches, there is a very high correlation between the final taxon abundance tables. The highest dissimilarities are introduced at genus level taxon annotation, which might differ in some samples (look at the family and genus level comparison between QIIME closed-ref und QIIME denovo).

Also when comparing the UniFrac distance matrices using “compare_distance_matrices.py” we will see, that all comparisons are highly significant:

$ compare_distance_matrices.py -i \
output/analysis/QIIME_closed_ref/bdiv_even5000/weighted_unifrac_dm.txt,\
output/analysis/QIIME_denovo/bdiv_even5000/weighted_unifrac_dm.txt \
--method mantel -o output/comp/mantel_qcr_qdn

$ compare_distance_matrices.py -i \
output/analysis/QIIME_open_ref/bdiv_even5000/weighted_unifrac_dm.txt,\
output/analysis/QIIME_denovo/bdiv_even5000/weighted_unifrac_dm.txt \
--method mantel -o output/comp/mantel_qor_qdn
$ compare_distance_matrices.py -i \
output/analysis/QIIME_dnVS/bdiv_even5000/weighted_unifrac_dm.txt,\
output/analysis/QIIME_denovo/bdiv_even5000/weighted_unifrac_dm.txt \
--method mantel -o output/comp/mantel_vdn_qdn

The comparison used here is the so-called “Mantel-Test”, which can be used to evaluate the correlation between two matrices.

6. Microbiome analysis using R
6.1 Data from DADA2

Unfortunately, DADA2 is not supported by the R version installed in the QIIME 1.9 VirtualBox. However, the abundance data we will now use has been created following the default DADA2 protocol that can be accessed here: http://benjjneb.github.io/dada2/tutorial.html, apart from the database used for sequence annotation, for which we used the same (GreenGenes 13.8) that is default in QIIME. As mentioned earlier, DADA2 does not rely OTU similarity threshold, but tries to estimate error frequencies and in the end to identify exact sequence variants. DADA2 provides a complete workflow, including quality control, error estimation, dereplication, chimara (bimera) removal and sequence annotation, and can thus be used as a complete replacement for QIIME, USEARCH, Mothur, etc.

The abundance and taxonomic annotation table can be found in the ‘Shared_Folder/dada2_out’ directory and we will need this a little later in the analysis.
6.2 Data conversion (PLEASE READ CARFEFULLY!!!)

While the R package ‘phyloseq’ (which we will install in the next step) actually provides a function to import data in biom-format, the versions of QIIME and phyloseq are somehow not compatible, throwing errors when trying to import. However, this is not a huge problem, as we can convert biom-data to text files for convenient import. For example:

$ cd $wosho
$ gunzip output/analysis/QIIME_open_ref/table_even5000.biom.gz
$ biom convert -i output/analysis/QIIME_open_ref/table_even5000.biom \
-o output/analysis/QIIME_open_ref/otu_table.txt --to-tsv \
--header-key taxonomy

When you go to the specified output folder and open the ‘otu_table.txt’ file in a text editor, you will notice, that the file starts with a header stating ‘# Constructed from biom file”. The second line again starts with the “#” symbol, followed by “OTU ID” and the sample names, ending with the column “taxonomy”, in which the taxonomy information for the given OTU (each line is one OTU) is stored. To import this into R, we first need to remove the “#” symbol from the second line, since otherwise R would skip this line while importing. Please do so and save the file.

6.3 Installing R packages

The beautiful thing about the R statistics software is, that it is very flexible and virtually anyone can contribute new packages or code pieces tailored to one’s own analysis workflows. A great place for packages and applications from the different fields of computation biology is Bioconductor (https://www.bioconductor.org). One package we will need today is “phyloseq” (https://joey711.github.io/phyloseq), which provides a large toolbox for the analysis and visualization of microbiome data, and is very extensively documented with a lot of examples and example datasets on the website. 

To install phyloseq, we first need to start R:

$ R

From now on, all commands that should be run in R will be indicated by ‘R>’ at the beginning of the line. In R, long commands can simply be continued in the next line, not needing any special characters. Now we have to install the phyloseq package and activate it (when the window asks you, if you want to use a local directory and if you want to update packages, type “y” for both and hit enter; this will take a few minutes).

R> source('http://bioconductor.org/biocLite.R')
	R> biocLite('phyloseq')
R> library(phyloseq)

6.4 Importing data into R

One thing before we go on: in R you can always get more information about what a command does and how it is used by simply typing a “?” followed by the command, e.g.:

	R> ?library

You can scroll by using the arrows on the keyboard and leave the page again by hitting “q”.

In R, you can import data in text format using the “read.table” command. We will now import the OTU table we just converted, telling R to use “tab” as field sepearator (sep=”\t”), check the dimensions of the imported table (dim(seqtab)) and then look at the first lines of the imported table (head() command):

R> seqtab<-read.table("output/analysis/QIIME_open_ref/otu_table.txt", 
head=T, sep="\t",stringsAsFactors=F,row.names=1)
R> dim(seqtab)
R> head(seqtab)

The table has roughly 1600 rows and 23 columns. You see that we now have our OTUs in rows and the samples in columns, with the last column being the taxonomic annotation. We need to split this now into a table containing only the abundance data and one containing the annotation:

R> otutab <- seqtab[,1:22]
	R> tax <- seqtab[,23,drop=F]
R> head(tax)

If you haven’t worked in R before, this needs some explanation: We can access rows and columns in R using the square brackets [ and ]. The comma separates  between rows (in front of the comma) and columns (after the comma). So the first command here is creating a new object called “otutab” and saving columns 1 to 22 of the data frame “seqtab” in it. The assignment is made by the arrow “<-”. In analogy, column 23 of “seqtab” is saved in the object “tax”. The tax object now has to be converted to a matrix with columns being the assignments from kingdom to species level:

R> taxtab <- t(sapply(tax$taxonomy, function(x) { vv<-rep(NA,7);
y<-strsplit(gsub("; ","\t",x),split="\t")[[1]];
vv[1:length(y)]<-y; return(vv)}))
R> taxtab[grep("[kpcofgs][_][_]$",taxtab)]<-NA
R> rownames(taxtab) <- rownames(otutab)
R> colnames(taxtab) <- c("Kingdom","Phylum","Class","Order","Family",
"Genus","Species")
R> head(taxtab)
R> dim(taxtab)

This command is a little harder to comprehend, and it is also not that important, how it actually works, in the end you just see, that there is now a table with one row of taxonomy information for each OTU. 

The last thing we need now is the sample metadata. We import it from the map file we also used in QIIME:

R> metadata<-read.table("resources/map.tsv",head=F, 
row.names=1)[,1:2]
R> colnames(metadata)<- c("BearID","Season")


These three tables can now be assembled to a single phyloseq object (putting brackets around the command will show you the content of the created object):

R> (ps<-phyloseq(otu_table(otutab,taxa_are_rows=T), tax_table(taxtab), 
sample_data(metadata)))





6.5 Data exploration with phyloseq

For the data exploration, we will not work with the OTU level abundance, but with the genus level annotation. In phyloseq we can easily create a phyloseq object at this taxonomic level using tax_glom:

R> (ps.gen <- tax_glom(ps, taxrank="Genus",NArm=F))

The output will immediately show you, that instead of ~ 1600 taxa we are now left with roughly 250. The option NArm=F tells R to keep also those taxonomic groups that do not have genus level annotations. Let’s get an overview of the samples. For this we first identify the ten most abundant genera in the dataset and then plot their distributions in the samples:

R> gen.sum <- tapply(taxa_sums(ps.gen), tax_table(ps.gen)[, "Genus"], 
sum, na.rm=TRUE)
R> top10gen = names(sort(gen.sum, TRUE))[1:10]
R> ps.gen.sub = prune_taxa((tax_table(ps.gen)[, "Genus"] 
%in% top10gen), ps.gen)
R> plot_bar(ps.gen.sub,fill="Genus",facet_grid=~Season)

Already by visual inspection we can see several things:
Winter and Summer microbiota differ a lot
While the Summer microbiota is dominated by Helicobacter and Streptococcus, we see the genera Bacteroides and Enterococcus as the most abundant tax in Winter.

Let’s proceed with an ordination plot to investigate these differences more. phyloseq provides a pretty simple function for that. We will create an ordination using Multidimensional Scaling (“MDS”; same as PCoA) and using the Bray-Curtis dissimilarity. In the plot we color the samples by Season:

R> ps.gen.ord_bc <- ordinate(ps.gen, "MDS", "bray")
R> plot_ordination(ps.gen, ps.gen.ord_bc, type="samples",
color="Season", title="Bray-Curtis")

We see, that the ordination using Bray-Curtis dissimilarity is very similar to what we saw in the UniFrac ordination plots created by QIIME. What we can now look at to gain more information is looking at the ‘biplot’. The biplot can be used to identify which taxa influence the structure in the ordination. We first store the output of the ordination in a data frame (ord_df) and then keep only those lines in the data frame that are either in the top 10 most abundant genera (as we defined before) or belong to the samples:

R> ord_df <- plot_ordination(ps.gen, ps.gen.ord_bc, type="biplot",
color="Season", title="Bray-Curtis",justDF=T)
R> ord_df.sub <-  ord_df[ord_df$Genus %in% top10gen | 
ord_df$id.type=="Samples",]

To plot this we need to use the package ggplot2 (https://ggplot2.tidyverse.org), which provides great and very flexible plotting functions. Its syntax is not completely intuitive from the beginning, however, once you got used to it, it really is a great and powerful tool for data visualization (especially when it comes to a little more complex figures). We plot:

R> library(ggplot2)
R> ggplot(ord_df.sub,aes(x=Axis.1,y=Axis.2, col=Season, 
shape=id.type))+geom_point(size=4)+geom_label(aes(label=Genus))

Nice huh? Okay, you probably now also want to know what this plot means: The dots represent you samples colored by season and the text is where the 10 selected taxa are placed in the plot. This placement relies on the abundance of the taxa in the samples, meaning that samples going lying in the same direction as the label (starting from the origin) have higher abundance values for this taxon. For example: Bacteroides, Fusobacteria and Porphyromonas are in the lower-left quadrant of the plot, thus the samples being also there have higher values of these taxa. Helicobacter goes straight to the right of the origin, so samples right of the origin have higher values in Helicobacter.

Okay, of course we need to put “a number” on that, to convince everyone, what the differences between summer- and winter-microbiota really are. But first, let us take a step back again, and return to the comparison of the processing pipeline.
6.6 Pipeline comparison

We prepared for you a file with all the pipeline outputs imported into R format on genus level (including the DADA2 output), because as you saw before this can be pretty tedious and time consuming, and also, so all can work with the same data now. Additionally, all sequences assigned to “Chloroplast” were now removed from the dataset.

Import the R-Object file with the complete dataset:

R> load("dada2_out/psall.Robj")
R> psall

The dataset includes all the metadata, with the used processing pipeline being encode in the “Pipe” column.  We re-plot our barplots of the most abundant genera.

R> gen.sum <- tapply(taxa_sums(psall), tax_table(psall)[, "Genus"], 
sum, na.rm=TRUE)
R> top10gen = names(sort(gen.sum, TRUE))[1:10]
R> ps.gen.sub = prune_taxa((tax_table(psall)[, "Genus"] 
%in% top10gen), psall)
R> plot_bar(ps.gen.sub,"Pipe",fill="Genus",facet_grid=~SeqID)+
scale_fill_brewer(palette="Paired")


We can see immediately, that the variation due to the used processing pipeline is absolutely minor. The largest differences are between DADA2 and everything run in QIIME, however, this is mostly the two genera Escherichia and Klebsiella. Theses are two genera belonging to the family Enterobacteriaceae, which are known to be problematic to be assigned with genus level specificity, because of the evolutionary close relationship and also the sometimes historically inaccurate naming, e.g. the genus Shigella is debated to actually be a member of the genus Escherichia, and they are hardly (if at all) separable by 16S sequencing. 

This becomes especially evident, when we plot families instead of genera:

R> fam.sum <- tapply(taxa_sums(psall), tax_table(psall)[, "Family"], 
sum, na.rm=TRUE)
R> top10fam = names(sort(fam.sum, TRUE))[1:10]
R> ps.fam.sub = prune_taxa((tax_table(psall)[, "Family"] 
%in% top10fam), psall)
R> plot_bar(ps.fam.sub,"Pipe",fill="Family",facet_grid=~SeqID)+
scale_fill_brewer(palette="Paired")

This is good news! It shows, that irrespective of what pipeline you use, you will get similar results! Let’s look at the biplot again (family and genus level):

R> gen.ord_bc <- ordinate(psall, "MDS", "bray")
R> ord_df.gen <- plot_ordination(psall, gen.ord_bc, type="biplot",
color="Season", title="Genus Bray-Curtis",justDF=T)
R> ord_df.gen.sub <-  ord_df.gen[ord_df.gen$Genus %in% top10gen | 
ord_df.gen$id.type=="Samples",]
R> plot1=ggplot(ord_df.gen.sub,aes(x=Axis.1,y=Axis.2, col=SeqID, 
shape=Pipe)) + geom_point(size=4)
	R> plot2=ggplot(ord_df.gen.sub,aes(x=Axis.1,y=Axis.2, col=Season, 
shape=id.type))+geom_point(size=4)+geom_label(aes(label=Genus))
R> psall.fam <- tax_glom(psall, taxrank="Family",NArm=F)
R> fam.ord_bc <- ordinate(psall.fam, "MDS", "bray")
R> ord_df.fam <- plot_ordination(psall.fam, 
fam.ord_bc, type="biplot",color="Season",
title="Family Bray-Curtis",justDF=T)
R> ord_df.fam.sub <-  ord_df.fam[ord_df.fam$Family %in% top10fam | 
ord_df.fam$id.type=="Samples",]
	R> plot3=ggplot(ord_df.fam.sub,aes(x=Axis.1,y=Axis.2, col=SeqID, 
shape=Pipe)) + geom_point(size=4)
	R> plot4=ggplot(ord_df.fam.sub,aes(x=Axis.1,y=Axis.2, col=Season, 
shape=id.type))+geom_point(size=4)+geom_label(aes(label=Family))
R> library(gridExtra)
R> grid.arrange(plot1,plot2,plot3,plot4)

You should now see four plots: two genus level on top and to family level on bottom. In the left plots, the data points are colored by sample and the shapes used for plotting are the used pipelines. As seen before, points from the same sample cluster closely together, only the DADA2 output is a little off in the genus level plot, however we see almost perfect concordance in the family level plots. Keep in mind, that we also downsampled to only 5000 sequences per sample, which might introduce some uncertainty (and you have seen in the rarefaction curves that not all possible diversity is covered by this downsampling).

On the right hand side, you see the plots similar to the ordination we did before, and I guess you can make sense of it yourself now. However, what we also see is, that there seem to be two clusters in the summer samples. We will have a look at this later.

One key point now is assessing, if the summer and winter samples differ significantly in their beta diversity. For this we can use the adonis function from the vegan package, which is designed for the analysis of ecological data. This analysis is based on a Multivariate Analysis of Variance (MANOVA) based on distance-matrices introduced by Marti J. Anderson  (“A new method for non-parametric multivariate analysis of variance “, Austral Ecol, 2001). Briefly, the method uses within- and between-group variability to estimate how much variance is explained using the factor under investigation, followed by permutation of the group labels to assess, if this amount of variation explained is greater than random. The adonis function makes use of the so-called “formula interface” to specify factors (covariates), i.e. temperature, age, altitude, etc., for which the variation should be excluded first, to assess only the remaining variation explained by the chosen explanatory variable. Due to the lack of covariates in our analysis, this is  straightforward:


  	R> library(vegan)
R> adonis(distance(psall, method="bray") ~ Season,
data=data.frame(sample_data(psall)), permutations=10000)


The resulting table summarizes how strongly winter- and summer-microbiota differ in terms of beta-diversity. The R2 value tells us, that roughly 32% of the total variation in the dataset are explained by the sampling Season (thus by hibernation). The column “Pr(>F)” gives the permutational P-value, indicating, that no random had better explanatory power than the true Season labels.

Alpha diversity can easily be accessed in phyloseq by the “plot_richness” function. Omitting the “measures” option will result in a plot of all alpha diversity measures available in phyloseq. However, we will no only look at the shannon diversity:

R> plot_richness(psall, x="Season",col="SeqID", 
measures=c("Shannon")) + facet_grid(~Pipe)

Here we actually see partially inconsistent results based on the used pipeline. While DADA2 and closed-ref OTU picking show no clear differences between summer and winter, the others show lower diversity in the winter samples. However, all identify a very-low diversity subgroup in the Summer samples.

We can also calculate alpha-diversity using the vegan package for R and the included “diversity” function:

	R> sha <- diversity(t(otu_table(psall)))

Using a Wilcoxon Rank Sum test with each Pipeline subset reveals no significant P-values:

R> t(sapply(unique(as.character(sample_data(psall)$Pipe)), 
function(x) wilcox.test(sha ~ sample_data(psall)$Season, 
subset=sample_data(psall)$Pipe==x)$p.value))

Finally, we are interested in, which microbial taxa are found to be differentially abundant in summer and winter. We will investigate this using the DADA2-based genus level abundance.

R> ps.dada<-prune_samples(sample_data(psall)$Pipe=="DADA", psall)

As microbiome data is usually highly dimensional, mostly with more microbial taxa present in the dataset than samples, one needs to reduce dimensions to have enough statistical power to detect differential abundance even after multiple testing correction. We do this by reducing the dataset to only those taxa, that show a mean abundance of at least 1% (=50 sequence counts):

R> (ps.dada.sub<-prune_taxa(rowMeans(otu_table(ps.dada))>=50,ps.dada))


This leaves us with 16 taxa we can test. Let us have a look a the distribution of the taxa first:
R> par(mfrow=c(4,4));apply(otu_table(ps.dada.sub),1,hist,nclass=10)


These histograms show us some very important properties of microbiome data (even though only based on 22 samples). What you can see easily, is that microbiome abundance data is rarely normally distributed, meaning following a nice bell-shaped curve. What you rather see is a highly skewed distribution, with a lot of zero or close to zero values and occasional outliers. This is one of the main reasons why finding a suitable model for the analysis is sometimes not that easy. 

One path you can basically always follow is using non-parametric test like the already mentioned Wilcoxon rank sum test (when comparing two groups) or the Kruskal-Wallis-test (when comparing more than two groups) to name only two.  
 
We can run a small loop to test our subset of taxa for differential abundance using the Wilcoxon test, afterwards we adjust the p-values for multiple testing:

R> outvec<-vector(length=16)
R> names(outvec)<-rownames(otu_table(ps.dada.sub))
R> for(i in 1:16){
	tax=rownames(otu_table(ps.dada.sub))[i]
	outvec[tax]<-wilcox.test(as.numeric(otu_table(ps.dada.sub)[i,])~ 
		sample_data(ps.dada.sub)$Season)$p.value
	}
R> (outtab=data.frame(p=outvec,p.adj=p.adjust(outvec)))
R> (outtab[outtab$p.adj<0.05,])

The last command writes the filtered table of only those taxa with differential abundance after correction for multiple testing. We find seven of the 16 candidate taxa to show seasonal variation in hibernating brown bears. Let’s visualize this with boxplots:

R> par(mfrow=c(2,4))
	R> for(tax in rownames(outtab[outtab$p.adj<0.05,])){
taxshort=strsplit(tax,split=";")[[1]][6]
boxplot(as.numeric(otu_table(ps.dada.sub)[tax,])~
sample_data(ps.dada.sub)$Season,
main=paste(taxshort, signif(outtab[tax,"p.adj"])))
}

These show us, that we can replicate three of the four OTUs described in the paper as being differentially abundant (Bacteroides, Wolbachia and Helicobacter). Interestingly, we can even add further candidates from the family Enterobacteriaceae and the genus Enterococcus to the list of seasonally changing microbial taxa in brown bears. The genus Streptococcus also showed a p-value below 0.05 in the Wilcoxon-test (p=0.044), however this did not survive the correction for multiple testing. 

(Depending on the time, if you like you can now try and see, if these results change when using a different pipeline output; SPOILER: it does!)

With this we have pretty extensively analyzed a rather small dataset consisting of only 22 samples. In the next step we will increase sample size a little, looking at a publicly available dataset of fecal samples from treatment-naive new-onset cases with chronic inflammatory bowel disease.

R> quit

7. Treatment-naive Crohn’s disease dataset

Now that we got our feet wet with a small dataset, we will peek at one that is a little larger. The dataset comes from a study published in 2014 (https://dx.doi.org/10.1016%2Fj.chom.2014.02.005) and can be downloaded from QIITA (https://qiita.ucsd.edu/study/description/1939) which is “an open-source microbial study management platform”, and provides access to a wide range of QIIME pre-processed datasets. We provided the needed files in the Shared_Folder. Get R running again:

$ R

And load the data:

R> library(phyloseq) 
R> load("resources/gevers.Robj")
R> ps

As you can see, the dataset includes 321 samples. These samples were already pre-filtered to include only samples from the terminal ileum of early-onset Crohn’s disease patients (diagnosis: CD) and healthy controls (diagnosis: no) between 10 and 25 years, that did not yet receive any medication. The available metadata is relatively extensive:

R> head(sample_data(ps))

We subsample our dataset to an even depth of 3,000 sequences per sample:

	R> ps.even <- rarefy_even_depth(ps,sample.size=3000,replace=F)

7.1 Exercise 

Now it’s your turn. Using the commands from the previous section, please perform the following tasks (reminder: you can always use the help/documentation in R using the ?, see 6.4):

Collapse the data on genus level annotation (command: tax_glom) 
Make a biplot highlighting the 10 most abundant genera and color using the ‘diagnosis’ column (commands: ordinate, plot_ordination, etc.). Don’t forget to load the ggplot2 package.
perform beta-diversity MANOVA using the ‘diagnosis’ column (command: adonis; package: vegan); BONUS: can you add ‘sex’ and ‘age’ as covariates before testing for ‘diagnosis’ effects?
plot shannon diversity of the ‘diagnosis’ groups (command: plot_richness)

Please try if you can already do this yourself. If you need help, the commands needed are on the next page.

Please (also when using the commands on the next page), summarize briefly the results of the steps performed!















7.2 Exercise solutions

Collapse the data on genus level annotation (command: tax_glom) 

R> ps.gen<-tax_glom(ps.even, taxrank="Genus", NArm=F)

Make a biplot highlighting the 10 most abundant genera and color using the ‘diagnosis’ column (commands: ordinate, plot_ordination, etc.)

R> library(ggplot2)
R> gen.sum <- tapply(taxa_sums(ps.gen), tax_table(ps.gen)[, "Genus"], 
sum, na.rm=TRUE)
R> top10gen = names(sort(gen.sum, TRUE))[1:10]
R> gen.ord_bc <- ordinate(ps.gen, "MDS", "bray")
R> ord_df.gen <- plot_ordination(ps.gen, gen.ord_bc, type="biplot",
color="diagnosis", title="Genus Bray-Curtis",justDF=T)
R> ord_df.gen.sub <-  ord_df.gen[ord_df.gen$Genus %in% top10gen | 
ord_df.gen$id.type=="Samples",]
	R> ggplot(ord_df.gen.sub,aes(x=Axis.1,y=Axis.2, col=diagnosis, 
shape=id.type))+geom_point(size=4)+geom_label(aes(label=Genus))

perform beta-diversity MANOVA using the ‘diagnosis’ column (command: adonis); BONUS: can you add ‘sex’ and ‘age’ as covariates before testing for ‘diagnosis’ effects?

R> adonis(distance(ps.gen, method="bray") ~ diagnosis,
data=data.frame(sample_data(ps.gen)), permutations=10000)
R> adonis(distance(ps.gen, method="bray") ~ sex + age + diagnosis,
data=data.frame(sample_data(ps.gen)), permutations=10000)

plot shannon diversity of the ‘diagnosis’ groups (command: plot_richness)

R> plot_richness(ps.gen, x="diagnosis",col="sex", 
measures=c("Shannon"))


7.3 Model based differential abundance analysis

In analogy to the bear dataset, we want to explore, which taxa are differentially abundant between the CD patients. We again first define a “core microbiota”, with a mean abundance of at least 1% (= 30 reads) in the dataset:

R> (ps.gen.sub<-prune_taxa(rowMeans(otu_table(ps.gen))>=30,ps.gen))


This leaves us (again) with 16 taxa we can test. Let us have a look a the distribution of the taxa first:
R> par(mfrow=c(4,4));apply(otu_table(ps.gen.sub),1,hist,nclass=10)

Again, we see the same pattern as we saw before in the bear dataset: usually a highly skewed distribution with a lot of zeroes or low abundances. However, now that we have more samples, it is possible to use a more sophisticated model-based approach. These approaches try to (as the name says) model the data regarding certain properties, to “split” the variation into variation that can be explained by explanatory variables and covariates, and the remaining (residual) variation, which is due to factors not included in the model, either because they were not assessed while sampling or because of random/natural noise in the data.

The ‘easiest’ of all models is probably the ‘linear model’. Basic linear regression models try to fit the data to the assumption, that the observed or dependent variable Y (e.g. our taxon under investigation) can be modelled by a linear combination of a set of explanatory variables X and the error ɛ as:

Y = Xβ + ɛ

In the easiest case, this just constant or a simple line with y = β0 + β1x1 ; however in general these models can also become very complex. The problem with count data or abundance data, is that linear regression assumes normally distributed error terms (which is rarely the case for such data) and it models also possible negative values (which does evener happen for such data). One thing that is then normally done with the data, is that it is transformed by log(x), log(x+1) or square-root transformation. What is also seen often in microbiome analyses, is the arcsine-square-root transformation, or different kinds of power-transformations (Box-Cox-Transformation, Yeo-Johnson-Transformation).
All these are trying to 1.) establish a linear relationship between predictor and dependent variables and/or 2.) to make the data more “normal” distributed, however this is definitely not something you should do, as the errors of the resulting distribution are likely to NOT following a normal distribution, violating the assumptions of linear regeression (https://doi.org/10.1111/j.2041-210X.2010.00021.x), and additionally, there are models that were introduced a long time ago to deal with exactly these properties: count data with skewed distributions and non-normal error terms. These models are called generalized linear models (GLMs), as they offer a general approach to data modelling, not restricted to normally distributed errors.

Before we proceed, let us get our OTU table into a form that is more easily accessible:

R> otu<-t(otu_table(ps.gen.sub))
R> colnames(otu)<-apply(tax_table(ps.gen.sub),1,
function(x) x[min(which(is.na(x)))-1])

And have a look at the distribution of the genus Faecalibacterium:

R> par(mfrow=c(1,1))
R> hist(otu[,"Faecalibacterium"],200)


What we see, is that the number of samples with zero counts is high, meaning this taxon shows “zero-inflation” in its distribution. In the model, we will exclude all zero values, as this might influence the model fitting:

R> library(MASS)
R> faec.abu<-as.numeric(otu[,"Faecalibacterium"])
R> (mod <- glm.nb(faec.abu ~ sex + age + diagnosis, 
data=data.frame(sample_data(ps.gen.sub)),
subset=faec.abu>0))

Let’s first have a look at what we really did: first thing is, we imported the MASS package, which provides the needed functions.  Then we fit our generalized linear model by using the glm.nb function. The ‘nb’ stands for ‘negative binomial’, which is a distribution, that fits skewed count data very well (more about common distributions here). The Faecalibacterium abundance is used as dependent variable (left of the ~ symbol), and sex, age and diagnosis are the predictors. The subset option at the end specifies to only includes samples with non-zero values. As a formula, our model looks like this:

AFaec=0 + 1*sex + 2*age + 3*diagnosis + 

The function automatically converts ‘sex’ and ‘diagnosis’ to dummy variables encoding for e.g. male=1 and female=0 for the ‘sex’ variable. Same is true for diagnosis (CD=0, no=1). The output looks (something) like this:

Call:  glm.nb(formula = faec.abu ~ sex + age + diagnosis, data = data.frame(sample_data(ps.gen.sub)), 
    subset = faec.abu > 0, init.theta = 0.9045866223, link = log)

Coefficients:
(Intercept)      sexmale          age  diagnosisno  
   5.614786     0.048041    -0.005399     0.090790  

Degrees of Freedom: 285 Total (i.e. Null);  282 Residual
Null Deviance:	    333.7 
Residual Deviance: 333 	AIC: 3788

The “coefficients” are the values that were fitted for the βs in our model formula, so to estimate the abundance of Faecalibacterium, we can put these values into the formula:

AFaec=5.614 + 0.048*sex -0.005*age + 0.091*diagnosis

So this means, a 20 year old (age=20) female person (sex=0) without CD (diagnosis= 1) has an estimated Feacalibacterium abundance of:

AFaec=5.614 + 0.048*0 -0.005*20 + 0.091*1= 5.614 +0-0.1 +0.091 = 5.605

As the glm uses log-transformation internally, this translates toe5.605=271.8 counts.

The same performed for a 40 year old (age=4) male (sex=1) with CD results in:

AFaec=5.614 + 0.048*1 -0.005*40 + 0.091*0= 5.462; e5.462=235.6

We see, there is a small difference between the values, however now we want to assess, if there is a statistically significant influence of our explanatory variables. For this we use the summary.glm() function. This function uses the values estimated for each individual and the real, observed values and calculates the error between these. From this standard error, using Student’s t-distribution, the function calculatest, if these factors have a significant influence: 

R> summary.glm(mod)

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  5.614786   0.403081  13.930   <2e-16 ***
sexmale      0.048041   0.114391   0.420    0.675    
age         -0.005399   0.028245  -0.191    0.849    
diagnosisno  0.090790   0.114417   0.794    0.428    

As we see, only the Intercept term is (highly) significant, indicating, that a sample where all parameters are 0 (age, sex, diagnosis) is highly likely to have a Faecalibacterium abundance that deviates from 0. 

The following command is a shortcut, that iterates over all genera under investigation and puts out the p-values for the estimated coefficients:

R> outtab <- data.frame(t(apply(otu, 2, function(x){abu<-x; 
summary.glm(glm.nb(abu ~  sex + age + diagnosis, 
data = data.frame(sample_data(ps.gen.sub)),
subset=abu>0))$coefficients[,4]})))
R> outtab$p.adj<-p.adjust(outtab[,4])
R> outtab[order(outtab$p.adj),]

A look through the columns of this table which stores the p-values of the used explanatory variables in columns quickly tells us that the influence of sex and age are rather small for all genera (only few p values < 0.05).  However, there seems to be a strong influence of disease on several of the genera. This will again be easier to see when we visualize it:

R> par(mfrow=c(2,5))
	R> for(tax in rownames(outtab[outtab$p.adj<0.05,])){
boxplot(as.numeric(otu[,tax])~
sample_data(ps.gen.sub)$diagnosis,
main=paste(tax, "- p=",signif(outtab[tax,"p.adj"])))
}




7.4 LEfSe

The last tool we will have a look on is LEfSe, which is short for LDA Effect Size. 

Image taken from the LEfSe Galaxy Website

LEfSe works by using data in categories and subcategories and a multi-step approach. In a first step, LEfSe excludes all features (OTUs, taxa, etc.) that do not show a significant difference (p<0.05) in a Kruskal-Wallis-Test explained by the category variable, followed by pairwise comparisons of the sub-categories using a Wilcoxon Test. All features surviving these analyses will be submitted to a Linear Discriminant Analysis (LDA) to estimate effect size, which is then plotted. 

To use LEfSe we need to reformat our data a little:

R> otu2<-otu
R> colnames(otu2)<-apply(tax_table(ps.gen.sub),1,
function(x) {y=x[is.na(x)==F]; paste(y,collapse="|")})
R> lefse.in<-data.frame(t(otu2))
R> lefse.in["id",]<-colnames(lefse.in)
R> lefse.in["diagnosis",] <- data.frame(sample_data(ps.gen.sub))[,"diagnosis"]
R> write.table(lefse.in,file="lefse_in.txt",sep="\t",col.names=F,quote=F)


Now we can navigate to the online server running LEfSe as a Galaxy module (https://huttenhower.sph.harvard.edu/galaxy/). Upload your data (right panel) and run the analysis and visualization steps in LEfSe. 

LEfSe is nice for visualization and data exploration, however, it does not give p-values, is quite a “black-box” and does not inform well if it does a correction for multiple testing, which you should keep in mind!



8. Conclusion

This is the end of the workshop tutorial which should give you a starting point for the analysis of your own microbiome survey data from 16S rRNA gene sequencing. However, this is by far not a complete or exhaustive manual how to do it! 

Here are some more things to consider:
Every dataset is different: don’t blindly rely on an analysis someone else did to be 100% fitting for your own data
Know your data: problems might (will) arise while analysing, leaving you troubled, if you don’t know where your data is coming from and what happened to it before sequencing etc. So also make yourself familiar with how the samples are processed and how the sequencing technology works, so you are prepared.
There are many, many more methods. We did not even touch on basic techniques like hierarchical clustering, (partially) constrained ordinations, and other really key methods in ecology. I (strongly) encourage you to read more about this! Reading suggestions: Numerical Ecology (a lot of theory) and Numerical Ecology in R (something like the applied working book).
Microbiome 16S data is NOT quantitative!!! Keep in mind that all you are working with is actually compositional data, meaning, it all adds up to one complete community in the end, but you do not know, how large the community actually is. For further reading about the concept of compositionality and specific analysis strategy tailor to the whole concept, I suggest: Analyzing Compositional Data with R, including also the published scientific articles by the authors of the book.
The field is moving fast: new methods are published almost daily, getting more and more specific and complex. Whole journals are dedicated to the analysis of microbial communities. So you are likely to find help about how to handle your own analysis and new perspectives from some of them.
The microbiome is not only bacteria. There are archaea, fungi, viruses/phages, small/unicellular eukaryotes. Don’t forget about these fellas!


