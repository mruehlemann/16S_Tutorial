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
$ filter_fasta.py -f $out/slout_q20/B.1104.S/seqs.fna -o $out/usearch61_chimera/B.1104.S/seqs_clean.fna -s $out/usearch61_chimera/B.1104.S/chimeras.txt
-o $out/usearch61_chimera/B.1104.S -s $out/usearch61_chimera/B.1104.S/chimeras.txt
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
