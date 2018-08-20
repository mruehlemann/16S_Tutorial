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

qiime@qiime-190-virtual-box:~$ 

and a ‚flashing’ box. This is where you can enter your commands. Every line with a command will now be written starting simply with a dollar sign $. For example:
```bash
$ pwd
```

Please do not include the $ when copying the command. When you type this command, the terminal will put out the folder you are currently working in (pwd is short for ‚print work directory’). You can learn more about commands if you enter man followed by the command:
```
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
```
$ cd /home/qiime/Desktop
$ ls
```
We can see that in this folder there is another folder called raw_data. We now change into this directory, and list the first 10 files (head –n 10):
```
$ cd raw_data
$ ls | head -n 10
```
Between the commands ls and head, we see a vertical bar. This bar is called the pipe and is used to pass the output of one command (also called stdout) to another command, in this case the head command. The head command shows only a specified number of lines of the input. The number of lines to be displayed can be adjusted using the -n parameter.

What you see as output is the typical data you will get from your sequencing facility. Let us dissect the filename of the first file to make sense of it:

B-1104-S_S70_L001_R1_001.fastq.gz

B-1104-S:	This is the sample name and can be anything you want it to be. In this case this name is composed of ‘B’ for ‘Bear’, ‘1104’ is the ID of the bear, and ‘S’ which means that the sample was taken in summer. You will also see “B-1104-W” a few lines further down, which is the data of the same individual but taken in winter.

S70:		This part of the filename simply tells you that it is sample 70 of this specific sequencing run. Among other things to expect at this position of the filename are for example the index-sequences used in multiplexing.

L001:	This stands for ‘Lane 1’. When working with MiSeq data, this will always be Lane 1, as a MiSeq sequencer only has a single lane. When working with data from e.g. a HiSeq sequencer, this might be something between L001 and L008.

R1:		This tells you that this file contains the forward read (Read 1) of the sequencing data. When you look at the second file in the list, you will see that the filename is exactly the same apart from this position, where it says ‘R2’. This means that these two files belong together and contain the read pairs of a single sample.

001:		This is always 001. Don’t ask why.

fastq.gz:	This is the filename ending, indicating the file format. fastq indicates that the file is in FastQ format, the gz at the end means, that the file is compressed using the Gzip compression. More about this later.

Okay, now that we know what all those filename mean, let’s have a look at what is in the file. For this we need to first decompress the file (gunzip), pass the output to stdout (-c), pipe it to the next command and display the first four lines in the file:

$ gunzip -c B-1104-S_S70_L001_R1_001.fastq.gz | head -n 4





What you get will look something like this:

@M02765:9:000000000-A8H6D:1:1101:15136:1407 1:N:0:70
TACAGAGGATGCAAGCGTTATCCGGAATGATTGGGCGTAAAGCGTCTGTAGG
+
ABBBBFFBFFFFGGGGFGGEGFHGGEFFHGGHHGHHGGGGHHHGGGEGHHHH

This output is the standard format of the FastQ sequencing data format. One sequence in FastQ format always consists of four lines:

Line 1: 	This line is the start of the entry and always starts with the @ symbol, followed by a unique sequence ID. This ID is the same in both read files (R1+R2) and defines which sequences belong together. When looking at Illumina sequencing data, this ID consists of the unique instrument ID (M02765), the run ID (9), the flowcell ID (A8H6D) and the lane (1), followed by coordinates of the sequence cluster on the flowcell. After the space there are more informations on the sequence, however these do not belong to the unique sequence ID.

Line 2:	This is the actual sequence obtained from basecalling of the raw data. In our case this is the 16S amplicon sequence.

Line 3:	this line can contain additional information about the sequence, however usually only contains a mandatory ‚+’ symbol. 

Line 4:	In this line the quality of every single base in the sequence is encoded. The first position of this line corresponds to the first position of the sequence, the second to the second, and so on. Each ASCII-character endcodes a special value between 0 and 41 (in the most widely used Illumina 1.8 Phred+33 format; +33 stands for an offset of 33 in the ASCII characters).  The highest Phred score a base can get is in this case 41, encoded by the letter J (ASCII: 74), the lowest is 0, encoded by ! (ASCII: 33). The Phred score (Q) is defined as the probability (P), that the base at this position is incorrect: P(Q)=10-Q10. This corresponds to 0.01% error probability at Q=40 and 10% at Q=10.

Phred+33 encoding of quality scores:

Symbol	 !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ  
 |                                        |                 
ASCII	33                                       73
Q Score	 0.2......................26...31........41                             


Of course it is not possible to check each single sequence by ‘hand’, however there are a lot of very handy tools to check if the sequencing in general yielded satisfactory data. One tool we want to introduce is the very easy to use and freely available FastQC. Please go the website (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) using the browser in the VirtualBox, download the tool and unpack it. For this, please create a folder on the Desktop called ‘software’, double-click on the downloaded archive and move the ‘FastQC’ folder into the newly created ‘software’ folder. Go to the terminal again and change to the FastQC software folder:

$ cd /home/qiime/Desktop/software/FastQC

To use fastqc, we first have to make the binary file executable. For this type this into the terminal:

$ chmod +x fastqc

The command chmod can be used to change the mode or access permission of a file. Here we add permission to execute (+x) to the fastqc binary. No we can start FastQC:

$ ./fastqc

This should open the FastQC software. In the Menu bar go to: “File > Open …” and navigate to the raw_data folder in Shared_Folder. By holding the Shift key you can select multiple files at once. Please select the R1 and R2 read of the Sample “B-1105-W” and click the “OK” button. You will see the files being processed and quickly the reports should open, each as a tab in the FastQC window. On the left side you will see several green, orange and red symbols. These are the different quality criteria assessed by FastQC. However, let’s start with the first one called “Basic statistic”. Here you can already see some basic properties of the opened files, for example that each of the files holds 77868 sequences (as they are the paired reads of the same sample this should be the same; if it is not, there is something wrong), the sequence lengths are between 248 and 251 nt and the GC content differs sligthly between R1 and R2, however both being a little above 50%.

Next please have a look at the “Per base sequence quality” tab for both files. What you see here is a pretty typical (though not perfect) example of how amplicon sequencing data normally looks like:

the R1 is relatively good quality all along, only dropping a little at the end of the read, however mostly staying in the green area until ~180-200 nt.
the R2 is usually of worse quality than the R1. This is  is basically what you can expect and mainly due to technical reasons, so don’t worry.  Only if you see really strange patterns of very drastic drops in quality (also maybe only at a single position across all samples) things might have to be investigated more deeply.

Move on to the tab called “Per base sequence content”. This tab will likely be marked with a red cross, indicating a problem. When you look at the image in this tab, you can see the per-position composition of the four nucleotides and notice, that there is only very little variation in the first 8-9 nt. This is just as it is supposed to be, as the amplification primers have to bind in a non- or at least less-variable region of the 16S gene to work. As FastQC is intended to be used for any kind of sequencing data, it assumes that this is an error, however in our case, again, this is what we are expecting.

The same rules apply for GC content (of course there is a bias in GC content when only looking at a specific sequence) and duplicated/overrepresented sequences: we expect to have many duplicated or similar sequences, as this is exactly the purpose of using an amplicon.

In general we can say: our data look just like we expect it to look!

One more very convenient thing about FastQC is, that you can also use it for batch processing of FastQ file. Close the FastQC window and type:

$ wosho=/home/qiime/Desktop/workshop
$ echo $wosho
$ raw=$wosho/raw_data
$ echo $raw
$ ./fastqc $raw/B-1105-W_* --outdir /home/qiime/Desktop

In the first line, we set a variable named “wosho”. This is very easy in the bash environment and might look a little different in other programming languages. In general, many things can be stored in variables, however, here we use it for a path to a folder. In the third line, we set another variable “raw”, which has the same content as “wosho”, plus “/raw_data” at the end. As you can see, the content of the variables can be show by using “echo” followed by the variable with a dollar sign prepended. The last line starts the fastqc binary, applied to all files in the folder that is stored in the “raw” variable that start with “B-1105-W_”. The asterisk “*” is a so-called “wildcard” and matches every character when searching for filenames.

Now look at your desktop. There should be two html files you can open with your browser. These file hold all the information we saw before in the FastQC window and this way of processing can be easily applied to a whole batch of sequencing files for fast processing to get an idea if your data looks usable.

Another software that is very nice for getting an overview of sequencing data in general is MultiQC (http://multiqc.info), which is able to aggregate the output of different QC tools to create reports. However, this is currently more aiming at RNAseq or RRBS experiments, but worth keeping an eye on. 
