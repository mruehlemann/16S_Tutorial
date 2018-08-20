# 16S Amplicon Sequencing Workshop Tutorial #

Hosted by: Malte Rühlemann & Felix Sommer

## 1. Introduction ##

This workshop is intended to give a broad overview of the processing and analysis of sequencing data from 16S rRNA gene amplicon-based experiments aiming to survey bacterial communities. We will start with command line basics and handling of sequencing data in FastQ format, to introduce what to expect when performing such experiments. After this we will move through the different data processing steps that are necessary when dealing with amplicon data. Finally, we will have a look at a final dataset and compare the outcomes of different processing possibilities, using QIIME output and also the statistical software R. The latter we will also use to discuss properties of microbiome data that need to be taken into account when applying statistical test. 

Please feel free to explore any output file in between as much as you like. We will also put this tutorial online so you can redo the analysis again and as a general resource for your own analyses. What we provide here is one or a few of the possible ways to process and analyze microbiome data, however there is not a single “right” way to do it and looking through the literature, you will easily find people doing something different from what you have been introduced to today. That is why we are being extensive in the explanations in this tutorial, so you can really understand what you are doing and why you are doing it. 

As we will be using the bash command line a lot, you can go here for a general resource about the commands used: https://ss64.com/bash/

The VirtualBox on the provided computers should already be configured and working. Please start it up by clicking the funny looking blue icon in the taskbar on the left.

The dataset we will use in the tutorial is a subset of the data from the article ‚The Gut Microbiota Modulates Energy Metabolism in the Hibernating Brown Bear Ursus arctos’ by Sommer et al. (Cell Reports, 2016, doi.org/10.1016/j.celrep.2016.01.026).

To get all the data needed open the terminal and type:
´´´
> cd Desktop
> git clone http://github.com/mruehlemann/16S_Tutorial
> 
´´´
