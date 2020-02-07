# epidome
Databases, scripts and example data pertaining to the epidome. Subspecies classification of S. epidermidis based on amplicon sequencing.

## Contents
* [What is the epidome?](#what-is-the-epidome)
* [Why study <i>S. epidermidis</i>?](#why-study-S.-epidermidis)
* [Why not whole genome sequencing?](#why-not-whole-genome-sequencing)
* [Manual](#anual)
  * [Basics](#basics)


## What is the epidome?
An approach to examining the polyclonal nature of <i>Staphylococcus epidermidis</i> bacterial communities based on targeted sequencing, similar to classical 16s rRNA gene sequencing, but using targets that are specific to <i>S. epidermidis</i> and which can be used to differentiate between clones of <i>S. epidermidis</i> on a subspecies level

## Why study S. epidermidis?
Although <i>S. epidermidis</i> is mostly found as a commensal of the human nasal and skin microbiome it is becomming increasingly clear that its role as an opportunistic pathogen should not be underestimated. <i>S. epidermidis</i> infections associated with indwelling medical devices in particular has become a concern, and as the use of medical implants increases with changing demographics the medical and economic burden of these infections are bound to increase further. <i>S. epidermidis</i> has a natural proficiency in biofilm formation that has been linked to its success in causing infections in implants, and with the emergence and global dissemination of multidrug-resistant lineages a narrowing of treatment options is also becomming a concern.

## Why not whole genome sequencing?
Multiple excellent epidimiological and genomic studies based on WGS data have been performed trying to identify the factors necessary for <i>S. epidermidis</i> to be succesfully invasive. These are mostly focused on comparing isolates from infections to those found as commensals in nasal and skin swabs, and have helped underline the key role of antibiotic resistance in infections. While most clinical infections are likely monoclonal as a result of extensive prophylactic treatment it has become evident that commensal <i>S. epidermidis</i> communities are almost always polyclonal. Examining single isolates will therefore provide a highly incomplete picture of the diversity found in commensal <i>S. epidermidis</i> communities.


## Manual

### Basics

### Key functions

#### Load dada2 output and metadata into R, make sure rownames in metadata match names of isolates
epi01_table = read.table("epi01_dada_output.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("epi02_dada_output.csv",sep = ";",header=TRUE,row.names=1)
metadata_table = read.table("metadata_table.txt")

#### Setup object for easy handling of data
epidome_object = setup_epidome_object(epi01_table,epi02_table,metadata_table)




## Read stuff

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2807625/
https://www.frontiersin.org/articles/10.3389/fcimb.2017.00090/full
