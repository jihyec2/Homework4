# EE282 Homework4 (Due December 10) by Jihye Choi(94323474)
## <Summarize partitions of a genome assembly>  
### We will be revisiting the Drosophila melanogaster genome. As with Homework 3, start at flybase.org. Go to the most current download genomes section and download the gzipped fasta file for all chromosomes.
  
### Hint: The partitioning can be accomplished in many different ways. In my opinion, the easiest way is by using bioawk and faSize. The bioawk tool can be found in the module jje/jjeutils and the fa* utilities can be found in the module jje/kent.

## Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:

First of all, we change the current directory to **Homework4** and use **wget** to download the data from the given website _flybase.org_. Then we verify the file integrity and unzip the file.  

    $ cd Homework4
    $ wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz    
    $ md5sum dmel-all-chromosome-r6.24.fasta.gz  
    $ gunzip dmel-all-chromosome-r6.24.fasta.gz 

Next, we load the modules for the **bioawk** and **faSize** tools, and use **bioawk** to complete the partitioning and to save the files with new names. Then we use **faSize** to calculate the followings.   

    $ module load jje/jjeutils jje/kent #load modules 
    $ bioawk -c fastx 'length($seq)  <= 100000{ print ">"$name; print $seq }' dmel-all-chromosome-r6.24.fasta | sort -rn > dmel_fasta_leq100kb.fasta #save it as a new file after partitioning with less than or equal to 100kb
    $ bioawk -c fastx 'length($seq)  > 100000{ print ">"$name; print $seq }' dmel-all-chromosome-r6.24.fasta | sort -rn > dmel_fasta_gre100kb.fasta #save it as a new file after partitioning with greater than 100kb
    $ faSize dmel_fasta_leq100kb.fasta #for the sequences less than or equal to 100kb
    
        6178042 bases (662593 N's 5515449 real 5515449 upper 0 lower) in 1863 sequences in 1 files
        Total size: mean 3316.2 sd 108053.8 min 0 (211000022278031) max 4245830 (mitochondrion_genome) median 0
        N count: mean 355.7 sd 11351.2
        U count: mean 2960.5 sd 96726.1
        L count: mean 0.0 sd 0.0
        %0.00 masked total, %0.00 masked real
        
    $ faSize dmel_fasta_gre100kb.fasta # for the sequences greater than 100kb
    
        137547960 bases (490385 N's 137057575 real 137057575 upper 0 lower) in 7 sequences in 1 files 
        Total size: mean 19649708.6 sd 51988242.2 min 0 (2L) max 137547960 (X) median 0
        N count: mean 70055.0 sd 185348.1
        U count: mean 19579653.6 sd 51802894.1
        L count: mean 0.0 sd 0.0
        %0.00 masked total, %0.00 masked real
        
### For all sequences ≤ 100kb:
1.	Total number of nucleotides : 6178042 nucleotides
2.	Total number of Ns : 662593 N's
3.	Total number of sequences : 1863 sequences
### For all sequences > 100kb:
1.	Total number of nucleotides : 137547960 nucleotides
2.	Total number of Ns : 490385 N's
3.	Total number of sequences : 7 sequences


## Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:

### Hint: bioawk has a function called gc(). Don't forget about the CDF plotting utility we used in class.
First of all, we begin this problem by loading modules. 

    $ module load perl
    $ module load jje/jjeutils
    $ module load rstudio/0.99.9.9

### For the whole genome:
1. Sequence length distribution

        $ bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > dmel_all_seq.length
        $ plotCDF2 dmel_all_seq.length all_seq.png  #plot by using CDF plotting utility


2. Sequence GC% distribution

        $ bioawk -c fastx '{ print $name, gc($seq) } ' dmel-all-chromosome-r6.24.fasta > GC_all.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        #On RStudio script, we write the following code to get a plot for the sequence GC distribution
        rm(list = ls())
        library(ggplot2) 
        library(dplyr)
        library(splitstackshape)

        setwd("~/Homework4") #set working directory
        dmel_GC.df <- read.table("~/Homework4/dmel_fasta_GC", header=FALSE, sep = "")
        colnames(dmel_GC.df) <-c("genomicLocation", "length")


        GC_plot <- ggplot(dmel_GC.df, aes(x=percentGC)) + geom_histogram(fill="deepskyblue4") + 
          ggtitle("Histogram: Sequences by GC content") +
          labs(x="GC content (%) ", y="number of sequences") +
          theme(plot.title = element_text(family = "Trebuchet MS", color="#000066", face="bold", size=24)) +
          theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16)) 

3. Cumulative genome size sorted from largest to smallest sequences

### For all sequences ≤ 100kb:
1. Sequence length distribution    

        $ bioawk -c fastx ' { print length($seq) } ' dmel_fasta_leq100kb.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > dmel_leq_seq.length
        $ plotCDF2 dmel_leq_seq.length leq_seq.png #plot by using CDF plotting utility
        $ ls *.png #list all png files to check if the plot named 'les_seq' exists   

2. Sequence GC% distribution
3. Cumulative genome size sorted from largest to smallest sequences

### For all sequences > 100kb:
1. Sequence length distribution

        $ bioawk -c fastx ' { print length($seq) } ' dmel_fasta_gre100kb.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > dmel_gre_seq.length
        $ plotCDF2 dmel_gre_seq.length gre_seq.png #plot by using CDF plotting utility  
        $ ls *.png #list all png files to check if the plot named 'gre_seq' exists 


2. Sequence GC% distribution
3. Cumulative genome size sorted from largest to smallest sequences


## <Genome assembly>
  
## Assemble a genome from MinION reads

### Hint: Read up on miniasm here. We're using one of the simplest assembly approaches possible. This assembly can literally be accomplished with three lines of code. This will literally take only 3 command lines.

### 1.Download the reads from here
    $ module load jje/jjeutils
    $ module load perl
    $ cd Homework4 #change directory 
    $ mkdir GNassembly #Make directory 
    $ cd GNassembly #change directory to GNassembly
    $ wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz #Download the reads
    $ gunzip iso1_onp_a2_1kb.fastq.gz
    $ ln -sf iso1_onp_a2_1kb.fastq reads.fq

### 2.Use minimap to overlap reads
    #I relogged into hpc and qrsh into a 32 core node
    $ qrsh -q epyc,abio128,free88i,free72i -pe openmp 32
    $ cd /data/users/jihyec2/Homework4/GNassembly
    $ minimap -t 32 -Sw5 -L100 -m0 reads.fq{,} | gzip -1 > onp.paf.gz

### 3.Use miniasm to construct an assembly
    $ miniasm -f reads.fq onp.paf.gz > reads.gfa
    
## Assembly assessment

### Hint: For MUMmer, you should run nucmer, delta-filter, and mummerplot.

### 1.Calculate the N50 of your assembly (this can be done with only faSize+awk+sort or with bioawk+awk+sort) and compare it to the Drosophila community reference's contig N50 (here)

    $ n50 () {
      bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
      | sort -rn \
      | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
      } 
      
    $ awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' reads.gfa \
      | tee >(n50 /dev/stdin > n50.txt) \
      | fold -w 60 \
      > unitigs.fa
      
    $ N50 dmell-contig.fa
      

### 2.Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase using a dotplot constructed with MUMmer (Hint: use faSplitByN as demonstrated in class)
      $ faSplitByN dmel-all-chromosome-r6.24.fasta dmel-contig.fasta #making contig assembly
      $

# Need to first make a contig assembly

module load jje/jjeutils perl
faSplitByN dmel-all-chromosome-r6.24.fasta dmel-all-chromosome-cntg-r6.24.fasta 10

# Will do mummer in another folder as a job
mkdir mummer
ln -s /pub/jje/ee282/bsorouri/nanopore_assembly1/nanopore_assembly1/data/processed/unitigs.fa
ln -s /pub/jje/ee282/bsorouir/hmwk4/dmell-all-chromosome-cntg-r6.24.fasta
ls
touch mummer.sh
nano mummer.sh # Copy and paste the content below into your shell script, afterwards save and exit out of shell script

#!/bin/bash
#
#$ -N mummer
#$ -q free128,free72i,free56i,free48i,free40i,free32i,free64
#$ -pe openmp 8
#$ -R Y

    ###Loading of binaries via module load or PATH reassignment
    source /pub/jje/ee282/bin/.qmbashrc
    module load gnuplot

    ###Query and Reference Assignment. State my prefix for output filenames
    REF="dmel-contig.fasta"
    PREFIX="flybase"
    SGE_TASK_ID=1
    QRY=$(ls u*.fa | head -n $SGE_TASK_ID | tail -n 1)
    PREFIX=${PREFIX}_$(basename ${QRY} .fa)

    ###please use a value between 75-150 for -c. The value of 1000 is too strict.
    nucmer -l 100 -c 125 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
    mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta \
      -R ${REF} -Q ${QRY} --postscript


#### This is after you saved and exited out ######

  qsub mummer.sh

### 3.Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot (Hint: use plotCDF2 as demonstrated in class and see this example)
    bioawk -c fastx ' { print length($seq) } ' dmel-contig.fasta \
    | sort -rn \
    | awk ' BEGIN { print "Assembly\tLength\nFB\t0" } { print "FB\t" $1 } ' \
    >  dmel-contig.length

    bioawk -c fastx ' { print length($seq) } ' unitigs.fa \
    | sort -rn \
    | awk ' BEGIN { print "Assembly\tLength\nMinimap_Ctg\t0" } { print "Minimap_Ctg\t" $1 } ' \
    > unitigs.length

    plotCDF2 {dmel-contig,unitigs}.length assembly.png

### 4.Calculate BUSCO scores of both assemblies and compare them
pwd # make sure you are in hmwk4 directory
touch busco_final8.sh
nano busco_final8.sh ##### Input and save the code below

#!/bin/bash
#
#$ -N busco8
#$ -q free128,free72i,free56i,free48i,free40i,free32i,free64
#$ -pe openmp 8
#$ -R Y 

    module load augustus/3.2.1
    module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
    source /pub/jje/ee282/bin/.buscorc

    INPUTTYPE="geno"
    MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
    MYLIB="diptera_odb9"
    OPTIONS="-l ${MYLIBDIR}${MYLIB}"
    ##OPTIONS="${OPTIONS} -sp 4577"
    QRY="unitigs.fa"
    ###Please change this based on your qry file. I.e. .fasta or .fa or .gfa
    MYEXT=".fa" 

    #my busco run
    #you can change the value after -c to tell busco how many cores to run on. Here we are using only 1 core.
    BUSCO.py -c 1 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}      


################# After Saving And Exiting Out #####################

qsub busco_final8.sh
# Note: I used 8 nodes, because it was taking forever in the queue. It took longer to complete, but it still got done.

#### Alternative, manual way of running busco:
# Need to be careful how we log on and off, but can do the above script manually by first putting the code below, then continuing with the remainder of the code:
qrsh -q free128,free72i,free56i,free48i,free40i,free32i,free64 -pe openmp 32


Output:

C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799

13 Complete BUSCOs (C)

13 Complete and single-copy BUSCOs (S)

0 Complete and duplicated BUSCOs (D)

32 Fragmented BUSCOs (F)

2754 Missing BUSCOs (M)

2799 Total BUSCO groups searched

