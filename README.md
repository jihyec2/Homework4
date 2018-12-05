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

### 1. Sequence length distribution
### 2. Sequence GC% distribution
### 3. Cumulative genome size sorted from largest to smallest sequences


### Because the calculations will be for the whole genome and two genome partitions, there will be 9 total plots.




## <Genome assembly>

### Note: This part of homework 4 is still being arranged. When this note is gone, it should be ready.

## Assemble a genome from MinION reads

### Hint: Read up on miniasm here. We're using one of the simplest assembly approaches possible. This assembly can literally be accomplished with three lines of code. This will literally take only 3 command lines.

### 1.Download the reads from here
### 2.Use minimap to overlap reads
### 3.Use miniasm to construct an assembly


## Assembly assessment

### Hint: For MUMmer, you should run nucmer, delta-filter, and mummerplot.

### 1.Calculate the N50 of your assembly (this can be done with only faSize+awk+sort or with bioawk+awk+sort) and compare it to the Drosophila community reference's contig N50 (here)
### 2.Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase using a dotplot constructed with MUMmer (Hint: use faSplitByN as demonstrated in class)
### 3.Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot (Hint: use plotCDF2 as demonstrated in class and see this example)
### 4.Calculate BUSCO scores of both assemblies and compare them
