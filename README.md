# EE282 Homework4 (Due December 10) by Jihye Choi(94323474)
## <Summarize partitions of a genome assembly>  
  
First of all, we change the current directory to **Homework4** and use **wget** to download the data from the given website _flybase.org_. Then we verify the file integrity and unzip the file.  

    $ cd Homework4
    $ wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz    
    $ md5sum dmel-all-chromosome-r6.24.fasta.gz  
    $ gunzip dmel-all-chromosome-r6.24.fasta.gz 
  
## Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:
### 1.Total number of nucleotides
### 2.Total number of Ns
### 3.Total number of sequences
### Because the calculations will be for two genome partitions, there will be 6 total responses.




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
