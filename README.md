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

        $ bioawk -c fastx '{ print length($seq) }' dmel-all-chromosome-r6.24.fasta > seqlen_all.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        
        #On RStudio script, we write the following code to get a plot for the sequence length distribution
        
        library(ggplot2)
        all_len <- read.table("seqlen_all.txt", header = FALSE)
        all_len$Percentcut <-cut(x=all_len[,1], breaks = 10)
        a <- ggplot(data = all_len)+ geom_bar(mapping = aes(Percentcut))
        a + labs(title="Sequence Length (Whole Genome)", x="Length", y="Count") 
        
![1](https://blogfiles.pstatic.net/MjAxODEyMTBfMjg2/MDAxNTQ0NDM4MDc5MDU1.vWRQLxcdHCUj4krI574Hw1L2CLIH66ARY5ipcwuJLPcg.zdya5r5EJMIWdJDoFWQd0m--i_pJxmDDgGmaqqmtpt4g.JPEG.nayeonkim93/1.jpeg)

2. Sequence GC% distribution

        $ bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta > GC_all.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        
        #On RStudio script, we write the following code to get a plot for the sequence GC distribution
        
        library(ggplot2)
        all_GC <- read.table("GC_all.txt", header = FALSE)
        all_GC$Percentcut <-cut(x=all_GC[,2], breaks = 20)
        a <- ggplot(data = all_GC)+ geom_bar(mapping = aes(Percentcut))
        a + labs(title="GC Distribution (Whole Genome)", x="Percentage", y="Count") 
 
![2](https://blogfiles.pstatic.net/MjAxODEyMTBfMjc5/MDAxNTQ0NDM4MDU0MDQ1.f-2BX4esbz75UIT5ihd6-Gae_3Nv110wWVOJYm4HvtMg.pDaX9GOEKJLz5kJH74GkpSpj-aV4ZtSOS0eW2moqwRQg.JPEG.nayeonkim93/4.jpeg) 

3. Cumulative genome size sorted from largest to smallest sequences

        $ bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > dmel_all_seq.length
        $ plotCDF2 dmel_all_seq.length all_seq.png  #plot by using CDF plotting utility

![all](https://blogfiles.pstatic.net/MjAxODEyMDZfMjMx/MDAxNTQ0MDk2NDA5NjIx.dYbVr_AdfEXVGTCAO-t4T8N7p1Zjqttdrzm9SptmL6Ig.-BELOSgyJzbYZtHz7KGJvXxK93OADlm_55ro2tnSTH0g.PNG.nayeonkim93/length-whole-genome.png)

### For all sequences ≤ 100kb:
1. Sequence length distribution    

        $ bioawk -c fastx '{ print length($seq) }' dmel_fasta_leq100kb.fasta > seqlen_leq.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        
        #On RStudio script, we write the following code to get a plot for the sequence length distribution
        
        library(ggplot2)
        leq_len <- read.table("seqlen_leq.txt", header = FALSE)
        leq_len$Percentcut <-cut(x=leq_len[,1], breaks = 10)
        a <- ggplot(data = leq_len)+ geom_bar(mapping = aes(Percentcut))
        a + labs(title="Sequence Length (<= 100kb)", x="Length", y="Count") 

![3](https://blogfiles.pstatic.net/MjAxODEyMTBfMjUy/MDAxNTQ0NDM4MDQ3Nzc2.7VvAw_o1USAowuvm6I_eEjuvCk3BIWAh9JzlnKJigX4g.WJ79hA3Yle31WEOBfnkdkOGuBy6XYsPpn9bnCOFm5XQg.JPEG.nayeonkim93/2.jpeg)

2. Sequence GC% distribution

        $ bioawk -c fastx '{ print $name, gc($seq) }' dmel_fasta_leq100kb.fasta  > GC_leq.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        
        #On RStudio script, we write the following code to get a plot for the sequence GC distribution
        
        library(ggplot2)
        leq_GC <- read.table("GC_leq.txt", header = FALSE)
        leq_GC$Percentcut <-cut(x=leq_GC[,2], breaks = 20)
        a <- ggplot(data = leq_GC)+ geom_bar(mapping = aes(Percentcut))
        a + labs(title="GC Distribution (<= 100kb)", x="Percentage", y="Count") 
        
![4](https://blogfiles.pstatic.net/MjAxODEyMTBfMjMy/MDAxNTQ0NDM4MDU3MTk1.sFNP4sDatd1fTNT2ZHIL_7CddRbEoQHrGb3599WN754g.7Q1PtoYY_o-is9hiwYdK9u3yzgP7g6pD3DkBrMi5kjwg.JPEG.nayeonkim93/5.jpeg)

3. Cumulative genome size sorted from largest to smallest sequences

        $ bioawk -c fastx ' { print length($seq) } ' dmel_fasta_leq100kb.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > dmel_leq_seq.length
        $ plotCDF2 dmel_leq_seq.length leq_seq.png #plot by using CDF plotting utility
        
![leq](https://blogfiles.pstatic.net/MjAxODEyMDZfMTIz/MDAxNTQ0MDk2NDE1OTc1.wgP88M7T4ml8YY9kOQf7Kdan9EAsRma3Ob_IVxyzO7kg.W6tvCuNmKciW2upqL0UeNG0tAERO_PIs-lxgjT4KTqcg.PNG.nayeonkim93/length-less-genome.png)

### For all sequences > 100kb:
1. Sequence length distribution

        $ bioawk -c fastx '{ print length($seq) }' dmel_fasta_gre100kb.fasta > seqlen_gre.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        
        #On RStudio script, we write the following code to get a plot for the sequence length distribution
        
        library(ggplot2)
        gre_len <- read.table("seqlen_gre.txt", header = FALSE)
        gre_len$Percentcut <-cut(x=gre_len[,1], breaks = 10)
        a <- ggplot(data = gre_len)+ geom_bar(mapping = aes(Percentcut))
        a + labs(title="Sequence Length (> 100kb)", x="Length", y="Count") 

![5](https://blogfiles.pstatic.net/MjAxODEyMTBfNTUg/MDAxNTQ0NDM4MDUwNDg4.nNgAkgZnqbhR7cak33_I0uaviX-xxFa3sQyU_VmQqPAg.8wcIJvG0hxkIC4cGl964qyn5ThdmDV7XYbiy4CZSbv4g.JPEG.nayeonkim93/3.jpeg)

2. Sequence GC% distribution

        $ bioawk -c fastx '{ print $name, gc($seq) }' dmel_fasta_gre100kb.fasta  > GC_gre.txt #prepping data files
        $ rstudio #open rstudio through the terminal
        
        #On RStudio script, we write the following code to get a plot for the sequence GC distribution
        
        library(ggplot2)
        gre_GC <- read.table("GC_gre.txt", header = FALSE)
        gre_GC$cut <-cut(x=gre_GC[,2], breaks = 20)
        a <- ggplot(data = gre_GC)+ geom_bar(mapping = aes(Percent_cut))
        a + labs(title="GC Distribution (> 100kb)", x="Percentage", y="Count") 

![6](https://blogfiles.pstatic.net/MjAxODEyMTBfMjcw/MDAxNTQ0NDM4MDYwMzUx.0uib88OQ9hcAK0t8zXfonPNL21VYHJWxq8GGEa4L4y0g.hT5FufL2cTMVjxXcbf8pvgcWiXD9yZkV5v54GyUBmqwg.JPEG.nayeonkim93/6.jpeg)

3. Cumulative genome size sorted from largest to smallest sequences

        $ bioawk -c fastx ' { print length($seq) } ' dmel_fasta_gre100kb.fasta | sort -rn | awk ' BEGIN { print           "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > dmel_gre_seq.length
        $ plotCDF2 dmel_gre_seq.length gre_seq.png #plot by using CDF plotting utility 

![gre](https://blogfiles.pstatic.net/MjAxODEyMDZfMTI0/MDAxNTQ0MDk2NDIyNzk5.Vo5bukEg4yEmMbUl4Nj-TbPIzpZjW8ODyWV-ontO3O0g.hTDhCAvZEppVJodRyZd9sg8o_ELKCaaKU99m4m3L3eIg.PNG.nayeonkim93/length-great-genome.png)

## <Genome assembly>
  
## Assemble a genome from MinION reads

### Hint: Read up on miniasm here. We're using one of the simplest assembly approaches possible. This assembly can literally be accomplished with three lines of code. This will literally take only 3 command lines.

### 1.Download the reads from here
    qrsh -q epyc,abio128,free88i,free72i -pe openmp 32    
    #I relogged into hpc and qrsh into a 32 core node
    cd /data/users/jihyec2/Homework4
    wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz
    gunzip *.gz
    ln -sf iso1_onp_a2_1kb.fastq reads.fq
    
### 2.Use minimap to overlap reads
    minimap -t 32 -Sw5 -L100 -m0 reads.fq{,} | gzip -1 > onp.paf.gz

### 3.Use miniasm to construct an assembly
    miniasm -f reads.fq onp.paf.gz > reads.gfa
    
## Assembly assessment

### Hint: For MUMmer, you should run nucmer, delta-filter, and mummerplot.

### 1.Calculate the N50 of your assembly (this can be done with only faSize+awk+sort or with bioawk+awk+sort) and compare it to the Drosophila community reference's contig N50 (here)
    module load jje/jjeutils perl

    $ n50 () {
      bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
      | sort -rn \
      | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
      } 
      
    $ awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' reads.gfa \
      | tee >(n50 /dev/stdin > n50.txt) \
      | fold -w 60 \
      > unitigs.fa
      
    # When 'ls', we can see the file 'n50.txt'. 
    $ less n50.txt
    
    4494246
      

### 2.Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase using a dotplot constructed with MUMmer (Hint: use faSplitByN as demonstrated in class)

    faSplitByN dmel-all-chromosome-r6.24.fasta dmel-contig.fasta 10 
    module unload jje/jjeutils perl
    
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
    nucmer -l 100 -c 150 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
    mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta -R ${REF} -Q ${QRY} --png

I uploaded the dotplot which I constructed with MUMmer on Github because it becomes blurry if I link my png file from my personal blog. 

### Hi, You should have still embedded the png here

### 3.Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot (Hint: use plotCDF2 as demonstrated in class and see this example)
    
    bioawk -c fastx '{ print length($seq) }' unitigs.fa \
    | sort -rn \
    | awk ' BEGIN { print "Assembly\tLength\nassembly\t0" } { print "assembly\t" $1 } ' \
    > mini_unitigs

    bioawk -c fastx ' { print length($seq) } ' dmel-contig.fasta \
    | sort -rn \
    | awk ' BEGIN { print "Assembly\tLength\ncontig\t0" } { print "contig\t" $1 } ' \
    >  contig

    bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
    | sort -rn \
    | awk ' BEGIN { print "Assembly\tLength\nscaffold\t0" } { print "scaffold\t" $1 } ' \
    >  scaffold

    plotCDF2 mini_unitigs contig scaffold Comparison.png

![8](https://blogfiles.pstatic.net/MjAxODEyMTBfMTAx/MDAxNTQ0NDQxNzg2NzEx.H1BustIvOVqPuzXj9CimHQrhsZBGFoGdk5thJEpPCYog.9hBoDZTXcKop1G8EyYYkHq1F7wGmSp7KZ7gvKnXX6fUg.PNG.nayeonkim93/contig_plot.png)
### 4.Calculate BUSCO scores of both assemblies and compare them

    #!/bin/bash  
    #  
    #$ -N HW4_busco  
    #$ -q free128,free72i,free56i,free48i,free40i,free32i,free64  
    #$ -pe openmp 32 
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
    BUSCO.py -c 125 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}      

Output:

C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799  
13 Complete BUSCOs (C)  
13 Complete and single-copy BUSCOs (S)  
0 Complete and duplicated BUSCOs (D)  
32 Fragmented BUSCOs (F)  
2754 Missing BUSCOs (M)  
2799 Total BUSCO groups searched

