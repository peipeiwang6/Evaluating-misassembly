# Evaluating-misassembly

teaser: "<sup>e<sub>: joint first/corresponding. <sup>u<sub>:undergrad. <sup>h<sub>:high school student. <sup>v<sub>:visiting scholar. Bold: people from the Shiu lab."

Factors influencing and machine learning models predicting read coverage in a genome assembled using short reads 

> ## Step 01: Map genome sequencing reads to tomato genome.

> mapping reads to tomato genome

 - bwa mem -R "@RG\tID:1\tSM:SRR404081\tLB:SRR404081\tPL:SRR404081\tPU:SRR404081" -t 8 Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa  Sly.1.trimmed.fastq Sly.2.trimmed.fastq > 01_SRR404081.sam

> reorder the mapped reads

 - java -jar $PICARD/ReorderSam.jar I=01_SRR404081.sam O=02_SRR404081_reorder.sam REFERENCE=Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa

> transfer sam file to bam file

 - samtools view -bS 02_SRR404081_reorder.sam -o 03_SRR404081_reorder.bam

> sort the bam file by coordinate

 - java -jar $PICARD/SortSam.jar INPUT=03_SRR404081_reorder.bam OUTPUT=04_SRR404081_sorted.bam SORT_ORDER=coordinate

> ## Step 02: Determine the optimal bin size. Note that bin sizes of 50, 100, 150, 200, 250, 300 bp have been tested

> Determine read depth and regions with variable read coverage using CNVnator (Abyzov et al., 2011) with different bin sizes. Note that in outputs of CNVnator, regions with significantly higher (HC) is referred to as "duplication", while regions with significanly lower (LC) read coverages is "deletion". The remain regions were taken as BG.

 - cnvnator -genome Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa -root Sly.root -tree 04_SRR404081_sorted.bam -unique

 - cnvnator -genome Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa -root Sly.root -his bin_size -d genome/

 - cnvnator -genome Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa -root Sly.root -stat bin_size

 - cnvnator -genome Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa -root Sly.root -partition bin_size

 - cnvnator -genome Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa -root Sly.root -call bin_size > Sly.cnv

> Determine the optimal bin size. The optimal bin size was the bin size leading to a ratio of RD average to RD standard deviation of ~4-5 as suggested (Abyzov et al., 2011)

 - cnvnator -root my.root -his bin_size -d genome/

 - cnvnator -root my.root -eval bin_size > out_put

> ## Step 03: save the read depth (RD) in each bin region to a file

 - awk '{print $2}END{print "exit"}' Sly.cnv | cnvnator -root Sly.root -genotype bin_size > out_put
 
> ## Step 04: resample reads based on different strategies

Reads were resampled based on simulated RDs: i) the only possible RD values were 0 (LC), 1 (BG), or 2 (HC) regions; ii) the analysisoriginal RD values were discretized (rounded) to their closest integers; iii) the analysis RD were used; iiii) the analysis RD but resampled at X-fold coverage.

- python 04_01_resample_reads_based_on_different_strategies.py Sly.cnv Solanum_lycopersicum_GCF_000188115.3_SL2.50_genomic.fa strategies_to_resample_read   ### Note that strategies_to_resample_read can be one of "012", "rounded", "analysis"

- python 04_02_resample_reads_with_different_coverages.py read_file number 

> ## Step 05: HC/LC designation filtering

Adjust p-values (q-value), then filter HC/LC regions based on q-value and q0 value. Before that, q-value threshold was determined by evaluating the consistency between results using original reads and resampled reads. 

 - R --vanilla --slave --args input_file output_file < 05_01_pvalue_adjustment.r  ### Note that the input_file is the Sly.cnv resulted from CNVnator
 
 - python 05_02_compare_two_CNVnator_calling.py original_HC/LC/BG_calling simulated_HC/LC/BG_calling
 
 - python 05_03_draw_F1_at_different_qvalues.py 

> ## Step 06: Determine HC/LC/BG regions with high confidence

Compare RD in HC, LC, and BG regions, then determine the threshold to call HC/LC/BG regions with high confidence.

 - R --vanilla --slave --args RD_from_CNVnator HC_LC_designation_from_CNVnator HC_LC_designation_from_CNVnator_after_filtering output_HC output_BG output_LC < 06_threshold_to_distinguishing_HC_LC_BG_regions.r
 
> ## Step 07: Accuracy and precision of HC/LC/BG calling

To evaluate the accuracy and precision of HC/LC/BG calling, resample reads from genome with simulated RDs, and rerun CNVnator again and compare the original RDs (anslysis RDs) with the new RD. The simulated RDs are: i) the only possible RD values were 0 (LC), 1 (BG), or 2 (HC) regions; ii) the analysisoriginal RD values were discretized (rounded) to their closest integers; iii) the analysis RD were used.
          
> ## Step 08: Impact of read coverage on RD determination

To assess the extent to which read coverage impacts the RD determination, reads at 5-fold, 10-fold, 20-fold and 30-fold coverages were resampled from dataset2, and were used to rerun CNVnator.

> ## Step 09: Get the features to build the machine learning models
** GC content **

** Density of genes **

** K-mer **

** Tandem repeats, or Simple sequence repeat **


> ## Step 10: Build prediction models using Random Forest

- python ML_classification.py -df your_dataset -alg RF -cm T -plots T -gs T -n 100 -cv 10 -p 8 -apply all -threshold_test accuracy -cl_train your_classes  ### ML_classification.py see https://github.com/ShiuLab/ML-Pipeline, contributed by Christina Azodi

> ## Step 11: Get importance value of each feature in the prediction models
  
 - python 11_tree-based_FS.py your_dataset your_classes_list

> ## Step 12: Select K-mer and SSR features based on p-value of Kruskal-Wallis H test

- python 12_KS_test.py path_to_your_files output_name
       
> ## Step 13: Determine the correlation between densities of HC/LC/BG regions with densities of genome features.

 - python 13_density_across_genome.py  ### input files see notes in the script
 
> ## Step 14: Estimate the prevalence of HC/LC/BG across the genome

To assess how the observed correlation between densities of HC/LC/BG regions with densities of genome features derived from random expectation, HC/LC/BG regions were reshuffled across the genome 1000 times, where HC/LC/BG regions don't overlap with each other.

 - python 14_01_reshuffle_HC_LC_BG.py path_to_save HC_regions LC_regions Length_of_chromosomes
 
 - python 14_02_density_of_reshuffled_HC_LC_BG_in_500Kb.py Length_of_chromosomes correlation_NCBI_SGN path_to_files 

> ## Step 15: Estimate gene proproties located in HC regions

To evaluate the properties of genes located in HC regions, first the number of each type of gene overlapping with HC regions were determined. Next, regions with same number and length of HC regions were randomly selected from BG regions, which was repeat 10,000 times. The observed numbers of genes were compared with the random null distribution to see which type of gene is over- or under-representative.



> ## Step 16: Gene set enrichment analysis (Fisher’s exact test)
 
  - python 16_Test_Fisher.py your_file 0/1    ## contributed by Sahra Uygun
 



 


