# Evaluating-misassembly
Factors influencing and machine learning models predicting read coverage in a genome assembled using short reads

> ## Step 01: Map genome sequencing reads to tomato genome.

> ## Step 02: Determine the optimal bin size.

> ## Step 03: Determine the read depth (RD) in each bin region.

> ## Step 04: HC/LC designation

Determine regions with significantly higher (HC) or lower (LC) read coverages than genome wide average (background, BG). Note that in outputs of CNVnator, HC is referred to as "duplication", while LC is "deletion". The remain regions were taken as BG.
          
> ## Step 05: HC/LC designation filtering

Adjust p-values (q-value), then filter HC/LC regions based on q-value and q0 value. Before that, q-value threshold was determined by evaluating the consistency between results using original reads and resampled reads. Reads were resampled based on simulated RDs: i) the only possible RD values were 0 (LC), 1 (BG), or 2 (HC) regions; ii) the analysisoriginal RD values were discretized (rounded) to their closest integers; iii) the analysis RD were used.

> ## Step 06: Determine HC/LC/BG regions with high confidence

Compare RD in HC, LC, and BG regions, then determine the threshold to call HC/LC/BG regions with high confidence.

 - R --vanilla --slave --args RD_from_CNVnator HC_LC_designation_from_CNVnator HC_LC_designation_from_CNVnator_after_filtering output_HC output_BG output_LC < 06_threshold_to_distinguishing_HC_LC_BG_regions.r
 
> ## Step 07: Accuracy and precision of HC/LC/BG calling

To evaluate the accuracy and precision of HC/LC/BG calling, resample reads from genome with simulated RDs, and rerun CNVnator again and compare the original RDs (anslysis RDs) with the new RD. The simulated RDs are: i) the only possible RD values were 0 (LC), 1 (BG), or 2 (HC) regions; ii) the analysisoriginal RD values were discretized (rounded) to their closest integers; iii) the analysis RD were used.
          
> ## Step 08: Impact of read coverage on RD determination

To assess the extent to which read coverage impacts the RD determination, reads at 5-fold, 10-fold, 20-fold and 30-fold coverages were resampled from dataset2, and were used to rerun CNVnator.

> ## Step 09: Get the features to build the machine learning models
** GC content **


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
 


