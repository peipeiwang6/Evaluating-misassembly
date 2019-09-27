# Evaluating-misassembly
# Step 1: map genome sequencing reads to tomato genome.

# Step 2: determine the optimal bin size.

# Step 3: determine the read depth (RD) in each bin region.

# Step 4: determine regions with significantly higher (HC) or lower (LC) read coverages than genome wide average (background, BG). 
          In outputs of CNVnator, HC is duplication, LC is deletion. The remain regions were taken as BG.
          
# Step 5: Adjust p-values (q-value), then filter HC/LC regions based on q-value and q0 value. Before that q-value threshold was determined
          by evaluating the consistency between results using original reads and resampled reads. Reads were resampled based on simulated 
          RDs: i) the only possible RD values were 0 (LC), 1 (BG), or 2 (HC) regions; ii) the analysisoriginal RD values were discretized 
          (rounded) to their closest integers; iii) the analysis RD were used.

# Step 6: Compare RD in HC, LC, and BG regions, then determine the threshold to call HC/LC/BG regions with high confidence.

# Step 7: To evaluate the accuracy and precision of HC/LC/BG calling, resample reads from genome with simulated RDs, and rerun CNVnator 
          again and compare the original RDs (anslysis RDs) with the new RD. The simulated RDs are: i) the only possible RD values 
          were 0 (LC), 1 (BG), or 2 (HC) regions; ii) the analysisoriginal RD values were discretized (rounded) to their closest integers;
          iii) the analysis RD were used.
          
# Step 8: To assess the extent to which read coverage impacts the detection of HC/LC/BG regions, reads at 5-fold, 10-fold, 20-fold 
          and 30-fold coverages were resampled from dataset2, and were used to rerun CNVnator.
          
# Step 9: Determine the correlation between densities of HC/LC/BG regions with densities of genome features.
       
         
