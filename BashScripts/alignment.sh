
# Align using STAR aligner and saving the aligned reference genome in the memory 
# Therefore we can load several runs with more than one thread
# Place in the path where the preprocessed fastq are available
STAR --genomeLoad LoadAndExit --runThreadN 12 --genomeDir /home/m/mz/041/q041mz/scratch/genome/mmu_prim_star
 
for i in *.fastq
do 
	echo $i
	# parameters
	STAR --genomeDir /home/m/mz/041/q041mz/scratch/genome/mmu_prim_star\
	     --genomeLoad LoadAndKeep \
	     --limitBAMsortRAM 10000000000 \
	     --readFilesIn $i \
	     --runThreadN 12 \
	     --outFilterMismatchNoverLmax 0.6 \
	     --outFilterMultimapNmax 200 \
	     --alignIntronMax 1000000 \
	     --outFilterType BySJout \
	     --limitOutSJcollapsed 5000000 \
	     --alignMatesGapMax 1000000 \
             --alignSJDBoverhangMin 1 --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate --chimOutType WithinBAM --outFileNamePrefix ${i}_aligned
	
		

done

# should remove the index from the memory
# Be aware if crashes than manually invoke this command to free the memory
STAR --genomeLoad Remove --genomeDir /home/m/mz/041/q041mz/scratch/genome/mmu_prim_star

