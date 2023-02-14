# this should deduplicate the UMI counts
# Place in path where the aligned bam files for deduplication are located
for i in *bam_deduplicated.bam
do
	echo $i
	umi_tools dedup -I $i -S ${i}_final.bam
done
