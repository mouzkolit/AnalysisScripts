# should be used for deduplications
for i in *.bam
do
	samtools index $i # first we need to create an index using samtools
	umi_tools dedup -I $i --output-stats=deduplicated -S ${i}_deduplicated.bam
done
