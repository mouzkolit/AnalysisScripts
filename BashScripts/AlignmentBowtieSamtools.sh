
# Should be opened in the folder where the fastq files are located
# index should beplace where bowtie index is allocated
for i in *.fastq
do 
	output="${i}.sam"
	echo $output
	bowtie -v 1 --best --strata --all -m 50 -S index $i $output # 50 multimapping places defined
	samtools view -b $output > "${output}.bam"
	samtools sort "${output}.bam" >  "sorted_${output}.bam"
	samtools index "sorted_${output}.bam"  
done
