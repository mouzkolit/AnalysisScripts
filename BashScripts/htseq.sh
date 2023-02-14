
# important for counting after STAR aligment -> used for single place alignment of miRNAs
for i in *sortedByCoord.out.bam

do 
	python -m HTSeq.scripts.count s- yes -r pos -t miRNA -f bam $i  mmu.gff3
done
