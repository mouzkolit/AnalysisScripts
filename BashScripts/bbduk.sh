# this should trim poly a tail and short sequences
# put in path of the fastq files
for sample in *uma*fastq
do
	echo $i
	bbduk.sh in=$sample out=${sample}_clean.fastq ref=../../polyA.fa.gz,../../truseq.fa.gz k=13 ktrim=r useshortkmers=t mink=5 qtrim=t trimq=10 minlength=20  
done
