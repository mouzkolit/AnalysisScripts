for i in *.fastq
#needs to be in the path variable
do 
	echo $(basename $i)
	flexbar  -r $i -a adapter_illumina.fa -t trimmed_$(basename $i)
done

