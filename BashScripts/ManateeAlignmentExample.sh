# place in folder of interest
# manatee needs to path variable
#index = place where bowite index is 
for i in *.fastq
do
        sample="${i::6}"
        echo $sample
        # files need to be allocated in the current path otherwise change on location!
        perl manatee -i $i  -o $sample  -index index  -genome GRCh38.primary_assembly.genome.fa  -annotation final_gtf_pirna.gtf
done

