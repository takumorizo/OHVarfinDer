ruby ../../../utils/simulate_sam2.rb 80 0.3 ./ 1
sh ../../../utils/to_bam.sh tumor;sh ../../../utils/to_bam.sh normal
samtools mpileup -B -f ../random_ref.fasta normal.bam > normal.pileup; samtools mpileup -B -f ../random_ref.fasta tumor.bam > tumor.pileup 
sh ../../../utils/make_windows.sh tumor.pileup normal.pileup ./ 
cat windows
