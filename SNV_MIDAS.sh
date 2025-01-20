mapfile -t files < sample.txt

for sample in "${files[@]}"; do

time run_midas.py species ./MIDAS_out/"$sample" -1 ./03decom/"$sample"_paired_1.fastq -2 ./03decom/"$sample"_paired_2.fastq -t 10 --remove_temp

time run_midas.py snps ./MIDAS_out/"$sample" -1 ./03decom/"$sample"_paired_1.fastq -2 ./03decom/"$sample"_paired_2.fastq -t 10 --remove_temp

time run_midas.py genes ./MIDAS_out/"$sample" -1 ./03decom/"$sample"_paired_1.fastq -2 ./03decom/"$sample"_paired_2.fastq -t 10 --remove_temp

done

merge_midas.py snps ./MIDAS_merge -i ./MIDAS_out  -t dir --all_sites --all_samples