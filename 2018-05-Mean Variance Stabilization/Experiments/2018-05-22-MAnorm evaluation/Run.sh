cat sample1_rep1.txt > sample1_rep1.bed
cat sample2_rep1.txt > sample2_rep1.bed
$ awk '{print $2,$3,$4}' *.txt > collapsed_output.txt
Bedtools intersect -a sample1_rep1.bed -b sample2_rep1.bed -sortout >sample1_rep1_signals.bed
Bedtools intersect -a sample2_rep1.bed -b sample1_rep1.bed -sortout> sample2_rep1_signals.bed
awk '{print $1,$2,$3}' sample1_rep1_signals.bed > MAnorm.bed