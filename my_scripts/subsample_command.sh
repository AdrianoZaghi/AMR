Sample_Gbp_size=0.01
Average_read_length=247
Number_of_reads_to_peek=$(echo "$Sample_Gbp_size * 10^9 / $Average_read_length" | bc)

for FILE in samples_trim/DCWGS10*; do 
	filename=$(basename "$FILE" .gz)
	seqtk sample -s111 $FILE $Number_of_reads_to_peek > "samples_trim_sub/$filename";
	echo $filename
done