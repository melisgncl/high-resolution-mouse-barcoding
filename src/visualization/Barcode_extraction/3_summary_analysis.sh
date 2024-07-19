#Adrian Serohijos 

# where $1 is the sample list

# where $2 is the sample type (ex: no_drug)

###### SUMMARIZE CONDITION ######
#	Produce a summary table of the bartender output by sample
#################################

printf "#Sample,Total Barcodes,Unique Barcodes\n" > $2_Summary.txt

for SAMPLE in `cat $1`;
do
    cd ~/Documents/exp1_well-C3_NEE/extracted_barcodes/$2/$SAMPLE/

    UNIQ=`sort -t, -k3 -n ~/Documents/exp1_well-C3_NEE/extracted_barcodes/$2/$SAMPLE/${SAMPLE}_barcode.csv| tail -n1 | awk 'BEGIN{FS = ","}; {print $3}'`

    ALLBARCODES=`awk '{print NR}' ~/Documents/exp1_well-C3_NEE/extracted_barcodes/$2/$SAMPLE/${SAMPLE}_barcode.csv | tail -n1`

    cd -

    printf "%s,%d,%d\n" ${SAMPLE} ${ALLBARCODES} ${UNIQ} >> $2_Summary.txt

done
