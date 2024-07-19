#Adrian Serohijos 
# edited DGL

# where $1 is a list of samples

# where $2 is the data directory

# where $3 is the data type (ex: no_drug)


###### EXTRACT BARCODES ######
#   For a given sequencing sample from the list:
#   1. unpack all fastq files
#   2. pool in a single fastq
#   3. prepend a base to each read (required for bartender)
#   4. extract barcodes with bartender
#   5. cluster the barcodes
##############################

for SAMPLE in `cat $1`;
do
    echo "Processing " $SAMPLE
    cd $2/raw/$3/$SAMPLE/
    gunzip -d *fastq.gz

    #pool reads
    echo "Pooling barcodes from replicates ... "
    cat *.fastq > pooled.fastq

    gzip SRR*.fastq

    # prepend N to read
    # prepend ? to quality
    echo "Preparing fastq."
    awk 'NR % 4 == 2 {sub(/^/,"N")} {print}' pooled.fastq | awk 'NR % 4 == 0 {sub(/^/,"?")} {print}' > pooled.fastq-v2

    # extract the barcodes
    echo "Extracting barcodes."
    bartender_extractor_com -f pooled.fastq-v2 -o ${SAMPLE} -q + -p N[15]TATC -m 3

    rm pooled.fastq-v2

    cd -

    mkdir $2/extracted_barcodes/$3/$SAMPLE/

    mv $2/raw/$3/$SAMPLE/*_barcode.txt $2/extracted_barcodes/$3/$SAMPLE/

    cd $2/extracted_barcodes/$3/$SAMPLE/
	
    echo "Clustering pooled barcodes... "
    # cluster the barcodes
    bartender_single_com -f ${SAMPLE}_barcode.txt -o ${SAMPLE} -d 2

    cd -
    
done
