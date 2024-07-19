#Adrian Serohijos 
# edited by DGL

# where $1 is the condition list

# where $2 is the sample type (ex: no_drug)

# where $3 is the data directory

###### POOL BARCODES ######
#	Since we have multiple time points per condition (mouse), we pool them to get a consistent list of clusters
#	For each condition (mouse) in list:
#	1. Append extracted barcodes from all time points to file (pooled-timepoints.txt)
#	2. cluster the pooled barcodes
###########################

for COND in `cat $1`;
do
    echo "Processing..." ${COND}
    mkdir $3/pooled_barcodes/${COND}-pooled/

    cat $3/extracted_barcodes/$2/pooled_initial_libraries/*barcode.txt > $3/pooled_barcodes/${COND}-pooled/pooled-timepoints.txt
    cd $3/pooled_barcodes/${COND}-pooled/

    echo "Clustering pooled barcodes... "
    # cluster the barcodes
    bartender_single_com -f pooled-timepoints.txt -o ${COND} -d 2

    cd -
done
