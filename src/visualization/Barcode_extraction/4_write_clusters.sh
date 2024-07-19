#!/bin/bash

# where $1 is a pooled barcode list

# where $2 is mouse number

###### WRITE CLUSTERS ########
#	Write barcode reads in long format using the pooled output from step #2
#	This uses the unique ID attributed to barcode sequences by bartender
##############################

echo "start"

echo "Time	Reads	ID" > ../Sample_M$2_clustering.txt

{
	#skip header
	read

	declare -a array

	#read first cluster
	IFS=, read -r barcode count id

	 
	#echo "Barcode: $barcode $count $id"

	# grep the barcode in the extracted files from step #1 
	# 'grep -e PATTERN' searches for PATTERN in the file, lines that match are returned
	result=($(grep -e $barcode ../data/extracted_barcodes/M/Sample_M$2*/*_barcode.csv))

	# store result in array
	arraylength=${#result[@]}
	for (( i=0; i<${arraylength}; i++ ));
	do
   		#sample=$(echo "${result[$i]}" | cut -d/ -f3 | grep -o .$)

   		# trim path; NOTE: these numbers are based on the path length (this depends on how you name the samples)
   		sample=${result[$i]:39:2}
		sample=${sample//[!0-9]/}

		reads=$(echo "${result[$i]}" | cut -d, -f 2 -)

		# if exists, add to value
		# else assign
		if [[ ${array[$sample]} ]]; then array[$sample]=$((array[$sample] + $reads)); else array[$sample]=$(($reads)); fi

		array[$sample]=$reads

	done

	LAST=$id

	while IFS=, read -r barcode count id
	do
		# if id is same as previous
	    # add to total
	    if [ $LAST != $id ]; then
	    	for i in "${!array[@]}"; do 
  				printf "%s\t%s\t%s\n" "$i" "${array[$i]}" "$LAST"
			done

			#echo $LAST

			#reinitialize array
			unset array
			declare -a array
		fi

	    #echo "Barcode: $barcode $count $id"
	    result=($(grep -e $barcode ../data/extracted_barcodes/M/Sample_M$2*/*_barcode.csv))
	    arraylength=${#result[@]}

	    # get counts per time point
	    for (( i=0; i<${arraylength}; i++ ));
		do
			sample=${result[$i]:39:2}
			sample=${sample//[!0-9]/}

			reads=$(echo "${result[$i]}" | cut -d, -f 2 -)
   			
   			if [[ ${array[$sample]} ]]; then array[$sample]=$((array[$sample] + $reads)); else array[$sample]=$(($reads)); fi

		done

	    LAST=$id
	done 
} < $1 >> ../Sample_M$2_clustering.txt

