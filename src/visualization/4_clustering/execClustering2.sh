# Melis Gencel & Louis Gauthier 

# where $1 is a list of samples

ROOTDIR="/Users/melisG/Desktop/mouse_barcoding_last"
PROC_SAMPLES=$ROOTDIR/data/processed_samples
DATA_CLUST=$ROOTDIR/data/clustering


echo "Creating directories for clustering data depending on choosen treshold"
for cond in `cat condition.list`;
do
    for SAMPLE in `cat ${cond}.list`;
    do
        echo $cond $SAMPLE
        mkdir -p $ROOTDIR/reports/figures/clustering/$SAMPLE
        
  done
done



################
#### STEP 5 ####
################
#### choose the best treshold that defines lineage and clone dynamics###
###you have to give the tresholds 
echo "Plotting clusters and loess averages for a given treshold..."

##im##
cp $ROOTDIR/reports/figures/clustering_control/im1/0.34/* $ROOTDIR/reports/figures/clustering/im1/
cp $ROOTDIR/reports/figures/clustering_control/im2/0.32/* $ROOTDIR/reports/figures/clustering/im2/
cp $ROOTDIR/reports/figures/clustering_control/im3/0.61/* $ROOTDIR/reports/figures/clustering/im3/
cp $ROOTDIR/reports/figures/clustering_control/im4/0.37/* $ROOTDIR/reports/figures/clustering/im4/

##rm##
cp $ROOTDIR/reports/figures/clustering_control/rm1/0.55/* $ROOTDIR/reports/figures/clustering/rm1/
cp $ROOTDIR/reports/figures/clustering_control/rm2/0.6/* $ROOTDIR/reports/figures/clustering/rm2/
cp $ROOTDIR/reports/figures/clustering_control/rm3/0.65/* $ROOTDIR/reports/figures/clustering/rm3/
cp $ROOTDIR/reports/figures/clustering_control/rm4/0.38/* $ROOTDIR/reports/figures/clustering/rm4/
  
##gf##  
cp $ROOTDIR/reports/figures/clustering_control/gf1/0.33/* $ROOTDIR/reports/figures/clustering/gf1/
cp $ROOTDIR/reports/figures/clustering_control/gf2/0.39/* $ROOTDIR/reports/figures/clustering/gf2/
cp $ROOTDIR/reports/figures/clustering_control/gf3/0.52/* $ROOTDIR/reports/figures/clustering/gf3/
cp $ROOTDIR/reports/figures/clustering_control/gf4/0.39/* $ROOTDIR/reports/figures/clustering/gf4/
  






################
#### STEP 6 ####
################
echo "Generate dendrograms of clusters of chosen treshold  for each group..."
Rscript 6_cluster_of_clusters.R

################
#### STEP 7 ####
################
echo "Plot overlap index and simulate z scores within the cohort..."

Rscript 7_cluster_similarity.R





