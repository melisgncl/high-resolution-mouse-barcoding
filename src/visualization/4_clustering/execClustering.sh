# Melis Gencel & Louis Gauthier 



ROOTDIR="/Users/melisG/Desktop/mouse_barcoding_last"
PROC_SAMPLES=$ROOTDIR/data/processed_samples
DATA_CLUST=$ROOTDIR/data/clustering


echo "Creating directories for clustering data..."
for cond in `cat condition.list`;
do
    for SAMPLE in `cat ${cond}.list`;
    do
          for cutoff in `cat  cutoff.list`; 
          do
              echo $cond $SAMPLE
              mkdir -p $DATA_CLUST/$SAMPLE
              mkdir -p $ROOTDIR/reports/figures/clustering_control/$SAMPLE
              mkdir -p $ROOTDIR/reports/figures/clustering_control/$SAMPLE/$cutoff
    done
  done
done

################
#### STEP 1 ####
################
echo "Filtering trajectories for clustering..."

###im##
Rscript 1_filter_data.R $PROC_SAMPLES/im1/Sample_M1_clustering.txt $DATA_CLUST/im1/im1_filtered.csv 0 3 5e-6
Rscript 1_filter_data.R $PROC_SAMPLES/im2/Sample_M2_clustering.txt $DATA_CLUST/im2/im2_filtered.csv 0 3 5e-6
Rscript 1_filter_data.R $PROC_SAMPLES/im3/Sample_M3_clustering.txt $DATA_CLUST/im3/im3_filtered.csv 0 5 5e-6
Rscript 1_filter_data.R $PROC_SAMPLES/im4/Sample_M4_clustering.txt $DATA_CLUST/im4/im4_filtered.csv 0 7 5e-6

## rm ##
Rscript 1_filter_data.R $PROC_SAMPLES/rm1/Sample_M1_clustering.txt $DATA_CLUST/rm1/rm1_filtered.csv 0 12 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/rm2/Sample_M2_clustering.txt $DATA_CLUST/rm2/rm2_filtered.csv 0 12 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/rm3/Sample_M3_clustering.txt $DATA_CLUST/rm3/rm3_filtered.csv 0 12 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/rm4/Sample_M4_clustering.txt $DATA_CLUST/rm4/rm4_filtered.csv 0 12 5e-5

## Germ-free ##
Rscript 1_filter_data.R $PROC_SAMPLES/gf1/Sample_GFM1_clustering.txt $DATA_CLUST/gf1/gf1_filtered.csv 0 13 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/gf2/Sample_GFM2_clustering.txt $DATA_CLUST/gf2/gf2_filtered.csv 0 13 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/gf3/Sample_GFM3_clustering.txt $DATA_CLUST/gf3/gf3_filtered.csv 0 14 5e-5
Rscript 1_filter_data.R $PROC_SAMPLES/gf4/Sample_GFM4_clustering.txt $DATA_CLUST/gf4/gf4_filtered.csv 0 13 5e-5



################
#### STEP 2 ####
################
echo "Computing correlation matrix..."
for cond in `cat condition.list`;
do
    for SAMPLE in `cat ${cond}.list`;
    do
        echo Sample $SAMPLE
        python 2_apply_clustering.py $DATA_CLUST/$SAMPLE/${SAMPLE}_filtered.csv $DATA_CLUST/$SAMPLE $SAMPLE
    done
done

################
#### STEP 3 ####
################
# See 3_dendrogram_upgma.py for details
# but remember that the last argument is the threshold for flattening the clusters
# this threshold is set manually by the user based on inspection of the dendrogram/heatmap
# to choose a treshold run a range between 0.1-0.9 with step size 0.02 and check smallest distance between loess then choose a cut off
# and plot clusters and loess averages
echo "Generating dendrogram and flattening clusters..."
for cond in `cat condition.list`;
do
    for SAMPLE in `cat ${cond}.list`;
    do
        for cutoff in `cat  cutoff.list`; 
        do 
              echo $cond $SAMPLE $cutoff
              python3 3_dendrogram_upgma.py $DATA_CLUST/$SAMPLE/${SAMPLE}_dist.csv $ROOTDIR/reports/figures/clustering_control/$SAMPLE/$cutoff $SAMPLE $cutoff
              echo  "Plotting clusters and loess averages..."
              Rscript 4_plot_clusters_loess.R $ROOTDIR/reports/figures/clustering_control/$SAMPLE/$cutoff/clusters_${SAMPLE}_average $DATA_CLUST/$SAMPLE/${SAMPLE}_filtered.csv $ROOTDIR/reports_article/figures/clustering_control/$SAMPLE/$cutoff/$SAMPLE $cond

        done
    done
done








################
#### STEP 4 ####
################
#### quantify the clustering to choose for best treshold to  define clone dynamics
echo "Generate plot of hierhical cluster quantifcation for each group..."
Rscript 5_h_clustering_quantification.R







