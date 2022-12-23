## CODE BLOCK
## Louis Gauthier, August 2020

#######################    READ ME!   #####################################

##### Merge barcode list with clusters to fetch the top N barcodes per sample

#######################    ^^^^^^^^   #####################################

n_intersect=1000

###############
# read pooled cluster tables

gf1_cluster <- read_csv("data/processed_samples/gf1/GFM1_cluster.csv")
gf2_cluster <- read_csv("data/processed_samples/gf2/GFM2_cluster.csv")
gf3_cluster <- read_csv("data/processed_samples/gf3/GFM3_cluster.csv")
gf4_cluster <- read_csv("data/processed_samples/gf4/GFM4_cluster.csv")

rm1_cluster <- read_csv("data/processed_samples/rm1/M1_cluster.csv")
rm2_cluster <- read_csv("data/processed_samples/rm2/M2_cluster.csv")
rm3_cluster <- read_csv("data/processed_samples/rm3/M3_cluster.csv")
rm4_cluster <- read_csv("data/processed_samples/rm4/M4_cluster.csv")

im1_cluster <- read_csv("data/processed_samples/im1/M1_cluster.csv")
im2_cluster <- read_csv("data/processed_samples/im2/M2_cluster.csv")
im3_cluster <- read_csv("data/processed_samples/im3/M3_cluster.csv")
im4_cluster <- read_csv("data/processed_samples/im4/M4_cluster.csv")





fetchTop <- function(reshaped_df, sample_cluster) {
	df.top = unique(reshaped_df[,1:4])
	df.top.final = df.top[order(-df.top$final),]
	df.top.max = df.top[order(-df.top$max),]

	df.top.final = df.top.final[1:n_intersect,]
	df.top.max = df.top.max[1:n_intersect,]

	# match back to the barcode cluster sequence
	df.top.final = merge(df.top.final,sample_cluster,by.x = "ID",by.y = "Cluster.ID")
	df.top.final$Cluster.Score=NULL
	df.top.final$time_point_1=NULL

	df.top.final = df.top.final[order(-df.top.final$final),]

	df.top.max = merge(df.top.max,sample_cluster,by.x = "ID",by.y = "Cluster.ID")
	df.top.max$Cluster.Score=NULL
	df.top.max$time_point_1=NULL

	df.top.max = df.top.max[order(-df.top.max$max),]

	return(list(df.top.final,df.top.max))
}

###############
# get top N IDs based on max frequency


###gf###
x = fetchTop(gf1.df,gf1_cluster)
gf1.top.final = as.data.frame(x[1])
gf1.top.max = as.data.frame(x[2])

x = fetchTop(gf2.df,gf2_cluster)
gf2.top.final = as.data.frame(x[1])
gf2.top.max = as.data.frame(x[2])

x = fetchTop(gf3.df,gf3_cluster)
gf3.top.final = as.data.frame(x[1])
gf3.top.max = as.data.frame(x[2])

x = fetchTop(gf4.df,gf4_cluster)
gf4.top.final = as.data.frame(x[1])
gf4.top.max = as.data.frame(x[2])

###rm###
x= fetchTop(rm1.df,rm1_cluster)
rm1.top.final = as.data.frame(x[1])
rm1.top.max = as.data.frame(x[2])

x = fetchTop(rm2.df,rm2_cluster)
rm2.top.final = as.data.frame(x[1])
rm2.top.max = as.data.frame(x[2])

x = fetchTop(rm3.df,rm3_cluster)
rm3.top.final = as.data.frame(x[1])
rm3.top.max = as.data.frame(x[2])

x = fetchTop(rm4.df,rm4_cluster)
rm4.top.final = as.data.frame(x[1])
rm4.top.max = as.data.frame(x[2])

###im###
x= fetchTop(im1.df,im1_cluster)
im1.top.final = as.data.frame(x[1])
im1.top.max = as.data.frame(x[2])

x = fetchTop(im2.df,im2_cluster)
im2.top.final = as.data.frame(x[1])
im2.top.max = as.data.frame(x[2])

x = fetchTop(im3.df,im3_cluster)
im3.top.final = as.data.frame(x[1])
im3.top.max = as.data.frame(x[2])

x = fetchTop(im4.df,im4_cluster)
im4.top.final = as.data.frame(x[1])
im4.top.max = as.data.frame(x[2])

