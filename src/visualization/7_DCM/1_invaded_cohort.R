###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### In here we only analyze gavaged cohorts with microbiota (im and rm)

#######################    ^^^^^^^^   #####################################

source("~/Desktop/mouse_barcoding_last/src/visualization/7_DCM/DCM.R")
#################################################im########################################################
#################################################im########################################################
#################################################im########################################################

im1_16S_taxa <- read_delim("data/16S/taxa_long_format/im1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
im2_16S_taxa <- read_delim("data/16S/taxa_long_format//im2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
im3_16S_taxa <- read_delim("data/16S/taxa_long_format//im3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
im4_16S_taxa <- read_delim("data/16S/taxa_long_format//im4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)


##calculate z dot##
d_im1=get_time_derivative_invaded(im1_16S_taxa,"im1",0.25)
d_im2=get_time_derivative_invaded(im2_16S_taxa,"im2",0.25)
d_im3=get_time_derivative_invaded(im3_16S_taxa,"im3",0.25)
d_im4=get_time_derivative_invaded(im4_16S_taxa,"im4",0.25)

##get the time evolving eigenspace ###

e_im1=get_time_varying_eigen_values(d_im1[[1]],d_im1[[2]],10,"im1",d_im1[[4]])
e_im2=get_time_varying_eigen_values(d_im2[[1]],d_im2[[2]],10,"im2",d_im2[[4]])
e_im3=get_time_varying_eigen_values(d_im3[[1]],d_im3[[2]],10,"im3",d_im3[[4]])
e_im4=get_time_varying_eigen_values(d_im4[[1]],d_im4[[2]],10,"im4",d_im4[[4]])

##plot eigenvalues##
plot_eigen_value_evolution(e_im1[[1]],"im1")
plot_eigen_value_evolution(e_im2[[1]],"im2")
plot_eigen_value_evolution(e_im3[[1]],"im3")
plot_eigen_value_evolution(e_im4[[1]],"im4")

##data frame for UMAP
im=rbind(e_im1[[2]],e_im2[[2]],e_im3[[2]],e_im4[[2]])%>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 

###check all possible UMAPs with 9 different clustering##
get_all_umap_with_cluster_by_one_by(im,"im")

##cluster number size quantification##
plot_umap_cluster_quantification("im")

####check all possible umap and choose the that shows all the patterns##
get_chosen_umap(im,9,"centroid","im")

###plot z dot with z###
plot_derivative(d_im1[[3]],d_im1[[2]],"im1")
plot_derivative(d_im2[[3]],d_im2[[2]],"im2")
plot_derivative(d_im3[[3]],d_im3[[2]],"im3")
plot_derivative(d_im4[[3]],d_im4[[2]],"im4")



#################################################rm########################################################
#################################################rm########################################################
#################################################rm########################################################




rm1_16S_taxa <- read_delim("data/16S/taxa_long_format/rm1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm2_16S_taxa <- read_delim("data/16S/taxa_long_format//rm2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm3_16S_taxa <- read_delim("data/16S/taxa_long_format//rm3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
rm4_16S_taxa <- read_delim("data/16S/taxa_long_format//rm4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)



##calculate z dot##
d_rm1=get_time_derivative_invaded(rm1_16S_taxa,"rm1",0.2)
d_rm2=get_time_derivative_invaded(rm2_16S_taxa,"rm2",0.2)
d_rm3=get_time_derivative_invaded(rm3_16S_taxa,"rm3",0.2)
d_rm4=get_time_derivative_invaded(rm4_16S_taxa,"rm4",0.2)


##get the time evolving eigenspace ###

e_rm1=get_time_varying_eigen_values(d_rm1[[1]],d_rm1[[2]],10,"rm1",d_rm1[[4]])
e_rm2=get_time_varying_eigen_values(d_rm2[[1]],d_rm2[[2]],10,"rm2",d_rm2[[4]])
e_rm3=get_time_varying_eigen_values(d_rm3[[1]],d_rm3[[2]],10,"rm3",d_rm3[[4]])
e_rm4=get_time_varying_eigen_values(d_rm4[[1]],d_rm4[[2]],10,"rm4",d_rm4[[4]])

##plot eigenvalues##
plot_eigen_value_evolution(e_rm1[[1]],"rm1")
plot_eigen_value_evolution(e_rm2[[1]],"rm2")
plot_eigen_value_evolution(e_rm3[[1]],"rm3")
plot_eigen_value_evolution(e_rm4[[1]],"rm4")

##data frame for UMAP
rm=rbind(e_rm1[[2]],e_rm2[[2]],e_rm3[[2]],e_rm4[[2]])%>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 

###check all possible UMAPs with 9 different clustering##
get_all_umap_with_cluster_by_one_by(rm,"rm")

##cluster number size quantification##
plot_umap_cluster_quantification("rm")

####check all possible umap and choose the that shows all the patterns##
get_chosen_umap(rm,65,"centroid","rm")

###plot z dot with z###
plot_derivative(d_rm1[[3]],d_rm1[[2]],"rm1")
plot_derivative(d_rm2[[3]],d_rm2[[2]],"rm2")
plot_derivative(d_rm3[[3]],d_rm3[[2]],"rm3")
plot_derivative(d_rm4[[3]],d_rm4[[2]],"rm4")





