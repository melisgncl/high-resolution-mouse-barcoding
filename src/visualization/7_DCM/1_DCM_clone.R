###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### In here we only analyze clone stability (gf and rm clones)

#######################    ^^^^^^^^   #####################################

source("~/Desktop/mouse_barcoding_last/src/visualization/7_DCM/DCM.R")
#################################################gf########################################################
#################################################gf########################################################
#################################################gf########################################################
set.seed(170)
d_gf1=get_time_derivative_clone("gf1")
d_gf2=get_time_derivative_clone("gf2")
d_gf3=get_time_derivative_clone("gf3")
d_gf4=get_time_derivative_clone("gf4")


e_gf1=get_time_varying_eigen_values(d_gf1[[1]],d_gf1[[2]],10,"gf1",d_gf1[[3]])
e_gf2=get_time_varying_eigen_values(d_gf2[[1]],d_gf2[[2]],10,"gf2",d_gf1[[3]])
e_gf3=get_time_varying_eigen_values(d_gf3[[1]],d_gf3[[2]],10,"gf3",d_gf1[[3]])
e_gf4=get_time_varying_eigen_values(d_gf4[[1]],d_gf4[[2]],10,"gf4",d_gf1[[3]])


plot_eigen_value_evolution(e_gf1[[1]],"gf1")
plot_eigen_value_evolution(e_gf2[[1]],"gf2")
plot_eigen_value_evolution(e_gf3[[1]],"gf3")
plot_eigen_value_evolution(e_gf4[[1]],"gf4")


gf=rbind(e_gf1[[2]],e_gf2[[2]],e_gf3[[2]],e_gf4[[2]]) %>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 


get_all_umap_with_cluster_by_one_by(gf,"gf")

##cluster number size quantification##
plot_umap_cluster_quantification("gf")


get_chosen_umap(gf,21,"centroid","gf")

###We also check the UMAPs without gf3 since it partions from the other mouse in the cohort###

gf_without_gf3=rbind(e_gf1[[2]],e_gf2[[2]],e_gf4[[2]]) %>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 
get_all_umap_with_cluster_by_one_by(gf_without_gf3,"gf_without_gf3")

get_chosen_umap(gf_without_gf3,36,"centroid","gf_without_gf3")

##cluster number size quantification##
plot_umap_cluster_quantification("gf_without_gf3")

#################################################rm########################################################
#################################################rm########################################################
#################################################rm########################################################

d_clone_rm1=get_time_derivative_clone("rm1")
d_clone_rm2=get_time_derivative_clone("rm2")
d_clone_rm3=get_time_derivative_clone("rm3")
d_clone_rm4=get_time_derivative_clone("rm4")


e_clone_rm1=get_time_varying_eigen_values(d_clone_rm1[[1]],d_clone_rm1[[2]],10,"rm1",d_clone_rm1[[3]])
e_clone_rm2=get_time_varying_eigen_values(d_clone_rm2[[1]],d_clone_rm2[[2]],10,"rm2",d_clone_rm2[[3]])
e_clone_rm3=get_time_varying_eigen_values(d_clone_rm3[[1]],d_clone_rm3[[2]],10,"rm3",d_clone_rm3[[3]])
e_clone_rm4=get_time_varying_eigen_values(d_clone_rm4[[1]],d_clone_rm4[[2]],10,"rm4",d_clone_rm4[[3]])


plot_eigen_value_evolution(e_clone_rm1[[1]],"e_clone_rm1")
plot_eigen_value_evolution(e_clone_rm2[[1]],"e_clone_rm2")
plot_eigen_value_evolution(e_clone_rm3[[1]],"e_clone_rm3")
plot_eigen_value_evolution(e_clone_rm4[[1]],"e_clone_rm4")

clone_rm=rbind(e_clone_rm1[[2]],e_clone_rm2[[2]],e_clone_rm3[[2]],e_clone_rm4[[2]])%>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 

get_all_umap_with_cluster_by_one_by(clone_rm,"clone_rm")
##cluster number size quantification##

plot_umap_cluster_quantification("clone_rm")

####check all possible umap and choose the that shows all the patterns##
get_chosen_umap(clone_rm,65,"centroid","clone_rm")



#########rm-clones with gf-clones######

clone_rm_gf=rbind(e_clone_rm1[[2]],e_clone_rm2[[2]],e_clone_rm3[[2]],e_clone_rm4[[2]],
                  e_gf1[[2]],e_gf2[[2]],e_gf3[[2]],e_gf4[[2]])%>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero()


get_all_umap_with_cluster_pairwise(clone_rm_gf,"clone_rm_gf")

plot_umap_cluster_quantification("clone_rm_gf")

