###all fucntuon with to work on terminal##
source("~/Desktop/0106_mice_barcode_data/src_article//visualization/0_config/0_config.R")
library(tidyr)
library(umap)
library(ggforce)
library(factoextra )
library(NbClust)



`%notin%` <- Negate(`%in%`)
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

countZeroes = function(col) {
  sum(col==0)<7
}

transformTaxa = function(taxa.data,cluster.data){
  # exclude taxa with too many zero-value entries
  taxa.data = taxa.data[,apply(taxa.data,2,countZeroes)]
  
  ##### log transform data
  # deal with zero values (0 = resolution limit)
  # log(x+0.000001)
  taxa.series = tslist(t(log10(taxa.data+0.000001)))
  
  # interpolate taxa series to length of barcode cluster series
  taxa.series <- reinterpolate(taxa.series, new.length = nrow(cluster.data))
  
  return(taxa.series)
}

get_time_derivative_gavaged = function(taxa,sample){
  clusters.loess=read_csv(file= paste("reports_article//figures/clustering/",sample,"/",sample,"_clustered_loess_log10.csv",sep=""))
  clusters=clusters.loess
  clusters$time=NULL
  taxa.series = transformTaxa(taxa,clusters)
  cluster.series = tslist(t(clusters))
  ###linearly interpolated series
  series_normal_inter = append(clusters,taxa.series)
  ###modeled taxa##polynomial interpolation
  taxa.data = taxa[,apply(taxa,2,countZeroes)]
  taxa.data=taxa.data[c((min(clusters.loess$time)/10):(max(clusters.loess$time)/10)), ]
  taxa.data = tslist(t(log10(taxa.data+0.000001)))
  Time=seq(from=min(clusters.loess$time)/10, to=max(clusters.loess$time)/10)
  xx <- seq(from=min(clusters.loess$time)/10, to=max(clusters.loess$time)/10, length.out=length(series_normal_inter[[1]]))
  list_models_loess_taxa  = sapply(taxa.data, function(x) predict(loess(x~Time,span=0.2),xx,se = FALSE))
  series_whole_loess= append(clusters,tslist(t(list_models_loess_taxa)))
  ####time_derivative of series modelled taxa###central numerical differences####
  Time=seq(1:length(series_normal_inter[[1]]))
  deriv_series_loess=sapply(series_whole_loess, function(y)  tslist(t((diff(y,lag=2)/diff(Time/10,lag=2)))))
  return(list(series_whole_loess,deriv_series_loess,series_normal_inter,clusters.loess$time))
}

get_time_derivative_GFM=function(sample){
  clusters.loess=read_csv(file= paste("reports_article//figures/clustering/",sample,"/",sample,"_clustered_loess_log10.csv",sep=""))
  #clusters.loess$time=NULL
  Time=seq(1:nrow(clusters.loess))
  series_whole_loess=clusters.loess
  series_whole_loess$time=NULL
  deriv_series_loess=sapply(series_whole_loess, function(y)  tslist(t((diff(y,lag=2)/diff(Time/10,lag=2)))))
  options(scipen=0)
  return(list(tslist(t(series_whole_loess)),deriv_series_loess,Time))
  
}

####to get the time-progressive eigen_values and vectors##
get_time_varying_eigen_values=function(series_whole_loess,deriv_series_loess,window,sample,time){
  options(scipen=0)
  total <- nrow(as.data.table(series_whole_loess))
  spots <- seq(from=10, to=(total), by=window)
  spots[length(spots)]=nrow(as.data.table(deriv_series_loess))
  ##it is a centered we need to reduce 2 steps from the spot to adjust the
  result=list()  
  for(k in 1:length(spots)){
    result[[k]] <- sapply(1:ncol(as.data.table(deriv_series_loess)), function(i) {
      sapply(1:ncol(as.data.table(series_whole_loess)), function(j){
        resTmp <- cov(deriv_series_loess[[i]][(1:spots[k])],series_whole_loess[[j]][(1:spots[k])]) 
        
      })
    })
  }
  eigen_values=sapply(result, function(y) eigen(y, symmetric= FALSE)$values)
  real_eigen_value=sapply(as.data.table(eigen_values), function(y) Re(y))
  im_eigen_value=sapply(as.data.table(eigen_values), function(y) Im(y))
  eigen_vector=lapply(result, function(y) eigen(y, symmetric= FALSE)$vector)
  
  
  ###order eigenvalues###
  conjugate= sqrt(real_eigen_value^2+im_eigen_value^2)
  real_eigen_value=melt(real_eigen_value)  %>% mutate(eigen="R")
  im_eigen_value=melt(im_eigen_value)  %>% mutate(eigen="I")
  conjugate_value=melt(conjugate)  %>% mutate(eigen="C")
  all_eigen_value=rbind(real_eigen_value,im_eigen_value,conjugate_value) %>% spread(eigen,value)  %>% 
    mutate(Var2=sub("V", "", Var2) ,Var1=as.factor(Var1),Var2=(as.numeric(Var2)))  %>% group_by(Var2) %>%
    arrange(desc(C),.by_group = TRUE) %>%
    mutate(Var3=as.factor(row_number()))  
  
  all_eigen_umap=rbind(real_eigen_value,im_eigen_value,conjugate_value) %>% spread(eigen,value)  %>% 
    mutate(Var2=sub("V", "", Var2) ,Var1=as.factor(Var1),Var2=(as.numeric(Var2)))  %>% mutate(Var2= case_when(min(time)/10 == 2 ~ Var2 +1,
                                                                                                              TRUE~ Var2)) %>%
    group_by(Var2) %>%
    arrange(desc(C),.by_group = TRUE) %>% mutate(sample=sample) %>% 
    mutate(Var3=as.factor(row_number())) %>% select(-C,-Var1)  %>% pivot_longer(cols = -c(Var2,Var3,sample), names_to = "variable", values_to = "value") %>% 
    group_by(Var2) %>% mutate(Var4=paste0(variable,"_",Var3)) %>%  select(-Var3,-variable) %>% spread(Var4,value) 
  
  return(list(all_eigen_value,all_eigen_umap))
  
  
}  

###plot the eigenspace evolution###
plot_eigen_value_evolution=function(all_eigen_value,sample){
  #names(all_eigen_value)[names(all_eigen_value) %in% c("X2","X3")]<-c("Var2","Var3")
  
  eigen_space= sort(as.numeric(unique(all_eigen_value$Var2)))
  plots_eigen_value=list()
  for(i in seq_along(eigen_space)){
    plots_eigen_value[[i]] = ggplot(all_eigen_value[all_eigen_value$Var2 %in% c(1:eigen_space[i]) ,],aes(R,I,color=Var3))+
      geom_point(size=5,alpha=0.3) + 
      geom_path(size=1.2, aes(color=Var3), arrow = arrow(angle = 15, type = "closed"), alpha=0.5) +
      geom_point(aes(R,I,color=Var3),size=4,data=all_eigen_value[all_eigen_value$Var2==eigen_space[i],]) +
      theme_Publication() +xlab(~ paste("Re", "(",lambda [i],")"))+ ylab(~ paste("lm", "(",lambda [i],")")) +
      geom_vline(xintercept = 0)+geom_hline(yintercept = 0) + xlim(min(all_eigen_value$R),max(all_eigen_value$R)) +
      ylim(min(all_eigen_value$I),max(all_eigen_value$I)) + guides(color=FALSE) +scale_color_manual(values = eigen.colors)
    ggsave(plots_eigen_value[[i]],
           filename = paste("reports_article//figures/third_figure/eigen_values/",sample, "_eigen_value_",eigen_space[i],"_time_interval.eps", sep=""),width =10.5,height =7.5, limitsize = FALSE,device = cairo_ps )
    ggsave(plots_eigen_value[[i]],
           filename = paste("reports_article//figures/third_figure/eigen_values/",sample, "_eigen_value_",eigen_space[i],"_time_interval.jpeg", sep=""),width =10.5,height =7.5, limitsize = FALSE )
  }
  cowplot::plot_grid(plotlist = plots_eigen_value,align = "hv") %>%  
    ggsave(filename = paste("reports_article//figures/third_figure/eigen_values_panel/",sample, "_all_eigen_value_time_interval.eps", sep=""),width =40,height =32, limitsize = FALSE,device = cairo_ps )
  cowplot::plot_grid(plotlist = plots_eigen_value,align = "hv") %>%  
    ggsave(filename = paste("reports_article//figures/third_figure/eigen_values_panel/",sample, "_all_eigen_value_time_interval.svg", sep=""),width =40,height =32, limitsize = FALSE)
  return(plots_eigen_value)
}

##search all possible neighbour paramters of  UMAP with diffrent 9 clusters algorithm ##

get_all_umap_with_cluster_by_one_by=function(df,sample){
  tryCatch({ 
    #names(df)[names(df) %in% c("X2")]<-c("Var2")
    my_palette <- colorRampPalette(c("purple","darkblue","#66a266","darkgreen","orange","darkorange"))(n =length(unique(df$Var2)))
    M_df=df[-c(1,2)]
    digits.labels = df[, c("Var2","sample")]
    dim=seq_along(2:nrow(df))+1
    clustering_methods=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans")
    plot=list()
    cluster_info2=list()
    for( c in seq_along(clustering_methods)){
      cluster=list()
      cluster_info=list()
      for(i in dim) { 
        set.seed(170)
        digits.umap = umap(M_df, n_components =2,n_neighbors=i,nn_method="annoy",n_trees="100",min_dist=0.1)
        #digits.umap = umap(M_spec_df, n_components =2,n_neighbours=15,min_dist=0.2) #15 is the standar
        layout <- digits.umap[["layout"]] 
        layout <- data.frame(layout) 
        final <- cbind(layout, df[, c("Var2","sample")])
        colnames(final) <- c('umap1', 'umap2', 'phase',"mouse") 
        plot[[i]]=final %>% ggplot(aes(x = umap1, y = umap2,color = phase,shape=mouse))+
          geom_point(size=4) +
          theme(legend.position="bottom") + theme_Publication() +scale_color_manual(values=my_palette,labels = c("1"= "3h-6h",
                                                                                                                 "2"="3h-12h",
                                                                                                                 "3"="3h-24h",
                                                                                                                 "4"="3h-2d",
                                                                                                                 "5"="3h-3d",
                                                                                                                 "6"="3h-4d",
                                                                                                                 "7"="3h-5d",
                                                                                                                 "8"="3h-6d",
                                                                                                                 "9"="3h-7d",
                                                                                                                 "10"="3h-8d",
                                                                                                                 "11"="3h-9d",
                                                                                                                 "12"="3h-10d",
                                                                                                                 "13"="3h-11d",
                                                                                                                 "14"="3h-12d",
                                                                                                                 "15"="3h-13d",
                                                                                                                 "16"="3h-14d",
                                                                                                                 "17"="3h-15d")) +
          scale_shape_manual(labels=c("M1","M2","M3","M4"),values=c(15,16,17,18)) 
        ggsave(plot[[i]],
               filename=paste0("reports_article/figures/third_figure/all_UMAP_with_all_neighbours/", sample,"_umap_",i,"_.eps",sep=""),width = 10.5,height=10.5) 
        
        p2=plot[[i]]+geom_path(aes(group=mouse))
        ggsave(p2,filename=paste0("reports_article/figures/third_figure/all_UMAP_with_all_neighbours/", sample,"_umap_path",i,"_.eps",sep=""),width = 10.5,height=10.5) 
        
        plot[[i]]=plot[[i]]+theme(legend.position = "none") +ggtitle(i)
        ##clustering with different methods##
        mapping.umap <- data.frame(
          id    = 1:NROW(digits.umap$layout),
          umap1  = digits.umap$layout[, 1],
          umap2  = digits.umap$layout[, 2],
          label = digits.labels)
        
        cluster_names<- as.data.frame((mapping.umap %>%  
                                         select(umap1, umap2) %>% NbClust(
                                           max.nc = 7, method =clustering_methods[c]))$Best.partition) ### determining the best number of clusters
        names(cluster_names)[1]="cluster"
        
        
        
        cluster[[i]]= mapping.umap %>%
          as_tibble() %>%
          mutate(cluster = cluster_names[1]) %>%
          ggplot(aes(umap1, umap2, color = factor(cluster$cluster),shape=label.sample)) +theme_Publication()+
          geom_point(size=4) +scale_shape_manual(labels=c("M1","M2","M3","M4"),values=c(15,16,17,18))
        
        ggsave(cluster[[i]],
               filename=paste0("reports_article/figures/third_figure/all_UMAP_with_all_neighbours/",sample,"_umap_cluster_",i,"_",clustering_methods[c],"_.eps",sep=""),width = 10.5,height=10.5) 
        
        cluster[[i]]=cluster[[i]]+theme(legend.position = "none") +ggtitle(i)
        
        n_neighbour=i
        cluster_type=clustering_methods[c]
        cluster_size=max(unique(cluster_names))
        sample_name=sample
        
        cluster_info[[i]]=data.frame(n_neighbour,cluster_type,cluster_size,sample_name)
        
      }  
      cluster_info2[[c]]=do.call(rbind,cluster_info)
      cowplot::plot_grid(plotlist = cluster,align = "hv") %>%  
        ggsave(filename =paste0("reports_article/figures/third_figure/UMAP_PANEL/",sample,"_",clustering_methods[c],"_all_umap_cluster_time_interval.eps",sep=""), width =80,height =50,limitsize = FALSE)
      
      
      
    }
    tf=do.call(rbind,cluster_info2)
    write_csv(tf,file =paste0("reports_article/figures/third_figure/UMAP_csv/",sample,"_cluster_size_umap.csv",sep="") )
    cowplot::plot_grid(plotlist = plot,align = "hv") %>%  
      ggsave(filename =paste0("reports_article/figures/third_figure/UMAP_PANEL/",sample,"_all_umap_time_interval.eps",sep=""), width =80,height =50,limitsize = FALSE)
    
  },error=function(e){})
  
}

get_all_umap_with_cluster_pairwise=function(df,sample){
  tryCatch({ 
    #names(df)[names(df) %in% c("X2")]<-c("Var2")
    my_palette <- colorRampPalette(c("purple","darkblue","#66a266","darkgreen","orange","darkorange"))(n =length(unique(df$Var2)))
    library(tidyverse)
    digits.labels = df[, c("Var2","sample","sample2","mouse2")]
    M_df=df[ , !names(df) %in% c("Var2","sample","sample2","mouse2")]
    clustering_methods=c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans")
    plot=list()
    cluster2=list()
    dim=seq_along(2:nrow(df))+1
    cluster_info2=list()
    for( c in seq_along(clustering_methods)){
      cluster=list()
      cluster_info=list()
      for(i in dim) { 
        set.seed(170)
        digits.umap = umap(M_df, n_components =2,n_neighbors=i,nn_method="annoy",n_trees="100",min_dist=0.1)
        #digits.umap = umap(M_spec_df, n_components =2,n_neighbours=15,min_dist=0.2) #15 is the standar
        layout <- digits.umap[["layout"]] 
        layout <- data.frame(layout) 
        final <- cbind(layout, df[, c("Var2","sample","sample2","mouse2")])
        colnames(final) <- c('umap1', 'umap2', 'phase',"mouse","sample2","mouse_number") 
        plot[[i]]=final %>% ggplot(aes(x = umap1, y = umap2,color = phase,shape=mouse_number))+
          geom_point(size=4) +
          theme(legend.position="bottom") + theme_Publication() +scale_color_manual(values=my_palette,labels = c("1"= "3h-6h",
                                                                                                                 "2"="3h-12h",
                                                                                                                 "3"="3h-24h",
                                                                                                                 "4"="3h-2d",
                                                                                                                 "5"="3h-3d",
                                                                                                                 "6"="3h-4d",
                                                                                                                 "7"="3h-5d",
                                                                                                                 "8"="3h-6d",
                                                                                                                 "9"="3h-7d",
                                                                                                                 "10"="3h-8d",
                                                                                                                 "11"="3h-9d",
                                                                                                                 "12"="3h-10d",
                                                                                                                 "13"="3h-11d",
                                                                                                                 "14"="3h-12d",
                                                                                                                 "15"="3h-13d",
                                                                                                                 "16"="3h-14d",
                                                                                                                 "17"="3h-15d")) +
          scale_shape_manual(labels=c("M1","M2","M3","M4"),values=c(15,16,17,18))+facet_wrap(~sample2)
        ggsave(plot[[i]],
               filename=paste0("reports_article/figures/third_figure/all_UMAP_with_all_neighbours_pairwise/", sample,"_umap_",i,"_.eps",sep=""),width = 20,height=10.5,device = cairo_ps) 
        
        p2=plot[[i]]+geom_path(aes(group=mouse))
        ggsave(p2,filename=paste0("reports_article/figures/third_figure/all_UMAP_with_all_neighbours_pairwise/", sample,"_umap_path",i,"_.eps",sep=""),width = 20,height=10.5,device = cairo_ps) 
        
        plot[[i]]=plot[[i]]+theme(legend.position = "none")+ggtitle(i)
        
        ####clustering on UMAP
        mapping.umap <- data.frame(
          id    = 1:NROW(digits.umap$layout),
          umap1  = digits.umap$layout[, 1],
          umap2  = digits.umap$layout[, 2],
          label = digits.labels)
        sample2=unique(mapping.umap$label.sample2)
        cluster_pam=list()
        for( s in seq_along(sample2)){
          cluster_pam[[s]]<- as.data.frame((mapping.umap[mapping.umap$label.sample2==sample2[s],]  %>%  
                                              select(umap1, umap2) %>% NbClust(max.nc = 7, method =clustering_methods[c]))$Best.partition)
        }  
        
        cluster_names=as.data.frame((do.call(rbind,cluster_pam)))
        names(cluster_names)[1]="cluster"
        
        cluster[[i]]= mapping.umap %>%
          as_tibble() %>%
          mutate(cluster = cluster_names[1]) %>%
          ggplot(aes(umap1, umap2, color = factor(cluster$cluster),shape=label.mouse2)) +theme_Publication()+
          geom_point(size=4)+facet_wrap(~label.sample2)+scale_shape_manual(labels=c("M1","M2","M3","M4"),values=c(15,16,17,18))
        
        ggsave(cluster[[i]],
               filename=paste0("reports_article/figures/third_figure/all_UMAP_with_all_neighbours_pairwise/",sample,"_umap_cluster_",i,"_.eps",sep=""),width = 20,height=10.5,device = cairo_ps) 
        
        cluster[[i]]=cluster[[i]]+theme(legend.position = "none")+ggtitle(i)
        
        n_neighbour=i
        cluster_type=clustering_methods[c]
        cluster_size=max(unique(cluster_names))
        sample_name=sample
        
        cluster_info[[i]]=data.frame(n_neighbour,cluster_type,cluster_size,sample_name)
        
        
      }
      
      cluster_info2[[c]]=do.call(rbind,cluster_info)
      cowplot::plot_grid(plotlist = cluster,align = "hv") %>%  
        ggsave(filename =paste0("reports_article/figures/third_figure/UMAP_PANEL/",sample,"_",clustering_methods[c],"_all_umap_cluster_time_interval.eps",sep=""), width =140,height =100,limitsize = FALSE)
    }
    tf=do.call(rbind,cluster_info2)
    write_csv(tf,file =paste0("reports_article/figures/third_figure/UMAP_csv/",sample,"_cluster_size_umap.csv",sep="") )
    cowplot::plot_grid(plotlist = plot,align = "hv") %>%  
      ggsave(filename =paste0("reports_article/figures/third_figure/UMAP_PANEL/",sample,"_all_umap_time_interval.eps",sep=""), width =140,height =100,limitsize = FALSE)
    
  },error=function(e){})
  
}



#######NM taxa data#####
NM1_16S_taxa <- read_delim("data/16S/taxa_long_format_old/NM1_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
NM2_16S_taxa <- read_delim("data/16S/taxa_long_format_old//NM2_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
NM3_16S_taxa <- read_delim("data/16S/taxa_long_format_old//NM3_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)
NM4_16S_taxa <- read_delim("data/16S/taxa_long_format_old//NM4_16S_taxa.tsv", 
                           "\t", escape_double = FALSE, col_types = cols(Time = col_skip()), trim_ws = TRUE)


##Remove the day befores###
NM1_16S_taxa=NM1_16S_taxa[-c(1,2,3), ]
##interpolate the day10 fro NM1
NM1_16S_taxa[13,]=(NM1_16S_taxa[12,]+NM1_16S_taxa[14,])/2


NM2_16S_taxa=NM2_16S_taxa[-c(1,2,3), ]
NM3_16S_taxa=NM3_16S_taxa[-c(1,2,3), ]
NM4_16S_taxa=NM4_16S_taxa[-c(1,2,3), ]



#get_all_umap_with_cluster_by_one_by(NM,"NM")



##NM_clones###

d_clone_NM1=get_time_derivative_GFM("NM1")
d_clone_NM2=get_time_derivative_GFM("NM2")
d_clone_NM3=get_time_derivative_GFM("NM3")
d_clone_NM4=get_time_derivative_GFM("NM4")


e_clone_NM1=get_time_varying_eigen_values(d_clone_NM1[[1]],d_clone_NM1[[2]],10,"NM1",d_clone_NM1[[3]])
e_clone_NM2=get_time_varying_eigen_values(d_clone_NM2[[1]],d_clone_NM2[[2]],10,"NM2",d_clone_NM2[[3]])
e_clone_NM3=get_time_varying_eigen_values(d_clone_NM3[[1]],d_clone_NM3[[2]],10,"NM3",d_clone_NM3[[3]])
e_clone_NM4=get_time_varying_eigen_values(d_clone_NM4[[1]],d_clone_NM4[[2]],10,"NM4",d_clone_NM4[[3]])


clone_NM=rbind(e_clone_NM1[[2]],e_clone_NM2[[2]],e_clone_NM3[[2]],e_clone_NM4[[2]])%>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 

####GFM###

d_GFM1=get_time_derivative_GFM("GFM1")
d_GFM2=get_time_derivative_GFM("GFM2")
d_GFM3=get_time_derivative_GFM("GFM3")
d_GFM4=get_time_derivative_GFM("GFM4")


e_GFM1=get_time_varying_eigen_values(d_GFM1[[1]],d_GFM1[[2]],10,"GFM1",d_GFM1[[3]])
e_GFM2=get_time_varying_eigen_values(d_GFM2[[1]],d_GFM2[[2]],10,"GFM2",d_GFM1[[3]])
e_GFM3=get_time_varying_eigen_values(d_GFM3[[1]],d_GFM3[[2]],10,"GFM3",d_GFM1[[3]])
e_GFM4=get_time_varying_eigen_values(d_GFM4[[1]],d_GFM4[[2]],10,"GFM4",d_GFM1[[3]])


GFM=rbind(e_GFM1[[2]],e_GFM2[[2]],e_GFM3[[2]],e_GFM4[[2]]) %>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 



GFM_without_GFM3=rbind(e_GFM1[[2]],e_GFM2[[2]],e_GFM4[[2]]) %>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() 





###pairwise_clone NM_GFM


NM_GFM=rbind(e_GFM1[[2]],e_GFM2[[2]],e_GFM3[[2]],e_GFM4[[2]],e_clone_NM1[[2]],e_clone_NM2[[2]],e_clone_NM3[[2]],e_clone_NM4[[2]]) %>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() %>% mutate(sample2=case_when(sample %in% c("NM1","NM2","NM3","NM4") ~ "NM",
                                                                                         sample %in% c("GFM1","GFM2","GFM3","GFM4") ~ "GFM"),
                                                                       mouse2=case_when(sample %in% c("NM1","GFM1")~ "1",
                                                                                        sample %in% c("NM2","GFM2")~ "2",
                                                                                        sample %in% c("NM3","GFM3") ~ "3",
                                                                                        sample %in% c("NM4","GFM4")~ "4"))


get_all_umap_with_cluster_pairwise(NM_GFM,"NM_GFM_clone")


###pairwise_clone NM_GFM without GFM3

NM_GFM=rbind(e_GFM1[[2]],e_GFM2[[2]],e_GFM4[[2]],e_clone_NM1[[2]],e_clone_NM2[[2]],e_clone_NM3[[2]],e_clone_NM4[[2]]) %>%
  ungroup %>%  mutate(Var2=as.factor(Var2) ) %>%  na.zero() %>% mutate(sample2=case_when(sample %in% c("NM1","NM2","NM3","NM4") ~ "NM",
                                                                                         sample %in% c("GFM1","GFM2","GFM3","GFM4") ~ "GFM"),
                                                                       mouse2=case_when(sample %in% c("NM1","GFM1")~ "1",
                                                                                        sample %in% c("NM2","GFM2")~ "2",
                                                                                        sample %in% c("NM3","GFM3") ~ "3",
                                                                                        sample %in% c("NM4","GFM4")~ "4"))


get_all_umap_with_cluster_pairwise(NM_GFM,"NM_GFM_clone_without_GFM3")





