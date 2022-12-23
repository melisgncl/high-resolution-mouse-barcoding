######according to choosen UMAP clustering calculate the fitness ###

source("~/Desktop/0106_mice_barcode_data/src_article//visualization/0_config/0_config.R")
c_list=c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
`%notin%` <- Negate(`%in%`)##opposite off in function

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
  Time=seq(1:nrow(clusters.loess))
  series_whole_loess=clusters.loess
  series_whole_loess$time=NULL
  deriv_series_loess=sapply(series_whole_loess, function(y)  tslist(t((diff(y,lag=2)/diff(Time/10,lag=2)))))
  options(scipen=0)
  return(list(tslist(t(series_whole_loess)),deriv_series_loess,Time))
  
}


##give the choosen neighbour and cluster number###

get_chosen_umap=function(df,n,method){
  M_df=df[-c(1,2)]
  digits.labels = df[, c("Var2","sample")]
  set.seed(170)
  digits.umap = umap(M_df, n_components =2,n_neighbors=n,nn_method="annoy",n_trees="100",min_dist=0.1)
  #digits.umap = umap(M_spec_df, n_components =2,n_neighbours=15,min_dist=0.2) #15 is the standar
  layout <- digits.umap[["layout"]] 
  layout <- data.frame(layout) 
  final <- cbind(layout, df[, c("Var2","sample")])
  colnames(final) <- c('umap1', 'umap2', 'phase',"mouse")
  
  cluster_names<- as.data.frame((final %>%  
                                   select(umap1, umap2) %>% NbClust(
                                     max.nc = 7, method =method))$Best.partition) ### determining the best number of clusters
  
  names(cluster_names)[1]="cluster"
  
  final=final%>%
    mutate(cluster = cluster_names[[1]])
  
  
  return(final)
  
}

####plot_choosen_umap##

#####plot_derivatives_with_normal_series###

plot_deriv_with_normal=function(series_whole_loess,series_whole_loess){
      plot_data_i=melt(do.call(cbind,series_whole_loess))
      plot_data_i$model="series"
      plot_data_d=melt(do.call(cbind,series_whole_loess))
      plot_data_d$model="derivative"
      plot_data=rbind(plot_data_i,plot_data_d)
      p=ggplot(plot_data)+ geom_line(aes(Var1/10,value,color=model,group=model),size=1) +facet_wrap(~Var2,nrow = 1) +
        theme_Publication() + theme(strip.background = element_blank())+scale_color_manual(values=c("purple","black"),label=c("dx/dt","Frequency")) +xlab("Time")+ylab("")
      ggsave(p,filename = paste0("reports_article//figures/fitness_graphs/", sample, "_series_fitting.eps"),height = 6,width = 45)
}

get_selection_coeff = function(umap,sample,deriv_series_loess){
              umap=umap[umap$mouse==sample,]
              clusters.loess=read_csv(file= paste("reports_article//figures/clustering/",sample,"/",sample,"_clustered_loess_log10.csv",sep=""))
              ###change first hours to day like
              cluster_timepoints2=list()
              cluster=unique(umap$cluster)
              ##match phase with its time###
              for (i in seq_along(cluster)){
                cluster_timepoints=list()
                phase=as.numeric(as.character(umap[umap$cluster == cluster[i], ]$phase))
                for (k in  seq_along(phase)){
                  cluster_timepoints[[k]]=seq(from=phase[k]*10, to=(phase[k]+1)*10-1, by=1)
                }
                cluster_timepoints2[[i]]=as.data.frame(sort(melt(do.call(rbind,cluster_timepoints))$value))
                cluster_timepoints2[[i]]$cluster=cluster[i]
              }
              
              cluster=do.call(rbind,cluster_timepoints2)
              names(cluster)[1]="value"
              
              hist= as.data.frame(deriv_series_loess) %>% mutate(Time=seq(from=min(clusters.loess$time), to=max(clusters.loess$time-2), by=1))%>% 
                mutate(phase = case_when(Time %in% (as.integer(cluster[cluster$cluster == 1, ]$value)) ~ "First phase",
                                         Time %in% (as.integer(cluster[cluster$cluster == 2, ]$value)) ~ "Second phase",
                                         Time %in% (as.integer(cluster[cluster$cluster == 3, ]$value))~ "Third phase"))%>%  
                pivot_longer(cols=(1:ncol( as.data.frame(deriv_series_loess))))%>% 
                mutate(Time = case_when(Time >= 31 ~ as.numeric(Time)-10,
                                        Time < 31  ~ as.numeric(Time/3))) %>%  filter(!grepl("\\.", Time))%>%  ###removes floats numbers
                mutate(sample=sample,type=case_when(name %in% c_list ~ "Clone",
                                                    name %notin% c_list ~ "Species"))
              
        return(hist)
}
              
              
              
####NM fitness graphs
umap=readRDS("reports_article/figures/fitness_graphs/Umap_gavaged_65n")

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


#####time derivatives###

d_NM1=get_time_derivative_gavaged(NM1_16S_taxa,"NM1")
d_NM2=get_time_derivative_gavaged(NM2_16S_taxa,"NM2")
d_NM3=get_time_derivative_gavaged(NM3_16S_taxa,"NM3")
d_NM4=get_time_derivative_gavaged(NM4_16S_taxa,"NM4")


####get clustered fitness##

NM1_coef=get_selection_coeff(umap,"NM1",d_NM1[[2]])
NM2_coef=get_selection_coeff(umap,"NM2",d_NM2[[2]])
NM3_coef=get_selection_coeff(umap,"NM3",d_NM3[[2]])
NM4_coef=get_selection_coeff(umap,"NM4",d_NM4[[2]]) 



####
hist_statistic= rbind(NM1_coef,NM2_coef,NM3_coef,NM4_coef) %>%
  mutate(definition = case_when(name =="C1" ~ "Dominant clone",
                                name %in% c("C3","C2")  ~ "C2 and C3",
                                name =="Lachnospiraceae"  ~ "Lachnospiraceae",
                                name %in% c("C4","C5","C6","C7","C8","C9","C10")  ~ "Low segregated clones",
                                name =="Lactobacillaceae"  ~ "Lactobacillaceae",
                                TRUE ~ "Other species")) %>%     
  group_by(phase,definition,Time) %>% 
  #mutate(mean=mean(value), sd=sd(value))
  summarise(value=value,sample=sample) %>% 
  mutate(definition=factor(definition,levels=c("Dominant clone","Low segregated clones","C2 and C3","Lachnospiraceae","Lactobacillaceae","Other species")))



  p2=ggplot(hist_statistic,aes(x=value,color=phase)) + 
  #geom_step(aes(y=value2),size=1.2,position = "identity") + 
  stat_ecdf(geom = "step",size=2,position = "identity") +
  #geom_ribbon(aes(ymin = lower ,ymax = upper),alpha = 0.2) +
  facet_wrap(~as.factor(definition), scales = "free_x",nrow =2) + 
  theme(strip.background = element_blank()) + theme(aspect.ratio = 1)+
  scale_color_manual(name = "", 
                     values = c("First phase"="darkblue",
                                "Second phase"="darkgreen",
                                "Third phase"="darkorange")) + theme_Publication() +
  theme(strip.background = element_blank()) +xlab("") +ylab("") +geom_hline(yintercept = 0.5,linetype="dashed", col = 'black',size=2) +
  geom_vline(xintercept = 0,linetype="dashed", col = 'black',size=2) +xlim(c(-1,1))


ggsave(p2,filename ="reports_article//figures/fitness_graphs//fitness_NM.eps",height = 10,width = 30,limitsize = FALSE)
 
p2=p2+ facet_grid(as.factor(definition)~sample, scales = "free_x")
ggsave(p2,filename ="reports_article//figures/fitness_graphs//fitness_NM2.eps",height = 20,width = 60,limitsize = FALSE)

data=hist_statistic[c(1,2,4)]

data.ecdf = hist_statistic %>%
  mutate(value = ceiling(value / 0.05) * 0.05) %>%
  group_by(phase,definition,value) %>%
  summarize(num.runs = n()) %>%
  ungroup() %>%
  group_by(phase, definition) %>%
  arrange(phase,definition,value) %>%
  mutate(prob = cumsum(num.runs / sum(num.runs)))


p2=ggplot(data.ecdf) +
  geom_line((aes(value, prob,color=phase,group=phase)),size=2) +facet_wrap(~as.factor(definition),nrow =2) + 
  theme(strip.background = element_blank()) + theme(aspect.ratio = 1)+
  scale_color_manual(name = "", 
                     values = c("First phase"="darkblue",
                                "Second phase"="darkgreen",
                                "Third phase"="darkorange")) + theme_Publication() +
  theme(strip.background = element_blank()) +xlab("") +ylab("") +geom_hline(yintercept = 0.5,linetype="dashed", col = 'black',size=2) +
  geom_vline(xintercept = 0,linetype="dashed", col = 'black',size=2) +xlim(c(-1,1))

ggsave(p2,filename ="reports_article//figures/fitness_graphs//fitness_NM2_0.05.eps",height = 20,width = 40,limitsize = FALSE)

###GFM fitness###

get_selection_coeff_GF = function(umap,sample,deriv_series_loess){
  clusters.loess=read_csv(file= paste("reports_article//figures/clustering/",sample,"/",sample,"_clustered_loess_log10.csv",sep=""))
  umap=umap[umap$mouse==sample,]
  ###change first hours to day like
  cluster_timepoints2=list()
  cluster=unique(umap$cluster)
  ##match phase with its time###
  for (i in seq_along(cluster)){
    cluster_timepoints=list()
    phase=as.numeric(as.character(umap[umap$cluster == cluster[i], ]$phase))
    for (k in  seq_along(phase)){
      cluster_timepoints[[k]]=seq(from=phase[k]*10, to=(phase[k]+1)*10-1, by=1)
    }
    cluster_timepoints2[[i]]=as.data.frame(sort(melt(do.call(rbind,cluster_timepoints))$value))
    cluster_timepoints2[[i]]$cluster=cluster[i]
  }
  
  cluster=do.call(rbind,cluster_timepoints2)
  names(cluster)[1]="value"
  hist= as.data.frame(deriv_series_loess) %>% mutate(Time=seq(from=min(clusters.loess$time), to=max(clusters.loess$time-2), by=1))%>% 
    mutate(phase = case_when(Time %in% (as.integer(cluster[cluster$cluster == 1, ]$value)) ~ "Second phase",
                             Time %in% (as.integer(cluster[cluster$cluster == 2, ]$value)) ~ "First phase", 
                             TRUE  ~ "Third mouse")) %>%  
    pivot_longer(cols=(1:ncol( as.data.frame(deriv_series_loess))))%>% 
    mutate(Time = case_when(Time >= 31 ~ as.numeric(Time)-10,
                            Time < 31  ~ as.numeric(Time/3))) %>%  filter(!grepl("\\.", Time))  ###removes floats numbers
  return(hist)
}



umap=readRDS("reports_article/figures/fitness_graphs/Umap_GFM_21n")

d_GFM1=get_time_derivative_GFM("GFM1")
d_GFM2=get_time_derivative_GFM("GFM2")
d_GFM3=get_time_derivative_GFM("GFM3")
d_GFM4=get_time_derivative_GFM("GFM4")

GF1=get_selection_coeff_GF(umap,"GFM1",d_GFM1[[2]])
GF2=get_selection_coeff_GF(umap,"GFM2",d_GFM2[[2]])
GF3=get_selection_coeff_GF(umap,"GFM3",d_GFM3[[2]])
GF4=get_selection_coeff_GF(umap,"GFM4",d_GFM4[[2]])


hist_statistic= rbind(GF1,GF2,GF4)


p=ggplot(hist_statistic,aes(x=value,color=phase)) + 
  #geom_step(aes(y=value2),size=1.2,position = "identity") + 
  stat_ecdf(geom = "step",size=2,position = "identity") +
  #geom_ribbon(aes(ymin = lower ,ymax = upper),alpha = 0.2) +
  facet_wrap(~as.factor(name), scales = "free_x",nrow =3) + 
  theme(strip.background = element_blank()) + theme(aspect.ratio = 1)+
  scale_color_manual(name = "", 
                     values = c("First phase"="darkblue",
                                "Second phase"="darkgreen"))+
  theme_Publication() +
  theme(strip.background = element_blank()) +xlab("") +ylab("") +geom_hline(yintercept = 0.5,linetype="dashed", col = 'black',size=2) +
  geom_vline(xintercept = 0,linetype="dashed", col = 'black',size=2) +xlim(c(-2,2))

ggsave(p,filename ="reports_article//figures/fitness_graphs//GF_fitness_suplementary.eps",height = 15,width = 20,limitsize = FALSE)


data=hist_statistic[c(3,2,4)]

data.ecdf = hist_statistic %>%
  mutate(value = ceiling(value / 0.05) * 0.05) %>%
  group_by(phase,name,value) %>%
  summarize(num.runs = n()) %>%
  ungroup() %>%
  group_by(phase, name) %>%
  arrange(phase,name,value) %>%
  mutate(prob = cumsum(num.runs / sum(num.runs)))


p2=ggplot(data.ecdf) +
  geom_line((aes(value, prob,color=phase,group=phase)),size=2) +
  facet_wrap(~as.factor(name),nrow =3) + 
  theme(strip.background = element_blank()) + theme(aspect.ratio = 1)+
  scale_color_manual(name = "", 
                     values = c("First phase"="darkblue",
                                "Second phase"="darkgreen"))+
  theme_Publication() +
  theme(strip.background = element_blank()) +xlab("") +ylab("") +geom_hline(yintercept = 0.5,linetype="dashed", col = 'black',size=2) +
  geom_vline(xintercept = 0,linetype="dashed", col = 'black',size=2) +xlim(c(-2,2))

ggsave(p2,filename ="reports_article//figures/fitness_graphs//GF_fitness_suplementary_2.eps",height = 15,width = 20,limitsize = FALSE)
