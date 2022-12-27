
## R code
##  Melis Gencel

#######################    READ ME!   #####################################
##here we are plotting all the mixing index values that is get from from all cut off with resuffled data##


source("~/Desktop/mouse_barcoding_last/src/visualization/0_config/0_config.R")
library(plyr)
library(ggbeeswarm)

#####load all the stat data####
#mydir = "reports/figures/mixing_index/stat_csv/"
#myfiles = list.files(path=mydir, pattern="*stat.tsv", full.names=TRUE)
#df_stat_all = ldply(myfiles, read_tsv)
#detach("package:plyr")

###the file tales too long load so I already saved them as a csv to use it
#saveRDS(df_stat_all,"reports/figures/mixing_index/stat_vall_stat_new.csv")
df_stat_all=readRDS("reports/figures/mixing_index/stat_vall_stat_new.csv")

df_rm=df_stat_all[df_stat_all$Condition %in% c("rm_rm","nc_rm","rm_gf","nc_gf"),]
df_rm=df_rm  %>% mutate(Condition=case_when(Condition =="rm_rm" ~ "1",
                                            Condition =="nc_rm" ~ "2",
                                            Condition =="rm_gf" ~ "3",
                                            Condition =="nc_gf" ~ "4"))

df_im=df_stat_all[df_stat_all$Condition %in% c("im_im","nc_im","im_gf","nc_gf"),]
df_im=df_im  %>% mutate(Condition=case_when(Condition =="im_im" ~ "5",
                                            Condition =="nc_im" ~ "6",
                                            Condition =="im_gf" ~ "7",
                                            Condition =="nc_gf" ~ "8"))



df=rbind(df_rm,df_im)


p=ggplot() +geom_violin(data=df,aes(x=Condition, y=1-ks_stat,color=Condition)) +theme_Publication() +
  scale_color_manual(name = "", 
                     values = c("1" = "darkorange", 
                                "3"="darkgreen",
                                "4" ="#31a354",
                                "2"="#6B298C",
                                "5" = "#00469C",
                                "7"="darkgreen",
                                "8" ="#31a354",
                                "6"="#6B298C") ,
                     guide="none") +
  ylab("Mixing Index") + xlab("") +
  geom_quasirandom(data=df,aes(x=Condition, y=1-ks_stat,color=Condition),size=0.2) +
  ylab("Mixing Index") + xlab("")+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.x=element_text(color = "black", size=14
                                 , angle=30, vjust=.8, hjust=0.8)) 



ggsave(p,filename = "reports/figures/mixing_index/Mixing_index_new.eps",width =12,height = 7.25)



