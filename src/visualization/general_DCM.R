##########general DCM
library(reshape2)
library(ggplot2)
library(dplyr)
library(Metrics) 
library(kernlab)
library(dplyr)
library(tidyr) # For data manipulation


apply_loess_to_all_species <- function(data, spans = seq(0.2, 0.9, by = 0.1), n) {
  # n being the interpolation length
  # Initialize an empty data frame to store results
  results <- data_frame(Species = character(), Time = numeric(), Smoothed_Population = numeric())
  
  # Create interpolation points
  xx <- seq(from = min(data$Time), to = max(data$Time), length.out = n)
  
  # Iterate over each species
  unique_species <- unique(data$Species)
  for (species_name in unique_species) {
    # Filter data for the current species
    species_data <- filter(data, Species == species_name)
    
    # Find best span for current species
    best_span <- find_best_loess_span(species_data, spans)
    
    # Perform LOESS smoothing with the best span
    loess_fit <- loess(Population ~ Time, data = species_data, span = best_span)
    smoothed_values <- predict(loess_fit, newdata = data.frame(Time = xx), se = FALSE)
    
    # Append results for the current species
    species_results <- data_frame(Species = species_name, 
                                  Time = xx, 
                                  Smoothed_Population = smoothed_values)
    results <- bind_rows(results, species_results)
  }
  
  return(results)
}

# Helper function to find best span for LOESS
find_best_loess_span <- function(data, spans) {
  # Initialize an object to store errors for each span
  errors <- sapply(spans, function(span) {
    fit <- loess(Population ~ Time, data = data, span = span)
    predictions <- predict(fit, newdata = data)
    mse <- mean((data$Population - predictions)^2) # Calculate mean squared error
    print(mse)
    return(mse)
  })
  
  # Return the span with the minimum error
  best_span <- spans[which.min(errors)]
  print(best_span)
  return(best_span)
}

### calculate the  ż for time series
get_derivative<-function(data,n_zero,n_loess,interpolate=TRUE,log=TRUE){
  "n_zero specifies the maximum allowable number of zeros in a column for it to be included in the analysis"
  data <- data %>%
      group_by(Species) %>%
      mutate(ZeroCount = sum(Population == 0)) %>%
      filter(ZeroCount < n_zero) %>%
      select(-ZeroCount)
  ## n is the loess range to interpolate you can use the given loees range or you can give your own 
  if (interpolate) {
     df=apply_loess_to_all_species(data,n=n_loess)}
    ##if it is log scale give a pseduo number accoridng to resulution
  if (log) {
    df$Smoothed_Population=log10(df$Smoothed_Population+0.000001)}
   
  Time=unique(df$Time)
  complete_loess_series <- tslist(t(df %>%
    pivot_wider(
      id_cols = Time,
      names_from = Species,
      values_from = Smoothed_Population
    ) %>%select(-Time)
    ))
  # Calculate time derivative using central differences
  derivative_series=sapply(complete_loess_series, function(y)  tslist(t((diff(y,lag=2)/diff(Time,lag=2)))))
  return(list(complete_loess_series,derivative_series,Time))
  }
  
#calculate time dependent jacobian 
##you can change the start point and the window size
calculate_time_dependent_jacobian=function(complete_loess_series,derivative_series,window,start=1){
  options(scipen=0)
  total <- nrow(as.data.table(derivative_series))
  # Generate sequence points from the first time point to the last, adjusting for the window size
  # Windows are defined to span from one day to two days
  spots <- seq(from=start, to=(total), by=window)
  if (spots[length(spots)] != total) {
    spots <- c(spots, total)
  }
  print(spots)
  # Adjust the window by reducing 2 steps to center it appropriately
  # Calculate the community interaction matrix for each window
  
  jacobian=list()  
  for(k in 2:length(spots)){
    jacobian[[k]] <- sapply(1:ncol(as.data.table(derivative_series)), function(i) {
      sapply(1:ncol(as.data.table(complete_loess_series)), function(j){
        resTmp <- cov(derivative_series[[i]][(spots[1]:spots[k])],complete_loess_series[[j]][(spots[1]:spots[k])]) 
        
      })
    })
  }
  jacobian <- jacobian[-1]
  # Assign row and column names from derivative_series to each matrix in results
  for(i in seq_along(jacobian)) {
    rownames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
    colnames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
  }
  
  
  return(jacobian)
  
}


# Calculate eigenvalues and eigenvectors for each jacobian matrix  
get_eigen_values=function(jacobian,sample,time_increase=FALSE){
  
  
  eigen_values=sapply(jacobian, function(y) eigen(y, symmetric= FALSE)$values)
  # Extract real and imaginary parts of eigenvalues
  real_eigen_value=sapply(as.data.table(eigen_values), function(y) Re(y))
  im_eigen_value=sapply(as.data.table(eigen_values), function(y) Im(y))
  
  eigen_vector=lapply(jacobian, function(y) eigen(y, symmetric= FALSE)$vector)
  
  #Calculate the magnitude of eigenvalues for sorting
  conjugate= sqrt(real_eigen_value^2+im_eigen_value^2)
  
  # Organize data into data frames for further manipulation
  real_eigen_value=melt(real_eigen_value)  %>% mutate(eigen="R")
  im_eigen_value=melt(im_eigen_value)  %>% mutate(eigen="I")
  conjugate_value=melt(conjugate)  %>% mutate(eigen="C")
  
  ordered_eigen_values =rbind(real_eigen_value,im_eigen_value,conjugate_value) %>% spread(eigen,value)  %>% 
    mutate(Var2=sub("V", "", Var2) ,Var1=as.factor(Var1),Var2=(as.numeric(Var2)))  %>% group_by(Var2) %>%
    arrange(desc(C),.by_group = TRUE) %>%
    mutate(Var3=as.factor(row_number()))  
  
  ###prepare the data for dimension reduction by using time dependent eigenvalue 
  time_dependent_eigenspace= ordered_eigen_values %>% mutate(sample=sample) %>% 
    select(-C,-Var1)  %>% pivot_longer(cols = -c(Var2,Var3,sample), names_to = "variable", values_to = "value") %>% 
    group_by(Var2) %>% mutate(Var4=paste0(variable,"_",Var3)) %>%  select(-Var3,-variable) %>% spread(Var4,value) 
  ##if the first time not exist
  if (time_increase) {
    time_dependent_eigenspace$Var2 <- time_dependent_eigenspace$Var2 +1
  }
  
  
  return(list(ordered_eigen_values,time_dependent_eigenspace))
  

}
    
plot_eigen_value_evolution=function(all_eigen_value,sample){
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
    ##save the all plots
    #ggsave(plots_eigen_value[[i]],filename = paste("reports//figures/DCM/eigen_values/",sample, "_eigen_value_",eigen_space[i],"_time_interval.eps", sep=""),width =10.5,height =7.5, limitsize = FALSE,device = cairo_ps )
  }
  ##save them as a panel
  cowplot::plot_grid(plotlist = plots_eigen_value,align = "hv",ncol = 10) %>%  
    ggsave(filename = paste("reports//figures/DCM2/eigen_values_panel/",sample, "_all_eigen_value_time_interval.eps", sep=""),width =50,height =20, limitsize = FALSE,device = cairo_ps)
}

plot_kpca=function(df,sigma,sample){
  set.seed(172)
  df_numeric <- df[, -which(names(df) %in% c("Var2", "sample"))]
  ##center the dataset for PCA but not scale
  df_numeric=as.data.frame(scale(df_numeric , center = TRUE, scale = FALSE))
  kpca_result <- kpca(~., data=df_numeric, kernel="rbfdot", kpar=list(sigma=sigma), features = 2)
  # Evaluate the model performance
  
  final <- cbind(kpca_result@rotated, df[, c("Var2","sample")])
  
  colnames(final) <- c('kpca1', 'kpca2', 'phase',"mouse")
  p=ggplot(final) +
    geom_point(size=2, aes(color = phase,x = phase, y = kpca1)) +
    geom_line(size=1.2, aes(x = phase, y = kpca1),color="red") +
    geom_point(size=4, aes(color = phase,x = phase, y = kpca2)) +
    geom_line(size=1.2, aes(x = phase, y = kpca2),color="blue") +
    theme_Publication() 
  
  return(p)
}

######################First run a 4 species glv##################
##or you can upload your time series###
set.seed(3) # for reproducibility
r <- rep(1, 4)
A <- -matrix(c(10, 6, 12,0,
               14, 10, 2,0,
               8, 18, 10,0,
               2,8,14,10), 4,4, byrow = TRUE)/100

##with this matrix you will get this data use src/7_DCM/4_plot_DCM_glv.R function
long_data <- reshape2::melt(integrate_GLV(r, A, x0,100,1)[[1]],id.vars=c("time"))
names(long_data)=c("Time","Species","Population")


###normalize data if necesary
long_data <- long_data %>%
  group_by(Time) %>%
  mutate(
    Total_Daily_Population = sum(Population),
    Population = Population / Total_Daily_Population
  ) %>% select(-Total_Daily_Population)%>%
  ungroup() 

##plot dynamics
ggplot(long_data, aes(x = Time, y =Population, color = Species)) +
  geom_line() +
  labs(title = "Population Trends of 10 Species Over 20 Days",
       x = "Date",
       y = "Population") +
  theme_minimal() +
  theme(legend.title = element_blank())


##if necessary interpolate data interpolate=TRUE,
##if you wanna perform ar log scale log=TRUE

####give the range of interpolation also you can gave your own range
max.range = max(long_data$Time)-min(long_data$Time)
loess.range = (max.range*10)+1

##you can give your loess range 
derivative_tf <- get_derivative(long_data,n_zero=,7,n_loess =loess.range,interpolate=TRUE,log=TRUE)

###calculate jacobian starting point is 1 but you can arrange that
jacobian_tf <-calculate_time_dependent_jacobian(derivative_tf[[1]],derivative_tf[[2]],10,start=1)

##### Calculate eigenvalues and eigenvectors for each jacobian matrix
eigen_tf=get_eigen_values(jacobian_tf,"sample_DCM")


###you can plot your eigen_values
plot_eigen_value_evolution(eigen_tf[[1]],"sample_DCM")

##you can apply pca and check your dyanmical shift of time depedent jacobain
plot_kpca(eigen_tf[[2]],2)




