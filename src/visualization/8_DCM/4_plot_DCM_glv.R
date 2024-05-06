###  Melis Gencel
### DCM analysis
#######################    READ ME!   #####################################

##### In here we analyzesimulated glv data

#######################    ^^^^^^^^   #####################################

source("~/Desktop/mouse_barcoding_last/src/visualization/8_DCM/1_DCM.R")
#################################################im########################################################
#################################################im########################################################
###functions for GLV simulation
library(deSolve) # integrate ODEs
library(tidyverse) # plotting and wrangling




GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-8] <- 0 # prevent numerical problems
    dxdt <- x * (r + A %*% x)
    return(list(dxdt))
  })
}
# function to plot output
plot_ODE_output <- function(out,n_species){
  out <- as.data.frame(out)
  #out <- out/ rowSums(out)
  colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
  out <- as_tibble(out) %>% gather(species, density, -time)  %>% group_by(time) #%>%
  #mutate(normalized_density = density / sum(density))
  pl <- ggplot(data = out) + 
    aes(x = time, y = density, colour = species) + 
    geom_line(size=1.2) +theme_Publication() 
  show(pl)
  #ggsave(pl,filename=paste0("reports/figures/GLV/",n_species,"_species_dynamics.eps"))
  return(list(out,pl))
}

integrate_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.01,n_species){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- as.data.frame(ode(y = x0, times = times, 
                           func = GLV, parms = parameters, maxsteps = 10^9,
                           method = "ode45"))
  print(out)
  out <- plot_ODE_output(out,n_species)
  out[[1]]=spread(out[[1]], species, density)
  return(list(out[[1]],out[[2]]))
}



get_time_derivative_glv=function(clusters.loess){
  clusters.loess$time=NULL
  Time=seq(1:nrow(clusters.loess))
  series_whole_loess=clusters.loess
  series_whole_loess$time=NULL
  deriv_series_loess=sapply(series_whole_loess, function(y)  tslist(t((diff(y,lag=2)/diff(Time,lag=2)))))
  options(scipen=0)
  return(list(tslist(t(series_whole_loess)),deriv_series_loess,Time))
  
}


calculate_time_dependent_jacobian_glv=function(complete_loess_series,derivative_series,window,sample){
  options(scipen=0)
  total <- nrow(as.data.table(complete_loess_series))
  spots <- seq(from=2, to=(total), by=window)
  num_rows <- nrow(as.data.table(derivative_series))
  # Set the last element of spots to be equal to num_rows
  spots[length(spots)] <- num_rows
  # Remove elements larger than num_rows
  spots <- spots[spots <= num_rows]
  print(spots)
  # Adjust the window by reducing 2 steps to center it appropriately
  # Calculate the community interaction matrix for each window
  
  jacobian=list()  
  for(k in 1:length(spots)){
    jacobian[[k]] <- sapply(1:ncol(as.data.table(derivative_series)), function(i) {
      sapply(1:ncol(as.data.table(complete_loess_series)), function(j){
        resTmp <- cov(derivative_series[[i]][(1:spots[k])],complete_loess_series[[j]][(1:spots[k])]) 
        
      })
    })
  }
  
  # Assign row and column names from derivative_series to each matrix in results
  for(i in seq_along(jacobian)) {
    rownames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
    colnames(jacobian[[i]]) <- colnames(as.data.table(derivative_series))
  }
  
  
  return(jacobian)
  
}



plot_correlation_and_pvalues=function(df,A,steptime){
  library(reshape2)
  melted_list <- lapply(df, melt)
  # Add a column for list number
  melted_list <- lapply(seq_along(melted_list), function(i) {
    melted_data <- melted_list[[i]]
    melted_data$list_number <- i
    melted_data
  })
  # Combine all melted dataframes into one
  combined_df <- do.call(rbind, melted_list)
  # Rename the columns for clarity
  colnames(combined_df) <- c("V1", "V2", "Value","List_Number")
  combined_df$sp=paste0(combined_df$V2,"_",combined_df$V1)
  #####M calculated by (Dx)A
  #M=melt(t(diag(res_2_2$average) %*%  A_2))
  M=reshape2::melt(t(A))
  M$sp=paste0(M$Var2,"_",M$Var1)
  M=M[c(3,4)]
  ####corellation
  corr <- combined_df %>%
    group_by(List_Number) %>%
    summarize(correlation = cor.test(M$value, Value)$estimate)
  corr$steptime=steptime
  pvalue= combined_df %>%
    group_by(List_Number) %>%
    summarize(pvalue = cor.test(M$value, Value)$p.value)
  pvalue$steptime=steptime
  ggplot()+geom_line(data=corr,aes(List_Number,correlation)) +theme_Publication() +xlab("Progressive Time interval")
  ggplot()+geom_line(data=pvalue,aes(List_Number,-log10(pvalue))) +theme_Publication() +xlab("Progressive Time interval")
  
  return(list(corr,pvalue,M,combined_df))
  
}

###############3 species ossilcation sytems ############################


##################DATA#################
set.seed(3) # for reproducibility
r_3 <- rep(1, 3)
A_3 <- -matrix(c(10, 6, 12, 
                 14, 10, 2, 
                 8, 18, 10), 3, 3, byrow = TRUE)/100
# check the existence of feasible equilibrium
print(solve(A_3, -r_3)) # feasible
x0_3=c(10,10,10)

##solve the sytem

res_3 <- integrate_GLV(r_3, A_3, x0_3,100,1)

d_3=get_time_derivative_glv(res_3[[1]])

j3=calculate_time_dependent_jacobian_glv(d_3[[1]],d_3[[2]],1,"three_sp")

e3=get_eigen_values(j3,"three_sp")

####correlation of jacobian with Aij

k=plot_correlation_and_pvalues(j3,A_3,1)
###plot_cor and pvalues

ylim.prim <- c(0, 1)   # in this example, precipitation
ylim.sec <- c(min(-log10(k[[2]]$pvalue)), max(-log10(k[[2]]$pvalue)))    

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


ggplot() +
  geom_line(data = k[[1]], aes(x = List_Number, y = correlation), size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dotdash", color = "red", size = 1) +
  geom_line(data = k[[2]], aes(x = List_Number, y = a+(-log10(pvalue)*b)), size = 1.2, color = "darkblue") +
  
  geom_hline(yintercept = 1.3*b , linetype = "dashed", color = "darkblue") +  # Example threshold for p-value
  
  scale_y_continuous(
    name = "Correlation",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - a)/b, name = "-log10(p-value)")
  ) +
  xlab("Time") +
  theme_Publication()

