#script to make mullerplots of output from gillespie_simulate_SC_compartment.R

#load required packages
library(MullerPlot)

iniNum = 200 #starting number of SC
p_mut = 0.005 #mutation rate per division
wd = "/Users/alyne/ownCloud/anne-marie/clonal_evo_paper_writing/code_for_paper/"
adv = 0.1 #fitness parameter. 0 for neutral mutations. 

setwd(wd)

#load output from gillespie simulation
load(paste0("results_gillespie_N_",iniNum,"_pmut_",gsub(".","",p_mut,fixed=T),
            "_selection_",gsub(".","",adv,fixed=T),".RData"))

#set directory for figure output
fig_dir=paste0(wd,"figures/")

#which figures to plot
to_plot = c(1,4) #plot reps 1 and 4

for(i in to_plot){
  
  #name figure appropriately
  png(paste(fig_dir,"mullerplot_","N_",iniNum,"_pmut_",gsub(".","",p_mut,fixed=T),
            "_selection_",gsub(".","",adv,fixed=T),
            "_rep_",i,".png",sep=""),width = 1000, height = 750)
  par(mar=c(7,8,4,2)+0.5,mgp=c(5, 2, 0))
  Attributes = ATtributes[[i]]
  if(T_save[[i]][length(T_save[[i]])] != T_mut[[i]][length(T_mut[[i]])]){
    Population.Data = cbind(Population.data[[i]],SpNum[[i]])
    colnames(Population.Data) = c(colnames(Population.data[[i]]),T_save[[i]][length(T_save[[i]])])
  }else{
    Population.Data = Population.data[[i]]
  }
  Muller.plot(attributes = Attributes, population.data = Population.Data,
              data.method = "table", time.interval.method = "linear",
              cex.main=3, cex.axis=3, cex.lab=3, ylab = "Frequency", xlab = "Generation",
              main = paste0(ncell[i]," cells, Index = ",format(round(ind_save2[i],3),nsmall = 3)))
  dev.off()
  
}


