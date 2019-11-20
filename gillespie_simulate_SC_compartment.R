
#code to simulate stem cell compartment using gillespie algorithm (Gillespie, 1977).
#Stem cells can undergo three reactions:
# SC -> Diff + Diff
# SC -> SC + Diff
# SC -> SC + SC

#in selection scenario, if simulated reaction is self renewal (SR, reaction 3), 
#cells with mutations more likely to be chosen.
#The more mutations a stem cell has, the more its SR probability increases

iniNum = 200 #starting number of SC
p_mut = 0.005 #mutation rate per division
tsk = 1 #used when running multiple instances of same code
wd = "/Users/alyne/ownCloud/anne-marie/clonal_evo_paper_writing/code_for_paper/"
simTime = 100 #number of time periods simulation runs for (cells divide on avergae once per time period)
probAS = 1/3 #probability of asymmetric division (prob SR = prob symm diff = (1-probAS)/2)
adv = 0.1 #fitness parameter. 0 for neutral mutations. 

source("compute_index_function.R")

set.seed(tsk) # set seed differently for each task

setwd(wd)
#parameters for Gillespie algorithm determining which cells divide and when
###########################################################################
#code SC = 1, Diff = 2
stoich_mat = matrix(c(-1,2,0,1,1,0),nrow=3,ncol=2,byrow=TRUE)
#define probabilities of reactions
probSC = c((1-probAS)/2,probAS,(1-probAS)/2)
rateSC = 1 #number of reactions on average per time period per cell
cSC = probSC * rateSC #rate of each reaction type
sum_cSC = sum(cSC) 
###########################################################################

nreps = 20 #number of independent repeats to do
nmut_tot = numeric(nreps) #total number of mutations occurring during simulation
nmut_ret = numeric(nreps) #number of mutations which end up in SC compartment (should be ~0.5*nmut_tot)
ind_save2 = numeric(nreps) #jaccard index
ncell = numeric(nreps) #number of remaining cells at end of simulation
nevent = numeric(nreps) #total number of reactions
nmut_final = numeric(nreps) #number of unique mutations in stem cell compartment at final time point
T_save = list() #times of reactions
T_mut = list() #times of mutations which end up in SC compartment
Mut_save = list() #clones defined by mutations accumulated in SC compartment
SpNum = list() #number of cells in the above populations
ATtributes = list() #required for Muller plots
Population.data = list() #ditto
for(k in 1:nreps){
  print(k)
  spNum = iniNum #species (number of cells per clone) count
  t = 0
  t_save = t
  t_mut = t
  
  mut_save = list()
  mut_save[[1]] = vector('numeric')
  mut = 1 #mutations which end up in SC compartment
  mut_lost = 0 #mutations which end up in differentiated compartment
  mut_mat = matrix(iniNum,nrow=1,ncol=1)
  while(t<simTime){
    #sim time tau and which reaction
    lambda = sum(spNum)*sum_cSC
    tau = rexp(1,rate=lambda)
    reac = which(rmultinom(1, 1, probSC)==1)
    
    #choose which population to update
    #if reaction is self renewal, then give higher probability to mutants
    if(length(spNum)==1){
      #only unmutated cells
      div_prob = 1
    }else{
      if(reac==3){
        #SR, count number of mutations in each clone
        rj = 1 + adv*sapply(mut_save,length) #fitness advantage depends on number of mutatins each clone has
        div_prob = rj*spNum/sum(rj*spNum) #probability depends on fitness advantage plus clone size
      }else{
        div_prob = spNum #probability proportional to clone size
      }
    }
    pop = which(rmultinom(1, 1, div_prob)==1)
    #update overall cell numbers
    rd = runif(1) #mutation?
    if(rd>=p_mut || reac==1){
      #either no mutation or reac 1 is SC -> diff + diff (i.e. mutation is lost anyway)
      #update population
      spNum[pop] = spNum[pop] + stoich_mat[reac,1]
      if(rd<p_mut){ 
        #must be reaction 1, mutation lost
        mut_lost = mut_lost + 1
      }
      if(reac==1){
        if(ncol(mut_mat)%%2){
          #update mutation matrix of cells with each mutation
          mut_mat = cbind(mut_mat,numeric(length=nrow(mut_mat)) )
          mut_mat[,ncol(mut_mat)] = spNum
          t_mut = append(t_mut,t+tau)
        }
      }
    }else{
      if(reac==2 && rd<p_mut){
        # SC -> SC + Diff and mutation has occurred
        if(runif(1)<0.5){
          #mutation stays in population
          spNum = c(spNum,1) #start new clone of size one
          spNum[pop] = spNum[pop] - 1 #number in original population decreases by 1
          mut_save[[mut+1]] = c(mut_save[[pop]],mut) #new cell has original mutations plus new one
          mut = mut + 1
          mut_mat = cbind(mut_mat,numeric(length=nrow(mut_mat)) )
          mut_mat = rbind(mut_mat,numeric(length=ncol(mut_mat)) )
          mut_mat[,ncol(mut_mat)] = spNum #update mutation matrix
          t_mut = append(t_mut,t+tau)
        }else{
          #mutation lost from population
          mut_lost = mut_lost + 1
        }
      }else{
        if(reac==3 && rd<p_mut){
          # SC -> SC + SC and mutation
          spNum = c(spNum,1) #start new clone of size 1
          mut_save[[mut+1]] = c(mut_save[[pop]],mut) #clone has previous mutations plus new one
          mut = mut + 1
          mut_mat = cbind(mut_mat,numeric(length=nrow(mut_mat)) )
          mut_mat = rbind(mut_mat,numeric(length=ncol(mut_mat)) )
          mut_mat[,ncol(mut_mat)] = spNum
          t_mut = append(t_mut,t+tau)
        }
      }
    }
    
    #update time
    t = t + tau
    t_save = append(t_save,t)
    #end if there are no SC or Diff cells left
    if(sum(spNum)==0){
      break
    }
  }
  #save output for rep k
  T_save[[k]] = t_save
  T_mut[[k]] = t_mut
  SpNum[[k]] = spNum
  Mut_save[[k]] = mut_save
  nmut_tot[k] = mut - 1 + mut_lost
  nmut_ret[k] = mut - 1
  ncell[k] = sum(spNum)
  nevent[k] = length(t_save)
  nmut_final[k] = length(unique(unlist(mut_save[spNum>0])))
  
  #compute index
  index = compute_index(spNum, mut_save)
  ind_save2[k] = index
  
  #make structures required for muller plots
  Attributes = matrix(ncol=2,nrow=length(mut_save))
  Attributes[,1] = as.character(seq(1,length(mut_save))-1)
  Attributes[,2] = as.character(rep(0,length(mut_save)))
  Attributes[1,2] = NA
  colnames(Attributes) = c("names","parents")
  lens = sapply(mut_save,length)
  if(max(lens)>1){
    for(i in 2:max(lens)){
      Attributes[cbind(sapply(mut_save[which(lens==i)],function(x) x[i])+1,rep(2,length(which(lens==i))))] = 
        as.character(sapply(mut_save[which(lens==i)],function(x) x[i-1]))
    }
  }
  population.data = mut_mat
  rownames(population.data) = Attributes[,1]
  colnames(population.data) = t_mut
  ATtributes[[k]] = Attributes
  Population.data[[k]] = population.data
}

save(nmut_tot,nmut_ret,ind_save2,ncell,nmut_final,nevent,iniNum,p_mut,
     T_save,T_mut,Mut_save,SpNum,ATtributes,Population.data,
     file=paste0("results_gillespie_N_",iniNum,"_pmut_",
                 gsub(".","",p_mut,fixed=T),"_selection_",
                 gsub(".","",adv,fixed=T),".RData"))

#nmut_tot is all mutations including those lost
