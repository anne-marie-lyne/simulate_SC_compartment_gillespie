compute_index = function(spNum, mut_save) {
  
  #function to compute linearity index for final time point
  
  #find all remaining mutations
  fin_mut = unique(unlist(mut_save[spNum>0]))
  
  if(length(fin_mut)>0){
    mut_save = lapply(mut_save,function(x) c(0,x)) #add marker for unmutated cells
    fin_mut = unique(unlist(mut_save[spNum>0]))
    #make intersection and union matrices
    int = matrix(0,nrow=length(fin_mut),ncol=length(fin_mut))
    uni = matrix(0,nrow=length(fin_mut),ncol=length(fin_mut))
    #iterate through all pairwise mutations
    for(i in 1:(length(fin_mut)-1)){
      for(j in (i+1):length(fin_mut)){
        #check if any clones contain mutation i and j
        logi = sapply(mut_save,function(x) (is.element(fin_mut[i],x) && is.element(fin_mut[j],x)))
        #intersection is sum of number of cells containing both i and j
        int[i,j] = sum(spNum[logi])
        #union is sum of cells containing i or j
        logi = sapply(mut_save,function(x) (is.element(fin_mut[i],x) || is.element(fin_mut[j],x)))
        uni[i,j] = sum(spNum[logi])
      }
    }
    index = sum(int)/sum(uni)
  }else{
    index = NA
  }
}