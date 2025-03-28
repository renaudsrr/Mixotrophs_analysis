
Shannon = function(data){
  H = numeric(length(data$Sample))
  
  for (lake in 1:length(data$Sample)){
    sum_H = 0  
    
    for (genus in 2:ncol(data)){  # Boucle sur all genus
      biomass_i = data[lake, genus]  # Biomasse sp i
      
      if (biomass_i > 0){  # Eviter log(0)
        p_i = biomass_i / summary_data_phyto$biomass_tot[lake]
        sum_H = sum_H + (p_i * log(p_i))}
    }
    H[lake] = -sum_H
  }
  return(H)
}

Simpson = function(data){
  D = numeric(length(data$Sample)) 
  
  for (lake in 1:length(data$Sample)){
    sum_D = 0
    
    for (genus in 2:ncol(data)){  # Boucle sur all genus
      biomass_i = data[lake, genus]  # Biomasse sp i
      
      if (biomass_i > 0){  # evite calcul + 0
        p_i = biomass_i / summary_data_phyto$biomass_tot[lake]
        sum_D = sum_D + p_i**2}
    }
    D[lake] = 1-sum_D # indide de Simpson complémentaire + intuitif
  }
  return(D)
}
