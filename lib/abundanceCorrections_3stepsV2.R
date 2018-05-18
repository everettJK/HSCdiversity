###########Correction of the data based on a particular cutoff of abundance, confusion matrix and PA=initial blood cell count##################
Correction_CutData_new <- function(data,cutoff, confusion, PA,PB){
  
  
  ###cut the data in order to consider only the high abundant IS
  Cut_data <- data[data$cellCount>cutoff,]
  Cut_data_bis <- Cut_data[,3:7]
  
  #########################Correction of the Cut-data for unbalanced sampling with the weight P(dataset)/P(postsort)##############################################
  PC_Cut_data = colSums(Cut_data_bis)/ sum(colSums(Cut_data_bis)) #P(dataset)
  weights_Ini_data = PC_Cut_data/PB
  
  WeightHT1 = rbind(PA,PC_Cut_data,weights_Ini_data)
  row.names(WeightHT1)=c("PA","PC","WeightHT1")
  
  M=matrix(0,ncol=5,nrow=5)
  diag(M)=weights_Ini_data^(-1) 
  
  Cut_data_corrected2=  round(as.matrix(Cut_data_bis)%*% M,2) # If we do not wan all the row to be identical , we can removed the round step
  
  #########################Correction2 of the Cut-data for unbalanced sampling with the weight + confusion matrix ###########
  
  Confusion_prof_row = prop.table(confusion,1)
  
  Cut_data_corrected2C = data.frame()
  Cut_data_corrected2C = Cut_data[,1:2]
  Cut_data_corrected2C[, 3:7] = as.matrix(Cut_data_corrected2) %*% Confusion_prof_row
  colnames(Cut_data_corrected2C) = colnames(Cut_data)[1:7]
  
  
  #########################Correction3 of the Cut-data for unbalanced sampling with the weight + confusion matrix + weight P(datacorr)/P(blood)###########
  
  
  
  Pdataset_corr=colSums( Cut_data_corrected2C[, 3:7])/ sum( Cut_data_corrected2C[, 3:7]) #calcul of the "new" pdataset
  # weights = PB/PA 
  weights = Pdataset_corr/PA 
  
  
  
  
  WeightHT2 = rbind(PA,PC_Cut_data,weights)
  row.names(WeightHT1)=c("PA","PC","WeightHT2")
  
  
  M=matrix(0,ncol=5,nrow=5)
  diag(M)=weights^(-1) 
  
  
  Cut_data_corrected3 = Cut_data_corrected2C[,1:2]
  Cut_data_corrected3[,3:7]=data.frame(as.matrix(Cut_data_corrected2C[,3:7])%*% M) 
  
  
  
  #####Generation of the dataframe in proportion from Cut_data_corrected3####################################
  DF1 = data.frame()
  DF1 = Cut_data_corrected3[,1:2]
  DF1[,3:7] = Cut_data_corrected3[,3:7]/rowSums(Cut_data_corrected3[,3:7])
  colnames(DF1)=colnames(Cut_data_corrected3)
  
  return( list(Cut_data =Cut_data, Cut_data_bis=Cut_data_bis,Cut_data_corrected2C=Cut_data_corrected2C,Cut_data_corrected3=Cut_data_corrected3,WeightHT1=WeightHT1,WeightHT2=WeightHT2 ))
}


