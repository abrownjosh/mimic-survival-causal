
#############################################################
#Functions to get pseudo observations or multiple imputations
#############################################################



pacman::p_load(gdata)

#Get MI
sum_function=function(X)
{
  temp=array(X,c(length(X),length(X)))
  lowerTriangle(temp,diag=F)=0
  apply(temp,2,sum,na.rm=TRUE)
}

getMI <- function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest){
  # @param n; number of subjects
  # @param M; numer of imputes desired
  # @param ID; vector of SSIDs
  # @param Z; covariates used for regression
  # @param X;  vector of observed time points for each subject (exact same data as that which is contained in T_t)
  # @param T_t; matrix of time points for each subject observed
  # @param delta; a vector indicating whether or not the person had an event.
  # @param tp.interest; the time point we are imputing to.
  
  
  # tau+little.t[length(little.t)] = time point of interest
  
  #Identify patients needing imputation
  impute_subjects=ID[round(X,digits = 1)<tp.interest & 
                       apply(delta,1,function(x){x[sum(!is.na(x))]})==0]  #Only impute those who loss of follow up
  n_impute_subjects=length(impute_subjects)
  imputed_data=T_t
  imputed_data[impute_subjects,]=imputed_data[impute_subjects,]*delta[impute_subjects,] #anywhere 0 needs to be imputed
  
  #this list will store our imputed datasets.
  imputed_data_list=rep(list(imputed_data),M)
  
  #Identify times when we completely observed subjects!
  for(k in 1:n_impute_subjects)
  {
    #Identify risk set for patient k
    patient_id=impute_subjects[k]
    largest_possible_t=0 # when we start follow up! The original code was written for a more complex case...
    risk_set=rep(0,n)
    j=0
    epsilon=0.01
    while(sum(risk_set)<5)
    {
      # iterate epsilon to keep a larger risk set.
      epsilon=epsilon+j*0.001
      # collect people into the risk set based on three criteria:
      risk_set[X>X[patient_id] &  # have to survive beyond censoring time for a person
                 T_t[,1]> T_t[patient_id,1]  # had to be at risk at time t
               
                # THE BELOW IS WHERE THE DISTANCE METRIC GETS PLACED!
               & abs(t(t(fitted.values)%*%t(Z)) -  #difference between fitted values for all subjects!
                       sum(fitted.values*Z[patient_id,]))<epsilon # imputed subject's fitted value.
               # the fitted covariates must be within epsilon for the person we are imputing
               ]=1
      # print(risk_set)
      j=j+1
      print(epsilon)
    }
    
    # gets risk set.
    X_riskset=T_t[risk_set==1,1]
    delta_riskset=delta[risk_set==1,1]
    
    n_riskset=length(X_riskset)
    time_riskset=sort(unique(X_riskset))
    n_time_riskset=length(time_riskset)
    
    #Estimate survival - NA estimator
    temp3=array(X_riskset,c(n_riskset,n_time_riskset)) # gets number of unique survival times at each failure time
    temp4=t(array(time_riskset,c(n_time_riskset,n_riskset)))
    delta_array_riskset=array(delta_riskset,c(n_riskset,n_time_riskset))
    dN_T_riskset=array(as.numeric(temp3==temp4 & delta_array_riskset==1),
                       c(n_riskset,n_time_riskset))
    Y_riskset=array(as.numeric(temp3>=temp4),c(n_riskset,n_time_riskset))
    
    denominator=apply(Y_riskset,2,sum,na.rm=TRUE)
    temp5=apply(t(dN_T_riskset)/denominator,1,sum,na.rm=TRUE)
    CH_Rk=sum_function(temp5)
    surv_Rk=exp(-CH_Rk)
    time_Rk=time_riskset
    ### this is another estimator of the survival curve! ###
    #plot(time_Rk,surv_Rk,type="s",col="red")
    
    
    #Sampling a valid impute
    l=1
    subject_imputes=NULL
    while(length(subject_imputes)<M & l<=100)
    {
      u=runif(1)
      if(u<=min(surv_Rk)){ # then we can consider that the subject has survived.
        subject_imputes_one=largest_possible_t+ tp.interest
        }
      if(u>min(surv_Rk)){ # if u is greater than the minimum survival probability, we
        # need to find the corresponding correct survival time.
        impute_index=(length(time_Rk)-apply(array(surv_Rk,c(length(surv_Rk),1))<=t(array(u,c(1,length(surv_Rk)))),2,sum)) + 1
        imputation_times=time_Rk[impute_index]
        
        temp=array(X_riskset,c(n_riskset,1))==t(array(time_Rk[impute_index],c(1,n_riskset)))
        risk_set_imputes_subjects=array(ID[risk_set==1],c(length(ID[risk_set==1]),1))[temp]
        # get the residuals for the imputed subjects. This creates random error for our imputed time. #### THIS NEEDS TO INCLUDE YOUR FITTED VALUES
        residuals=(apply(cbind(T_t[risk_set_imputes_subjects,1],rep( tp.interest,1)),1,min))-sum(fitted.values*Z[risk_set_imputes_subjects,]) 
        
        ###                     NEEDS FITTED VALUES!
        subject_imputes_one=(sum(fitted.values*Z[patient_id,]) + residuals) + largest_possible_t
      }
      if(subject_imputes_one[1]>X[patient_id]) # checks to make sure that we have a valid impute.
      {
        subject_imputes=c(subject_imputes,subject_imputes_one[1])
      }
      print(paste("Imputing",l))
      l=l+1
    }
    #If don't have enough imputes after 100 sampling use the largest censoring time
    # subject_imputes <- subject_imputes[1:M]
    # subject_imputes <- c(subject_imputes, rep(largest_possible_t+ tp.interest,M-length(subject_imputes)))
    
    #Put imputes in each data
    for (m in 1:M){
      # for (ti in 1:which(little.t==largest_possible_t)){
        if (imputed_data_list[[m]][patient_id,1]==0){
          imputed_data_list[[m]][patient_id,1]=subject_imputes[m]-0
        }

    }
    print(paste("Patient",k))
  }
  return(imputed_data_list)
}








