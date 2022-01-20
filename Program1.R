##############################################################################################
# Simulation I of "Bayesian Semi-parametric Inference for Dynamic Marginal Structural Models"
# Authors: Daniel Rodriguez Duque, David A. Stephens, Erica E.M. Moodie, Marina Klein
##############################################################################################

rm(list=ls(all=TRUE))
library(gtools) #this library allows us to draw from the Dirichlet distribution

#-------------------------------------------------------------------------------
# Producing results Table 1
#
#  (a) Recommended settings for testing row 4 and 5 in Table 1 :
        B=1
        R=2
        #For better results set: R=100 (Runtime is around 2h with this setting) 
#                                      Note: 
#                                          1)Settings in original simulation:B=1; R=500  ~
#                                          2)With these settings, the program should finish running in several hours
#  (b) Recommended setting for testing row 9 and 10 in Table 1: 
         #B=100 
         #R=2
          #For better results set: R=100 (Runtime is several hours with this setting)
#                                      Note:
#                                           1) This takes several hours to run and gives a good approximation to the Bayesian results
#                                           2) Settings in original simulation: B=1000; R=500
#                                           3) These settings were originally run in parallel.It would take 
#                                              more than one week for this to finish running serially. 
#-------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------    
# Note that the data-generating mechanism is 
# Y=X1 -(-Th1+X1)*(A1_Opt-A1)-(-Th2+X2)*(A2_Opt-A2)+rand_comp
# Ignoring constants,the second product can be re-expressed as:
#                                                               A2 + A2_opt + X2*A2 + X2*A2_opt
# Marginalized versions of the terms above are referred to as:  T1 + T2     + T3    + T4 
#----------------------------------------------------------------------------------------------                                                         
  

expit=function(x) 1 / (1 + exp(-x))

#functions "T3" and "extra_term_compute" help perform marginalization
#over stage 2 covariates, so that outcome model conditional on only stage 1 covariates can be fit. 
T3=function(x,epsilon,Th2){
  #-----------------------------------------------------------------------------
  #This function helps marginalize the X2*A2 term in the data-generating mechanism.
  #when A2 is assigned according to a specific regime treat of X2>Th2.
  #-----------------------------------------------------------------------------
  A1=x[1]
  X1=x[2]
  to_return=(A1+0.5*X1+epsilon)* (Th2<(A1+0.5*X1+epsilon))
  return(mean(to_return))
}

extra_term_compute=function(A1,X1, Th1, Th2,extra=0){
  #-----------------------------------------------------------------------------
  #This function helps marginalize over 
  #A2 and X2*A2 when A2 is assigned according to a specific regime treat when X2>Th2
  #-----------------------------------------------------------------------------
  to=cbind(A1,X1)
  epsilon=rnorm(10000)
  T1_ev=pnorm(Th2-A1-0.5*X1,lower.tail=FALSE)
  T3_ev=apply(to,1,FUN=T3,epsilon=epsilon,Th2=Th2)
  if(extra==0){
    return(data.frame(T1_ev,T3_ev))
  }else{
    return(data.frame("T1_regime"=T1_ev, "T3_regime"=T3_ev))
  }
}


to_produce=function(Th1,Th2){
  #----------------------------
  # This function produces replicates datasets for analysis
  #----------------------------
  dati=array(0,dim=c(n,9,R))
  dimnames(dati)[[2]]=c("X1","A1","A1_Opt","X2","A2","A2_Opt","Y","T1_ev","T3_ev")
  for(i in 1:R){
    set.seed(i)
    rand_comp=rnorm(n, 0, 0.5)  # Random error in Y
    X1 = rnorm(n)  
    A1 = rbinom(n, 1, expit(1.5*X1))  
    A1_Opt = as.numeric(X1 >Th1) 
    X2=rnorm(n,A1+0.5*X1,sd=1)
    A2=rbinom(n, 1, expit(2*X2-0.5*A1))
    A2_Opt=as.numeric(X2>Th2)
    
    Y=X1 -(-Th1+X1)*(A1_Opt-A1)-(-Th2+X2)*(A2_Opt-A2)+rand_comp
    
    extract=extra_term_compute(A1,X1, Th1, Th2) #We are right away able to compute over the terms that marginalize over the optimal regime

    attach(extract)
    dati[,,i]=cbind(X1,A1,A1_Opt,X2,A2,A2_Opt,Y,T1_ev,T3_ev)
    detach(extract)
  }
  return(dati)
}


val_eval=function(index,theta){
  #-----------------------------------------------------------------------------
  #For a fixed regime of interest, with thresholds Theta1, Theta2
  #This function computes the value of the regime using both the double robust and IPW estimators
  #-----------------------------------------------------------------------------
  
  Theta1=as.numeric(theta[index,1])
  Theta2=as.numeric(theta[index,2])
  
  track_DR=matrix(0,R,B)
  track_IPW=matrix(0,R,B)
  
  for(i in 1:R){
    set.seed(i)
    attach(as.data.frame(dati[,,i])) #attaching the dataset corresponding to the ith replicate
    
    A1_regime = as.numeric(X1> Theta1)
    A2_regime = as.numeric(X2> Theta2)
    
    epsilon=rnorm(10000)
    to=cbind(A1,X1)
    #T2_ev & T4_ev marginalize terms with patients following regime X2>Theta2 and A1 varying as in the study population 
    T2_ev=pnorm(Theta2-A1-0.5*X1,lower.tail=FALSE)
    T4_ev=apply(to,1,FUN=T3,epsilon=epsilon,Th2=Theta2)

    to=cbind(A1_regime,X1)
    #T2_regime & T4_regime marginalize terms with patients following regime X2>Theta2 and A1 varying as in the study population 
    T2_regime=pnorm(Theta2-A1_regime-0.5*X1,lower.tail=FALSE)
    T4_regime=apply(to,1,FUN=T3,epsilon=epsilon,Th2=Theta2)
    extract_g=extra_term_compute(A1=A1_regime,X1, Th1, Th2,extra=1)

    DR_Vect4=vector()
    IPW=vector()

    for (o in 1:B){ ## posterior samples
      set.seed(o)
      richi=n*t(rdirichlet(1, rep(1,n)))
      if(B==1){
        richi=rep(1,n)
      }
      #################################################
      #Outcome Modeling with Correctly specified models
      #-----------------------------------------------#
      om2=lm(Y~-1+X1*A1_Opt + X1*A1+A2_Opt + A2 + A2_Opt:X2+ X2:A2, weights = richi)
      pseudo_outcome2=predict(om2, newdata=data.frame(X1=X1, A1_Opt=A1_Opt, A1=A1, A2_Opt=A2_Opt, A2=A2_regime,X2=X2))

      if(abs(Theta2-Th2)>0.001){ #can't use direct equals on floats
        om1=lm(pseudo_outcome2~X1*A1_Opt+X1*A1+T1_ev+T2_ev+T3_ev+T4_ev, weights=richi)
        attach(extract_g) # this dataset has T1_ev evaluated at the regime of interest
        pseudo_outcome1=predict(om1, newdata=data.frame(X1=X1, A1_Opt=A1_Opt, A1=A1_regime,T1_ev=T1_regime,T2_ev=T2_regime,T3_ev=T3_regime,T4_ev=T4_regime))
        detach(extract_g)
      }else{ #When the regime of interest is the same as the optimal regime, the fit needs to be different to prevent co-linearity. 
        om1=lm(pseudo_outcome2~X1*A1_Opt+X1*A1+T1_ev+T3_ev, weights=richi)
        attach(extract_g) # this dataset has T1_ev evaluated at the regime
        pseudo_outcome1=predict(om1, newdata=data.frame(X1=X1, A1_Opt=A1_Opt, A1=A1_regime,T1_ev=T1_regime,T3_ev=T3_regime))
        detach(extract_g)
      }

      phi2_34=predict(om2, newdata=data.frame(X1=X1, A1_Opt=A1_Opt, A1=A1_regime, A2_Opt=A2_Opt, A2=A2_regime,X2=X2))
      phi1_34=pseudo_outcome1

      ##################
      # Treatment Modeling
      ##################
      trt_model1_24 = glm(A1 ~-1+X1, family=quasibinomial,weights=richi)# Fit stage 1 propensity score model
      trt_model2_24 = glm(A2 ~-1+X2+A1, family=quasibinomial,weights=richi) # Fit stage 2 propensity score model
  
      C1_d = as.numeric(A1==A1_regime)  # Indicator of treatment matching regime recommendation
      C2_d = as.numeric(A2==A2_regime)    # Indicator of treatment matching regime recommendation
      
      #--------------
      p1_24 = predict(trt_model1_24, type="response")  # Estimated propensity score
      p2_24 = predict(trt_model2_24, type="response")  # Estimated propensity score
      
      p1_d_24 = A1*p1_24 + (1 - A1)*(1 - p1_24) 
      p2_d_24 = A2*p2_24 + (1 - A2)*(1 - p2_24) 
      
      AugY4=phi1_34+(C1_d)*(1/p1_d_24)*(phi2_34-phi1_34)+(C1_d)*(C2_d)*(1/(p2_d_24*p1_d_24))*(Y-phi2_34)
     
      #returning final values
      IPW[o]=mean(richi*C1_d*C2_d*Y/(p1_d_24*p2_d_24))
      DR_Vect4[o]=mean(richi*AugY4);
    }
    track_DR[i,]=DR_Vect4;DR_Vect4
    track_IPW[i,]=IPW;IPW
    detach(as.data.frame(dati[,,i]))
  }
  return(list("track_DR"=track_DR,"track_IPW"=track_IPW))
}


#sample size
n=500
#Thresholds corresponding to optimal regimes
Th1 = 0.4
Th2 = 0.8

#Generating Datasets for each replication
dati=to_produce(Th1,Th2)


#producing list of all regimes of interest
sequence1=seq(0,1,0.1)
sequence2=seq(0,1,0.1)
theta=as.matrix(expand.grid(sequence1,sequence2))


DR_Est=array(0,dim=c(R,B,nrow(theta)))
IPW_Est=array(0,dim=c(R,B,nrow(theta)))

#iterating through each regime of interest
for (i in 1:nrow(theta)){
     start_time <- Sys.time()
    outi=val_eval(seq(1,nrow(theta))[i],theta)
    DR_Est[,,i]=outi$track_DR
    IPW_Est[,,i]=outi$track_IPW
    end_time <- Sys.time()
    print(end_time - start_time)
}



#-------------------------------------------------------------------------------
# Producing row 4 and 5 of Table 1 (Frequentist)
# Note: Only run this when producing results for row 4 and 5, else go to next section
#-------------------------------------------------------------------------------
if(B==1){ #if frequentist case
      freq_max_val=function(x){
        index=which(x==max(x))
        return(c(x[index],index))
      }
      
      set.seed(1)
      popn=10000
      norm1=rnorm(popn,0.6,1);norm2=rnorm(popn,0.1,sd=1);rand_comp=rnorm(popn)
      rando=matrix(c(norm1,norm2,rand_comp),popn,3)
      
      freq_func=function(x){
        O1=0+norm1; A1=(O1>x[1])+0
        O2=A1+0.5*O1+norm2; A2=(O2>x[2])+0
        A1_Opt=(O1>Th1)+0;A2_Opt=(O2>Th2)+0
        Y=mean(O1 -(-Th1+O1)*(A1_Opt-A1)-(-Th2+O2)*(A2_Opt-A2) +rand_comp)
        return(Y)
      }
      
      IPW_max=apply(IPW_Est[,1,],MARGIN=1,FUN=freq_max_val)
      DR_max=apply(DR_Est[,1,],MARGIN=1,FUN=freq_max_val)
      
      
      IPW_Test_Pop=apply(theta[IPW_max[2,],],MARGIN=1,FUN=freq_func)
      DR_Test_Pop=apply(theta[IPW_max[2,],],MARGIN=1,FUN=freq_func)
      
    #Results Frequentist
      print(c(mean(theta[DR_max[2,],1]) ,sd(theta[DR_max[2,],1])
              ,mean(theta[DR_max[2,],2]),sd(theta[DR_max[2,],2])
              ,mean(DR_max[1,])         ,sd(DR_max[1,])
              ,mean(DR_Test_Pop)        ,sd(DR_Test_Pop)))
    #Results for test setting with R=2:   0.45000000 0.07071068 0.65000000 0.21213203 0.04940815 0.06435320 0.57008044 0.00140487
    #Results for test setting with R=100: 0.42000000 0.21272641 0.77200000 0.16334879 0.01496356 0.05899287 0.58425398 0.01575183
    #Reproducing Exact output from Original Simulation : 0.41540000 0.18245388 0.79300000 0.16230739 0.01832975 0.05637198 0.58654344 0.01429267
          
      print(c(mean(theta[IPW_max[2,],1]),sd(theta[IPW_max[2,],1])
          ,mean(theta[IPW_max[2,],2]) ,sd(theta[IPW_max[2,],2])
          ,mean(IPW_max[1,])          ,sd(IPW_max[1,])
          ,mean(IPW_Test_Pop)         ,sd(IPW_Test_Pop)))
     #Results for test setting with R=2: 0.60000000 0.28284271 0.55000000 0.35355339 0.01890263 0.04955105 0.57008044 0.00140487
     #Results for test setting with R=100: 0.45200000 0.22851298 0.73900000 0.21028840 0.02731124 0.06465542 0.58425398 0.01575183
     #Reproducing Exact output from Original Simulation : 0.44080000 0.20508408 0.74720000 0.20904440 0.03458791 0.06426116 0.58654344 0.01429267
 
}
#-------------------------------------------------------------------------------
# Producing row 9 and 10 of Table 1 (Bayesian) 
# Note: Only run the following when looking to produce results for row 9 and 10
#-------------------------------------------------------------------------------
if(B>1){ #if Bayesian case
        getmode <- function(v) {
          uniqv <- unique(v)
          uniqv[which.max(tabulate(match(v, uniqv)))]
        }
       maxy_compute=function(ary){
         #This function computes the optimal threshold and the optimal value in the posterior distribution
          maxy=matrix(NA,R,B)
          maxy_val=matrix(NA,R,B)
          for(i in 1:R){
            for (j in 1:B){
              maximizers=which(ary[i,j,]==max(ary[i,j,]))
              if(length(maximizers)>1){
                to_take=sample(maximizers,1)
                maxy[i,j]=to_take
                maxy_val[i,j]=ary[i,j,to_take]
              }else{
                maxy[i,j]=maximizers
                maxy_val[i,j]=ary[i,j,maximizers]
              }
            }
          }
          mode_theta_opt=apply(maxy,1,FUN=getmode)
          mean_value_opt=apply(maxy_val,1,FUN=mean)
          return(list("To_Return"=cbind(mode_theta_opt,mean_value_opt),"maxy"=maxy))
        }    
        
        DR_Result=maxy_compute(DR_Est)
        IPW_Result=maxy_compute(IPW_Est)
        
      
      ############################################
      # Evaluating Coverage
      #-------------------------------------------
      Coverage_Func=function(x){
          thetavect=theta[x,]
          q1=quantile(thetavect[,1], c(0.025,0.975))
          q2=quantile(thetavect[,2], c(0.025,0.975))
          return(c(q1,q2))
      }
      CoverageDR=t(apply(DR_Result$maxy,MARGIN=1,FUN=Coverage_Func))
      CoverageIPW=t(apply(IPW_Result$maxy,MARGIN=1,FUN=Coverage_Func))
      
      CTO1=mean(round(CoverageDR[,1],5)<=0.4 & 0.4<=round(CoverageDR[,2],5))       
      CTO2=mean(round(CoverageDR[,3],5)<=0.8 & 0.8<=round(CoverageDR[,4],5))
      
      CIPW1=mean(round(CoverageIPW[,1],5)<=0.4 & 0.4<=round(CoverageIPW[,2],5))         
      CIPW2=mean(round(CoverageIPW[,3],5)<=0.8 & 0.8<=round(CoverageIPW[,4],5))
      
      ##############################################
      # Evaluating Value in Test Population
      #---------------------------------------------
      set.seed(1)
      popn=10000
      norm1=rnorm(popn,0.6,1);norm2=rnorm(popn,0.1,sd=1);rand_comp=rnorm(popn)
      rando=matrix(c(norm1,norm2,rand_comp),popn,3)
      per_person=function(x, thetavect){
        O1=0+x[1]
        A1=(O1>thetavect[,1])+0
        O2=A1+0.5*O1+x[2]
        A2=(O2>thetavect[,2])+0
        A1_Opt=(O1>Th1)+0
        A2_Opt=(O2>Th2)+0
        Y=O1 -(-Th1+O1)*(A1_Opt-A1) - (-Th2+O2)*(A2_Opt-A2) + x[3]
        return(Y)
      }
      Test_Pop_Func=function(x){
        thetavect=theta[x,]
        Y=apply(rando,MARGIN=1,FUN=per_person,thetavect=thetavect)
        return(rowMeans(Y))
      }
      TestPopDR=apply(DR_Result$maxy,MARGIN=1,FUN=Test_Pop_Func)
      TestPopIPW=apply(IPW_Result$maxy,MARGIN=1,FUN=Test_Pop_Func)
      
      TO=colMeans(TestPopDR)
      IPW=colMeans(TestPopIPW)
      
      
#Results Bayesian:      
      print(c(mean(theta[DR_Result$To_Return[,1],1]) ,sd(theta[DR_Result$To_Return[,1],1])
        ,mean(theta[DR_Result$To_Return[,1],2]),sd(theta[DR_Result$To_Return[,1],2])
        ,mean(DR_Result$To_Return[,2])         ,sd(DR_Result$To_Return[,2])
        ,mean(TO)                               ,sd(TO)))
        # Output for test settings with R=2:    0.1500000000 0.2121320344 0.8000000000 0.0000000000 0.0667681108 0.0593558435 0.5840619731 0.0006471301
        # Output from test settings with R=100: 0.418000000 0.226247384 0.783000000 0.168807463 0.024845061 0.058314006 0.586234696 0.008620513
      
      print(c(mean(theta[IPW_Result$To_Return[,1],1]) ,sd(theta[IPW_Result$To_Return[,1],1])
        ,mean(theta[IPW_Result$To_Return[,1],2]),sd(theta[IPW_Result$To_Return[,1],2])
        ,mean(IPW_Result$To_Return[,2])         ,sd(IPW_Result$To_Return[,2])
        ,mean(IPW)                              ,sd(IPW)))
       # Output for test settings with R=2:    0.400000000 0.565685425 0.900000000 0.141421356 0.049597558 0.051587701 0.567557786 0.004665416
       # Output from test settings with R=100: 0.48300000 0.24290124 0.76400000 0.21391965 0.04830347 0.06433437 0.57973910 0.00959519     
}     


