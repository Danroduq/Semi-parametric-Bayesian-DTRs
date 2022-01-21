rm(list=ls(all=TRUE))
library(gtools)

##############################################################################################
# Simulation II of "Bayesian Semi-parametric Inference for Dynamic Marginal Structural Models"
# Authors: Daniel Rodriguez Duque, David A. Stephens, Erica E.M. Moodie, Marina B. Klein
##############################################################################################

#---------------------------------------------------------------------------------
# Produce rows 4 and 5, columns 3,4,5 of Table 2
# Recommended settings for testing:
        n=500                                 #Sample size
        B=1                                   #Number of posterior draws; setting it to 1 gives the frequentist fit
        R=2                                   #Number of replicates
          #For better results set R=100       #This takes a few hours to run
        MM=1000                               #Note: this should not be too small
#       Note: 
#           1) MM=50000; R=500 is what was used in original simulations.
#           2) This program can also run Bayesian with B=500 simulation but this requires powerful computing system.
#---------------------------------------------------------------------------------

  expit=function(x) 1 / (1 + exp(-x))
        
  All_c3=function(X,epsilon1_1,epsilon2_1,Theta1,Theta2,Theta3,c3_data){
    #---------------------------------------------------------------------------
    # This is a helper function for C3_compute(); C2_compute(); c1_compute()
    # This computes all  the stage-specific terms that are marginalized 
    # conditional on information up the previous stage, with treatment being assigned
    # according to theta1, theta2, theta3
    #---------------------------------------------------------------------------
    
    A3=c3_data[X,1]
    X31=c3_data[X,2]
    X32=c3_data[X,3]
    XX=c3_data[X,4]
    
    A4_regime=(Theta1*(0.2*A3+0.1*X31+epsilon1_1)+
                 Theta2*(0.5*A3+0.1*X32+epsilon2_1)+Theta3*XX)>0.5
    #term T9_c2
    T12_c3=-0.5*mean(A4_regime)
    #term T7B_c2
    T10B_c3=0.1*mean(XX*A4_regime)
    #term T7_c2
    T10_c3=0.5*mean((0.2*A3+0.1*X31+epsilon1_1)*A4_regime)
    #term T8_c2
    T11_c3=0.5*mean((0.5*A3+0.1*X32+epsilon2_1)*A4_regime)
    return(c(T10B_c3,T10_c3,T11_c3,T12_c3))
  }
  All_c2=function(X,epsilon1_2, epsilon2_2, epsilon1_1,epsilon2_1, Theta1, Theta2, Theta3, c2_data){
    #---------------------------------------------------------------------------
    # This is a helper function for C2_compute(); c1_compute()
    # This computes all  the stage-specific terms that are marginalized 
    # conditional on information from two stages prior, with treatment being assigned
    # according to theta1, theta2, theta3
    #---------------------------------------------------------------------------
    A2=c2_data[X,1]
    X21=c2_data[X,2]
    X22=c2_data[X,3]
    XX=c2_data[X,4]
    
    X31_c2=0.2*A2+0.1*X21+epsilon1_2
    X32_c2=0.5*A2+0.1*X22+epsilon2_2
    A3_regime=(Theta1*X31_c2+Theta2*X32_c2+Theta3*XX)>0.5
    X41_c2=(0.2*A3_regime+0.1*X31_c2+epsilon1_1)
    X42_c2=(0.5*A3_regime+0.1*X32_c2+epsilon2_1)
    A4_regime=(Theta1*X41_c2+Theta2*X42_c2+Theta3*XX)>0.5
    
    #term T9_c1
    T9_c2=-0.5*mean(A4_regime)
    #term T7B_c1
    T7B_c2=0.1*mean(XX*A4_regime)
    #term T7_c1
    T7_c2=0.5*mean(X41_c2*A4_regime)
    #term T8_c1
    T8_c2=0.5*mean(X42_c2*A4_regime)
    
    return(c(T7B_c2,T7_c2,T8_c2,T9_c2))
  }

  All_c1=function(X,  epsilon1_3, epsilon2_3, 
                      epsilon1_2, epsilon2_2, 
                      epsilon1_1, epsilon2_1, Theta1, Theta2, Theta3, c1_data){
    #---------------------------------------------------------------------------
    # # This is a helper function for C1_compute()
    # This computes all  the stage-specific terms that are marginalized 
    # conditional on information from three stages prior , with treatment being assigned
    # according to theta1, theta2, theta3
    #---------------------------------------------------------------------------
    A1=c1_data[X,1]
    X11=c1_data[X,2]
    X12=c1_data[X,3]
    XX=c1_data[X,4]
    
    X21_c1=(0.2*A1+0.1*X11+epsilon1_3)
    X22_c1=(0.5*A1+0.1*X12+epsilon2_3)
    A2_regime=(Theta1*X21_c1+Theta2*X22_c1+Theta3*XX)>0.5
    X31_c1=(0.2*A2_regime+0.1*X21_c1+epsilon1_2)
    X32_c1=(0.5*A2_regime+0.1*X22_c1+epsilon2_2)
    A3_regime=(Theta1*X31_c1+Theta2*X32_c1+Theta3*XX)>0.5
    X41_c1=(0.2*A3_regime+0.1*X31_c1+epsilon1_1)
    X42_c1=(0.5*A3_regime+0.1*X32_c1+epsilon2_1)
    A4_regime=(Theta1*X41_c1+Theta2*X42_c1+Theta3*XX)>0.5
    
    #term T9_c1
    T6_c1=-0.5*mean(A4_regime)
    #term T7B_c1
    T4B_c1=0.1*mean(XX*A4_regime)
    #term T7_c1
    T4_c1=0.5*mean(X41_c1*A4_regime)
    #term T8_c1
    T5_c1=0.5*mean(X42_c1*A4_regime)
    
    return(c(T4B_c1,T4_c1,T5_c1,T6_c1))
  }

#----------------------------------
  c3_compute=function(c3_data,regime=0, 
                      Theta1,Theta2,Theta3,
                      Th1,Th2, Th3,
                      epsilon31, epsilon32, 
                      All_c3r_opt_dtr){
   #------------------------------------------------------------------------------------------
   #This function computes all relevant marginalized terms conditional on stage III information
   #------------------------------------------------------------------------------------------
    All_c3r=t(sapply(X=seq(1,n),FUN=All_c3,epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                     Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,c3_data=c3_data))
    colnames(All_c3r)=c("T10B_c3r","T10_c3r","T11_c3r","T12_c3r")
    
    if(regime==1){ #when it's a term that involves the optimal regime, but where current stage treatment is assigned
                   #according to another regime.
      All_c3r_opt=t(sapply(X=seq(1,n),FUN=All_c3,epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                           Theta1=Th1,Theta2=Th2,Theta3=Th3,c3_data=c3_data))
      colnames(All_c3r_opt)=c("T10B_c3r_opt","T10_c3r_opt","T11_c3r_opt","T12_c3r_opt")
    }else{
      All_c3r_opt=All_c3r_opt_dtr  #These terms have already been computed with the Marginals_At_Opt() function
                                   #We return them again just for symmetry in the function
    }
    #there are 8 terms to return
    return(cbind(All_c3r,All_c3r_opt))
  }

  c2_compute=function(c2_data,regime=0, 
                      Theta1,Theta2,Theta3,
                      Th1,Th2, Th3,
                      epsilon21, epsilon22,
                      epsilon31, epsilon32, 
                      All_c2r_opt_dtr,
                      All9_c2r_opt_dtr){
    #------------------------------------------------------------------------------------------
    #This function computes all relevant marginalized terms conditional on stage II information
    #------------------------------------------------------------------------------------------
    
      All_c2r=t(sapply(X=seq(1,n),FUN=All_c2,
                       epsilon1_2=epsilon21,epsilon2_2=epsilon22,
                       epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                       Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,c2_data=c2_data))
      colnames(All_c2r)=c("T10B_c2r","T10_c2r","T11_c2r","T12_c2r")
      
      All9_c2r=t(sapply(X=seq(1,n),FUN=All_c3,
                        epsilon1_1=epsilon21,epsilon2_1=epsilon22,
                        Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,c3_data=c2_data))
      colnames(All9_c2r)=c("T7B_c2r","T7_c2r","T8_c2r","T9_c2r")
    
    
    if(regime==1){ #when it's a term that involves the optimal regime, but where current stage treatment is assigned
                   #according to another regime.
      All_c2r_opt=t(sapply(X=seq(1,n),FUN=All_c2,
                           epsilon1_2=epsilon21,epsilon2_2=epsilon22,
                           epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                           Theta1=Th1,Theta2=Th2,Theta3=Th3,c2_data=c2_data))
      colnames(All_c2r_opt)=c("T10B_c2r_opt","T10_c2r_opt","T11_c2r_opt","T12_c2r_opt")
      
      All9_c2r_opt=t(sapply(X=seq(1,n),FUN=All_c3,
                            epsilon1_1=epsilon21,epsilon2_1=epsilon22,
                            Theta1=Th1,Theta2=Th2,Theta3=Th3,c3_data=c2_data))
      colnames(All9_c2r_opt)=c("T7B_c2r_opt","T7_c2r_opt","T8_c2r_opt","T9_c2r_opt")
    }else{
      All_c2r_opt=All_c2r_opt_dtr    #These terms have already been computed with the Marginals_At_Opt() function
                                     #We return them again just for symmetry in the function
      All9_c2r_opt=All9_c2r_opt_dtr  
    }
    #there are 16 terms to return
    return(cbind(All_c2r,All9_c2r,All_c2r_opt,All9_c2r_opt))
  }

  c1_compute=function(c1_data,regime=0,
                      Theta1,Theta2,Theta3,
                      Th1,Th2,Th3,
                      epsilon11, epsilon12, 
                      epsilon21, epsilon22,
                      epsilon31, epsilon32, 
                      All_c1r_opt_dtr,
                      All9_c1r_opt_dtr,
                      All6_c1r_opt_dtr){
    #------------------------------------------------------------------------------------------
    #This function computes all relevant marginalized terms conditional on stage I information
    #------------------------------------------------------------------------------------------
    All_c1r=t(sapply(X=seq(1,n),FUN=All_c1,
                     epsilon1_3=epsilon11,epsilon2_3=epsilon12,
                     epsilon1_2=epsilon21,epsilon2_2=epsilon22,
                     epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                     Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,c1_data=c1_data))
    colnames(All_c1r)=c("T10B_c1r","T10_c1r","T11_c1r","T12_c1r")
    #------------
    All9_c1r=t(sapply(X=seq(1,n),FUN=All_c2,
                      epsilon1_2=epsilon11,epsilon2_2=epsilon12,
                      epsilon1_1=epsilon21,epsilon2_1=epsilon22,
                      Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,c2_data=c1_data))
    colnames(All9_c1r)=c("T7B_c1r","T7_c1r","T8_c1r","T9_c1r")
    #------------
    All6_c1r=t(sapply(X=seq(1,n),FUN=All_c3,
                      epsilon1_1=epsilon11,epsilon2_1=epsilon12,
                      Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,c3_data=c1_data))
    colnames(All6_c1r)=c("T4B_c1r","T4_c1r","T5_c1r","T6_c1r")
    
    if(regime==1){ #when it's a term that involves the optimal regime, but where current stage treatment is assigned
                   #according to another regime.
      All_c1r_opt=t(sapply(X=seq(1,n),FUN=All_c1,
                           epsilon1_3=epsilon11,epsilon2_3=epsilon12,
                           epsilon1_2=epsilon21,epsilon2_2=epsilon22,
                           epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                           Theta1=Th1,Theta2=Th2,Theta3=Th3,c1_data=c1_data))
      colnames(All_c1r_opt)=c("T10B_c1r_opt","T10_c1r_opt","T11_c1r_opt","T12_c1r_opt")
      
      #------------
      All9_c1r_opt=t(sapply(X=seq(1,n),FUN=All_c2,
                            epsilon1_2=epsilon11,epsilon2_2=epsilon12,
                            epsilon1_1=epsilon21,epsilon2_1=epsilon22,
                            Theta1=Th1,Theta2=Th2,Theta3=Th3,c2_data=c1_data))
      colnames(All9_c1r_opt)=c("T7B_c1r_opt","T7_c1r_opt","T8_c1r_opt","T9_c1r_opt")
      #------------
      All6_c1r_opt=t(sapply(X=seq(1,n),FUN=All_c3,
                            epsilon1_1=epsilon11,epsilon2_1=epsilon12,
                            Theta1=Th1,Theta2=Th2,Theta3=Th3,c3_data=c1_data))
      colnames(All6_c1r_opt)=c("T4B_c1r_opt","T4_c1r_opt","T5_c1r_opt","T6_c1r_opt")
    }else{
      All_c1r_opt=All_c1r_opt_dtr           #These terms have already been computed with the Marginals_At_Opt() function
                                            #We return them again just for symmetry in the function
      All9_c1r_opt=All9_c1r_opt_dtr
      All6_c1r_opt=All6_c1r_opt_dtr}
    
    #there are 24 terms to return
    return(cbind(All_c1r,All9_c1r,All6_c1r, All_c1r_opt, All9_c1r_opt, All6_c1r_opt))
  }

  Marginals_At_Opt=function(XX, X11,X12, A1, A1_opt,
                                X21,X22, A2, A2_opt,
                                X31,X32, A3, A3_opt,
                                X41,X42, A4, A4_opt, 
                                Theta1,Theta2, Theta3, 
                                Th1, Th2, Th3,
                                epsilon11, epsilon12,
                                epsilon21, epsilon22, 
                                epsilon31, epsilon32){
    # This function computes the expected value of terms in the data-generating process
    # conditional on stage 3, stage 2, and stage 1 information. 
    # when treatment is assigned according to the optimal regime at the stage of marginalization.
    # And when previous stage treatment varies as in the study population 
    # When previous stage treatment varies according to a specific regime, then these terms must be computed separately.
    
    c3_data=cbind(A3,X31,X32,XX)
    c2_data=cbind(A2,X21,X22,XX)
    c1_data=cbind(A1,X11,X12,XX)
    #---------------Conditional on 3
        #there are a total of four terms that need to be marginalized 
    All_c3r_opt_dtr=t(sapply(X=seq(1,n),FUN=All_c3,epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                         Theta1=Th1,Theta2=Th2,Theta3=Th3,c3_data=c3_data))
    colnames(All_c3r_opt_dtr)=c("T10B_c3r_opt","T10_c3r_opt","T11_c3r_opt","T12_c3r_opt")
    
    #---------------Conditional on 2
       #there are a total of eight terms that need to be marginalized
    All_c2r_opt_dtr=t(sapply(X=seq(1,n),FUN=All_c2,
                         epsilon1_2=epsilon21,epsilon2_2=epsilon22,
                         epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                         Theta1=Th1,Theta2=Th2,Theta3=Th3,c2_data=c2_data))
    colnames(All_c2r_opt_dtr)=c("T10B_c2r_opt","T10_c2r_opt","T11_c2r_opt","T12_c2r_opt")
    
    All9_c2r_opt_dtr=t(sapply(X=seq(1,n),FUN=All_c3,
                          epsilon1_1=epsilon21,epsilon2_1=epsilon22,
                          Theta1=Th1,Theta2=Th2,Theta3=Th3,c3_data=c2_data))
    colnames(All9_c2r_opt_dtr)=c("T7B_c2r_opt","T7_c2r_opt","T8_c2r_opt","T9_c2r_opt")
    
    #--------Conditional on 1---------
      #there are a total of 12 terms that need to be marginalized
    All_c1r_opt_dtr=t(sapply(X=seq(1,n),FUN=All_c1,
                             epsilon1_3=epsilon11,epsilon2_3=epsilon12,
                             epsilon1_2=epsilon21,epsilon2_2=epsilon22,
                             epsilon1_1=epsilon31,epsilon2_1=epsilon32,
                         Theta1=Th1,Theta2=Th2,Theta3=Th3,c1_data=c1_data))
    colnames(All_c1r_opt_dtr)=c("T10B_c1r_opt","T10_c1r_opt","T11_c1r_opt","T12_c1r_opt")
    #------------
    All9_c1r_opt_dtr=t(sapply(X=seq(1,n),FUN=All_c2,
                          epsilon1_2=epsilon11,epsilon2_2=epsilon12,
                          epsilon1_1=epsilon21,epsilon2_1=epsilon22,
                          Theta1=Th1,Theta2=Th2,Theta3=Th3,c2_data=c1_data))
    colnames(All9_c1r_opt_dtr)=c("T7B_c1r_opt","T7_c1r_opt","T8_c1r_opt","T9_c1r_opt")
    #------------
    All6_c1r_opt_dtr=t(sapply(X=seq(1,n),FUN=All_c3,
                          epsilon1_1=epsilon11,epsilon2_1=epsilon12,
                          Theta1=Th1,Theta2=Th2,Theta3=Th3,c3_data=c1_data))
    colnames(All6_c1r_opt_dtr)=c("T4B_c1r_opt","T4_c1r_opt","T5_c1r_opt","T6_c1r_opt")
    
    #--------------
    
    return(
      list( "All_c3r_opt_dtr"=All_c3r_opt_dtr,
            "All_c2r_opt_dtr"=All_c2r_opt_dtr, "All9_c2r_opt_dtr"=All9_c2r_opt_dtr,
            "All_c1r_opt_dtr"=All_c1r_opt_dtr, "All9_c1r_opt_dtr"=All9_c1r_opt_dtr,  "All6_c1r_opt_dtr"=All6_c1r_opt_dtr))
  }

Value_Each=function(XX,  X11,X12, A1, A1_Opt,
                           X21,X22, A2, A2_Opt,
                           X31,X32, A3, A3_Opt,
                           X41,X42, A4, A4_Opt, 
                           Theta1,Theta2, Theta3, 
                           Th1, Th2, Th3,
                           epsilon11, epsilon12,
                           epsilon21, epsilon22, 
                           epsilon31, epsilon32,
                           All_c3r_opt_dtr,
                           All_c2r_opt_dtr, All9_c2r_opt_dtr,
                           All_c1r_opt_dtr, All9_c1r_opt_dtr, All6_c1r_opt_dtr){
  

    A1_regime=((Theta1*X11+Theta2*X12+Theta3*XX)>0.5)+0
    A2_regime=((Theta1*X21+Theta2*X22+Theta3*XX)>0.5)+0
    A3_regime=((Theta1*X31+Theta2*X32+Theta3*XX)>0.5)+0
    A4_regime=((Theta1*X41+Theta2*X42+Theta3*XX)>0.5)+0
    
    #-----------------------------------------------------
    #  Marginal terms conditional on stage 3 covariates
    #-----------------------------------------------------
    c3_data=cbind(A3,X31,X32,XX)
    All_c3r=c3_compute(c3_data=c3_data,
                       Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,
                       Th1=Th1,Th2=Th2, Th3=Th3,
                       epsilon31=epsilon31, epsilon32=epsilon32, 
                       All_c3r_opt_dtr=All_c3r_opt_dtr)
    
    #producing terms needed to compute pseudo outcome
    c3_regime_data=cbind(A3_regime,X31,X32,XX)
    All_c3r_regime=c3_compute(c3_data=c3_regime_data,regime=1,
                              Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,
                              Th1=Th1,Th2=Th2, Th3=Th3,
                              epsilon31=epsilon31, epsilon32=epsilon32)
    
    Dat_c3=data.frame(cbind(XX,
                            X11, X12, A1_Opt, A1, 
                            X21, X22, A2_Opt, A2, 
                            X31, X32, A3_Opt, A3,
                            All_c3r))
    
    #-----------------------------------------------------
    #  Marginal terms conditional on stage 2 covariates
    #-----------------------------------------------------
    c2_data=cbind(A2,X21,X22,XX)
    All_c2r=c2_compute(c2_data=c2_data,
                       Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,
                       Th1=Th1,Th2=Th2, Th3=Th3,
                       epsilon21=epsilon21, epsilon22=epsilon22,
                       epsilon31=epsilon31, epsilon32=epsilon32, 
                       All_c2r_opt_dtr=All_c2r_opt_dtr, All9_c2r_opt_dtr=All9_c2r_opt_dtr)
    
    #producing terms needed to compute pseudo outcome and phis
    c2_regime_data=cbind(A2_regime,X21,X22,XX)
    All_c2r_regime=c2_compute(c2_data=c2_regime_data,regime=1,
                              Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,
                              Th1=Th1,Th2=Th2, Th3=Th3,
                              epsilon21=epsilon21, epsilon22=epsilon22,
                              epsilon31=epsilon31, epsilon32=epsilon32)
    
    Dat_c2=data.frame(cbind(XX,
                            X11, X12, A1_Opt, A1, 
                            X21, X22, A2_Opt, A2, 
                            All_c2r))
    
    
    #-----------------------------------------------------
    #  Marginal terms conditional on stage 1 covariates
    #-----------------------------------------------------
    c1_data=cbind(A1,X11,X12,XX)
    All_c1r=c1_compute(c1_data, 
                       Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,
                       Th1=Th1,Th2=Th2, Th3=Th3,
                       epsilon11=epsilon11, epsilon12=epsilon12, 
                       epsilon21=epsilon21, epsilon22=epsilon22,
                       epsilon31=epsilon31, epsilon32=epsilon32, 
                       All_c1r_opt_dtr=All_c1r_opt_dtr, All9_c1r_opt_dtr=All9_c1r_opt_dtr, All6_c1r_opt_dtr=All6_c1r_opt_dtr)
    
    #producing terms needed to compute phis
    c1_regime_data=cbind(A1_regime,X11,X12,XX)
    All_c1r_regime=c1_compute(c1_regime_data,regime=1,
                              Theta1=Theta1,Theta2=Theta2,Theta3=Theta3,
                              Th1=Th1,Th2=Th2, Th3=Th3,
                              epsilon11=epsilon11, epsilon12=epsilon12, 
                              epsilon21=epsilon21, epsilon22=epsilon22,
                              epsilon31=epsilon31, epsilon32=epsilon32)
    
    Dat_c1=data.frame(cbind(XX,X11, X12, A1_Opt, A1, All_c1r))
   #---------------------------------------------------------
    
    DR4v=vector()
    IPWv=vector()
    
    for(b in 1:B){
      set.seed(b)
      richi=n*t(rdirichlet(1, rep(1,n)))
      if(B==1){
        richi=rep(1,n)
      }
    #####################################################
    # Fitting Stage 4 model and computing pseudo-outcome
    #####################################################
    Dat_c4=data.frame(cbind(XX,X11,X12,A1_Opt,A1,
                            X21,X22,A2_Opt,A2,
                            X31,X32,A3_Opt,A3,
                            X41,X42,A4_Opt,A4))
    om4=lm(Y~X11+X12+
             X11:A1_Opt +X12:A1_Opt +XX:A1_Opt+A1_Opt+
             X11:A1     +X12:A1     +XX:A1    +A1+ 
             X21:A2_Opt +X22:A2_Opt +XX:A2_Opt+A2_Opt+
             X21:A2     +X22:A2     +XX:A2    +A2+ 
             X31:A3_Opt +X32:A3_Opt +XX:A3_Opt+A3_Opt+
             X31:A3     +X32:A3     +XX:A3    +A3+
             X41:A4_Opt +X42:A4_Opt +XX:A4_Opt+A4_Opt+
             X41:A4     +X42:A4     +XX:A4    +A4,data=Dat_c4,weights=richi)
    
    pseudo_outcome4=predict(om4, newdata=data.frame(XX=XX,
                                                    A1_Opt=A1_Opt, A1=A1,        X11=X11,X12=X12,
                                                    A2_Opt=A2_Opt, A2=A2,        X21=X21,X22=X22,
                                                    A3_Opt=A3_Opt, A3=A3,        X31=X31,X32=X32,
                                                    A4_Opt=A4_Opt, A4=A4_regime, X41=X41,X42=X42)) 
    ############################################################################
    # Fitting Stage 3 model and computing pseudo-outcome
    ############################################################################
      if(abs(Theta1-Th1)>0.001 & abs(Theta3-Th2)>0.001){ #when regime of interest is not the optimal regime
        om3=lm(pseudo_outcome4~X11+X12+
                 X11:A1_Opt +X12:A1_Opt +XX:A1_Opt+A1_Opt+
                 X11:A1     +X12:A1     +XX:A1+A1+ 
                 X21:A2_Opt +X22:A2_Opt +XX:A2_Opt+A2_Opt+
                 X21:A2     +X22:A2     +XX:A2+A2+ 
                 X31:A3_Opt +X32:A3_Opt +XX:A3_Opt+A3_Opt+
                 X31:A3     +X32:A3     +XX:A3+A3+ 
                 T10B_c3r+T10_c3r+T11_c3r+T12_c3r+
                 T10B_c3r_opt+T10_c3r_opt+ T11_c3r_opt+ T12_c3r_opt,data=Dat_c3,weights=richi)
      
      }else{ #When regime of interest is the optimal regime
        om3=lm(pseudo_outcome4~X11+X12+
                 X11:A1_Opt +X12:A1_Opt +XX:A1_Opt+A1_Opt+
                 X11:A1     +X12:A1     +XX:A1+A1+ 
                 X21:A2_Opt +X22:A2_Opt +XX:A2_Opt+A2_Opt+
                 X21:A2     +X22:A2     +XX:A2+A2+ 
                 X31:A3_Opt +X32:A3_Opt +XX:A3_Opt+A3_Opt+
                 X31:A3     +X32:A3     +XX:A3+A3+ 
                 T10B_c3r_opt+T10_c3r_opt+ T11_c3r_opt+ T12_c3r_opt,data=Dat_c3,weights=richi)
      }
      
     pseudo_outcome3=predict(om3, newdata=data.frame(XX=XX, 
                                                    X11=X11,X12=X12,A1_Opt=A1_Opt, A1=A1,
                                                    X21=X21,X22=X22,A2_Opt=A2_Opt, A2=A2,
                                                    X31=X31,X32=X32,A3_Opt=A3_Opt, A3=A3_regime,
                                                    All_c3r_regime))
    
    ############################################################################
    # Fitting Stage 2 model and computing pseudo-outcome
    ############################################################################
      if(abs(Theta1-Th1)>0.001 & abs(Theta3-Th2)>0.001){ #can't use direct equals on floats
      om2=lm(pseudo_outcome3~X11+X12+
               X11:A1_Opt +X12:A1_Opt +XX:A1_Opt +A1_Opt+
               X11:A1     +X12:A1     +XX:A1     +A1+ 
               X21:A2_Opt +X22:A2_Opt +XX:A2_Opt +A2_Opt+
               X21:A2     +X22:A2     +XX:A2     +A2+ 
               T7B_c2r    +T7_c2r     +T8_c2r     +T9_c2r+
               T7B_c2r_opt+T7_c2r_opt +T8_c2r_opt +T9_c2r_opt+
               T10B_c2r     +T10_c2r     +T11_c2r     +T12_c2r+
               T10B_c2r_opt +T10_c2r_opt +T11_c2r_opt +T12_c2r_opt,data=Dat_c2,weights=richi)
      }else{
        om2=lm(pseudo_outcome3~X11+X12+
                 X11:A1_Opt +X12:A1_Opt +XX:A1_Opt +A1_Opt+
                 X11:A1     +X12:A1     +XX:A1     +A1+ 
                 X21:A2_Opt +X22:A2_Opt +XX:A2_Opt +A2_Opt+
                 X21:A2     +X22:A2     +XX:A2     +A2+ 
                 T7B_c2r_opt+T7_c2r_opt +T8_c2r_opt +T9_c2r_opt+
                 T10B_c2r_opt +T10_c2r_opt +T11_c2r_opt +T12_c2r_opt,data=Dat_c2,weights=richi)
      }
      pseudo_outcome2=predict(om2, newdata=data.frame(XX=XX, 
                                                      X11=X11,X12=X12,A1_Opt=A1_Opt, A1=A1,
                                                      X21=X21,X22=X22,A2_Opt=A2_Opt, A2=A2_regime,
                                                      All_c2r_regime)) 
    ############################################################################
    # Fitting Stage 1 model 
    ############################################################################
    if(abs(Theta1-Th1)>0.001 & abs(Theta3-Th2)>0.001){ #can't use direct equals on floats
     om1=lm(pseudo_outcome2~X11+X12+
             X11:A1_Opt  +X12:A1_Opt  +XX:A1_Opt   +A1_Opt+
             X11:A1      +X12:A1      +XX:A1       +A1+ 
             T4B_c1r     +T4_c1r      +T5_c1r      +T6_c1r+
             T4B_c1r_opt +T4_c1r_opt  +T5_c1r_opt  +T6_c1r_opt+
             T7B_c1r     +T7_c1r      +T8_c1r      +T9_c1r+
             T7B_c1r_opt +T7_c1r_opt  +T8_c1r_opt  +T9_c1r_opt+
             T10B_c1r    +T10_c1r     +T11_c1r     +T12_c1r+
             T10B_c1r_opt+T10_c1r_opt +T11_c1r_opt +T12_c1r_opt,data=Dat_c1,weights=richi)
    }else{
      om1=lm(pseudo_outcome2~X11+X12+
               X11:A1_Opt   +X12:A1_Opt  +XX:A1_Opt   +A1_Opt+
               X11:A1       +X12:A1      +XX:A1       +A1+ 
               T4B_c1r_opt  +T4_c1r_opt  +T5_c1r_opt  +T6_c1r_opt+
               T7B_c1r_opt  +T7_c1r_opt  +T8_c1r_opt  +T9_c1r_opt+
               T10B_c1r_opt +T10_c1r_opt +T11_c1r_opt +T12_c1r_opt,data=Dat_c1,weights=richi)
    }
    ############################
    # Finally computing the phis
    ############################
    #fully conditional
    phi4_34=predict(om4, newdata=data.frame(XX=XX,  
                                            A1_Opt=A1_Opt, A1=A1_regime,X11=X11,X12=X12,
                                            A2_Opt=A2_Opt, A2=A2_regime,X21=X21,X22=X22,
                                            A3_Opt=A3_Opt, A3=A3_regime,X31=X31,X32=X32,
                                            A4_Opt=A4_Opt, A4=A4_regime,X41=X41,X42=X42)) 
    #conditional stage III and before
    phi3_34=predict(om3, newdata=data.frame(XX=XX,
                                            A1_Opt=A1_Opt, A1=A1_regime,X11=X11,X12=X12,
                                            A2_Opt=A2_Opt, A2=A2_regime,X21=X21,X22=X22,
                                            A3_Opt=A3_Opt, A3=A3_regime,X31=X31,X32=X32,
                                            All_c3r_regime))
    #conditional stage II and before
    phi2_34=predict(om2, newdata=data.frame(XX=XX,
                                            A1_Opt=A1_Opt, A1=A1_regime,X11=X11,X12=X12,
                                            A2_Opt=A2_Opt, A2=A2_regime,X21=X21,X22=X22,
                                            All_c2r_regime))
    #conditional on stage I
    phi1_34=predict(om1, newdata=data.frame(XX=XX,
                                            A1_Opt=A1_Opt, A1=A1_regime,X11=X11,X12=X12,
                                            All_c1r_regime))
    
    
    ##################
    # Treatment Modeling
    ##################
    trt_model1_24 = glm(A1 ~-1+X11+X12,    family=quasibinomial,weights=richi)# Fit propensity score model
    trt_model2_24 = glm(A2 ~-1+X21+X22+A1, family=quasibinomial,weights=richi) # Fit propensity score model
    trt_model3_24 = glm(A3 ~-1+X31+X32+A2, family=quasibinomial,weights=richi) # Fit propensity score model
    trt_model4_24 = glm(A4 ~-1+X41+X42+A3, family=quasibinomial,weights=richi) # Fit propensity score model
    
    
    C1_d = as.numeric(A1==A1_regime)    # Indicator of treatment matching regime recommendation
    C2_d = as.numeric(A2==A2_regime)    # Indicator of treatment matching regime recommendation
    C3_d = as.numeric(A3==A3_regime)    # Indicator of treatment matching regime recommendation
    C4_d = as.numeric(A4==A4_regime)    # Indicator of treatment matching regime recommendation
    
    #--------------
    p1_24 = predict(trt_model1_24, type="response")  # Estimated propensity score
    p2_24 = predict(trt_model2_24, type="response")  # Estimated propensity score
    p3_24 = predict(trt_model3_24, type="response")  # Estimated propensity score
    p4_24 = predict(trt_model4_24, type="response")  # Estimated propensity score
    
    
    p1_d_24 = A1*p1_24 + (1 - A1)*(1 - p1_24) 
    p2_d_24 = A2*p2_24 + (1 - A2)*(1 - p2_24) 
    p3_d_24 = A3*p3_24 + (1 - A3)*(1 - p3_24) 
    p4_d_24 = A4*p4_24 + (1 - A4)*(1 - p4_24) 
    
    
    AugY4=phi1_34+
      (C1_d)*                                              (1/p1_d_24)* (phi2_34-phi1_34)+
      (C1_d)*(C2_d)*                              (1/(p2_d_24*p1_d_24))*(phi3_34-phi2_34)+
      (C1_d)*(C2_d)*(C3_d)*               (1/(p3_d_24*p2_d_24*p1_d_24))*(phi4_34-phi3_34)+
      (C1_d)*(C2_d)*(C3_d)*(C4_d)*(1/(p4_d_24*p3_d_24*p2_d_24*p1_d_24))*(Y      -phi4_34)
    
    #returning final values
    IPWv[b]=mean(richi*C1_d*C2_d*C3_d*C4_d*Y/(p1_d_24*p2_d_24*p3_d_24*p4_d_24)) #note that richi have been multiplied by n in the Bayesian case to cancel out the division by n in this mean
    DR4v[b]=mean(richi*AugY4)
  }
    if(B==1){
      Final_Results=matrix(c("DR4v"=DR4v,"IPWv"=IPWv),ncol=2)
    }else{
      Final_Results=as.matrix(cbind("DR4v"=DR4v,"IPWv"=IPWv))
    }
    return(Final_Results)
}

  
#generating list of regimes of interest
sequence1=seq(0.2,0.8,0.05)
sequence2=1-sequence1
sequence3=seq(-0.3,0.3,0.1)
theta=as.matrix(expand.grid(sequence1,sequence3))


DR_Est=array(0,dim=c(R,B,nrow(theta)))
IPW_Est=array(0,dim=c(R,B,nrow(theta)))

for(repi in 1:R){
  set.seed(repi)
  #Note: the implemented code is tailored to the data-generating mechanism presented below. 
  start_time <- Sys.time()
  
  #Generating Data for each replicate
  XX=3*rbinom(n,1,0.5)
  rand_comp=rnorm(n, 0, 0.1)  # Random error in Y
  X11=rnorm(n,mean=1,sd=1)
  X12=rnorm(n,sd=1)
  A1=rbinom(n, 1, expit(.5*X12+1*X11))
  A1_Opt=((0.5*X11+0.5*X12+0.1*XX)>0.5)+0
  
  X21=rnorm(n,0+0.2*A1+0.1*X11,sd=1)
  X22=rnorm(n,0+0.5*A1+0.1*X12,sd=1)
  A2=rbinom(n, 1, expit(.5*X22-0.6*A1+ 1*X21))
  A2_Opt=((0.5*X21+0.5*X22+0.1*XX)>0.5)+0
  
  X31=rnorm(n,0+0.2*A2+0.1*X21,sd=1)
  X32=rnorm(n,0+0.5*A2+0.1*X22,sd=1)
  A3=rbinom(n, 1, expit(.5*X32-0.6*A2+ 1*X31))
  A3_Opt=((0.5*X31+0.5*X32+0.1*XX)>0.5)+0
  
  X41=rnorm(n,0+0.2*A3+0.1*X31,sd=1)
  X42=rnorm(n,0+0.5*A3+0.1*X32,sd=1)
  A4=rbinom(n, 1, expit(.5*X42-0.6*A3+ 1*X41))
  A4_Opt=((0.5*X41+0.5*X42+0.1*XX)>0.5)+0
  
  Y=X11+X12 - (0.5*X11+0.5*X12+0.1*XX-0.5)*(A1_Opt-A1) -
    (0.5*X21+0.5*X22+0.1*XX-0.5)*(A2_Opt-A2)-
    (0.5*X31+0.5*X32+0.1*XX-0.5)*(A3_Opt-A3)-
    (0.5*X41+0.5*X42+0.1*XX-0.5)*(A4_Opt-A4)+rand_comp
  
  epsilon11=rnorm(MM,sd=1)  #randomo components of stage I terms
  epsilon12=rnorm(MM,sd=1)
  
  epsilon21=rnorm(MM,sd=1)  #random components of stage II terms
  epsilon22=rnorm(MM,sd=1)
  
  epsilon31=rnorm(MM,sd=1)  #random components of stage III terms
  epsilon32=rnorm(MM,sd=1)
  
  
  #for each replication, we only need to compute the marginals
  #conditional on stage i information once
  Mar_opt_dtr=Marginals_At_Opt(XX=XX,X11=X11,X12=X12, A1=A1, A1_opt=A1_opt, 
                 X21=X21,X22=X22, A2=A2, A2_opt=A2_opt,
                 X31=X31,X32=X32, A3=A3, A3_opt=A3_opt,
                 X41=X41,X42=X42, A4=A4, A4_opt=A4_opt, 
                 Th1=0.5, Th2=0.5, Th3=.1,
                 epsilon11=epsilon11, epsilon12=epsilon12,
                 epsilon21=epsilon21, epsilon22=epsilon22, 
                 epsilon31=epsilon31, epsilon32=epsilon32)

  
  for(i in 1:nrow(theta)){
    #estimating the value for each DTR of interest
    outi=Value_Each(XX=XX,X11=X11,X12=X12, A1=A1, A1_Opt=A1_Opt, 
                          X21=X21,X22=X22, A2=A2, A2_Opt=A2_Opt,
                          X31=X31,X32=X32, A3=A3, A3_Opt=A3_Opt,
                          X41=X41,X42=X42, A4=A4, A4_Opt=A4_Opt, 
                          Theta1=theta[i,1],Theta2=(1-theta[i,1]), Theta3=theta[i,2], 
                          Th1=0.5, Th2=0.5, Th3=.1,
                          epsilon11=epsilon11, epsilon12=epsilon12,
                          epsilon21=epsilon21, epsilon22=epsilon22, 
                          epsilon31=epsilon31, epsilon32=epsilon32,
                          All_c3r_opt_dtr=Mar_opt_dtr$All_c3r_opt_dtr,
                          All_c2r_opt_dtr=Mar_opt_dtr$All_c2r_opt_dtr, All9_c2r_opt_dtr=Mar_opt_dtr$All9_c2r_opt_dtr,
                          All_c1r_opt_dtr=Mar_opt_dtr$All_c1r_opt_dtr, All9_c1r_opt_dtr=Mar_opt_dtr$All9_c1r_opt_dtr,  All6_c1r_opt_dtr=Mar_opt_dtr$All6_c1r_opt_dtr)
          
    DR_Est[repi,,i]=outi[,1]
    IPW_Est[repi,,i]=outi[,2]
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  
}

#-------------------------------------------------------------------------------
# Producing row 4 and 5, columns 3,4,5 of Table 2 (Frequentist)
# Note: based on this program above and from Program 1, it is straightforward to generate results for the Bayesian 
#       setting, though it is much more computationally complex
#-------------------------------------------------------------------------------
if(B==1){
  freq_max_val=function(x){
    index=which(x==max(x))
    return(c(x[index],index))
  }
  
  IPW_max=apply(IPW_Est[,1,],MARGIN=1,FUN=freq_max_val)
  DR_max=apply(DR_Est[,1,],MARGIN=1,FUN=freq_max_val)
  
  
  
  print(c(mean(theta[DR_max[2,],1]) ,sd(theta[DR_max[2,],1])
         ,mean(theta[DR_max[2,],2]) ,sd(theta[DR_max[2,],2])
         ,mean(DR_max[1,])          ,sd(DR_max[1,])))
  # Output from test settings with R=2:   0.50000000 0.00000000 0.10000000 0.00000000 0.99953437 0.02553246
  # Output from test settings with R=100: 0.49250000 0.03579896 0.09900000 0.01000000 0.99943209 0.06338475
  print(c(mean(theta[IPW_max[2,],1]) ,sd(theta[IPW_max[2,],1])
         ,mean(theta[IPW_max[2,],2]) ,sd(theta[IPW_max[2,],2])
         ,mean(IPW_max[1,])          ,sd(IPW_max[1,])))
  # Output from test settings with R=2:   0.37500000 0.24748737 0.15000000 0.21213203 1.26884918 0.03705943
  # Output from test settings with R=100: 0.4835000 0.1681833 0.0790000 0.1387498 1.1685376 0.1442115
}



