require('tidyverse')
require('gtsummary')
require('sjmisc')
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(lmtp)


# CJP ----
if (ClaCos_anal==FALSE) {
  ITT_CC_CJP <- list();ITT_CC_CJP_pooled <- list();
  ITT_MI_CJP <- list();ITT_MI_CJP_pooled <- list();
  ITT_IPCW_CJP <- list();ITT_IPCW_CJP_pooled <- list();
  PP_IPCW_CJP <- list();PP_IPCW_CJP_pooled <- list();
  PP_STD_CJP <- list();PP_STD_CJP_pooled <- list();
  i <- 1
  for (i in 1:length(all_outcomes_CJP)) {
    
    if (MI==FALSE) {
      
      ITT_CC_CJP[[i]] <- run_tmle(df=df.full_imputed , outcome=all_outcomes_CJP[i], baseline=baseline_CJP,
                                  cens=NULL, MI = FALSE, anal="ITT", CompCovar=TRUE)
      
      #run anal with imputed outcomes and covariates
      if(run_ITT_MI) {
        
      ITT_MI_CJP[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed_MI[[m_idx]] , outcome=all_outcomes_CJP[i], baseline=baseline_CJP_MI,
                 cens=NULL, MI = TRUE, anal="ITT")
      } )
      ITT_MI_CJP_pooled[[i]] <- extract_Rubin(ITT_MI_CJP[[i]])
      }
      
      ITT_IPCW_CJP[[i]] <- run_tmle(df=df.full_imputed , outcome=all_outcomes_CJP[i], baseline=baseline_CJP,
                                    cens="cens", MI = FALSE, anal="ITT", CompCovar=TRUE)
      
      PP_IPCW_CJP[[i]] <- run_tmle(df=df.full_imputed , outcome=all_outcomes_CJP[i], baseline=baseline_CJP,
                                   cens="cens", MI = FALSE, anal="PP", CompCovar=TRUE)
      
      PP_STD_CJP[[i]] <- run_tmle(df=df.full_imputed , outcome=all_outcomes_CJP[i], baseline=baseline_CJP,
                                  cens=NULL, MI = FALSE, anal="PP", CompCovar=TRUE)
      
    } else if (MI==TRUE) {
      
      ##ITT remove missing outcome ----
      ITT_CC_CJP[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=all_outcomes_CJP[i], baseline=baseline_CJP, cens=NULL, MI = FALSE, anal="ITT")
      } )
      ITT_CC_CJP_pooled[[i]] <- extract_Rubin(ITT_CC_CJP[[i]])
      
      ##ITT impute missing outcome ----
      # Multiple Imputation is the most common partner for ITT. It is generally recommended because it maintains the randomized group sizes.
      ITT_MI_CJP[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=all_outcomes_CJP[i], baseline=baseline_CJP, cens=NULL, MI = TRUE, anal="ITT")
      } )
      ITT_MI_CJP_pooled[[i]] <- extract_Rubin(ITT_MI_CJP[[i]]$tmle_cont)
      
      ##ITT IPCW  ----
      #same as impute and IPCW
      ITT_IPCW_CJP[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=all_outcomes_CJP[i], baseline=baseline_CJP, cens="cens", MI = FALSE, anal="ITT")
      } )
      ITT_IPCW_CJP_pooled[[i]] <- extract_Rubin(ITT_IPCW_CJP[[i]])
      
      
      ##PP IPCW ----
      PP_IPCW_CJP[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=all_outcomes_CJP[i], baseline=baseline_CJP, cens="cens", MI = FALSE, anal = "PP")
      } )
      PP_IPCW_CJP_pooled[[i]] <- extract_Rubin(PP_IPCW_CJP[[i]])
      
      ##PP w/o IPCW ----
      PP_STD_CJP[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=all_outcomes_CJP[i], baseline=baseline_CJP, cens=NULL, MI = FALSE, anal = "PP")
      } )
      PP_STD_CJP_pooled[[i]] <- extract_Rubin(PP_STD_CJP[[i]])
      
    }
    
    print(all_outcomes_CJP[i])
    
  }
  
  save(ITT_CC_CJP,ITT_MI_CJP,ITT_IPCW_CJP,PP_IPCW_CJP,
       ITT_CC_CJP_pooled,ITT_MI_CJP_pooled,ITT_IPCW_CJP_pooled,PP_IPCW_CJP_pooled,PP_STD_CJP_pooled,
       file="CJP.RData")
}


#CJS ----
if (ClaCos_anal==FALSE) {
  ITT_CC_CJS_t1 <- list();ITT_CC_CJS_t1_pooled <- list()
  ITT_MI_CJS_t1 <- list();ITT_MI_CJS_t1_pooled <- list()
  ITT_IPCW_CJS_t1 <- list();ITT_IPCW_CJS_t1_pooled <- list()
  PP_IPCW_CJS_t1 <- list();PP_IPCW_CJS_t1_pooled <- list()
  PP_STD_CJS_t1 <- list();PP_STD_CJS_t1_pooled <- list()
  
  ITT_CC_CJS_t2 <- list();ITT_CC_CJS_t2_pooled <- list()
  ITT_MI_CJS_t2 <- list();ITT_MI_CJS_t2_pooled <- list();
  ITT_IPCW_CJS_t2 <- list();ITT_IPCW_CJS_t2_pooled <- list();
  PP_IPCW_CJS_t2 <- list();PP_IPCW_CJS_t2_pooled <- list();
  PP_STD_CJS_t2 <- list();PP_STD_CJS_t2_pooled <- list();
  for (i in 1:length(all_outcomes_CJS)) {
    
    if (MI==FALSE) {
      ITT_CC_CJS_t1[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"),
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                     cens=NULL, MI = FALSE, anal="ITT")
      ITT_CC_CJS_t2[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"),
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                     cens=NULL, MI = FALSE, anal="ITT")
      
      #run anal with imputed outcomes and covariates
      if(run_ITT_MI) {
        ITT_MI_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
          run_tmle(df=df.full_imputed_MI[[m_idx]],
                   outcome=paste0(all_outcomes_CJS[i], "_t1"),
                   baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")),
                   cens=NULL, MI = TRUE, anal="ITT")
        } )
        ITT_MI_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t1[[i]])
        ITT_MI_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
          run_tmle(df=df.full_imputed_MI[[m_idx]],
                   outcome=paste0(all_outcomes_CJS[i], "_t2"),
                   baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")),
                   cens=NULL, MI = TRUE, anal="ITT")
        } )
        ITT_MI_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t2[[i]])
      }
      
      ITT_IPCW_CJS_t1[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"),
                                       baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                       cens="cens", MI = FALSE, anal="ITT")
      ITT_IPCW_CJS_t2[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"),
                                       baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                       cens="cens", MI = FALSE, anal="ITT")
      
      PP_IPCW_CJS_t1[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"),
                                      baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                      cens="cens", MI = FALSE, anal = "PP")
      PP_IPCW_CJS_t2[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"),
                                      baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                      cens="cens", MI = FALSE, anal = "PP")
      
      PP_STD_CJS_t1[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"),
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                     cens=NULL, MI = FALSE, anal = "PP")
      PP_STD_CJS_t2[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"),
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")),
                                     cens=NULL, MI = FALSE, anal = "PP")
      
    } else if (MI==TRUE) {
      
      ## ITT remove missing outcome ----
      ITT_CC_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal="ITT")
      } )
      ITT_CC_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_CC_CJS_t1[[i]])
      
      ITT_CC_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal="ITT")
      } )
      ITT_CC_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_CC_CJS_t2[[i]])
      
      ##ITT impute missing outcome ----
      # Multiple Imputation is the most common partner for ITT. It is generally recommended because it maintains the randomized group sizes.
      ITT_MI_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = TRUE, anal="ITT")
      } )
      ITT_MI_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t1[[i]])
      
      ITT_MI_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = TRUE, anal="ITT")
      } )
      ITT_MI_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t2[[i]])
      
      ##ITT IPCW  ----
      #same as impute and IPCW
      ITT_IPCW_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal="ITT")
      } )
      ITT_IPCW_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_IPCW_CJS_t1[[i]])
      
      ITT_IPCW_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal="ITT")
      } )
      ITT_IPCW_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_IPCW_CJS_t2[[i]])
      
      ##PP IPCW ----
      PP_IPCW_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal = "PP")
      } )
      PP_IPCW_CJS_t1_pooled[[i]] <- extract_Rubin(PP_IPCW_CJS_t1[[i]])
      
      PP_IPCW_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal = "PP")
      } )
      PP_IPCW_CJS_t2_pooled[[i]] <- extract_Rubin(PP_IPCW_CJS_t2[[i]])
      
      ##PP w/o IPCW ----
      PP_STD_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal = "PP")
      } )
      PP_STD_CJS_t1_pooled[[i]] <- extract_Rubin(PP_STD_CJS_t1[[i]])
      
      PP_STD_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"),
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal = "PP")
      } )
      PP_STD_CJS_t2_pooled[[i]] <- extract_Rubin(PP_STD_CJS_t2[[i]])
    }
    
    print(all_outcomes_CJS[i])
    
  }
  
  save(ITT_CC_CJS_t1,ITT_MI_CJS_t1,ITT_IPCW_CJS_t1,PP_IPCW_CJS_t1,PP_STD_CJS_t1,
       ITT_CC_CJS_t1_pooled,ITT_MI_CJS_t1_pooled,ITT_IPCW_CJS_t1_pooled,PP_IPCW_CJS_t1_pooled,PP_STD_CJS_t1_pooled,
       file="CJS_t1.RData")
  save(ITT_CC_CJS_t2,ITT_MI_CJS_t2,ITT_IPCW_CJS_t2,PP_IPCW_CJS_t2,PP_STD_CJS_t2,
       ITT_CC_CJS_t2_pooled,ITT_MI_CJS_t2_pooled,ITT_IPCW_CJS_t2_pooled,PP_IPCW_CJS_t2_pooled,PP_STD_CJS_t2_pooled,
       file="CJS_t2.RData")
}


#CJS: ClaCos ----
if (ClaCos_anal) {
  ITT_CC_CJS_t1 <- list();ITT_CC_CJS_t1_pooled <- list()
  ITT_MI_CJS_t1 <- list();ITT_MI_CJS_t1_pooled <- list()
  ITT_IPCW_CJS_t1 <- list();ITT_IPCW_CJS_t1_pooled <- list()
  PP_IPCW_CJS_t1 <- list();PP_IPCW_CJS_t1_pooled <- list()
  PP_STD_CJS_t1 <- list();PP_STD_CJS_t1_pooled <- list()
  
  ITT_CC_CJS_t2 <- list();ITT_CC_CJS_t2_pooled <- list()
  ITT_MI_CJS_t2 <- list();ITT_MI_CJS_t2_pooled <- list();
  ITT_IPCW_CJS_t2 <- list();ITT_IPCW_CJS_t2_pooled <- list();
  PP_IPCW_CJS_t2 <- list();PP_IPCW_CJS_t2_pooled <- list();
  PP_STD_CJS_t2 <- list();PP_STD_CJS_t2_pooled <- list();
  for (i in 1:length(all_outcomes_CJS)) {
    
    if (MI==FALSE) {
      ITT_CC_CJS_t1[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                     cens=NULL, MI = FALSE, anal="ITT")
      ITT_CC_CJS_t2[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                     cens=NULL, MI = FALSE, anal="ITT")
      
      #run anal with imputed outcomes and covariates
      if(run_ITT_MI) {
        ITT_MI_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
          run_tmle(df=df.full_imputed_MI[[m_idx]],
                   outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                   baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), 
                   cens=NULL, MI = TRUE, anal="ITT")
        } )
        ITT_MI_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t1[[i]])
        ITT_MI_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
          run_tmle(df=df.full_imputed_MI[[m_idx]],
                   outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                   baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), 
                   cens=NULL, MI = TRUE, anal="ITT")
        } )
        ITT_MI_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t2[[i]])
      }
      
      ITT_IPCW_CJS_t1[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                                       baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                       cens="cens", MI = FALSE, anal="ITT")
      ITT_IPCW_CJS_t2[[i]] <- run_tmle(df=df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                                       baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                       cens="cens", MI = FALSE, anal="ITT")
      
      PP_IPCW_CJS_t1[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                                      baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                      cens="cens", MI = FALSE, anal = "PP")
      PP_IPCW_CJS_t2[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                                      baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                      cens="cens", MI = FALSE, anal = "PP")
      
      PP_STD_CJS_t1[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                     cens=NULL, MI = FALSE, anal = "PP")
      PP_STD_CJS_t2[[i]] <- run_tmle(df= df.full_imputed, outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                                     baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0"), paste0(all_outcomes_CJS[i], "_t0_NA")), 
                                     cens=NULL, MI = FALSE, anal = "PP")
      
    } else if (MI==TRUE) {
      
      ## ITT remove missing outcome ----
      ITT_CC_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal="ITT")
      } )
      ITT_CC_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_CC_CJS_t1[[i]])
      
      ITT_CC_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal="ITT")
      } )
      ITT_CC_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_CC_CJS_t2[[i]])
      
      ##ITT impute missing outcome ----
      # Multiple Imputation is the most common partner for ITT. It is generally recommended because it maintains the randomized group sizes.
      ITT_MI_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = TRUE, anal="ITT")
      } )
      ITT_MI_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t1[[i]])
      
      ITT_MI_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = TRUE, anal="ITT")
      } )
      ITT_MI_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_MI_CJS_t2[[i]])
      
      ##ITT IPCW  ----
      #same as impute and IPCW
      ITT_IPCW_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal="ITT")
      } )
      ITT_IPCW_CJS_t1_pooled[[i]] <- extract_Rubin(ITT_IPCW_CJS_t1[[i]])
      
      ITT_IPCW_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df=df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal="ITT")
      } )
      ITT_IPCW_CJS_t2_pooled[[i]] <- extract_Rubin(ITT_IPCW_CJS_t2[[i]])
      
      ##PP IPCW ----
      PP_IPCW_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal = "PP")
      } )
      PP_IPCW_CJS_t1_pooled[[i]] <- extract_Rubin(PP_IPCW_CJS_t1[[i]])
      
      PP_IPCW_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens="cens", MI = FALSE, anal = "PP")
      } )
      PP_IPCW_CJS_t2_pooled[[i]] <- extract_Rubin(PP_IPCW_CJS_t2[[i]])
      
      ##PP w/o IPCW ----
      PP_STD_CJS_t1[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t1"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal = "PP")
      } )
      PP_STD_CJS_t1_pooled[[i]] <- extract_Rubin(PP_STD_CJS_t1[[i]])
      
      PP_STD_CJS_t2[[i]] <- lapply(1:m, function(m_idx) {
        run_tmle(df= df.full_imputed %>% complete(m_idx) , outcome=paste0(all_outcomes_CJS[i], "_t2"), 
                 baseline=c(baseline_CJS, paste0(all_outcomes_CJS[i], "_t0")), cens=NULL, MI = FALSE, anal = "PP")
      } )
      PP_STD_CJS_t2_pooled[[i]] <- extract_Rubin(PP_STD_CJS_t2[[i]])
    }
    
    print(all_outcomes_CJS[i])
    
  }
  
  save(ITT_CC_CJS_t1,ITT_MI_CJS_t1,ITT_IPCW_CJS_t1,PP_IPCW_CJS_t1,PP_STD_CJS_t1,
       ITT_CC_CJS_t1_pooled,ITT_MI_CJS_t1_pooled,ITT_IPCW_CJS_t1_pooled,PP_IPCW_CJS_t1_pooled,PP_STD_CJS_t1_pooled,
       file="CJS_t1_ClaCos.RData")
  save(ITT_CC_CJS_t2,ITT_MI_CJS_t2,ITT_IPCW_CJS_t2,PP_IPCW_CJS_t2,PP_STD_CJS_t2,
       ITT_CC_CJS_t2_pooled,ITT_MI_CJS_t2_pooled,ITT_IPCW_CJS_t2_pooled,PP_IPCW_CJS_t2_pooled,PP_STD_CJS_t2_pooled,
       file="CJS_t2_ClaCos.RData")
}


stop("next step: reporting results")



