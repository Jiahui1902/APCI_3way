rm(list = ls())
# load packages: if packages is not successfully loaded, install the corresponding
# packages
pkgs <- c("haven","survey","tidyverse","magrittr",'data.table')
installpackages <- lapply(pkgs,function(x){
  if(x %in% rownames(installed.packages()) == FALSE) {install.packages(`x`)}
})
loadpackages <- lapply(pkgs,function(x){
  library(`x`,character.only = T)
})
rm(list = c("installpackages","loadpackages","pkgs"))


int_model <- function(fm,
                      data,
                      outcome = "inlfc",
                      age = "acc",
                      period = "pcc",
                      cohort = NULL,
                      weight = NULL,
                      covariate = NULL,
                      family = "quasibinomial",
                      threeway = FALSE,
                      ...){
  ###############
  #no missing data
  ###############
  if(family=="binomial"){
    family <- "quasibinomial"
  }
  
  data2 = na.omit( data[ , unique(c(outcome,age,period,cohort,weight,covariate)[!is.null(c(outcome,age,period,cohort,weight,covariate))]) ] )
  
  if(is.null(weight)){
    weight <- 1
  }
  
  A = length(unique(data2[,age]))
  P = length(unique(data2[,period]))
  C = A + P - 1
  
  data2$acc <- data2[,age]
  data2$pcc <- data2[,period]
  
  data2$acc <- as.factor(data2$acc)
  data2$pcc <- as.factor(data2$pcc)
  
  options(contrasts=c("contr.sum","contr.poly"), na.action = na.omit)
  wtdata2 = survey::svydesign(id=~1, strata=NULL, data=data2,
                              weights=as.formula(paste0("~",weight)))
  
  temp6_formula <- as.formula(fm)
  print(temp6_formula)
  temp6 = survey::svyglm(temp6_formula,
                         wtdata2,
                         family = get(family))
  
  # design matrix ---
  model <- temp6
  m <- model.matrix(model)
  # 2 way interactions
  m2 <- m[,str_detect( colnames(m), "^acc[0-9]*:pcc[0-9]$")]
  m2 <- m2[!duplicated(m2),]
  # 3 way interactions
  m3 <- m[,str_detect( colnames(m), "^acc[0-9]*:pcc[0-9]*:region[0-9]*$")]
  m3 <- m3[!duplicated(m3),]
  
  # age-period-cohort effects -----
  r6 = model$coefficients[stringr::str_detect(names(model$coefficients) , "acc|pcc|(Intercept)")]
  r6se = summary(model)$coef[,"Std. Error"]
  r6p = summary(model)$coef[,"Pr(>|t|)"]
  # age ----
  fullae = array(rep(0, A), dim=c(A, 1))
  fullas = array(rep(0, A), dim=c(A, 1))
  S1 = array(rep(0, A*(A-1)), dim=c(A, (A-1)))
  ind = A*1:(A-1)
  newind = 1:(A*(A-1))
  newind = newind[-ind]
  S1[newind]  = diag(A-1)
  S1[ind]    = rep(-1,(A-1))
  
  fullae = as.vector(S1%*%model$coef[stringr::str_detect(names(model$coef) , "^acc([0-9])*$" )])
  
  row_ind <- stringr::str_detect(rownames(vcov(model)) , "^acc([0-9])*$")
  col_ind <- stringr::str_detect(colnames(vcov(model)) , "^acc([0-9])*$")
  
  fullas = sqrt(diag(S1%*%vcov(model)[row_ind, col_ind]%*%t(S1)))
  fullat = fullae/fullas
  fullap = pt(-abs(fullat),df.residual(model))*2
  
  sig = rep('   ', A)
  sig[fullap<.05] = '*  '
  sig[fullap<.01] = '** '
  sig[fullap<.001] = '***'
  fullasig = sig
  fulla=cbind(fullae, fullas, fullap, fullasig)
  age_results <- fulla %>% as.data.frame
  colnames(age_results) = c("age_estimate", "age_se", "age_p", "sig")
  rownames(age_results) = c()
  
  # period ----
  fullpe = array(rep(0, P), dim=c(P, 1))
  fullps = array(rep(0, P), dim=c(P, 1))
  S2 = array(rep(0, P*(P-1)), dim=c(P, (P-1)))
  ind = P*1:(P-1)
  newind = 1:(P*(P-1))
  newind = newind[-ind]
  S2[newind]  = diag(P-1)
  S2[ind]    = rep(-1,(P-1))
  
  fullpe = as.vector(S2%*%model$coef[stringr::str_detect(names(model$coef) , "^pcc([0-9])*$" )])
  row_ind <- stringr::str_detect(rownames(vcov(model)) , "^pcc([0-9])*$")
  col_ind <- stringr::str_detect(colnames(vcov(model)) , "^pcc([0-9])*$")
  
  fullps = sqrt(diag(S2%*%vcov(model)[row_ind,col_ind]%*%t(S2)))
  fullpt = fullpe/fullps
  fullpp = pt(-abs(fullpt),df.residual(model))*2
  
  sig = rep('   ', P)
  sig[fullpp<.05] = '*  '
  sig[fullpp<.01] = '** '
  sig[fullpp<.001] = '***'
  fullpsig = sig
  fullp=cbind(fullpe, fullps, fullpp, fullpsig)
  period_results <- fullp %>% as.data.frame
  colnames(period_results) = c("period_estimate", "period_se", "period_p", "sig")
  rownames(period_results) = c()
  
  # region ----
  if(threeway==TRUE){
    data$region <- factor(data$region)
    R <- nlevels(data$region)
    fullRe = array(rep(0, R), dim=c(R, 1))
    fullRs = array(rep(0, R), dim=c(R, 1))
    S2 = array(rep(0, R*(R-1)), dim=c(R, (R-1)))
    ind = R*1:(R-1)
    newind = 1:(R*(R-1))
    newind = newind[-ind]
    S2[newind]  = diag(R-1)
    S2[ind]    = rep(-1,(R-1))
    fullRe = as.vector(S2%*%model$coef[stringr::str_detect(names(model$coef) , "^region([0-9])*$" )])
    row_ind <- stringr::str_detect(rownames(vcov(model)) , "^region([0-9])*$")
    col_ind <- stringr::str_detect(colnames(vcov(model)) , "^region([0-9])*$")
    
    fullRs = sqrt(diag(S2%*%vcov(model)[row_ind,col_ind]%*%t(S2)))
    fullRt = fullRe/fullRs
    fullRp = pt(-abs(fullRt),df.residual(model))*2
    
    sig = rep('   ', R)
    sig[fullRp<.05] = '*  '
    sig[fullRp<.01] = '** '
    sig[fullRp<.001] = '***'
    fullRsig = sig
    fullR=cbind(fullRe, fullRs, fullRp, fullRsig)
    region_results <- fullR %>% as.data.frame
    colnames(region_results) = c("region_estimate", "region_se", "region_p", "sig")
    rownames(region_results) = c()
  }else{
    region_results <- NULL
  }
  
  # intercept ----
  inte = as.vector(r6[1])
  intse = r6se[1]
  intp = r6p[1]
  intsig = rep('   ', 1)
  intsig[r6p[1]<.05] = '*  '
  intsig[r6p[1]<.01] = '** '
  intsig[r6p[1]<.001] = '***'
  fullint = cbind(inte,intse,intp,intsig)
  
  # interaction effect of age and period -----
  T <- m2%>%as.numeric%>%matrix(.,A*P)
  row_ind <- stringr::str_detect(rownames(vcov(model)) , "^acc([0-9])*:pcc([0-9])*$")
  col_ind <- stringr::str_detect(colnames(vcov(model)) , "^acc([0-9])*:pcc([0-9])*$")
  
  row_ind_r6 <- stringr::str_detect(names(r6) , "^acc([0-9])*:pcc([0-9])*$")
  col_ind_r6 <- stringr::str_detect(names(r6) , "^acc([0-9])*:pcc([0-9])*$")
  
  iatemp = vcov(model)[row_ind,col_ind]
  iavcov = T%*%iatemp%*%t(T)
  df = model$df.residual
  
  
  iaesti = as.vector(T%*%r6[row_ind_r6])
  iase   = sqrt(diag(iavcov))
  iap    = pt(-abs(iaesti/iase), df)*2
  
  cindex <- sapply(1:P,function(j){
    seq((A+j-1),j, -1)
  })
  
  sig = rep('   ', (A*P))
  sig[iap<.05]  = '*  '
  sig[iap<.01]  = '** '
  sig[iap<.001] = '***'
  iasig = sig
  
  cohortindex = as.vector(cindex)
  ia          = as.data.frame(cbind(iaesti,iase,iap,iasig, cohortindex))
  
  ####################### inter-cohort changes 
  cohortint <- sapply(1:C,function(k){
    O = sum(cindex == k)
    k1 = rep(1/O, O)
    k2 = rep(0, A*P)
    k2[cindex == k] = k1
    
    contresti = k2%*%iaesti
    contrse = sqrt(t(k2)%*%iavcov%*%k2)
    
    t = contresti/contrse
    
    if (t > 0){
      p = 2*pt(t, df, lower.tail=F)
    } else {
      p = 2*pt(t, df, lower.tail=T)
    }
    
    sig <- '   '
    if (p<.05){
      sig <- '*  '
    }
    if(p<.01){
      sig <- '** '
    }
    if(p<.001){
      sig <- '***'
    }
    
    
    c(contresti,contrse,t,p,sig)
  })%>%t%>%
    as.data.frame%>%
    `colnames<-`(c("cohort_average","cohort_average_se",
                   "cohort_average_t","cohort_average_p","sig"))
  
  cohortint$cohort_group = seq(1, C)
  cohortint = cohortint[,c("cohort_group", "cohort_average","cohort_average_se",
                           "cohort_average_t","cohort_average_p", "sig")]
  
  
  # interaction effect of age and period and region -----
  if(threeway==TRUE){
    data$region <- as.factor(data$region)
    R <- nlevels(data$region)
    m3_1 <- m3[,(colnames(m3) %in% colnames(vcov(model)))]
    T <- m3_1%>%as.numeric%>%matrix(.,A*P*R)
    
    row_ind <- stringr::str_detect(rownames(vcov(model)) , "^acc([0-9])*:pcc([0-9])*:region([0-9])*$")
    col_ind <- stringr::str_detect(colnames(vcov(model)) , "^acc([0-9])*:pcc([0-9])*:region([0-9])*$")
    
    
    row_ind_r6 <- stringr::str_detect(names(r6)[(names(r6) %in% rownames(vcov(model)))] , "^acc([0-9])*:pcc([0-9])*:region([0-9])*$")
    col_ind_r6 <- stringr::str_detect(names(r6)[(names(r6) %in% colnames(vcov(model)))] , "^acc([0-9])*:pcc([0-9])*:region([0-9])*$")
    
    iatemp = vcov(model)[row_ind,col_ind]
    iavcov = T%*%iatemp%*%t(T)
    df = model$df.residual
    
    # T <- m3%>%as.numeric%>%matrix(.,A*P*R)
    iaesti = as.vector(T%*%r6[(names(r6) %in% colnames(vcov(model)))][row_ind_r6])
    iase   = sqrt(diag(iavcov))
    iap    = pt(-abs(iaesti/iase), df)*2
    
    cindex <- sapply(1:P,function(j){
      seq((A+j-1),j, -1)
    })
    cindex <- rep(cindex,R)
    rindex <- rep(1:4,each = A*P)
    
    sig = rep('   ', (A*P*R))
    sig[iap<.05]  = '*  '
    sig[iap<.01]  = '** '
    sig[iap<.001] = '***'
    iasig = sig
    
    cohortindex = as.vector(cindex)
    regionindex <- as.vector(rindex)
    ia          = as.data.frame(cbind(iaesti,iase,iap,iasig, cohortindex,regionindex))
    
    ####################### inter-cohort changes 
    cohortint3 <- lapply(1:R, function(r){
      sapply(1:C,function(k){
        O = sum(cindex == k&rindex==r)
        k1 = rep(1/O, O)
        k2 = rep(0, A*P)
        k2[cindex[rindex==r] == k] = k1
        
        contresti = k2%*%(iaesti[rindex==r])
        contrse = sqrt(t(k2)%*%(iavcov[rindex==r,rindex==r])%*%k2)
        
        t = contresti/contrse
        
        if (t > 0){
          p = 2*pt(t, df, lower.tail=F)
        } else {
          p = 2*pt(t, df, lower.tail=T)
        }
        
        sig <- '   '
        if (p<.05){
          sig <- '*  '
        }
        if(p<.01){
          sig <- '** '
        }
        if(p<.001){
          sig <- '***'
        }
        c(contresti,contrse,t,p,sig,r)
      })%>%t%>%
        as.data.frame%>%
        `colnames<-`(c("cohort_average","cohort_average_se",
                       "cohort_average_t","cohort_average_p","sig","region"))
    })%>%bind_rows()
    
    cohortint3$cohort_group = rep(1:C,R)
    cohortint3 = cohortint3[,c("cohort_group", "cohort_average","cohort_average_se",
                               "cohort_average_t","cohort_average_p", "sig","region")]
  }else{
    cohortint3 <- NA
  }
  
  
  # output
  list(A=A,P=P,C=C,model=temp6,age_main = age_results,
       period_main = period_results,region_main = region_results,
       age_period = cohortint,
       age_period_region = cohortint3)
}

# readfiles
allfiles <- list.files("/Users/xujiahui/OneDrive - The Pennsylvania State University/nhis bmi data",
                       full.names = T,pattern="*.txt")
# file names in the environment ultimately
nm <- allfiles%>%str_split(.,"/")%>%
  unlist%>%
  matrix(.,ncol = 6,byrow = T)%>%
  .[,6]%>%
  str_replace_all(.,".txt| ","")

# dataset <- lapply(setNames(allfiles[c(1:18,20:31)], make.names(
dataset <- lapply(setNames(allfiles[c(1,2,5,6)], make.names(
  # nm[c(1:18,20:31)])),function(file_name){
  nm[c(1:2,5,6)])),function(file_name){
    df <- read.table(file_name)
    #col_types = cols(Zip4 = col_character()))
    df
  })

# data <- all_bm

run_model <- function(data){
  data$acc <- factor(data$acc)
  data$pcc <- factor(data$pcc)
  data$obesity <- ifelse(data$bmicalc>=30,1,0)
  data$overweight <- ifelse(data$bmicalc>=25,1,0)
  data$region <- factor(data$region)
  xxx <- int_model(fm = "obesity ~ acc + pcc + region + acc:pcc + acc:pcc:region",
                   outcome = "obesity",
                   weight = "wt",
                   age = "acc",
                   period = "pcc",
                   covariate = "region",
                   data = data[order(data$region,data$pcc,data$acc),],
                   threeway = TRUE)
  xxx
}

bm <- run_model(data = dataset$all_bm)
bw <- run_model(data = dataset$all_bw)
wm <- run_model(data = dataset$all_wm)
ww <- run_model(data = dataset$all_ww)

# make plots ------
viz <- function(data){
  
  df <- data.frame(odds = (rep(data$region_main$region_estimate, each = (data$A+data$P-1))%>%as.character%>%as.numeric +
                             (data$age_period_region$cohort_average)%>%as.character%>%as.numeric) %>% exp,
                   region = data$age_period_region$region %>%factor(.,labels = c("Northeast","Midwest","South","West")),
                   cohort_group = data$age_period_region$cohort_group %>% factor(.,labels = seq(1910,1995,5)%>%as.character),
                   sig = data$age_period_region$sig)
  
  df <- df%>%
    group_by(cohort_group)%>%
    mutate(sig_n = sum(sig!="   "))
  df$odds <- ifelse(df$sig_n>0,df$odds,NA)
  df$star <- ifelse(df$sig!="   ",df$odds,NA)
  
  p <- ggplot(df,aes(group=region,x=cohort_group,y = odds))+
    geom_bar(aes(x=cohort_group,fill=region),stat="identity",position=position_dodge(),col="black")+
    scale_fill_brewer(palette="Greys",name = "Region",direction = 1)+    
    theme_minimal()+
    labs(x = "",names = "",y = "")+    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
    scale_x_discrete(breaks = seq(1910,1995,10)%>%as.character)+  
    guides(color = FALSE, size = FALSE)+
    geom_text(aes(x=cohort_group,y=star*1),label = "*",
              position = position_dodge(width = .9),
              color="red",size = 7)
  p
}

g9 <- viz(data = bm)+labs(title = "Black Men")
g10 <- viz(data = bw)+labs(title = "Black Women")
g11 <- viz(data = wm)+labs(title = "White Men")
g12 <- viz(data = ww)+labs(title = "White Women")
text <- "Cohort Group"
text.p1 <- ggpubr::ggparagraph(text = text, face = "italic", size = 12, color = "black")
ggpubr::ggarrange(g9,g10,g11,g12,NULL,text.p1,
                  ncol = 2,nrow = 3,heights = c(1, 1,0.15),
                  common.legend = T,legend = "bottom")

