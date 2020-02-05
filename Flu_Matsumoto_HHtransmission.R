library(Rcpp)
library(RcppArmadillo)
library(LaplacesDemon)

# data import
HHdata<-as.matrix(read.csv("HHdata_FluMatsumoto2014-15.csv"))
HHdata[,1:5]=pmax(HHdata[,1:5],HHdata[,6:10]) # If HH member < infected member, assume HH member is underreported and increase household member


# data with the "single parent" category
index_SP=HHdata[,"Father.num"]+HHdata[,"Mother.num"]==1
HHdata_SP=cbind(HHdata[,1:5],SP.num=0,HHdata[,6:10],SP.inf=0,HHdata[,11])
HHdata_SP[index_SP,c("SP.num","SP.inf")]=HHdata_SP[index_SP,c("Father.num","Father.inf")]+HHdata_SP[index_SP,c("Mother.num","Mother.inf")]
HHdata_SP[index_SP,c("Father.num","Father.inf","Mother.num","Mother.inf")]=0

## Rcpp functions
# Longini-Koopman model likelihood
code_cppLK=c('
  int vectok(arma::vec vcases,arma::vec vmaxcase){
    unsigned int vlen=vcases.n_elem;
    arma::vec vbase=arma::ones<arma::vec>(vlen);
    if(vlen>1){
      arma::vec base_temp=arma::cumprod(vmaxcase+1);
      vbase.tail(vlen-1)=base_temp.head(vlen-1);
    }
    return arma::dot(vcases,vbase);
  }
  arma::vec ktovec(int k,arma::vec vcases){
    if(vcases.n_elem==1){
      arma::vec ret=arma::vec(1);
      ret(0)=k;
      return ret;
    }else{
      arma::vec vhead=arma::vec(1);
      vhead(0)=k%(int)(vcases(0)+1.5);
      return join_cols(vhead,ktovec(floor(k/(vcases(0)+1)),vcases.tail(vcases.n_elem-1)));
    }
  }
  double LK_nosurvivor(arma::vec vbeta_s,arma::mat cmat,arma::vec vlambda,arma::vec vcases,arma::vec nfactor){
    double sum_p=0;
    for(int k=0;k<prod(vcases+1)-1;++k){
      arma::vec tempcase=ktovec(k,vcases);
      double lch=0;
      for(unsigned int i=0;i<vcases.n_elem;++i){lch+=Rf_lchoose(vcases(i),tempcase(i));}
      sum_p+=exp(lch+(-arma::dot(vcases-tempcase,(arma::diagmat(vbeta_s%nfactor)*cmat*tempcase)+vlambda)))*LK_nosurvivor(vbeta_s,cmat,vlambda,tempcase,nfactor);
    }
    return 1.0-sum_p;
  }
  double logLK(arma::vec vbeta_s,arma::mat cmat,arma::vec vlambda,arma::vec vcases,arma::vec vfam,double npower){
    arma::vec N_eff=(cmat-arma::diagmat(arma::diagvec(cmat)/arma::max(vfam,vfam*0.0+1.0)))*vfam/mean(cmat*arma::ones<arma::vec>(vfam.n_elem)); // calculate C_k (normalised)
    if(arma::dot(N_eff,vfam)==0)N_eff=vfam*0.0+1.0;  // avoid zero division caused by a household of size is 1  
    arma::vec nfactor=1.0/pow(N_eff,npower);
    arma::vec vsurvs=vfam-vcases;
    double lcsum=0;
    for(unsigned int c=0;c<vfam.n_elem;++c){lcsum+=Rf_lchoose(vfam(c),vcases(c));}// combination (N-n,n)
    return lcsum - arma::dot(vsurvs,arma::diagmat(vbeta_s%nfactor)*cmat*vcases+vlambda) /* <- p of survival*/ + log(LK_nosurvivor(vbeta_s,cmat,vlambda,vcases,nfactor));
  }
  double ldbinom(arma::vec x,arma::vec size,arma::vec p){
    int lgth=size.n_elem;
    double res=0;
    for(int i=0;i<lgth;++i){
      res+=R::dbinom(x(i),size(i),p(i),true);
    }
    return res;
  }
  ','
  double logLK_HH(arma::vec vbeta_s,arma::mat cmat,arma::vec vlambda,arma::mat HHdata, double npower,arma::vec ascertain){
    int nclass=(HHdata.n_cols-1)/2;
    arma::vec freq=HHdata.tail_cols(1);
    arma::mat mfam=HHdata.head_cols(nclass);
    arma::mat mcases=HHdata.cols(nclass,2*nclass-1);

    double ret=0;
    if(sum(ascertain)<0.0){
      for(unsigned int r=0;r<mcases.n_rows;++r){
        ret+=freq(r)*logLK(vbeta_s,cmat,vlambda,mcases.row(r).t(),mfam.row(r).t(),npower);
      }
    }else{
      for(unsigned int r=0;r<mcases.n_rows;++r){
        double marg_lik=0; //marginal likelihood given observed cases
        
        for(int k=0;k<prod(mfam.row(r)-mcases.row(r)+1);++k){ //sum only for trucase>mcases
          arma::vec truecase=ktovec(k,mfam.row(r).t()-mcases.row(r).t())+mcases.row(r).t();
          marg_lik+=exp(logLK(vbeta_s,cmat,vlambda,truecase,mfam.row(r).t(),npower)+ldbinom(mcases.row(r).t(),truecase,ascertain));
        }
        ret+=freq(r)*log(marg_lik);
      }
    }
    return ret;
  }')
cppFunction(depends='RcppArmadillo',includes=code_cppLK[1],code=code_cppLK[2])

## Code testing
source("codetest.r")
maxfam=c(1,4,1,1,3) # Upper limit for each category
HHcomp_ptn=expand.grid(0:maxfam[1],0:maxfam[2],0:maxfam[3],0:maxfam[4],0:maxfam[5])[-1,] # all possible HH compositions within upper limit
test_sum1(runif(5),matrix(runif(25),5),runif(5),HHcomp_ptn,runif(1),-1) # Checks if the sum of all possible infection patterns in a specified HH composition is 1
test_nocase(runif(5),matrix(runif(25),5),runif(5),HHcomp_ptn,runif(1),-1) # Checks if the loglikelihood of a household with no infection matches the mathematical derivation

## MCMC function
runMCMC<-function(model,HHdataset=list(HHdata,HHdata_SP),WBIC=F){
  ntype<-5+model$single_parent*1
  parlist=list(risk_community=rep(0.3,(ntype)^as.numeric(model$het_rcom)),HH_transmission=rep(0.5,ntype^as.numeric(model$prop)))
  logical_powfree<-!(model$pow==0||model$pow==1)
  HHdata_MCMC=HHdataset[[1]]
  if(logical_powfree)parlist$gamma=0.5 else gamma=as.numeric(model$pow)
  Data=list(N=sum(HHdata_MCMC[,11]),mon.names=c("nlogl","beta"),parm.names=as.parm.names(parlist),pos=as.list(sapply(names(parlist),function(x)grep(x,as.parm.names(parlist)),simplify = F)))
  inf_index<-c("Father.inf","Mother.inf","Other.inf","Sibling.inf")
  fam_index<-c("Father.num","Mother.num","nonsiblings.num","siblings.num")
  if(model$single_parent){
    inf_index=c(inf_index,"SParent.inf")
    fam_index=c(fam_index,"SParent.fam")
    HHdata_MCMC=HHdataset[[2]]
  }
  Model<-function(parm,Data){
    # parameters
    gm<-ifelse(logical_powfree,parm[Data$pos$gamma],gamma)
    rcom<-interval(parm[Data$pos$risk_community],0)
    c_HH<-interval(parm[Data$pos$HH_transmission],0)
  
    vbeta<-rep(0,ntype)+c_HH
    vrcom<-rep(0,ntype)+rcom
    
    # Construct contact intensity matrix c_kl
    ntype=length(vbeta)
    cmat<-matrix(vbeta[ntype],ntype,ntype)
    cmat[,1:2]=c(vbeta[1],vbeta[-ntype])
    cmat[1:2,]=rep(c(vbeta[1],vbeta[-ntype]),each=2)
    
    
    ## Model variants with different c_kl
    
    ## Model 12a: intense contact within couples
    # cmat[5,1:2]=cmat[1:2,5]=vbeta[2]
    # cmat[3:4,3:4]=vbeta[4]
    
    ## Model 12b: mother acting as a hub
    # cmat[5,1:2]=cmat[1:2,5]=vbeta[2]
    # cmat[3,4]=cmat[4,3]=cmat[4,5]=cmat[5,4]=vbeta[4]
    
    ## Model 12c: generation assortative
    # cmat[c(3,5),1:2]=cmat[1:2,c(3,5)]=vbeta[5]
    # cmat[2:3,2:3]=vbeta[2]
    # cmat[5,5]=vbeta[4]
    
    beta<-mean(cmat%*%rep(1,ntype))^gm*cmat[2,3]^(1-gm) # compute beta
    
    # ascertainment bias
    ascertain<- -1 # ascertainment bias not considered when ascertain = -1
    if(model$ascertain){
      ascertain<-numeric(ntype)+0.8
      ascertain[c(1,5)]=0.8
    }
    
    ll<-logLK_HH(numeric(ntype)+1,cmat,vrcom,HHdata_MCMC,gm,ascertain)
    
    parm=interval(parm,0);parm[Data$pos$gamma]=gm[1]
    return(list(LP=ll/log(Data$N)^WBIC,Dev=-2*ll,Monitor=c(ll,beta),yhat=NULL,parm=parm))
  }
  fit<-LaplacesDemon(Model=Model,Data=Data,Initial.Values=unlist(parlist),Covar=NULL,Iterations=10000,Status=10,Thinning=10,Algorithm='HARM',Specs=NULL)
  fit<-LaplacesDemon(Model=Model,Data=Data,Initial.Values=as.initial.values(fit),Covar=fit$Covar,Iterations=30000,Status=10000,Thinning=10,Algorithm='AM',Specs=list(Adaptive=300,Periodicity=100))
  fit<-LaplacesDemon(Model=Model,Data=Data,Initial.Values=as.initial.values(fit),Covar=fit$Covar,Iterations=40000,Status=1000,Thinning=4,Algorithm='RWM')
  fit<-Combine(fit,Data,0)
  fit$Data=Data
  return(fit)  
}


## WBIC computation
modelpattern<-expand.grid(prop=c(F,T),pow=c(F,T,0.5),het_rcom=c(F,T),single_parent=c(F,T),ascertain=F)
fit_WBIC<-list()
for(mID in 1:nrow(modelpattern)){fit_WBIC[[mID]]<-runMCMC(model=modelpattern[mID,],list(HHdata,HHdata_SP),WBIC=T)}

## MCMC implementation
mID=12 # best model
fit_MCMC<-runMCMC(model=modelpattern[mID,],list(HHdata,HHdata_SP))

# Posterior samples
epsilon_k<-fit_MCMC$Posterior2[,1:5]
c_kl<-fit_MCMC$Posterior2[,5+1:5]
gamma<-fit_MCMC$Posterior2[,11]
beta<-fit_MCMC$Monitor[,2]

