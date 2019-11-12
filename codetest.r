library(matrixStats)
# Test that sum of all possible pattrens is 1
test_sum1<-function(vbeta_s,cmat,vlambda,HHcomposition,npower,ascertain=-1){
    nclass=ncol(HHcomposition)
    assert<-all(apply(HHcomposition,1,function(vfam){
        #mcases<-[1:nclass+nclass]
        expandmcase<-expand.grid(0:vfam[1],0:vfam[2],0:vfam[3],0:vfam[4],0:vfam[5]) # all possible infection pattern in family composition vfam
        line_iter<-c(vfam,vfam*0,1)
        
        issum1<-all.equal(0, logSumExp(
                sapply(1:nrow(expandmcase),function(r){
                    line_iter[1:nclass+nclass]=unlist(expandmcase[r,])
                    return(logLK_HH(vbeta_s,cmat,vlambda,matrix(line_iter,1),npower,ascertain))
                })
            )
        )

        return(issum1)
    }))
    return(assert)
}

test_nocase<-function(vbeta_s,cmat,vlambda,HHcomposition,npower,ascertain=-1){
    nclass=ncol(HHcomposition)
    assert<-all(apply(HHcomposition,1,function(HHcomp_line){
        vfam=HHcomp_line
        #mcases<-[1:nclass+nclass]
        vcases<-numeric(nclass) # all possible infection pattern in family composition vfam
        data_line<-c(vfam,vcases,1)
        
        iscorrect<-all.equal(sum(-vlambda*vfam),
                logLK_HH(vbeta_s,cmat,vlambda,matrix(data_line,1),npower,ascertain)
        )
        

        return(iscorrect)
    }))
    return(assert)
}
