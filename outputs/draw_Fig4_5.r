library(abind)
load("data_Fig4_5.RData")
effcontact<-data_Fig4_5[["effcontact"]]
breakdown<-data_Fig4_5[["breakdown"]]
totalcontact<-data_Fig4_5[["totalcontact"]]
beta=0.201956234147969

draw_effcontact<-function(effcontact){
    bindEC<-abind(effcontact[[1]],effcontact[[2]],effcontact[[3]],effcontact[[4]],along=3)
    par(mfrow=c(2,2),cex=1,xaxt="s",xaxs="r",las=2)
    for(typ in 1:4){
        barplot(t(bindEC[bindEC[,1,typ]>=0,,typ][order(totalcontact[typ,bindEC[,1,typ]>=0]),]),ylim=c(0,2.5),col=rev(c("mediumpurple","tomato","dodgerblue","yellowgreen")),ylab="Effective household contact",xlab="",axes=F)
        axis(side=2,at=0:5/10/beta,labels=format(0:5/10,nsmall=2))
        legend(x=0,y=2.6, legend=c("Child","Mother","Father","Other"),ncol=2,lty=0,bty="n",pch=19,col=rev(c("mediumpurple","dodgerblue","tomato","yellowgreen")),y.intersp=0.6,x.intersp=0.5,cex=0.9)
    }
}
draw_riskbreakdown<-function(breakdown){
    bindRB<-abind(breakdown[[1]],breakdown[[2]],breakdown[[3]],breakdown[[4]],along=3)
    par(mfrow=c(2,2),cex=1,xaxt="s",xaxs="r",las=2)
    for(typ in 1:4){
        barplot(t(1-exp(-bindRB[bindRB[,2,typ]>=0,,typ][order(totalcontact[typ,bindRB[,2,typ]>=0]),])),ylim=c(0,c(0.3,rep(0.15,3))[typ]),col=rev(c("mediumpurple","tomato","dodgerblue","yellowgreen","grey")),ylab="Risk of infection",xlab="")
        legend(x=0,y=c(0.32,rep(0.155,3))[typ], legend=c("From Outside","From Child","From Mother","From Father","From Other"),ncol=2,lty=0,bty="n",pch=19,col=rev(c("mediumpurple","dodgerblue","tomato","yellowgreen","grey")),y.intersp=0.75,x.intersp=0.5,cex=0.9)
    }
}

draw_effcontact(effcontact)
draw_riskbreakdown(breakdown)