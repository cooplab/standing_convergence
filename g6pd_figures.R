
source('~/Dropbox/postdocs/Peter/Spatial_adaptation/Spatial_adaptation/standing-variation-fns.R')
#install.packages("wesanderson")
library("wesanderson")

my.col<-wes.palette(5, "Darjeeling")  #wes.palette(2, "Royal1")

 mu=1e-6; rhos=c(0.2,2); sb=0.05; sdvals=seq(0.0001,0.3,length=800); 

sigmas=c(100,50,10)

layout(t(1:2))
par(mar=c(4,4.1,2,0.3))
for(j in 1:2){
	rho=rhos[j]
	if(j == 1) plot(sdvals,sdvals,type="n",ylim=c(100,800),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab=expression(chi),cex.lab=1.5,cex.axis=1.5,main=expression(rho==0.2))
	if(j == 2) plot(sdvals,sdvals,type="n",ylim=c(100,800),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="",cex.lab=1.5,cex.axis=1.5,main=expression(rho==2))
	
	for(i in 1:3){
	
	sigma=sigmas[i]; charlength.sd<-sapply(sdvals, function (sd) { charLength(mu,rho,sb,sd,sigma)$value } )
	 lines(sdvals,charlength.sd,lwd=2,col=my.col[i])
	 abline(h= newCharLength(mu, rho, sb, sd=NA, sigma)$value,col=my.col[i],lwd=2,lty=2)  ##new mut. only
	}
	
	 ##contribution of standing var only, i.e. old variation
	 charlength.standing.only<-sapply(sdvals, function (sd) {  oldCharLength(mu, rho, sb, sd, sigma)$value})
	 lines(sdvals,charlength.standing.only ,lwd=2,col="black",lty=3)
	  if(j == 2) legend("topright",col=c("grey","grey","black",my.col[1:3]),lty=c(2,3,1,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("No Standing var.","No new mut.","Both",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
 	}
 	
dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_charlengths.pdf")
 
###############################################################
####Mean time standing
###############################################################
 layout(t(1:2))
par(mar=c(4,4.1,2,0.3))
for(j in 1:2){
	rho=rhos[j]
	if(j == 1) plot(sdvals,sdvals,type="n",ylim=c(0,250),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="Generations",cex.lab=1.5,cex.axis=1.5,main=expression(rho==0.2))
	if(j == 2) plot(sdvals,sdvals,type="n",ylim=c(0,250),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="",cex.lab=1.5,cex.axis=1.5,main=expression(rho==2))
	
	for(i in 1:3){
	
	sigma=sigmas[i]; 
	charlength.sd<-sapply(sdvals, function (sd) {  meanTime(mu, rho, sb, sd, sigma)$value } )
	 lines(sdvals,charlength.sd,lwd=2,col=my.col[i])
	 abline(h= meanTime(mu, rho, sb, sd=1, sigma)$value,col=my.col[i],lwd=2,lty=2)  ##new mut. only
 	 charlength.sd<-sapply(sdvals, function (sd) {  meanTime(mu, rho, sb, sd, sigma,include.new=FALSE)$value } )
	 lines(sdvals,charlength.sd,lwd=2,lty=3,col=my.col[i])

	}
	

#	 charlength.standing.only<-sapply(sdvals, function (sd) {  oldCharLength(mu, rho, sb, sd, sigma)$value})
#	 lines(sdvals,charlength.standing.only ,lwd=2,col="black",lty=3)
	  if(j == 2) legend("topright",col=c("grey","grey","grey",my.col[1:3]),lty=c(2,3,1,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("No Standing var.","No new mut.","Both",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
 	}
 	
dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_chartimes.pdf")

 
###############################################################
####Proportion of range & alleles that are standing
###############################################################

layout(t(1:2))
par(mar=c(4,4.1,2,0.3))
for(j in 1:2){
		rho=rhos[j]
		plot(sdvals,charlength.sd,type="n",ylim=c(0,1),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="",cex.lab=1.5,cex.axis=1.5)

		if(j==1){ mtext(side=3,expression(rho==0.2),cex=1.5);mtext(side=2,line=3,"Proportion",cex=1.5)}
		if(j==2){ mtext(side=3,expression(rho==2),cex=1.5)}
		
			  if(j == 1) legend("topright",col=c(NA,"grey","grey",my.col[1:3]),lty=c(NA,1,2,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("Contribution of Standing Variation as:","Proportion of Area,","Proportion of Alleles.",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
		
	for(i in 1:3){
		sigma=sigmas[i]
		prop.area.sd<-sapply(sdvals, function (sd){  standingProportionArea(mu, rho, sb, sd=sd, sigma)$value })
		prop.number.sd<-sapply(sdvals, function(sd){  standingProportionNumbers(mu, rho, sb, sd=sd, sigma)$value })
		lines(sdvals,prop.area.sd,lty=1,col=my.col[i],lwd=2)
	 	lines(sdvals,prop.number.sd,lty=2,col=my.col[i],lwd=2)
	}
}
dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_standing_var_proportion.pdf")
 
 
 
# as.expression( c(lapply(rate , function(x) substitute(2*N*nu==myS, list(myS=x) ) )
 as.expression( c(lapply(rate , function(x) substitute(2*N*nu==myS, list(myS=x) ) ),"Top-hat","Step","Theory"))
 ##areas with Malaria
 cia<-read.csv("~/Downloads/rawdata_2147_new.csv",as.is=TRUE,head=FALSE)

malaria<-read.csv("~/Dropbox/postdocs/Peter/Spatial_adaptation/malaria_present.csv",as.is=TRUE,head=FALSE)
malaria<-malaria[-1,]
malaria.eurasia<-malaria[malaria$V2 %in% c("Middle-East","SoutheastAsia","Euroasia","EasternAsia","CentralAsia","CentralAsia","Transcaucasia"),]

malaria.eurasia[!(malaria.eurasia$V1 %in% cia$V2),]
missing.these<-c("Laos","Korea North", "Korea South","Burma","Syria","Vietnam","Saudi Arabia")

cia$V2[ cia$V2 %in% malaria.eurasia$V1 |   cia$V2 %in% missing.these]