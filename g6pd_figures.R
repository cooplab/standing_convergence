source("~/Dropbox/postdocs/Peter/standing_parallelism/standing-variation-fns.R")
#install.packages("wesanderson")
library("wesanderson")

my.col<-wes_palette(name="Darjeeling",n=5)

######  Figure 4
###################################################
###################Char length as function of sd
###################################################
 mu=1e-6; rhos=c(0.2,2); sb=0.05; sdvals=seq(0.0001,0.3,length=800); 

sigmas=c(100,50,10)

layout(t(1:2))
par(mar=c(4,4.3,2,0.5))

for(j in 1:2){
	rho=rhos[j]
	if(j == 1) plot(sdvals,sdvals,type="n",ylim=c(100,800),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab=expression(list(chi,km)),cex.lab=1.5,cex.axis=1.5,main=expression(rho==0.2))
	if(j == 2) plot(sdvals,sdvals,type="n",ylim=c(100,800),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="",cex.lab=1.5,cex.axis=1.5,main=expression(rho==2))
	
	for(i in 1:3){
	
	sigma=sigmas[i]; charlength.sd<-sapply(sdvals, function (sd) { charLength(mu,rho,sb,sd,sigma)$value } )
	 lines(sdvals,charlength.sd,lwd=2,col=my.col[i])
	 abline(h= newCharLength(mu, rho, sb, sd=NA, sigma)$value,col=my.col[i],lwd=2,lty=2)  ##new mut. only
	}
	
	 ##contribution of standing var only, i.e. old variation
	 charlength.standing.only<-sapply(sdvals, function (sd) {  oldCharLength(mu, rho, sb, sd, sigma)$value})
	 lines(sdvals,charlength.standing.only ,lwd=2,col="black",lty=3)
	  if(j == 2) legend("topright",col=c("grey","grey","grey",my.col[1:3]),lwd=2,lty=c(2,3,1,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("No Standing var.","No new mut.","Both",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
 	}
 	
dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_charlengths.pdf",useDingbat=FALSE)

######  Figure 5
###################################################
###################Char length as function of sb
###################################################

 mu=1e-6; rhos=c(0.2,2); sd=0.05; sbvals=seq(0.00001,0.3,length=1600); 

sigmas=c(100,50,10)

layout(t(1:2))
par(mar=c(4,4.3,2,0.5))
for(j in 1:2){
	rho=rhos[j]
	if(j == 1) plot(sbvals,sbvals,type="n",ylim=c(100,1500),lwd=2,col=my.col[1],xlab=expression(s[b]),ylab=expression(list(chi, km)),cex.lab=1.5,cex.axis=1.4,main=expression(rho==0.2))
	if(j == 2) plot(sbvals,sbvals,type="n",ylim=c(100,800),lwd=2,col=my.col[1],xlab=expression(s[b]),ylab="",cex.lab=1.5,cex.axis=1.4,main=expression(rho==2))
	for(i in 1:3){
	sigma=sigmas[i]; charlength.sb<-sapply(sbvals, function (sb) { charLength(mu,rho,sb,sd=sd,sigma)$value } )
	 lines(sbvals,charlength.sb,lwd=2,col=my.col[i]); 
	 
	 charlength.new.only<-sapply(sbvals, function (sb) {  newCharLength(mu, rho, sb, sd=NA, sigma)$value})
	 lines(sbvals,charlength.new.only,col=my.col[i],lwd=2,lty=2)  ##new mut. only
	}
	 ##contribution of standing var only, i.e. old variation
	 charlength.standing.only<-sapply(sbvals, function (sb) {  oldCharLength(mu, rho, sb, sd=sd, sigma)$value})
	 lines(sbvals,charlength.standing.only ,lwd=2,col="black",lty=3)
	  if(j == 2) legend("topright",col=c("grey","grey","grey",my.col[1:3]),lwd=2,lty=c(2,3,1,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("No Standing var.","No new mut.","Both",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
 	}

dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_charlengths_sb.pdf",useDingbat=FALSE) 





#####Figure 6
###############################################################
####Mean time standing as func. of sd & sb
###############################################################
 layout(t(1:2))
par(mar=c(4,4.3,2,0.5))
 mu=1e-6; rhos=c(0.2,2); sd=0.05; sbvals=seq(0.001,0.3,length=800); 

for(j in 1:2){
	rho=rhos[j]
	if(j == 1) {	plot(sdvals,sdvals,type="n",ylim=c(0,250),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="Time, generations",cex.lab=1.5,cex.axis=1.5,main=expression(rho==0.2))
		}
	if(j == 2){ plot(sdvals,sdvals,type="n",ylim=c(0,250),lwd=2,col=my.col[1],xlab=expression(s[b]),ylab="",cex.lab=1.5,cex.axis=1.5,main=expression(rho==0.2))
}	

	for(i in 1:3){
		if(j == 1){ 
			sigma=sigmas[i]; 
			charlength.sd<-sapply(sdvals, function (sd) {  meanTime(mu, rho, sb, sd, sigma)$value } )
			 lines(sdvals,charlength.sd,lwd=2,col=my.col[i])
		 	abline(h= meanTime(mu, rho, sb, sd=1, sigma)$value,col=my.col[i],lwd=2,lty=2)  ##new mut. only
 		 	charlength.sd<-sapply(sdvals, function (sd) {  meanTime(mu, rho, sb, sd, sigma,include.new=FALSE)$value } )
	 		lines(sdvals,charlength.sd,lwd=2,lty=3,col=my.col[i])
		}

		if(j == 2){ 
			sigma=sigmas[i]; 
			charlength.sb<-sapply(sbvals, function (sb) {  meanTime(mu, rho, sb, sd=sd, sigma)$value } )
			 lines(sbvals,charlength.sb,lwd=2,col=my.col[i])

			charlength.sb<-sapply(sbvals, function (sb) {  meanTime(mu, rho, sb, sd=1, sigma)$value } )
			lines(sbvals,charlength.sb,col=my.col[i],lwd=2,lty=2)  ##new mut. only

 		 	charlength.sb<-sapply(sbvals, function (sb) {  meanTime(mu, rho, sb, sd, sigma,include.new=FALSE)$value } )
	 		lines(sbvals,charlength.sb,lwd=2,lty=3,col=my.col[i])
		}
	
	}
	

#	 charlength.standing.only<-sapply(sdvals, function (sd) {  oldCharLength(mu, rho, sb, sd, sigma)$value})
#	 lines(sdvals,charlength.standing.only ,lwd=2,col="black",lty=3)
	  if(j == 2) legend("topright",col=c("grey","grey","grey",my.col[1:3]),lwd=2,lty=c(2,3,1,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("No Standing var.","No new mut.","Both",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
 	}
 	
dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_chartimes_sd_sb.pdf",,useDingbat=FALSE)

##### Figure 7
###############################################################
####Proportion of range & alleles that are standing
###############################################################

layout(t(1:2))
par(mar=c(4,4.1,2,0.5))
for(j in 1:2){
		rho=rhos[j]
		plot(sdvals,charlength.sd,type="n",ylim=c(0,1),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="",cex.lab=1.5,cex.axis=1.5)

		if(j==1){ mtext(side=3,expression(rho==0.2),cex=1.5);mtext(side=2,line=3,"Proportion",cex=1.5)}
		if(j==2){ mtext(side=3,expression(rho==2),cex=1.5)}
		
			  if(j == 1) legend("topright",col=c(NA,"grey","grey",my.col[1:3]),lwd=2,lty=c(NA,1,2,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("Contribution of Standing Var. as:","Proportion of Area,","Proportion of Alleles.",
	 	lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
		
	for(i in 1:3){
		sigma=sigmas[i]
		prop.area.sd<-sapply(sdvals, function (sd){  standingProportionArea(mu, rho, sb, sd=sd, sigma)$value })
		prop.number.sd<-sapply(sdvals, function(sd){  standingProportionNumbers(mu, rho, sb, sd=sd, sigma)$value })
		lines(sdvals,prop.area.sd,lty=1,col=my.col[i],lwd=2)
	 	lines(sdvals,prop.number.sd,lty=2,col=my.col[i],lwd=2)
	}
}
dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_standing_var_proportion.pdf",useDingbat=FALSE)
 
 #######################
 ###  Figure 8
 ####effect of pleiotropy on standing variation
  more.cols<-wes_palette("FantasticFox",n=3)
  sd.1<-10^(-seq(-log10(0.05),5,length=100));
  plot(x=range(sd.1),y=c(0,1),type="n",log="x",xlab=expression(s[d1]),ylab=expression(p[1]),cex.lab=1.5,cex.axis=1.5)
 mu.2s<-c(1e-7,1e-6,1e-5)
 
 
 for(i in 1:length(mu.2s)){
  	prob<-sapply(sd.1,function(x){MultipleTypes(mu=c(1e-8,mu.2s[i]), rho=2, sb=c(0.05,0.05), sd=c(x,0.05), sigma=10)$value})
 	lines(sd.1,prob,col=more.cols[i],lty=1,lwd=2)
  	prob<-sapply(sd.1,function(x){MultipleTypes(mu=c(1e-8,mu.2s[i]), rho=2, sb=c(0.05,0.05), sd=c(x,0.05), sigma=100)$value})
 	lines(sd.1,prob,col=more.cols[i],lty=2,lwd=2)

  	prob<-sapply(sd.1,function(x){MultipleTypes(mu=c(1e-8,mu.2s[i]), rho=0.2, sb=c(0.05,0.05), sd=c(x,0.05), sigma=10)$value})
 	lines(sd.1,prob,col=more.cols[i],lty=3,lwd=2)
  	prob<-sapply(sd.1,function(x){MultipleTypes(mu=c(1e-8,mu.2s[i]), rho=0.2, sb=c(0.05,0.05), sd=c(x,0.05), sigma=100)$value})
 	lines(sd.1,prob,col=more.cols[i],lty=4,lwd=2)
 
	prob<-sapply(sd.1,function(x){MultipleTypes(mu=c(1e-8,mu.2s[i]), rho=2, sb=c(0.05,0.05), sd=c(x,0.05), 	sigma=20,standing.var.calc=TRUE)$value})
 	lines(sd.1,prob,lty=3,col=more.cols[i],type="b",lwd=2) 
 }
# legend(x="topright",legend=rep("blah",4),lty=1:4)

other.params<-as.expression( apply(cbind(c(2,2,0.2,0.2),c(10,100,10,100) ),1, function(x) {substitute(list( rho== myrho, sigma==mysigma), list(myrho=x[1],mysigma=x[2]))})) 
mut.rates<-as.expression( lapply(c(-7,-6,-5) , function(x) {substitute( mu[2]== 10^mymu, list(mymu=x))}))
my.exp<- as.expression(c(other.params,"Standing Var. Only",mut.rates))
 
 legend("topright",col=c(rep("grey",5),more.cols[1:3]),lwd=2,lty=c(1:4,rep(NA,4)),pch=c(rep(NA,4),1,rep(19,3)),legend=my.exp )
 dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/pleiotropy_calc.pdf",useDingbat=FALSE)


 
# ##areas with Malaria
 # cia<-read.csv("~/Downloads/rawdata_2147_new.csv",as.is=TRUE,head=FALSE)

# malaria<-read.csv("~/Dropbox/postdocs/Peter/Spatial_adaptation/malaria_present.csv",as.is=TRUE,head=FALSE)
# malaria<-malaria[-1,]
# malaria.eurasia<-malaria[malaria$V2 %in% c("Middle-East","SoutheastAsia","Euroasia","EasternAsia","CentralAsia","CentralAsia","Transcaucasia"),]

# malaria.eurasia[!(malaria.eurasia$V1 %in% cia$V2),]
# missing.these<-c("Laos","Korea North", "Korea South","Burma","Syria","Vietnam","Saudi Arabia")

# cia$V2[ cia$V2 %in% malaria.eurasia$V1 |   cia$V2 %in% missing.these]



# ######  Figure 
# ###############################################################
# ####Mean time standing as func. of sd
# ###############################################################
 # layout(t(1:2))
# par(mar=c(4,4.1,2,0.3))
# for(j in 1:2){
	# rho=rhos[j]
	# if(j == 1) plot(sdvals,sdvals,type="n",ylim=c(0,250),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="Time, generations",cex.lab=1.5,cex.axis=1.5,main=expression(rho==0.2))
	# if(j == 2) plot(sdvals,sdvals,type="n",ylim=c(0,250),lwd=2,col=my.col[1],xlab=expression(s[d]),ylab="",cex.lab=1.5,cex.axis=1.5,main=expression(rho==2))
	
	# for(i in 1:3){
	
	# sigma=sigmas[i]; 
	# charlength.sd<-sapply(sdvals, function (sd) {  meanTime(mu, rho, sb, sd, sigma)$value } )
	 # lines(sdvals,charlength.sd,lwd=2,col=my.col[i])
	 # abline(h= meanTime(mu, rho, sb, sd=1, sigma)$value,col=my.col[i],lwd=2,lty=2)  ##new mut. only
 	 # charlength.sd<-sapply(sdvals, function (sd) {  meanTime(mu, rho, sb, sd, sigma,include.new=FALSE)$value } )
	 # lines(sdvals,charlength.sd,lwd=2,lty=3,col=my.col[i])

	# }
	

# #	 charlength.standing.only<-sapply(sdvals, function (sd) {  oldCharLength(mu, rho, sb, sd, sigma)$value})
# #	 lines(sdvals,charlength.standing.only ,lwd=2,col="black",lty=3)
	  # if(j == 2) legend("topright",col=c("grey","grey","grey",my.col[1:3]),lwd=2,lty=c(2,3,1,rep(NA,3)),pch=c(NA,NA,NA,rep(19,3)),legend= as.expression( c("No Standing var.","No new mut.","Both",
	 	# lapply(sigmas , function(x){  {substitute( sigma==myS, list(myS=x) )} }))))
 	# }
 	
# dev.copy2pdf(file="~/Dropbox/postdocs/Peter/standing_parallelism/G6PD_chartimes.pdf")
