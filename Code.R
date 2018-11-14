############################################################################
# ASL data Analyses R-code
############################################################################

############################################################################
#  0.0  Initialize working Environment                                          
############################################################################
rm(list=ls(all=TRUE))
# Include mgcv library
library(mgcv)
library(gplots)
# Specify working directory: Where your file is located 
Data_dir <- 'C:/Projects/Yukon_River/JTC/Size Committee/Border correction/Data Analyses/'
setwd(Data_dir)

# Enter the name of data file  
data_file1 <- 'Canada_fishwheel_Elizabeth.csv'
data_file2 <- 'Eagle_out.csv'
data_file3 <- 'Current_border_age.csv' 

#Bromaghin Correction factor
brcf <- c(0.049967,0.015221,0.209598,0.636551,1,0.702612)

############################################################################
#   1.0 Read source file 
############################################################################
windows(width=20,height=20,record=TRUE)

dataclean<-function(data_dir,data_file){
size <- read.csv(paste0(data_dir, data_file), na.strings ='', stringsAsFactors = FALSE,header = TRUE)
# Rename data
names(size) <- c('Year','Sex','Gear','Mesh','Species','Run','Age','Length','N')
# Remove data mesh size is less than equal to 1.
#size <- size[size$Mesh>1,]
#size <- size[which(size$Length > 300 & !is.na(size$Age)),]
# Fall chum is Chum salmon of Fall run 
size$Species[with(size, Species=='Chum' & Run =='Fall Run')] <- 'Fall Chum'
# Summer chum is Chum salmon of Summer  run 
size$Species[with(size, Species=='Chum' & Run =='Summer Run')] <- 'Summer Chum'
# Flatten data 
size <- with(size, data.frame(rep(Species,N),rep(Year,N),rep(Gear,N),rep(Mesh,N),rep(Sex,N),rep(Age,N),rep(Length,N)))
names(size) <- c('Species','Year','Gear','Mesh','Sex','Age','Length')
size 
}

############################################################################
#   2.0 Length Selectivity function 
############################################################################  
select <- function(par,lpr){
    sigma <- par[1]
    theta <- par[2]
    lamda <- par[3]
	tau <- par[4]
  dm1 <- lamda/(2.0 * theta)
  dm2 <- (1 + dm1^2)^theta	
  dm3 <- (lpr - sigma * dm1 - tau) / sigma
  dm4 <- (1 + dm3^2)^(-theta)
  dm5 <- exp(-lamda * (atan(dm3)+ atan(dm1)))
  g <- dm2 * dm4 * dm5     
  return(g)
  }  
############################################################################

############################################################################
#   2.1 Length Selectivity likelihood function 
############################################################################  
likelihood.s <- function(par,likedat){
  likedat1 <- likedat[which(likedat$loc=='fishwheel'),c('Year','Length','n','lpr')]
  likedat2 <- likedat[which(likedat$loc=='Eagle'),c('Year','Length','n')]
  sigma <- par[1]
  theta <- par[2]
  lamda <- par[3]
  tau <- 0.0001	
  par.1 <- c(sigma,theta,lamda,tau)     
  likedat1$g <- select(par.1,likedat1[ ,'lpr'] )  
# Add tangle parameter
#    likedat1$f <- pmax((tau<likedat1[,'lpr'])*tangle,g)
# Calculate selectivity adjusted nunmber
    likedat1$nf <- with(likedat1,n/g)
# Calculate sum
     sums <- aggregate(nf~Year,sum,data=likedat1)
# merge 
	likedat1 <- merge(likedat1,sums,by=c('Year'))
# merge with Eagle data (Assume Eagle length comp is correct)	
	likedat2 <- merge(likedat2,likedat1,by=c('Year','Length'))
	likelihood <- with(likedat2,-n.x*log(nf.x/nf.y))
#### Likelihood calculation  ###############################################    
loglink  <- sum(likelihood)
return(loglink) 
}

############################################################################
#   2.2 Sex separated Length Selectivity likelihood function 
############################################################################  
likelihood.as <- function(par,likedat){
  likedat1 <- likedat[which(likedat$loc=='fishwheel'),c('Year','Sex','Length','n','lpr')]
  likedat2 <- likedat[which(likedat$loc=='Eagle'),c('Year','Sex','Length','n')]
  sigma.m <- par[1]
  theta.m <- par[2]
  lamda.m <- par[3]
  sigma.f <- par[4]
  theta.f <- par[5]
  lamda.f <- par[6]
  tau <- 0.001	
  par.1 <- c(sigma.m,theta.m,lamda.m,tau)  
  par.2 <- c(sigma.f,theta.f,lamda.f,tau) 
  likedat1$g <- with(likedat1,ifelse(Sex=='male',
			   select(par.1,likedat1[,'lpr']),
			   select(par.2,likedat1[,'lpr'])))
# Add tangle parameter
#    likedat1$f <- pmax((tau<likedat1[,'lpr'])*tangle,g)
# Calculate selectivity adjusted nunmber
    likedat1$nf <- with(likedat1,n/g)
# Calculate sum
     sums <- aggregate(nf ~ Year, sum, data = likedat1)
# merge 
	likedat1 <- merge(likedat1, sums, by = c('Year'))
# merge with Eagle data (Assume Eagle length comp is correct)	
	likedat2 <- merge(likedat2,likedat1,by=c('Year','Length'))
	likelihood <- with(likedat2,-n.x*log(nf.x/nf.y))
#### Likelihood calculation  ###############################################    
loglink  <- sum(likelihood)
return(loglink) 
}

############################################################################
#   2.3 AS Selectivity  likelihood function 
############################################################################  
   likelihood.b <- function(par,likedat){
    likedat1 <- likedat[which(likedat$loc=='fishwheel'),c('Year','Sex','Age','p')]
    likedat2 <- likedat[which(likedat$loc=='Eagle'),c('Year','Sex','Age','Length.x','p')]
	m.4 <-  1
	m.5 <-  par[1]
	m.6 <-  par[2]
	m.7 <-  par[3]
	f.4 <-  par[4]
	f.5 <-  par[5]
	f.6 <-  par[6]
	f.7 <-  par[7]

# Calculate selectivity adjusted nunmber
   likedat1$nf <- with(likedat1,ifelse(Sex=='male'& Age==4,p/m.4,
	ifelse(Sex=='male'& Age==5,p/m.5,
	ifelse(Sex=='male'& Age==6,p/m.6,
	ifelse(Sex=='male'& Age==7,p/m.7,
	ifelse(Sex=='female'& Age==4,p/f.4,
	ifelse(Sex=='female'& Age==5,p/f.5,
	ifelse(Sex=='female'& Age==6,p/f.6,p/f.7))))))))
# merge 
   sums <-aggregate(nf~Year,sum,data=likedat1)
	likedat1 <- merge(likedat1,sums,by=c('Year'))
	likedat1$p <- with(likedat1,nf.x/nf.y)
	likedat1 <- merge(likedat2,likedat1,by=c('Year','Sex','Age'))
# merge with Eagle data (Assume Eagle length comp is correct)	
	likelihood <- with(likedat1,-Length.x*log(p.y))
#### Likelihood calculation  ###############################################    
loglink  <- sum(likelihood)
return(loglink) 
}


############################################################################
# 3.0 read ASL table and clean up for analyses. 
############################################################################
# Read fish wheel data 
cafw <- read.csv(paste0(Data_dir,data_file1),na='')
cafw$loc <- 'fishwheel'
# Convert fork length to MEF length
cafw$Length <- 1.446+0.898*cafw$Length
# Extract only data needed 
cafw <- cafw[,c('Year','Sex','Age','Length','loc')]
aggregate(Length~Year,min, data=cafw)
aggregate(Length~Year,max, data=cafw)
aggregate(Length~Year,median, data=cafw)
aggregate(Length~Age+Sex,mean, data=cafw[cafw$Year>2006,])


# Read Eagle Sonar data
eagl <- dataclean(Data_dir,data_file2)
# Get only Chinook data 
eagl <- eagl[eagl$Species =='Chinook',]
# Change unknown Sex to NA
eagl$Sex[eagl$Sex=='unknown'] <- NA
# Add location name 
eagl$loc <- 'Eagle'
# Restrict mesh size
#eagl <- eagl[eagl$Mesh == 5.25|eagl$Mesh == 7.5,]
#temp <-aggregate(Length~Year+Mesh,length, data=eagl)
# Extract only data needed 
eagl <- eagl[,c('Year','Sex','Age','Length','loc')]
# Combine fish wheel and test fishery data 
esc.u <- rbind(cafw,eagl)


######################################################################%#####
# 4.0 Prepare data for likelihood analyes:  
############################################################################  
# Get data when both fishwheel and test-fish were operated 
yr <- c(2005,2006,2007,2008,2010,2011,2012)
likedata.o <- esc.u[which(esc.u$Year %in% yr),]
# Combine Age 3 and 4:  Age 3 is minimal  
likedata.o$Age[which(likedata.o$Age ==3)]<- 4

######################################################################%#####
# 4.1 Obtain summarized data  
############################################################################  
# Get number of samples by location, year. age. sex
as <-aggregate(Length~loc+Year+Age+Sex,length, data=likedata.o)
# Get number of samples by location, year
as.t <-aggregate(Length~loc+Year,sum, data=as)
# Get number of samples by location, year
as.a <-aggregate(Length~loc+Year+Age,length, data=likedata.o)
# Get totoal number of samples by location, year
as.p <- merge(as,as.a,by=c('loc','Year','Age'))
# Merge the by location, year
as.tp <- merge(as,as.t,by =c('loc','Year'))
# Calculated Age Sex proporion by year 
as.tp$p <- with(as.tp,Length.x/Length.y)
as.p$p <- with(as.p,Length.x/Length.y)


######################################################################%#####
# 4.2 Diagnostic Graohics between Fishwheel and Gillnet   
############################################################################  
######################################################################%#####
# 4.2.1 Compare Age proportion   
############################################################################  
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(i in 1:length(yr)){
dat1 <- as[which(as$Year==yr[i]),]
tbl <- xtabs(Length~ loc+Age, data=dat1)
barplot(t(prop.table(tbl,1)),main=yr[i],ylim=c(0,1),names.arg=c('Test fishery','Fish wheel'))
#print(paste(yr[i]))
#print(tbl)
#print(chisq.test(tbl))
}
plot.new()
legend('topleft',legend=c(4,5,6,7),fill=gray.colors(4), bg='white', horiz=FALSE,bty='n')
###########  Add Texts  ####################################################
#mtext('Border Eagle Test Fishery vs Fish wheel',line = 1.2 ,side = 3, outer = TRUE)
mtext('Age Proportion (4,5,6,7)', side = 2, line = 1.8, outer = TRUE)
#mtext("Age", side = 1, line = 1.8, outer = TRUE)
############################################################################

######################################################################%#####
# 4.2.2 Compare Sex-proportion  
############################################################################  
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(i in 1:length(yr)){
dat1 <- as[which(as$Year==yr[i]),]
tbl <- xtabs(Length~ loc+Sex, data=dat1)
barplot(t(prop.table(tbl,1)),main=yr[i],ylim=c(0,1),names.arg=c('Test fishery','Fish wheel'))
#print(paste(yr[i]))
#print(tbl)
#print(chisq.test(tbl))
}
plot.new()
legend('topleft',legend=c('Female','Male'),fill=gray.colors(2), bg='white', horiz=FALSE,bty='n')
#
###########  Add Texts  ####################################################
#mtext('Border Eagle vs Fish wheel',line = 1.2 ,side = 3, outer = TRUE)
mtext('Female Proportion', side = 2, line = 1.8, outer = TRUE)
#mtext("Age", side = 1, line = 1.8, outer = TRUE)
############################################################################


######################################################################%#####
# 4.2.2 Compare Sex-proportion by Age 
############################################################################  
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(i in 1:length(yr)){
dat1 <- as.p[which(as.p$Year==yr[i]),]
dat1 <- dat1[which(dat1$Sex=='female'),]
tbl <- xtabs(p~ loc+Age, data=dat1)
barplot((tbl),main=yr[i],beside = TRUE,ylim=c(0,1),legend =rownames(tbl),args.legend = list(x='topleft',bty='n'))
print(paste(yr[i]))
print(tbl)
}
###########  Add Texts  ####################################################
mtext('Border Eagle vs Fishwheel',line = 1.2 ,side = 3, outer = TRUE)
mtext('Female Proportion', side = 2, line = 1.8, outer = TRUE)
mtext("Age", side = 1, line = 1.8, outer = TRUE)
############################################################################


######################################################################%#####
# 4.2.3 Compare Mean length by age   
############################################################################  
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(j in 1:length(yr)){
dat1 <- likedata.o[which(likedata.o$Year==yr[j]),]
dat1 <- dat1[,c('Length','Age','loc')]
boxplot(Length~loc+Age,data=dat1,main=yr[j],col=c('white','grey'),ylim=c(450,1000),xaxt='n')
axis(1,at=seq(1,8,by=2),labels=c(4,5,6,7))
}
plot.new()
legend("topleft", c("Test fishery","Fish wheel"), fill=c('white','grey'),bty='n')
###########  Add Texts  ####################################################
#mtext('Border Eagle vs Fishwheel',line = 1.2 ,side = 3, outer = TRUE)
mtext('Length (mm)', side = 2, line = 1.8, outer = TRUE)
mtext("Age", side = 1, line = 1.8, outer = TRUE)
############################################################################

######################################################################%#####
# 4.2.4 Compare Mean length by age: Male    
############################################################################  
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(j in 1:length(yr)){
dat1 <- likedata.o[which(likedata.o$Year==yr[j]),]
dat1 <- dat1[which(dat1$Sex=='male'),]
dat1 <- dat1[,c('Length','Age','loc')]
boxplot(Length~loc+Age,data=dat1,main=yr[j],col=c('white','grey'),ylim=c(450,1000),xaxt='n')
axis(1,at=seq(1,8,by=2),labels=c(4,5,6,7))
legend("topleft", c("Test fishery","Fish wheel"), fill=c('white','grey'),bty='n')
}
###########  Add Texts  ####################################################
mtext('Border Eagle vs Fishwheel',line = 1.2 ,side = 3, outer = TRUE)
mtext('Male Length (mm)', side = 2, line = 1.8, outer = TRUE)
mtext("Age", side = 1, line = 1.8, outer = TRUE)
############################################################################

######################################################################%#####
# 4.2.5 Compare Mean length by age: Female
############################################################################  
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(j in 1:length(yr)){
dat1 <- likedata.o[which(likedata.o$Year==yr[j]),]
dat1 <- dat1[which(dat1$Sex=='female'),]
dat1 <- dat1[,c('Length','Age','loc')]
boxplot(Length~loc+Age,data=dat1,main=yr[j],col=c('white','grey'),ylim=c(450,1000),xaxt='n')
axis(1,at=seq(1,8,by=2),labels=c(4,5,6,7))
}
plot.new()
legend("topleft", c("Eagle","Fishwheel"), fill=c('white','grey'),bty='n')

###########  Add Texts  ####################################################
mtext('Border Eagle vs Fishwheel',line = 1.2 ,side = 3, outer = TRUE)
mtext('female Length (mm)', side = 2, line = 1.8, outer = TRUE)
mtext("Age", side = 1, line = 1.8, outer = TRUE)
############################################################################

######################################################################%#####
# 5.0  Likelihood Analyses:  Estimate fishwheel selectivity   
############################################################################  
likedat <- likedata.o
######################################################################%#####
# 5.1  Data preparation   
############################################################################  
# Change fish size by 20mm
likedat$Length <- round(likedat$Length/20+0.0001)*20
# Calculate LPR
likedat$lpr <- likedat$Length/(5.25*25.4*2)
# Calculate freqiemcy by length
likedata.s <-aggregate(Age~loc+Year+Length+lpr, length, data=likedat,na.action=na.pass)
# Calculate < total frequency
sums <-aggregate(Age~loc+Year, sum, data=likedata.s)
# Get proportion by length 
likedata.s <- merge(likedata.s,sums, by=c('loc','Year'))
# p is a proportion by length 
likedata.s$p <- with(likedata.s, Age.x/Age.y)
# Calculate LPR 
#likedata.s$lpr <- likedata.s$Length/(5.25*25.4*2)
# Rename
names(likedata.s) <- c('loc','Year','Length','lpr','n','N','p')

######################################################################%#####
# 5.2  Likelihood calculation  
############################################################################  
### set  initial values and bounds #########################################
init <- c(0.8,0.5,2.4)
lb <-  c(0.0001,0.001,-4.0)
ub <-  c(1.5,10.0,5.0)

########## Calculate Likelihood ############################################
ptm <- proc.time()
nll.fl <- optim(par=init,fn=likelihood.s,method="L-BFGS-B",lower=lb, upper = ub, 
    likedat=likedata.s[likedata.s$Year>2006, ], hessian = TRUE)
proc.time() - ptm
Rprof()
nll.fl
#1: Hessian Matrix
hessian <- nll.fl$hessian
est <- nll.fl$par

# Create a variance-covariance matrix
var_covar_mat <- solve(hessian)
# Pull out diagonal nu
var_obs <- diag(var_covar_mat)
# Calculate standard error
std_err <- sqrt(var_obs)
############################################################################
 
############################################################################
#   5.3  Diagnostic Graphics
############################################################################  
########## Calculate selectivity ###########################################
# Separate data between Eagle and fishwheel
   likedat.fl <- likedata.s[which(likedata.s$loc=='fishwheel'),c('Year','Length','p','n','lpr')]
   likedat2 <- likedata.s[which(likedata.s$loc=='Eagle'),c('Year','Length','p','n','lpr')]
# Add selectivity to fishwheel data 
   likedat.fl$sel <- select(c(nll.fl$par,0.001),likedat.fl$lpr)
   likedat.fl$nf <- with(likedat.fl,n/sel)
   sums <-aggregate(nf~Year,sum,data=likedat.fl)
# merge data 
	likedat.fl <- merge(likedat.fl,sums,by=c('Year'))
# Calculate p.sel: selectivity adjusted length proportion	
	likedat.fl$p.sel <- with(likedat.fl,nf.x/nf.y)
# 
temp <- merge(likedat2,likedat.fl, by = c('Year', 'Length'), all.y=TRUE)
temp[is.na(temp)] <- 0
temp$gf <- with(temp,abs(p.x-p.sel))
gf.sums <-aggregate(gf~Year,mean,data=temp)

############################################################################
# 5.4  Graphics 
############################################################################
windows(width=10,height=10, record = TRUE)
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.7)
	lpr <- seq(0,5,by=0.05)
    plot(lpr*(5.25*25.4*2),select(c(nll.fl$par,0.001),lpr),type='l',xlim=c(300,1100),main='Fish wheel Selectivity')
for(i in 1:length(yr)){
	dat1 <- likedat.fl[which(likedat.fl$Year==yr[i]),]
	dat1 <- dat1[order(dat1$Length),]
	dat2 <- likedat2[which(likedat2$Year==yr[i]),]
	dat2 <- dat2[order(dat2$Length),]
	plot(dat2$Length,dat2$p,type='l',xlim=c(300,1100),ylim=c(0,0.15),col='gray40',lwd=2,main=yr[i])
	lines(dat1$Length,dat1$p,col=1,lty=2)
	lines(dat1$Length,dat1$p.sel,col=1,lwd=2)
	}    
    plot.new()
	legend("topleft",bty='n',xjust=0,lty =c(1,1,2), lwd = c(2,2,1), col=c('gray40',1,1), legend = c('Test fishery','Fish wheel adjusted','Fish wheel'))

###########  Add Texts  ####################################################
#mtext('Border Fishwheel Selectivity',line = 1.2 ,side = 3, outer = TRUE)
mtext('Length proportion', side = 2, line = 1.8, outer = TRUE)
mtext("METF Length (mm)", side = 1, line = 1.8, outer = TRUE)
############################################################################

############################################################################
# 6.0  Sex Separated Likelihood Analyses:  Estimate Correctuon factor    
############################################################################  
### set  initial values and bounds #########################################
init <- c(0.3,0.8,0.1,0.4,0.1,0.1,0.05)
lb <-  c(0.00001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001)
ub <- c(1,1,1,1,1,1,1)

########## Calculate Likelihood ############################################
ptm <- proc.time()
nll <- optim(par=init,fn=likelihood.b,method="L-BFGS-B",lower=lb, upper = ub, likedat=as.tp[as.tp$Year>2006,], hessian = TRUE)
min_NLL <- nll$value
proc.time() - ptm
Rprof()
nll
#1: Hessian Matrix
hessian <- nll$hessian
est <- nll$par

# Create a variance-covariance matrix
var_covar_mat <- solve(hessian)
# Pull out diagonal nu
var_obs <- diag(var_covar_mat)
# Calculate standard error
std_err <- sqrt(var_obs)


############################################################################
#   6.0 Estimate Corretion facotor 
############################################################################  
# Estimate fishwheel adjusted 

likedat1 <- as.tp[which(as.tp$loc=='fishwheel'),c('loc','Year','Sex','Age','p')]
as.e <- as.tp[which(as.tp$loc=='Eagle'),c('Year','Sex','Age','Length.x','p')]
names(as.e) <- c('Year','Sex','Age','n','p')
as.sel <- c(1,nll$par)

likedat1$nf <- with(likedat1,ifelse(Sex=='male'& Age==4,p/as.sel[1],
  ifelse(Sex=='male'& Age==5,p/as.sel[2],ifelse(Sex=='male'& Age==6,p/as.sel[3],
  ifelse(Sex=='male'& Age==7,p/as.sel[4],ifelse(Sex=='female'& Age==4,p/as.sel[5],
  ifelse(Sex=='female'& Age==5,p/as.sel[6],ifelse(Sex=='female'& Age==6,p/as.sel[7],p/as.sel[8]))))))))
# merge 
	sums <-aggregate(nf~loc+Year,sum,data=likedat1)
	likedat1 <- merge(likedat1,sums,by=c('loc','Year'))
	likedat1$p.b <- with(likedat1,nf.x/nf.y)
    as.b <- likedat1[,c('Year','Sex','Age','p','p.b')]

############################################################################
#   6.1 Compare performace between net-selectivty and correction factors 
############################################################################  
# Separate data between Eagle and fishwheel
#fl <- likedata.o[which(likedata.o$loc=='fishwheel'),]
fl <- likedat[which(likedat$loc=='fishwheel'), ]
# Calculate LPR 
fl$lpr <- fl$Length/(5.25*25.4*2)
# Calculate selectivity
# Add selectivity to fishwheel data 
#fl$sel <- with(fl,ifelse(Sex=='male', 
#               select(c(nll.s[1:3],0.001),lpr),
#			   select(c(nll.s[4:6],0.001),lpr)))	
fl$sel <- select(c(nll.fl$par,0.001),fl$lpr)
# Calculate selction correciton factor  
fl$nsel <- 1/fl$sel
# Calculate selection corrected n by Year, Age, Sex
as.s <-aggregate(nsel~Year+Age+Sex, sum, data=fl)
# Calculate selection corrected n by Year
ass.s <-aggregate(nsel~Year, sum,data=as.s)
# mearge two data 
as.s <- merge(as.s,ass.s,by=c('Year'))
# Calculate selectivity adjusted AS proportion
as.s$p.s <- with(as.s,nsel.x/nsel.y)
# Keep only data needed 
as.s <- as.s[,c('Year','Age','Sex','p.s')]
# Get correction factor corrected AS proportion
# Merge two data 
as.c <- merge(as.s,as.b,by=c('Year','Sex','Age'))
names(as.c) <- c('Year','Sex','Age','p.s','p.f','p.b')
# Merge Eagle as data 
as.c <- merge(as.e,as.c,by=c('Year','Sex','Age'),all=TRUE)
#as.c <- merge(as.c,as.t[as.t$loc=='Eagle',2:3],by='Year')
as.c[is.na(as.c) ]<- 0
as.c$p <- ifelse(is.na(as.c$p),0,as.c$p)
as.c$ls <- with(as.c,abs(p-p.s))
as.c$lb <- with(as.c,abs(p-p.b))
as.c$lf <- with(as.c,abs(p-p.f))
as.c$sb <- with(as.c,abs(p.s-p.b))
as.like <- aggregate(cbind(p,p.f,p.b,p.s)~Age+Sex,FUN=mean, data=as.c[as.c$Year > 2006, ])
as.like <- aggregate(cbind(ls,lb,sb)~Age+Sex,FUN=mean, data=as.c[as.c$Year > 2006, ])
as.like <- aggregate(cbind(lf,lb,ls,sb)~Year,FUN=mean, data=as.c)


age.comp <-aggregate(cbind(p,p.f,p.b,p.s) ~ Year+Age, sum, data=as.c)
age.comp$ls <- with(age.comp,(abs(p-p.s)))
age.comp$lb <- with(age.comp,(abs(p-p.b)))
age.comp$lf <- with(age.comp,(abs(p-p.f)))
age.comp$sb <- with(age.comp,(abs(p.s-p.b)))
age.like <- aggregate(cbind(ls,lb,lf,sb)~Year,FUN=sum, data=age.comp)

sex.comp <-aggregate(cbind(p,p.f,p.b,p.s) ~ Year+Sex, sum, data=as.c)
sex.comp$ls <- with(sex.comp,(abs(p-p.s)))
sex.comp$lb <- with(sex.comp,(abs(p-p.b)))
sex.comp$lf <- with(sex.comp,(abs(p-p.f)))


############################################################################
par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(i in 1:length(yr)){
dat1 <- age.comp[which(age.comp$Year==yr[i]),]
barplot(as.matrix(dat1[3:6]),main=yr[i],ylim=c(0,1),col=rainbow(4),names.arg=c('TF','FW','AS','LS'))
# use rainbow(4) for color 
}
plot.new()
legend('topleft',legend=c(4,5,6,7),fill=rainbow(4), bg='white', horiz=FALSE, bty='n')
###########  Add Texts  ####################################################
mtext('Age proportion', side = 2, line = 1.8, outer = TRUE)

par(mfrow=c(3,3), mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
for(i in 1:length(yr)){
dat1 <- sex.comp[which(sex.comp$Year==yr[i]),]
barplot(as.matrix(dat1[3:6]),main=yr[i],ylim=c(0,1),col=rainbow(2),names.arg=c('TF','FW','AS','LS'))
}
plot.new()
legend('topleft',legend=c('female','male'),fill=rainbow(2), bg='white', horiz=FALSE,bty='n')
###########  Add Texts  ####################################################
mtext('Sex proportion', side = 2, line = 1.8, outer = TRUE)


############################################################################
#   7.0 Create historical AS comp 
############################################################################  
# Extract data fishwheel before 2005
fs <- cafw[cafw$Year<2007,]
# Change fish size by 20mm
fs$Length <- round(fs$Length/20+0.0001)*20
fs$lpr <- fs$Length/(5.25*25.4*2)
# Calculate selectivity
# Add selectivity to fishwheel data 
fs$sel <- select(c(nll.fl$par,0.001),fs$lpr)
# Calculate selction correciton factor  
fs$nsel <- 1/fs$sel
fs.f <- fs

# Get fish wheel full ages 
as <-aggregate(nsel~Year+Sex+Age,length,data=fs.f)
# Get number of samples by location, year
as.t <-aggregate(nsel~Year,sum, data=as)
as.p <- merge(as,as.t,by=c('Year'))
# Calculate selectivity adjusted AS proportion
as.p$p.f <- with(as.p,nsel.x/nsel.y)
as.p <- as.p[,c('Year','Sex','Age','p.f')]

# Get selectivity method full ages 
as.f <-aggregate(nsel~loc+Year+Sex+Age,sum,data=fs.f)
# Get number of samples by location, year
as.f.t <-aggregate(nsel~loc+Year,sum, data=as.f)
as.f.p <- merge(as.f,as.f.t,by=c('Year'))
# Calculate selectivity adjusted AS proportion
as.f.p$p.s <- with(as.f.p,nsel.x/nsel.y)
as.f.p <- as.f.p[,c('Year','Sex','Age','p.s')]
as.f.p.w <- reshape(as.f.p,idvar = c('Sex','Year'),timevar = c(('Age')),direction='wide')
# merge 
as.f.p <- merge(as.p,as.f.p,by=c('Year','Sex','Age'))
	

age.comp <-aggregate(cbind(p.f,p.s) ~ Year+Age, sum, data=as.f.p)
sex.comp <-aggregate(cbind(p.f,p.s) ~ Year+Sex, sum, data=as.f.p)


# Age-sex selectivity method  
# Combine Age 3 and 4:  Age 3 is minimal  
fs$Age[which(fs$Age ==3)]<- 4
fs$Age[which(fs$Age ==8)]<- 7
as <-aggregate(loc~Year+Age+Sex,length, data=fs)
# Get number of samples by location, year
as.t <-aggregate(loc~Year,sum, data=as)
# Merge the by location, year
fs.tp <- merge(as,as.t,by =c('Year'))
# Calculated Age Sex proporion by year 
fs.tp$p <- with(fs.tp,loc.x/loc.y)
fs.tp$nf <- with(fs.tp,ifelse(Sex=='male'& Age==4,p/nll$par[1],
  ifelse(Sex=='male'& Age==5,p/nll$par[2],ifelse(Sex=='male'& Age==6,p/nll$par[3],
  ifelse(Sex=='male'& Age==7,p/nll$par[4],ifelse(Sex=='female'& Age==4,p/nll$par[5],
  ifelse(Sex=='female'& Age==5,p/nll$par[6],ifelse(Sex=='female'& Age==6,p/nll$par[7],p))))))))
# merge data
	sums <-aggregate(nf~Year,sum,data=fs.tp)
	fs.tp <- merge(fs.tp,sums,by=c('Year'))
	fs.tp$p.b <- with(fs.tp,nf.x/nf.y)
    fs.tp <- fs.tp[,c('Year','Sex','Age','p','p.b')]
age.comp.b <-aggregate(cbind(p.b,p) ~ Year+Age, sum, data=fs.tp)
age.comp.b.w <- reshape(age.comp.b,idvar = c('Year'),timevar = 'Age',direction='wide')
sex.comp.b <-aggregate(cbind(p.b,p) ~ Year+Sex, sum, data=fs.tp)
sex.comp.b.w <- reshape(sex.comp.b,idvar = c('Year'),timevar = 'Sex',direction='wide')

# Get original data 
as.o <- read.csv(data_file3,na='')
#as.o2 <- as.o 
#as.o2$Age.4 <- as.o2$Age.3+as.o2$Age.4
#as.o2$Age.7 <- as.o2$Age.7+as.o2$Age.8
#as.o2 <- as.o2[,c(1,3:6,8)]
as.o$male <- 1-as.o$Female

windows(width=10,height=10, record = TRUE)
m <- c(1,1,1,2,2,2,3,3,3,4)
layout(m)
par(mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
dat1 <- age.comp[,c('Year','Age','p.f')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Age',direction='wide')
dat1.w[is.na(dat1.w)]<- 0
names(dat1.w) <- c('Year', 'Age.3','Age.4','Age.5','Age.6','Age.7','Age.8')
barplot(as.matrix(t(rbind(dat1.w[order(dat1.w$Year)[1:22],2:7],as.o[23:34,2:7]))),main='Raw Data',col=rainbow(6),ylim=c(0,1),names.arg=as.o$Year)
barplot(as.matrix(t(as.o[ ,2:7])),main='Original',ylim=c(0,1),col=rainbow(6),names.arg=as.o$Year)
dat1 <- age.comp[,c('Year','Age','p.s')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Age',direction='wide')
dat1.w[is.na(dat1.w)]<- 0
names(dat1.w) <- c('Year', 'Age.3','Age.4','Age.5','Age.6','Age.7','Age.8')
barplot(as.matrix(t(rbind(dat1.w[order(dat1.w$Year),2:7],as.o[25:34,2:7]))),main='Length Selectivity',col=rainbow(6),ylim=c(0,1),names.arg=as.o$Year)
plot.new()
legend('top',legend=c(3,4,5,6,7,8),fill=rainbow(6), bg='white', horiz=TRUE)


m <- c(1,1,1,2,2,2,3,3,3,4)
layout(m)
par(mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
dat1 <- sex.comp[,c('Year','Sex','p.f')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Sex',direction='wide')
names(dat1.w) <- c('Year','Female','male')
barplot(as.matrix(t(rbind(dat1.w[1:22,2:3],as.o[23:34,8:9]))),col=rainbow(2),main='Raw Data',ylim=c(0,1),names.arg=as.o$Year)
barplot(as.matrix(t(as.o[,8:9])),main='Original',ylim=c(0,1),col=rainbow(2),names.arg=as.o$Year)
dat1 <- sex.comp[,c('Year','Sex','p.s')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Sex',direction='wide')
names(dat1.w) <- c('Year','Female','male')
barplot(as.matrix(t(rbind(dat1.w[,2:3],as.o[25:34,8:9]))),col=rainbow(2),main='Length Selectivity',ylim=c(0,1),names.arg=as.o$Year)
plot.new()
legend('top',legend=c('female','male'),fill=rainbow(2), bg='white', horiz=TRUE)


windows(width=10,height=10, record = TRUE)
m <- c(1,1,1,2,2,2,3,3,3,4)
layout(m)
par(mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
dat1 <- age.comp.b[,c('Year','Age','p')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Age',direction='wide')
dat1.w[is.na(dat1.w)]<- 0
names(dat1.w) <- c('Year','Age.4','Age.5','Age.6','Age.7')
barplot(as.matrix(t(dat1.w[order(dat1.w$Year),2:5])),main='Raw Data',col=rainbow(6),ylim=c(0,1),names.arg=as.o$Year[1:24])
dat1 <- age.comp.b[,c('Year','Age','p.b')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Age',direction='wide')
dat1.w[is.na(dat1.w)]<- 0
names(dat1.w) <- c('Year', 'Age.4','Age.5','Age.6','Age.7')
barplot(as.matrix(t(dat1.w[order(dat1.w$Year),2:5])),main='Age-sex Selectivity',col=rainbow(6),ylim=c(0,1),names.arg=as.o$Year[1:24])
plot.new()
legend('top',legend=c(4,5,6,7),fill=rainbow(6), bg='white', horiz=TRUE)

windows(width=10,height=10, record = TRUE)
m <- c(1,1,1,2,2,2,3,3,3,4)
layout(m)
par(mar=c(2,2,2,2),oma = c(3,4,3,3),cex=0.6)
dat1 <- sex.comp.b[,c('Year','Sex','p')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Sex',direction='wide')
dat1.w[is.na(dat1.w)]<- 0
names(dat1.w) <- c('Year','Female','male')
barplot(as.matrix(t(dat1.w[order(dat1.w$Year),2:3])),main='Raw Data',col=rainbow(2),ylim=c(0,1),names.arg=as.o$Year[1:24])
dat1 <- sex.comp.b[,c('Year','Sex','p.b')]
dat1.w <- reshape(dat1,idvar = c('Year'),timevar = 'Sex',direction='wide')
dat1.w[is.na(dat1.w)]<- 0
names(dat1.w) <- c('Year','Female','male')
barplot(as.matrix(t(dat1.w[order(dat1.w$Year),2:3])),main='Age-sex Selectivity',col=rainbow(2),ylim=c(0,1),names.arg=as.o$Year[1:24])
plot.new()
legend('top',legend=c('female','male'),fill=rainbow(2), bg='white', horiz=TRUE)

