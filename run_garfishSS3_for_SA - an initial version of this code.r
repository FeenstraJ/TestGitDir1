
# rm(list=ls(all=TRUE))   # clears all objects from the primary environment
# Get latest (master) version of r4ss from Github:

# devtools::install_github("r4ss/r4ss")

# Set up for analyses --------------------------------------------------------

library(rforSS3)
library(r4ss)
library(codeutils)
library(hplot)

ddir <- getDBdir()
source(pathtopath(ddir,"projects/SA_SS3_Garfish/sagarfish/R/ss3_utilities.R"))
#library(r4maps)

options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)

# wkdir is the directory within which the directory structure is setup.
# save a copy of this R file in the working directory to secure any customizations
# that are made
wkdir <- pathtopath(getDBdir(),"projects/SA_SS3_Garfish/sagarfish/")

#source("C:/A_CSIRO/Rcode/SESSF/ss3/oro2017/oro2017utils.R")
# The directory structure inside the wkdir will include 'store' and 'calc'
# 'store' contains a directory for each step in the move from the previous
# balanced basecase to the new balanced basecase.
# 'calc' is the directory in which all the SS3 files are put ready for analysis.
# SS3.exe needs to be put in 'calc'.
store <- pathtopath(wkdir,"basecase/")  # define directories
calc <- pathtopath(wkdir,"calc/")
profile <- pathtopath(wkdir,"profile/")
# generate these two directories
dirExists(store)
dirExists(calc)
dirExists(profile)
# what sub-directories will be used within 'basecase'?
# Change these names to suit those that will actually be used.
# The order represents the order in which data will be added to the earliest of
# simplest model leading to the final balanced basecase
basecase <- c("SGbasecase",
              "SGbasecase_domed",
              "SGbasecase_h7",
              "SGbasecase_M4",
              "SGbasecase_M45",
              "SGbasecase_M47",
              "GSVbasecase",
              "GSVbasecase6")         ## 8 ages have small delta         
numdirs <- length(basecase)
# now safely generate the directories
for (direct in 1:numdirs) dirExists(paste0(store,basecase[direct]))
printV(basecase)
{
print("Now ensure that you have X.ctl, X.dat, X.par, X.for, and X.sta file",quote=FALSE)
print("in each directory; the .for will become forcast.ss and the .sta file",quote=FALSE)
print("will become the starter.ss file",quote=FALSE)
}


# sscopyto(origin=store,fromdir="SGbasecase",todir="SGbasecase_M47")


# now ensure that you have X.ctl, X.dat, X.par, X.for, and X.sta file in each
# directory where X is the name of each directory. The .ctl, .dat, and .par
# files relate to files of the same name but the .for file will become the
# forcast.ss file and the .sta file will become the starter.ss file. 
# The idea is to conduct a bridging analysis from an old assessment
# to a new assessment that includes all the new data and is balanced.
#
# loop through the directories of 'basecase' if you have prepared the five
# required files in each directory, that is the  X.ctl, X.dat, X.par, X.for,
# and X.sta file. That would entail including extra data in the dat file and
# adjusting the .ctl file accordingly.
#
# remember 'basecase' is the directory in which these sub-directories are kept.
#L -----
# Fit a model ------
#numiter <- length(basecase)
item <- 1
getCase(index=item,basecase)   # this lists the basecase indices to the screen

#executable <- c("SS","SS","SS","SS","SS","SS","SS3","SS3")
starttime <- Sys.time()
    analysis <- getCase(index=item,basecase)  # 
    cat("\n\n")
    print("New Analysis")
    print(analysis)
    destination <- pathtopath(store,analysis)
    print(destination)
    
    copyfiles(analysis,store,calc) # copy and rename the needed files into calc
    # now call SS3 -nohess twice, second time using the pars from the
    # first, then the third time to calculate the hessian
     # run(dir=calc,exe="ss3",extras="-stopph 0 -nohess",show_in_console = TRUE,
     #       skipfinished = FALSE)


    run(dir=calc,exe="ss3",extras="-nohess -maxfn 700",show_in_console = TRUE,
        skipfinished = FALSE)   

    fixstarter(calc,findtext="init_values_src")
 # estimate the inverse Hessian
    run(dir=calc,exe="ss3",extras="-maxfn 700",show_in_console = TRUE,
        skipfinished = FALSE)

    # Now return the required files back into the source directory defined
    # in store and basecase
    storeresults(calc,destination)
    # now print and plot the results using r4ss
    # first go to the directory where the results have been stored
    print(store)
    setwd(store)
    print(paste0("Running from ",destination))
    fileout <- pathtopath(destination,paste0(analysis,".txt"))
    sink(fileout)
    Btarget <- 0.50  # 
    Blimit  <- 0.20  # default SA Limit Reference Point

    plotreport  <- SS_output(dir=destination, 
                             repfile = "Report.sso",
                             compfile = "CompReport.sso",
                             covarfile = "covar.sso",
                             forefile = "Forecast-report.sso",
                             wtfile = "wtatage.ss_new",
                             warnfile = "warning.sso",
                             forecast=TRUE,
                             covar=TRUE,
                             readwt=FALSE,
                             verbose=TRUE)

    finalplot <- SS_plots(replist=plotreport, plot=1:26,
                          btarg=Btarget,minbthresh=Blimit,
                          uncertainty=TRUE, datplot=TRUE,
                          forecastplot=TRUE, 
                          png=TRUE,html=FALSE)
    SS_html(replist=plotreport,plotdir=pathtopath(destination,"plots"),
            title=analysis)
    
    
    sink()   # close off fileout containing screen dump from SS_output and SS_plots
    
   filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
    save(plotreport,file=filename)
    cat("\n\nplotreport saved to ",filename)
    # load(filename)
    cat("SS_output and SS_plots txt sent to ",fileout,"\n")
    # Tidy up; return to wkdir
    setwd(wkdir)
 

endtime <- Sys.time()
print(endtime - starttime)

round(printV(summarizeSS3(plotreport)$answer),6)

round(getStatus(plotreport),7)

# End model fit --------------------------------------

# Balance Tuning variances --------------------------------------

{
fltnames <- plotreport$FleetNames  
pickfleets <- c(2)
nfleet <- length(pickfleets)
tuneinfo <- r4ss:::get_tuning_table(plotreport,fleets=pickfleets,
                                    option="Francis")[,c(1,2,3,5,6,7)]
i <- 1  # for loop required if more than one fleet with comp data
cat("   ",tuneinfo[1,1],tuneinfo[1,2],tuneinfo[1,3],
    "  # Variance adjustment_",fltnames[pickfleets[i]],"_Length  \n")
cat("   ",tuneinfo[2,1],tuneinfo[2,2],tuneinfo[2,3],
    "  # Variance adjustment_",fltnames[pickfleets[i]],"_Age \n")
}




# Adjust Recruit bias Ramp-----------------------

biasramp <- SS_fitbiasramp(plotreport,
                           verbose = FALSE,
                           startvalues = NULL,
                           method = "BFGS",
                           twoplots = TRUE,
                           transform = FALSE,
                           plot = TRUE,
                           print = FALSE, # print to png files?
                           plotdir = "default",
                           shownew = TRUE,
                           oldctl = NULL,
                           newctl = NULL,
                           altmethod = "nlminb",
                           exclude_forecast = FALSE,
                           pwidth = 6.5,pheight = 5,punits = "in",
                           ptsize = 10,res = 300,cex.main = 1)
newpars <- biasramp$df
{
for (i in 1:4) {
  cat("  ",newpars[i,"value"],newpars[i,"label"],"\n")
}
cat("  ",newpars[5,"value"],
    "#_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0) \n")
}


# Profile -----------------------------------
item <- 1
getCase(index=item,basecase)   # this lists the basecase indices to the screen

analysis <- getCase(index=item,basecase)  # 
print(analysis)
destination <- pathtopath(store,analysis)
print(destination)

copyfiles(analysis,store,calc) 

run(dir=calc,exe="ss3",extras="-nohess -maxfn 700",show_in_console = TRUE,
    skipfinished = FALSE)   
# all that should ensure that the control.ss.new is correct
profilestarter(calc=calc,findtext=c("init_values_src","prior_like"),
                           verbose=TRUE) 

profvec <- c(0.41,0.43,0.45,0.47,0.49,0.51)  # NatM
profvec <- c(8.84,8.845,8.85,8.855,8.86,8.865,8.87,8.875,8.88,8.885)       # LnR0

profans <- r4ss::profile(dir=calc,oldctlfile="control.ss_new",
                         newctlfile="ss.ctl",
                         linenum=123,profilevec=profvec,
                         usepar=TRUE,parlinenum=49,exe="ss3",
                         extras="-nohess -maxfn 600",
                         show_in_console = TRUE,
                         prior_check=TRUE)



if (any(profans[,"converged[whichruns]"]) == FALSE) {
  cat("Some parameter values failed to converge \n")
} else {
  cat("all trials succeeded inconverging  \n")
}
 
invar <- "LnR0"
minlike <- 0.01
keepcols <- which(abs(colSums(profans)) > minlike)
limprofans <- profans[,keepcols]
colans <- ncol(limprofans)
numcol <- colans - 1
columns <- colnames(limprofans)[2:colans]
numrow <- length(profvec)
likes <- matrix(0,nrow=numrow,ncol=numcol,dimnames=list(profvec,columns))
mins <- apply(limprofans[,2:colans],2,min,na.rm=TRUE) 
for (i in 1:numcol) {
  likes[,i] <- limprofans[,(i+1)]-mins[i]
}
round(limprofans[,,5])
round(likes,5)

pickfinal <- which(apply(likes,2,max,na.rm=TRUE) > minlike)
limlikes <- likes[,pickfinal]
round(limlikes,5)
numfinal <- length(pickfinal)
label <- colnames(limlikes)

maxy <- getmax(limlikes)
plotprep(width=10, height=6)
parset()
plot(profvec,limlikes[,"TOTAL"],type="l",lwd=1,col=1,xlab=invar,
     ylab="Likelihood Difference",ylim=c(0,maxy),panel.first = grid())
for (i in 2:numfinal) lines(profvec,limlikes[,i],col=i,lwd=3)
abline(h=c(0,1.92),lwd=1,col=1,lty=c(1,3))
legend("topright",legend=label,lwd=4,col=1:numfinal,bty="n",cex=1.5)

# [3:37 PM] Fay Helidoniotis
# https://connect.fisheries.noaa.gov/ss3-helper/SS3-Helper


# Tabulate results---------------------------------------------
# generates a Tables directory in destination
analysis <- "SGbasecase"
analysis <- "GSVbasecase6"
analysis <- "SGbasecase_h7"
destination <- pathtopath(store,analysis)

filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)

sort(names(plotreport))

SSexecutivesummary(replist=plotreport,
  plotfolder = "default",
  ci_value = 0.95,
  es_only = FALSE,
  fleetnames = NULL,
  add_text = "SG",
  so_units = "millions of eggs",
  divide_by_2 = FALSE,
  endyr = NULL,
  adopted_ofl = NULL,
  adopted_abc = NULL,
  adopted_acl = NULL,
  forecast_ofl = NULL,
  forecast_abc = NULL,
  verbose = TRUE
)
# end tabulate results


# Get All PARAMETERS---------------------------------------

library(rforSS3)
library(r4ss)
library(codeutils)
library(hplot)
wkdir <- pathtopath(getDBdir(),"projects/SA_SS3_Garfish/sagarfish/")
store <- pathtopath(wkdir,"basecase/")  # define directories
analysis <-  "SGbasecase"      # chaneg to GSVbasecase to get GSV results
tmpdest <- paste0(store,analysis,"/")
filename <- pathtopath(tmpdest,paste0("plotreport_",analysis,".Rdata"))
load(filename)


pars <- plotreport$parameters
pickP <- which((pars[,"Phase"] > 0) & (pars[,"Value"] != 0))
fittedpars <- pars[pickP,]
usecols <- c(3,5,12)
fittedpars[,usecols]

pickrec <- grep("RecrDev",fittedpars[,"Label"])
nrec <- length(pickrec)

mainpars <- fittedpars[-pickrec,usecols]  # 

# Estimated parameters
knitr::kable(mainpars[,usecols],row.names=TRUE,digits=c(5,0,12))

pars <- plotreport$parameters
pickP <- which((pars[,"Phase"] < 0) & (pars[,"Value"] != 0))
fixedpars <- pars[pickP,c(3,5)]
# Fixed parameters
knitr::kable(fixedpars,digits=c(6,0))

# Quick Report-----------------------------------

analysis <-  "SGbasecase"
store <- "c:/Users/Malcolm/DropBox/projects/SA_SS3_Garfish/sagarfish/basecase/"
destination <- paste0(store,analysis,"/")
filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)
printV(round(summarizeSS3(plotreport)$answer,6))

# SELECTIVITY--------------------------------------
analysis <-  "SGbasecase"
store <- "c:/Users/Malcolm/DropBox/projects/SA_SS3_Garfish/sagarfish/basecase/"
destination <- paste0(store,analysis,"/")
filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)


plotselex(plotreport,sex="Female",yrs=c(1984,2004,2016),upbound=365)

plotselex(plotreport,sex="Male",yrs=c(1984,2004,2016))

# BIOMASS _at_AGE ----------------------------------------------
batage <- plotreport$batage
printV(colnames(batage))

pickcol <- c(3,8,10,11,12:19)
batage <- batage[,pickcol]
batage$totB <- NA
batage$totB <- rowSums(batage[,7:12])

pickB <- which((batage$Sex == 1) & (batage[,"Beg/Mid"] == "B"))
bfem <- batage[pickB,]
bfem$depl <- NA
bfem$depl <- bfem$totB/bfem$totB[1]
bfem

pickB <- which((batage$Sex == 2) & (batage[,"Beg/Mid"] == "B"))
bmal <- batage[pickB,]
bmal$depl <- NA
bmal$depl <- bmal$totB/bmal$totB[1]
bmal

yrs <- 1984:2027
picky <- match(yrs,bmal[,"Yr"])
maxy <- getmax(bmal[picky,"totB"])
plotprep(width=9,height=4.5)
parset()
plot(yrs,bmal[picky,"totB"],type="l",lwd=2,ylim=c(0,maxy),xlab="",
     ylab="Total Biomass",panel.first=grid())
lines(yrs,bfem[picky,"totB"],lwd=2,col=2)


# NUMBERES-AT-AGE-----------------------------------------------------

natage <- plotreport$natage_annual_2_with_fishery
properties(natage)

yrs <- sort(unique(natage[,"Yr"]))
nyr <- length(yrs)
nagef <- natage[which(natage[,"Sex"] == 1),]
nagem <- natage[which(natage[,"Sex"] == 2),]

sumfreqf <- rowSums(nagef[,5:10],na.rm=TRUE)
sumfreqm <- rowSums(nagem[,5:10],na.rm=TRUE)

sumfreqcntf <- sumfreqcntm <- sumfreqf
for (i in 1:nyr) {
  sumfreqcntf[i] <- sum(nagef[i,5:10] * 1:6) 
  sumfreqcntm[i] <- sum(nagem[i,5:10] * 1:6) 
}
avagef <- sumfreqcntf/sumfreqf
avagem <- sumfreqcntm/sumfreqm
cbind(yrs,avagef,avagem)


# end Numbers-at-age

# Exploitabke Biomass Depletion -------------------
library(rforSS3)
library(r4ss)
library(codeutils)
library(hplot)
ddir <- getDBdir()
source(pathtopath(ddir,"projects/SA_SS3_Garfish/sagarfish/R/ss3_utilities.R"))
options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)
wkdir <- pathtopath(getDBdir(),"projects/SA_SS3_Garfish/sagarfish/")
store <- pathtopath(wkdir,"basecase/")  # define directories
calc <- pathtopath(wkdir,"calc/")
profile <- pathtopath(wkdir,"profile/")
dirExists(store)
dirExists(calc)
dirExists(profile)
basecase <- c("SGbasecase",
              "GSVbasecase",
              "GSVbasecase6")         ## 8 ages have small delta         
numdirs <- length(basecase)
for (direct in 1:numdirs) dirExists(paste0(store,basecase[direct]))

analysis <- "SGbasecase/tables/"
destination <- pathtopath(store,analysis)

refpts <- read.table(file=pathtopath(destination,"e_ReferencePoints_ES.csv"),
                     sep=",",row.names=1)

B0 <- as.numeric(refpts[3,1])

timeser <- read.table(file=pathtopath(destination,"TimeSeries.csv"),
                      sep=",",row.names=1)
dimen <- dim(timeser)
columns <- timeser[1,]
rows <- rownames(timeser)[2:dimen[1]]
yrs <- as.numeric(rows)
times <- matrix(0,nrow=(dimen[1]-1),ncol=dimen[2],
                dimnames=list(rows,columns))
for (i in 1:8) times[,i] <- as.numeric(timeser[2:dimen[1],i])

depleB <- times[,3]/B0
printV(depleB)


# CPUE alternative plots-------------------------------------------

#analysis <- "SGbasecase"
analysis <- "GSVbasecase"
destination <- pathtopath(store,analysis)
filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)
sort(names(plotreport))

cpue <- plotreport$cpue
colnames(cpue) <- tolower(colnames(cpue))
fleets <- unique(cpue[,"fleet_name"])
nfleet <- length(fleets)
plotprep(width=8, height=8)
parset(plots=c(nfleet,2),margin=c(0.3,0.45,0.05,0.1))
for (fl in 1:nfleet) { # fl = 1
  dat <- cpue[cpue[,"fleet_name"]==fleets[fl],]
  rown <- nrow(dat)
  maxy <- getmax(dat[,c("obs","exp")])
  plot(dat[,"yr"],dat[,"exp"],type="l",lwd=2,xlab="",ylab=fleets[fl],
       ylim=c(0,maxy),panel.first=grid())
  points(dat[,"yr"],dat[,"obs"],pch=16,cex=1.0,col=2)
  plot(dat[,"yr"],dat[,"dev"],type="p",pch=16,cex=1.0,xlab="",ylab="Deviate")
  abline(h=0.0,lwd=1.0,col=1)
  for (i in 1:rown) lines(c(dat[i,"yr"],dat[i,"yr"]),c(0.0,dat[i,"dev"]),lwd=1)
}


# Spawning Biomass depletion-----------------------------
analysis <- "SGbasecase"
analysis <- "SGbasecase_dog"
analysis <- "GSVbasecase"
analysis <- "GSVbasecase6"

destination <- pathtopath(store,analysis)
filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)
sort(names(plotreport))
plotreport$current_depletion

sert <- plotreport$timeseries

pickcols <- c(2,3,5,7,8,grep("obs_cat",colnames(sert)),grep("dead\\(B\\)",colnames(sert)))

usedat <- sert[,pickcols]
columns <- c("year","era","expB","matB","recruit","RECcat","HNTcat","HNNTcat",
             "DPNcat","chartcat","RECpred","HNTpred","HNNTpred","DPNpred",
             "chartpred")
colnames(usedat) <- columns
usedat
rown <- nrow(usedat)
deplsB <- usedat[,"matB"]/usedat[1,"matB"]

cbind(usedat[,"year"],deplsB)
yrs <- usedat[2:rown,"year"]

totC <- rowSums(usedat[,c(6:10)],na.rm=TRUE)
predC <- rowSums(usedat[,c(11:15)],na.rm=TRUE)
pickr <- 2:rown

allC <- c(totC[2:41],predC[42:46])
cbind(yrs,allC)

plotprep(width=9,height=7)
parset(plots=c(2,1),cex=1.0)
maxy <- getmax(c(deplsB[pickr], 0.5))
plot(yrs,deplsB[pickr],type="l",lwd=2,ylim=c(0,maxy),xlab="",
     ylab="Spawning Biomass Depletion",panel.first=grid(),yaxs="i")
points(yrs,deplsB[pickr],pch=16,cex=1)
abline(h=c(0.2,0.5),lwd=1,lty=2,col=2)
abline(v=c(2002.5,2014.5,2022.5),lwd=1,lty=2,col=1)
text(x=1992,y=0.025,"LML=210",cex=1.0)
text(x=2009,y=0.025,"LML=230",cex=1.0)
text(x=2019,y=0.025,"LML=250",cex=1.0)
text(x=2023,y=0.025,"Projection",cex=1,pos=4)
maxy <- getmax(totC[pickr])
plot(yrs,allC,type="l",lwd=2,ylim=c(0,maxy),xlab="",ylab="Total Catch (t)",
      panel.first=grid(),yaxs="i")
abline(v=c(2002.5,2014.5,2022.5),lwd=1,lty=2,col=1)

# compare two scenarios
# scenario 1
analysis <- "SGbasecase"
destination <- pathtopath(store,analysis)
filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)
sort(names(plotreport))
plotreport$current_depletion
sert <- plotreport$timeseries
pickcols <- c(2,3,5,7,8,grep("obs_cat",colnames(sert)),
              grep("dead\\(B\\)",colnames(sert)))
usedat <- sert[,pickcols]
columns <- c("year","era","expB","matB","recruit","RECcat","HNTcat","HNNTcat",
             "DPNcat","chartcat","RECpred","HNTpred","HNNTpred","DPNpred",
             "chartpred")
colnames(usedat) <- columns
usedat

rown <- nrow(usedat)
deplsB <- usedat[,"matB"]/usedat[1,"matB"]
#scenario 2
analysis <- "SGbasecase_dog"
destination <- pathtopath(store,analysis)
filename <- pathtopath(destination,paste0("plotreport_",analysis,".Rdata"))
load(filename)
sort(names(plotreport))
plotreport$current_depletion
sert2 <- plotreport$timeseries
pickcols <- c(2,3,5,7,8,grep("obs_cat",colnames(sert2)),
              grep("dead\\(B\\)",colnames(sert2)))
usedat2 <- sert2[,pickcols]
columns <- c("year","era","expB","matB","recruit","RECcat","HNTcat","HNNTcat",
             "DPNcat","chartcat","RECpred","HNTpred","HNNTpred","DPNpred",
             "chartpred")
colnames(usedat2) <- columns

rown <- nrow(usedat2)
deplsB2 <- usedat2[,"matB"]/usedat2[1,"matB"]

pickr <- 2:rown
yrs <- usedat[pickr,"year"]

plotprep(width=9,height=5)
parset(plots=c(1,1),cex=1.0)
maxy <- getmax(c(deplsB[pickr], 0.5))
plot(yrs,deplsB[pickr],type="l",lwd=2,ylim=c(0,maxy),xlab="",
     ylab="Spawning Biomass Depletion",panel.first=grid(),yaxs="i")
points(yrs,deplsB[pickr],pch=16,cex=1)
lines(yrs,deplsB2[pickr],lwd=2,col=2)
points(yrs[40:45],deplsB2[pickr[40:45]],pch=16,cex=1,col=2)
abline(h=c(0.2,0.5),lwd=1,lty=2,col=2)
abline(v=c(2002.5,2014.5,2022.5),lwd=1,lty=2,col=1)
text(x=1992,y=0.025,"LML=210",cex=1.0)
text(x=2009,y=0.025,"LML=230",cex=1.0)
text(x=2019,y=0.025,"LML=250",cex=1.0)
text(x=2023,y=0.025,"Projection",cex=1,pos=4)
legend("top",c("SGbasecase","SGbasecase_dog"),lwd=3,col=c(1,2),bty="n",
       cex=1.25)









# end of spawning biomass plots






# end-of-file--------
