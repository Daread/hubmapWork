
# Commands to set up R 4.1.2 to use 
# module unload R/4.0.0
# module unload pcre2/10.35
# module load pcre2/10.39
# module load R/4.1.2

# Working off of a vignette at http://rstudio-pubs-static.s3.amazonaws.com/371411_b572e753d5e94c959cee7f05bc969839.html
#   for fitting a beta binomial mixed model using glmmTMB

# library("devtools")
# devtools::install_github("jhmaindonald/qra")

library("qra")
library(glmmTMB)

# Getting and formatting data


HawCon <- as.data.frame(qra::HawCon)
# HawCon <- within(HawCon, {
#   trtGp <- paste0(CN,LifestageTrt, sep=":")
#   trtGp <- gsub("Fly","",trtGp)
#   trtGpRep <- paste0(CN,LifestageTrt,RepNumber)
#   scTime <- scale(TrtTime)
# })


# # Set up model
# glmmFit <- glmmTMB(cbind(Dead,Live)~0+trtGp/TrtTime+(scTime|trtGpRep),
#                       family=list(family="betabinomial",link="logit"),
#             data=HawCon)


HawCon <- as.data.frame(qra::HawCon)

HawCon <- within(as.data.frame(qra::HawCon),{
	trtGp <- paste0(CN,LifestageTrt)
	trtGp <- gsub("Fly","",trtGp)
	trtGp <- factor(trtGp,levels=sort(unique(trtGp))[c(1,5,2,6,3,7,4,8)])
	gp <- paste0(CN,LifestageTrt,":",RepNumber)
	scTime <- scale(TrtTime)
	##Modelmayfitmorereadilywithacenteredandscaledvariable
})


# From vignetta at https://maths-people.anu.edu.au/~johnm/r-book/4edn/ch7-BetaBinomial.pdf
# bbFit <- glmmTMB ( cbind ( Dead , Live ) ~0+ trtGp / TrtTime +(TrtTime | gp ) ,
# 	dispformula=~ trtGp+poly ( TrtTime ,2) ,
# 	family=betabinomial ( link="logit"), data=HawCon )



bbFit <- glmmTMB ( cbind ( Dead , Live ) ~ trtGp + TrtTime +(1 | gp ) ,
	#dispformula=~ trtGp+poly ( TrtTime ,2) ,
	family=betabinomial ( link="logit"), data=HawCon )


bbFitDisp <- glmmTMB ( cbind ( Dead , Live ) ~ trtGp + TrtTime +(1 | gp ) ,
	dispformula=~ trtGp + TrtTime  ,
	family=betabinomial ( link="logit"), data=HawCon )







