## functions

## takes: matrix of age X years
## and a start year and end year
## and age at start year
## returns: someone's (diagonal) yearly hazard of death


##' Get diagonal hazard
##'
##' TODO
##' @title Projected death risks
##' @param yr0 initial year
##' @param age0 initial age
##' @param yr1 final year
##' @param M matrix to extrapolate
##' @return vector hazard for years 
##' @author Pete Dodd
##' @export
getHzhist <- function(yr0,age0,yr1,M){
  ans <- akima::bilinear(x=magps,y=ymps,z=M,
                  x0=age0 + 0:(yr1-yr0),y0=c(yr0:yr1))
  ans <- ans$z
  end <- -Inf
  if(!all(ans[-1]>0)){
    end <- which(ans[-1]>0)
    if(length(end)>0)
      end <- max(end)+1           #last-one-carried fwd for old
    else
      end <- -Inf
  }
  if(end > -Inf & end<length(ans))      #-Inf is if all >0
    ans[end:length(ans)] <- ans[end]
  ans
}

##' Get matrix from data in right format
##'
##' @title Get input data
##' @param cn iso3 code for country
##' @param sx sex = Female/Male/Total
##' @return matrix ages, years
##' @author Pete Dodd
##' @export
getCNdata <- function(cn,sx){
  ZM <- matrix(ncol=length(ymps),nrow=length(magps),
               data=TT[iso3==cn & Sex==sx,mx])
  ZM
}

##' Discounted life years
##'
##' TODO
##' @title discounted life-years
##' @param iso3 country
##' @param sex Female/Male/Total
##' @param age age now
##' @param yearnow year now
##' @param endyear end year
##' @param HR hazard ratio of death if relevant
##' @param dr discount rate (exponential version)
##' @return discounted life years
##' @author Pete Dodd
discly <- function(iso3,                #country
                   age,                 #initial age
                   yearnow,             #year now
                   sex='Total',         #sex:Female/Male/Total
                   endyear=2098,        #year end 2098 max
                   HR=1,                #hazard ratio for mortality
                   dr=0                 #discount rate (exponential version)
                   ){
  zm <- getCNdata(iso3,sex)
  tmp <- getHzhist(yearnow,age,endyear,zm)
  sum(exp(-cumsum(HR*tmp+dr)))
}
