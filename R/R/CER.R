#'Calculate the cumulated excess risk due to radiation exposure
#'
#'@param agex age(s) at exposure (a single value or a vector)
#'@param doseGy dose in Gy or Sv (a single value or a vector)
#'@param sex sex 1:male 2:female
#'@param maxage maximum age to follow up (an integer value)
#'@param riskmodel risk model (a list object, which contains two list objects for excess relative risk model (err) and excess absolute risk model (ear), each of which contains a vector of parameter values (para), a matrix of variance covariance matrix (var), and a function to compute the risk given a parameter vector, a dose value, an age at exposure, an attained age and sex.
#'@param wgt weights for ERR vs EAR transfer (e.g. c(1,0) (default), which indicates ERR transfer)
#'@param baseline age- and sex-specific baseline rate (data.frame with columns for age, male, female); default: all solid cancer mortality rate of Japan in 2018
#'@param mortality age- and sex-specific all cause mortality rate (data.frame with columns for age, male, female); default: all cause mortality rate of Japan in 2018
#'@param alpha significance level for the confidence interval (default value = 0.05, which corresponds to 95\% confidence interval)
#'@param n_mcsamp number of Monte Carlo sample size (default:10000)
#'@param seed random number seed (a single value; default is no seed)
#'
#'@return risk information(vector)
#'
#'@examples
#'    # The following examples use default data provided in the JHPSrisk package
#'    # for riskmodel (LSS R14 all solid cancer model),
#'    #     baseline (all solid cancer mortality rates in Japan 2018) and
#'    #     mortality  (all cause mortality rates in Japan 2018)
#'
#'    # Cumulated excess risk for male exposed to 0.1 Gy at age 10
#'    #  followed up to age 90 with err transfer
#'    CER( agex=10, doseGy=0.1, sex=1, maxage=90 )
#'
#'    # Cumulated excess risk for female exposed to 0.01 Gy/year at ages 10-19
#'    #  followed up to age 100 with 7:3 weights for ERR and EAR transfers
#'    CER( agex=10:19, doseGy=rep(0.01,10), sex=2, maxage=100, wgt=c(.7,.3) )
#'
#'@importFrom MASS mvrnorm
#'#'@export
CER <- function( agex, doseGy, sex, riskmodel=LSS_allsolid_mortality, wgt=c(1,0),
                 baseline=allsolid_mortality_Japan2018, mortality=mortality_Japan2018,
                 maxage=max(mortality$age),
                 alpha=0.05, n_mcsamp=10000, seed=NULL )
{
  if( length(agex)!=length(doseGy)) stop(message="lengths of agex and doseGy must be equal.")
  if( length(wgt)!=2 | sum(wgt)!=1 )stop(message="Invalid weights (wgt) for ERR and EAR transfers.")
  if( sex!=1 & sex!=2 ) stop(message="sex must be specified by 1 (male) or 2 (female).")
  if( min(agex) < min(baseline$age) |  max(agex) > max(baseline$age) ) stop(message="agex must be in baseline$age.")
  if( maxage < min(baseline$age) |  maxage > max(baseline$age) ) stop(message="maxage must be in baseline$age.")
  if( maxage < max(agex) ) stop(message="maxage must be equal to or greater than max(agex).")
  if( sum( names(riskmodel) != c("err","ear") ) ) stop(message="riskmodel must have list objects; err and ear.")
  if( sum( names(riskmodel$err) != c("para","var","f") ) ) stop(message="riskmodel$err must have objects; para, var, f.")
  if( sum( names(riskmodel$ear) != c("para","var","f") ) ) stop(message="riskmodel$ear must have objects; para, var, f.")

  if( length(unique(diff(mortality$age))) != 1 | ( diff(mortality$age)[1] != 5 & diff(mortality$age)[1] != 1 ) ) stop(message="mortality$age must be equally spaced by 5 or 1.")
  if( length(unique(diff(baseline$age))) != 1 | ( diff(baseline$age)[1] != 5 & diff(baseline$age)[1] != 1 ) ) stop(message="baseline$age must be equally spaced by 5 or 1.")

  interpolate_rates <- function( rates_by5 ){ # interpolate rates in 5-year intervals to those in 1-year intervals
    interpolate0 <- function( mr0 ){
      ( mr0 <- c(0,mr0) )
      ( res <- (mr0[2]-mr0[1])/2.5 * (0:2+.5) )
      for( i in 2:(length(mr0)-1) )
        res <- c( res, mr0[i] + (mr0[i+1]-mr0[i])/5 * (1:5) )
      res
    }
    mr <- interpolate0( rates_by5$male )
    fr <- interpolate0( rates_by5$female )
    data.frame( age=1:length(mr), male=mr, female=fr )[1:pmin(100,length(mr)),]
  }
  if( diff(mortality$age)[1] == 5 ) mortality <- interpolate_rates( mortality )
  if( diff(baseline$age)[1] == 5 ) baseline <- interpolate_rates( baseline )


  set.seed(seed)
  ages <- baseline$age
  nexp <- length(agex)
  sexlab <- c("male","female")[sex]
  mrate <- mortality[[sexlab]]
  brate <- baseline[[sexlab]]

  survp <- c( 1, exp(-cumsum(mrate))[-length(mrate)] )

  mc_paras_err <- MASS::mvrnorm(n=n_mcsamp, mu=riskmodel$err$para, Sigma=riskmodel$err$var)
  mc_paras_ear <- MASS::mvrnorm(n=n_mcsamp, mu=riskmodel$ear$para, Sigma=riskmodel$ear$var)

  sim_err <- apply( mc_paras_err, 1, function(bet){
    a <- sapply( 1:nexp, function(i)
      riskmodel$err$f( bet, dose=doseGy[i], agex=agex[i], age=ages, sex=sex ) )
    apply(a, 1, sum) } )

  sim_ear <- apply( mc_paras_ear, 1, function(bet){
    a <- sapply( 1:nexp, function(i)
      riskmodel$ear$f( bet, dose=doseGy[i], agex=agex[i], age=ages, sex=sex ) )
    apply(a, 1, sum) } )

  ok <- min(agex) < ages  &   ages <= maxage

  CERs_err <- apply( sim_err, 2, function(err){ ( brate*err*survp )[ok]  / survp[min(agex)+1] } )
  CERs_ear <- apply( sim_ear, 2, function(ear){ ( ear*survp )[ok] / survp[min(agex)+1] } )


  res_byage <- wgt[1] * CERs_err + wgt[2] * CERs_ear
  res_total <- apply( res_byage, 2, sum )
  ci_total  <- quantile( res_total, prob=c(alpha/2, 1-alpha/2)  )
  ci_byage <- apply( res_byage, 1, quantile, prob=c(alpha/2, 1-alpha/2) )

  CERbyAGE <- data.frame( age=ages[ok], mean=apply( res_byage, 1, mean ), CI=t(ci_byage) )
  CER <- c( mean=mean(res_total), median=median(res_total), CI=ci_total )
  print(CER)
  #  invisible( structure(c(CER=CER, CERbyAGE=CERbyAGE), class="CER") )
  invisible( structure(list(CER=CER, CERbyAGE=CERbyAGE), class="CER") )
}
