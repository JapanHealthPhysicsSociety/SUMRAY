#' All solid cancer risk mortality models from Life Span Study Report 14.
#'
#' All solid cancer risk mortality models derived from LSS Report 14 1950-2003
#'  (Ozasa et al., Radiat Res, 2012) for use in LAR function
#'
#'
#' @format A list object with 2 sub list objects for risk model information
#'          for excess relative risk (err) and excess absolute risk (ear) models.
#' \describe{
#'   \item{$err}{
#'       \describe{
#'               \item{$para}{a vector of parameter values}
#'               \item{$var}{ a variance covariance matrix}
#'               \item{$f}{a function to compute the excess relative risk}
#'        }
#'      }
#'   \item{$ear}{
#'       \describe{
#'               \item{$para}{a vector of parameter values}
#'               \item{$var}{ a variance covariance matrix}
#'               \item{$f}{a function to compute the excess absolute risk}
#'               }
#'        }
#' }
#'
LSS_allsolid_mortality <- list(
  err=list(
    para=c( 0.422650, -0.345890, -0.857428,  0.344079 ),
    var=matrix( c(
      0.00253235,  0.00140692,  0.00692967, -0.00028615,
      0.00140692,  0.00661873, -0.01588490, -0.00101646,
      0.00692967, -0.01588490,  0.17912000,  0.00354537,
      -0.00028615, -0.00101646,  0.00354537,  0.00771445), ncol=4, byrow=TRUE),
    f=function (beta, dose, agex, age, sex) {
      beta[1] * dose * exp(beta[2] * (agex - 30)/10 + beta[3] *
                             log(age/70)) *
        (1 + c(-1, 1)[sex] * beta[4])/(1 + exp(-(age - agex - 7.5)))
    }
  ),
  ear=list(
    para=c( 0.00263965, -0.21346100,  3.38439000,  0.06769100 ),
    var=matrix( c(
      9.87492e-08,  6.41078e-06,  4.17608e-05, -1.25134e-05,
      6.41078e-06,  5.12184e-03, -1.22614e-02, -2.25000e-05,
      4.17608e-05, -1.22614e-02,  1.33872e-01, -5.84032e-03,
      -1.25134e-05, -2.25000e-05, -5.84032e-03,  9.31031e-03), ncol=4, byrow=TRUE),
    f=function (beta, dose, agex, age, sex) {
      beta[1] * dose * exp(beta[2] * (agex - 30)/10 + beta[3] *
                             log(age/70)) *
        (1 + c(-1, 1)[sex] * beta[4])/(1 + exp(-(age - agex - 7.5)))
    }
  ) )


