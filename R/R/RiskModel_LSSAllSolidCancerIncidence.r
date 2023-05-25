#' All solid cancer incidence risk models from Life Span Study
#'
#' All solid cancer incidence risk models derived from LSS 1958-2009
#'  (Grant et al., Radiat Res, 2017) for use in LAR function
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
LSS_allsolid_incidence <- list(
  err=list(para = c(d10col = 0.452246303147881, e30 = -0.189230676147737,
                    lage70 = -1.74186812049818, msex = 0.250755761008251),
           var = structure(c(0.00195131904979306,0.000413004853015443, 0.00671838485961983, -0.000156670327057745,
                             0.000413004853015443, 0.00444394708519854, -0.0127818040921709, -0.000421103579208728,
                             0.00671838485961983, -0.0127818040921709, 0.103802938359842, 0.00284023304347374,
                             -0.000156670327057745, -0.000421103579208728, 0.00284023304347374, 0.00491548662758437
           ), .Dim = c(4L, 4L), .Dimnames = list(c("d10col", "e30", "lage70","msex"),
                                                 c("d10col", "e30", "lage70", "msex"))),
           f = function (beta,dose, agex, age, sex) {
             beta[1] * dose * exp(beta[2] * (agex - 30)/10 + beta[3] *log(age/70)) *
               (1 + c(-1, 1)[sex] * beta[4])/(1 + exp(-(age -agex - 7.5)))
           } ),
  ear=list(
    para=c(d10col = 0.00505108100884472, e30 = -0.287599669064682,
           lage70 = 2.37403004629745, msex = 0.171952011757523),
    var=structure(c(2.44198984049119e-07, 5.36502503488343e-06, 6.21799394209355e-05, -1.29885990096665e-05,
                    5.36502503488343e-06, 0.00400540779279526, -0.0105161313375575, -0.000145731631208299,
                    6.21799394209355e-05, -0.0105161313375575, 0.0747308867068616, -0.00275992238361282,
                    -1.29885990096665e-05, -0.000145731631208299, -0.00275992238361282, 0.00505962128686142),
                  .Dim = c(4L, 4L), .Dimnames = list(c("d10col", "e30", "lage70", "msex"),
                                                     c("d10col", "e30", "lage70", "msex"))),
    f=function (beta, dose, agex, age, sex) {
      beta[1] * dose * exp(beta[2] * (agex - 30)/10 + beta[3] *
                             log(age/70)) *
        (1 + c(-1, 1)[sex] * beta[4])/(1 + exp(-(age - agex - 7.5)))
    }
  )
)
