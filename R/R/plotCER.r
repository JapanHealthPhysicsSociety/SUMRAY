#'Plot the age-specific risk (attributable probability rate) due to radiation exposure
#'
#'@param cer a CER class object (returned from function CER)
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
#'    cer <- CER( agex=10, doseGy=0.1, sex=1, maxage=90 )
#'    plot(cer)
#'

#'#'@export
plotCER <- function( cer, xlab="Attained age", ylab="Attributable probability rate", main="Attributable probability rate", col=1 )
{
  plot( range(cer$CERbyAGE$age), range(cer$CERbyAGE[,3:4]), type="n",
        xlab=xlab, ylab=ylab, main=main )
  lines(cer$CERbyAGE$age, cer$CERbyAGE$mean, lwd=2, col=col)
  lines(cer$CERbyAGE$age, cer$CERbyAGE[,3], lty=2, lwd=1, col=col)
  lines(cer$CERbyAGE$age, cer$CERbyAGE[,4], lty=2, lwd=1, col=col)
  invisible(NULL)
}
