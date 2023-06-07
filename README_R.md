![JHPS_elogo](https://github.com/JapanHealthPhysicsSociety/SUMRAY/assets/100466085/da86d36f-5f97-4b75-bf70-f92a17e744f0)

# **The R library provided in the SUMRAY code**
This is the R library to compute the attributable probability rates and the Cumulative Excess Risks (CERs) of all solid cancers from radiation exposures.

The library were written by the Research Group on the development of cancer risk estimation code associated with radiation exposure (FY2020 - FY2021) of Japan Health Physics Society (JHPS).

*Version 1.0. Released June 30, 2023.*

## Preparation
This library is installed and used in the R-language execution environment. To install the SUMRAY package of R from the Github repository, execute the following set of commands in the R environment.
  
* Installing the *remotes* package to download and install R packages from remote repositories.

```sh
install.packages("remotes")
```

* Installing the *SUMRAY* package stored in the Github repository.  

```sh
remotes::install_github("JapanHealthPhysicsSociety/SUMRAY", subdir="/R")
```

* Loading the *SUMRAY* package in your R environment.

```sh
library(SUMRAY)
```

You are now ready to use the library.
  
## Usage
For instructions on how to use this library, please refer to the [PDF file](https://github.com/JapanHealthPhysicsSociety/SUMRAY/files/11660915/R-Usage.pdf) or the supplement of the paper.

## Example
As an example calculation using this SUMRAY library, the cumulative excess risks (CERs) at attained age 90 are calculated for a woman chronically exposed to 0.01 Gy annually from age 18 to 65 with the risk transfers of 0.5 between EAR and ERR models. Here, all items in *Preparation* have been completed.

```sh
res <- CER( agex=18:65, doseGy=rep(0.01,48), sex=2, maxage=90, wgt=c(.5,.5))
```  

![Screenshot1](https://github.com/JapanHealthPhysicsSociety/SUMRAY/assets/100466085/db423b55-5f00-4cee-b8ee-316bb91da513)

In the previous command, the values of the CERs and the attributable probability rates are stored in the variable *res*, which can be referenced in the following command.

```sh
res
```  

![Screenshot2](https://github.com/JapanHealthPhysicsSociety/SUMRAY/assets/100466085/708312cc-44b4-4ba5-ab1b-672ee1432de1)

In additon, the following command draws the attributable probability rates stored in the variable *res*.

```sh
plotCER(res)
```  

![Rplot](https://github.com/JapanHealthPhysicsSociety/SUMRAY/assets/100466085/08020d1b-ff47-4520-bc35-67a8d9b85d1b)

## Citation
M Sasaki, K Furukawa, D Satoh, K Shimada, S Kudo, S Takagi, S Takahara, M Kai;  
SUMRAY: R and Python codes for calculating cancer risk due to radiation exposure of a population,  
Journal of Radiation Protection and Research, 2023.  
https://***

## License
SUMRAY has been released under the [MIT license](https://github.com/JapanHealthPhysicsSociety/SUMRAY/blob/main/LICENSE.md)

## Copyright
2020-2022 Japan Health Physics Society.  
E-mail: exec.off@jhps.or.jp
