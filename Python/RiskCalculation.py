#==============================================================================
#    This is a Python-3 script to compute the Cumulative Excess Risks (CERs) 
#    of all solid cancers from radiation exposures.
#------------------------------------------------------------------------------
#    Written by the Research Group on the development of cancer risk estimation 
#    code associated with radiation exposure (FY2020 - FY2021) of Japan Health 
#    Physics Society (JHPS)
#------------------------------------------------------------------------------
#    Version 1.1
#    November 07, 2022
#------------------------------------------------------------------------------
#    Copyright (c) 2020-2022 Japan Health Physics Society
#    Released under the MIT license
#    https://opensource.org/licenses/mit-license.php
#==============================================================================

#------------------------------------------------------------------------------
#    List of the CSV-files for the baseline data located in the data directory
#------------------------------------------------------------------------------

# File name of the all cause mortality data for the population.
file_all_cause_mortality = 'all_cause_mortality(2009-2018)E.csv'

# File name of the cancer-incident baseline data for the population.
file_cancer_incidence = 'cancer_incidenceNCR(2016-2018)E.csv'

# File name of the cancer-mortality baseline data for the population.
file_cancer_mortality = 'cancer_mortality(1958-2019)E.csv'

#------------------------------------------------------------------------------
#    List of the files for the EAR- and ERR-model parameters located in the data directory
#------------------------------------------------------------------------------

# Name of the file containing the EAR-model parameters and covariance matrix for cancer-incidence rates.
file_EAR_incidence = 'EAR_incidence.data'

# Name of the file containing the EAR-model parameters and covariance matrix for cancer-mortality rates.
file_EAR_mortality = 'EAR_mortality.data'

# Name of the file containing the ERR-model parameters and covariance matrix for cancer-incidence rates.
file_ERR_incidence = 'ERR_incidence.data'

# Name of the file containing the ERR-model parameters and covariance matrix for cancer-mortality rates.
file_ERR_mortality = 'ERR_mortality.data'

#------------------------------------------------------------------------------
#    Loading the modules
#------------------------------------------------------------------------------

import os
import sys
from turtle import color
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy import interpolate

os.chdir(os.path.dirname(os.path.abspath('__file__')))

#------------------------------------------------------------------------------
#    Reading the input parameters from input file
#------------------------------------------------------------------------------

inp = []
with open('Inputs.data', 'r', encoding='UTF-8') as fin:
    for line in fin.readlines():
        line = line.replace(' ','')
        line = line.replace('\t','')
        
        if line.find('\n') == 0:
            continue

        line = line.replace('\n','')
        ps = line.find('#')

        if ps == 0:
            continue
        elif ps > 0:
            line = line[0:ps]

        inp.append(line)

# Target risk (tgt0)
tgt0 = int(inp[0])

if tgt0 == 1:
    starget = 'All solid-cancer mortality'
elif tgt0 == 2:
    starget = 'All solid-cancer incidence'
else:
    print("ERROR: Invarid value (tgt0)")
    print(tgt0)
    sys.exit(1)

# Latency represented by a Sigmoid curve with a parameter (lpy0) in the center of the Sigmoid function
lpy0 = float(inp[1])

if lpy0 < 0.0:
    print("ERROR: Invarid value (lpy0)")
    print(lpy0)
    sys.exit(1)

# Step of age for baseline data (stp0)
stp0 = int(inp[2])

if stp0 == 0:
    sstep = 'User-provided values at each one year'
elif stp0 == 1:
    sstep = 'Original class values of five years'
elif stp0 == 2:
    sstep = 'Linearly interporated values at each one year'
elif stp0 == 3:
    sstep = 'PCHIP values at each one year'
else:    
    print("ERROR: Invarid value (stp0)")
    print(stp0)
    sys.exit(1)

interp = stp0

# Gender (sex0)
sex0 = int(inp[3])

if sex0 == 1:
    sgender = 'Male'
elif sex0 == 2:
    sgender = 'Female'
else:
    print("ERROR: Invarid value (sex0)")
    print(sex0)
    sys.exit(1)

# Attained age (maxage0)
maxage0 = int(inp[4])

if maxage0 <= 0 or 100 < maxage0:
    print("ERROR: Invarid value (maxage0)")
    print(maxage0)
    sys.exit(1)

maxage = maxage0
ages = np.arange(1,maxage+1)

# Risk-transfer weight; EAR(1-wgt0) ERR(wgt0)
wgt0 = float(inp[5])

if wgt0 < 0.0 or 1.0 < wgt0:
    print("ERROR: Invarid value (wgt0)")
    print(wgt0)
    sys.exit(1)

# Calculation mode (mode0)
mode0 = int(inp[6])

if mode0 < 0:
    print("ERROR: Invarid value (mode0)")
    print(mode0)
    sys.exit(1)

# Exposure events; (dose0, agex0)
datalist = []
for i in range(7,len(inp)):
    row = []
    toks = inp[i].split(',')
    for tok in toks:
        try:
            dat = float(tok)
        except ValueError as e:
            print(e, file=sys.stderr)
            sys.exit(1)
        row.append(dat)

    datalist.append(row)

df_input = pd.DataFrame(datalist,columns=['dose0','agex0'])
df_input['maxage0'] = maxage0
df_input['wgt0'] = wgt0
num_event = len(df_input)
if num_event <= 0:
    print("ERROR: There is no input data for exposure events.")
    sys.exit(1)
minage = min(df_input['agex0'])

#------------------------------------------------------------------------------
#    Input echo
#------------------------------------------------------------------------------

print("================================")
print("Inputs:")
print("--------------------------------")
print(starget)
print(sgender)
print("--------------------------------")
print("Baseline data")
print(sstep)
print("--------------------------------")
print("Exposure event")
print(df_input)
print("================================")

#------------------------------------------------------------------------------
#    Reading the baseline rates
#------------------------------------------------------------------------------

# This is a function to extract the target data from the population data.
def f_solidc(rates,sex,title,year,istart,Const):
    crates_all = rates[(rates["Sex"] == sex) & (rates[title] == year) & (rates["Site"] == "All cancers")]
    crates_lym = rates[(rates["Sex"] == sex) & (rates[title] == year) & (rates["Site"] == "Malignant lymphoma")]
    crates_mye = rates[(rates["Sex"] == sex) & (rates[title] == year) & (rates["Site"] == "Multiple myeloma")]
    crates_leu = rates[(rates["Sex"] == sex) & (rates[title] == year) & (rates["Site"] == "Leukemia")]

    solidc = (
             crates_all.iloc[:,istart:].to_numpy()
           - crates_lym.iloc[:,istart:].to_numpy()
           - crates_mye.iloc[:,istart:].to_numpy()
           - crates_leu.iloc[:,istart:].to_numpy()
             ) / Const

    return solidc

# The population data are normalized to 100,000 Person-Year.
Const100k = 100000

# The baseline data of the all solid cancer mortality for Japanese are read from the CSV file put in the *data* directory.
cmor0 = pd.read_csv(filepath_or_buffer='data/'+file_cancer_mortality, encoding='utf-8', sep=',')

cmor_m = f_solidc(rates=cmor0,sex="Male",title="Year at death",year=2018,istart=7,Const=Const100k)
cmor_f = f_solidc(rates=cmor0,sex="Female",title="Year at death",year=2018,istart=7,Const=Const100k)

cmor_m = np.append(cmor_m,np.repeat(cmor_m[:,-1],3))
cmor_f = np.append(cmor_f,np.repeat(cmor_f[:,-1],3))
cmor_m = np.append(0,cmor_m)
cmor_f = np.append(0,cmor_f)

# The baseline data of the all solid cancer incidence for Japanese are read from the CSV file put in the *data* directory.
cinc0 = pd.read_csv(filepath_or_buffer='data/'+file_cancer_incidence, encoding='utf-8', sep=',')

cinc_m = f_solidc(rates=cinc0,sex="Male",title="Year at diagnosis",year=2018,istart=6,Const=Const100k)
cinc_f = f_solidc(rates=cinc0,sex="Female",title="Year at diagnosis",year=2018,istart=6,Const=Const100k)

cinc_m = np.append(0,cinc_m)
cinc_f = np.append(0,cinc_f)

# The baseline data of the cancer mortality/incidence rates for Japanese are given as the class values of 5 years.
age_base = np.arange(2.5,105,5)
age_base = np.append(0,age_base)

#------------------------------------------------------------------------------
# The step of age for the baseline data is set according to the input parameter stp0  
#     stp0 = interp =;  
#     0 : for user-provided baseline data at each one year,  
#     1 : original class values of 5 years,   
#     2 : values at each one year with linear intepolation,  
#     3 : values at each one year with Piecewise Cubic Hermite Interpolating Polynominal (PCHIP)
#------------------------------------------------------------------------------

age_grid = np.arange(0.5,100)

if interp == 0:
    age_base = np.arange(0.5,100)
    age_base = np.append(0,age_base)
    fit_m = interpolate.interp1d(age_base,cmor_m)
    fit_f = interpolate.interp1d(age_base,cmor_f)
elif interp == 1:
    ip = lambda x, y: interpolate.interp1d(x,y, kind='nearest')
    fit_m = ip(age_base[:-3],cmor_m[:-3])
    fit_f = ip(age_base[:-3],cmor_f[:-3])
elif interp == 2:
    fit_m = interpolate.interp1d(age_base[:-3],cmor_m[:-3])
    fit_f = interpolate.interp1d(age_base[:-3],cmor_f[:-3])
elif interp == 3:
    fit_m = interpolate.PchipInterpolator(age_base[:-3],cmor_m[:-3])
    fit_f = interpolate.PchipInterpolator(age_base[:-3],cmor_f[:-3])

# We assume that the cancer mortality rates above 87.5 years old are constant.
rates_m = pd.DataFrame(age_grid+0.5,columns=['age'])
rates_f = pd.DataFrame(age_grid+0.5,columns=['age'])

if interp != 0:
    rates_m['mortality'] = fit_m(age_grid[-13])
    rates_m['mortality'].iloc[:-12] = fit_m(age_grid[:-12])
    rates_f['mortality'] = fit_f(age_grid[-13])
    rates_f['mortality'].iloc[:-12] = fit_f(age_grid[:-12])
else:
    rates_m['mortality'] = fit_m(age_grid)
    rates_f['mortality'] = fit_f(age_grid)

# The data for the cancer incidence rates.
if interp == 0:
    fit_m = interpolate.interp1d(age_base,cinc_m)
    fit_f = interpolate.interp1d(age_base,cinc_f)
elif interp == 1:
    ip = lambda x, y: interpolate.interp1d(x,y, kind='nearest')
    fit_m = ip(age_base,cinc_m)
    fit_f = ip(age_base,cinc_f)
elif interp == 2:
    fit_m = interpolate.interp1d(age_base,cinc_m)
    fit_f = interpolate.interp1d(age_base,cinc_f)
elif interp == 3:
    fit_m = interpolate.PchipInterpolator(age_base,cinc_m)
    fit_f = interpolate.PchipInterpolator(age_base,cinc_f)

rates_m['incidence'] = fit_m(age_grid)
rates_f['incidence'] = fit_f(age_grid)

if tgt0 == 1:
    rates_m['crate'] = rates_m['mortality']
    rates_f['crate'] = rates_f['mortality']
elif tgt0 == 2:
    rates_m['crate'] = rates_m['incidence']
    rates_f['crate'] = rates_f['incidence']

# Next is setting the baseline data of the all cause mortality for Japanese. The data are read from the CSV file put in the *data* directory.
mrates0 = pd.read_csv(filepath_or_buffer='data/'+file_all_cause_mortality, encoding='utf-8', sep=',')

mrates0_m = mrates0[(mrates0["Sex"] == "Male") & (mrates0["Year at death"] == 2018)]
mrates0_f = mrates0[(mrates0["Sex"] == "Female") & (mrates0["Year at death"] == 2018)]

mrates_m = mrates0_m.iloc[:,2:].to_numpy() / Const100k
mrates_f = mrates0_f.iloc[:,2:].to_numpy() / Const100k
mrates_m = np.append(0,mrates_m.ravel())
mrates_f = np.append(0,mrates_f.ravel())

# The data are fitted, and the results are stored in the PANDAS DataFrame.
if interp == 0:
    fit_m = interpolate.interp1d(age_base,mrates_m)
    fit_f = interpolate.interp1d(age_base,mrates_f)
elif interp == 1:
    ip = lambda x, y: interpolate.interp1d(x,y, kind='nearest')
    fit_m = ip(age_base,mrates_m)
    fit_f = ip(age_base,mrates_f)
elif interp == 2:
    fit_m = interpolate.interp1d(age_base,mrates_m)
    fit_f = interpolate.interp1d(age_base,mrates_f)
elif interp == 3:
    fit_m = interpolate.PchipInterpolator(age_base,mrates_m)
    fit_f = interpolate.PchipInterpolator(age_base,mrates_f)

rates_m['mrate'] = fit_m(age_grid)
rates_f['mrate'] = fit_f(age_grid)

#------------------------------------------------------------------------------
# The Exponential survival functions for all solid cancer mortality/incidece.  
# For all solid cancer mortality;  
# $$
# S(a,s)=exp[-\int^{a}_{0}\mu_{\mathrm{total}}(\tau,s)d\tau],
# $$
# for all solid cancer incidence;  
# $$
# S_{i}(a,s)=exp[-\int^{a}_{0}[\mu_{\mathrm{total}}(\tau,s)-\mu_{i}(\tau,s)+m_{i}(\tau,s)]d\tau],
# $$
# where $a$ is the attained age; $\tau$ is the current age; $s$ is the sex; $\mu_{\mathrm{total}}$ is the baseline all-cause mortality rate; $\mu_{i}$ is the baseline cancer mortality rate for the cancer site $i$; $m_{i}$ is the baseline cancer incidence rate for the cancer site ${i}$.
#------------------------------------------------------------------------------

if tgt0 == 1:
    rates_m['survp'] = np.append(1,np.exp(-1*np.cumsum(rates_m['mrate'].iloc[:-1])))
    rates_f['survp'] = np.append(1,np.exp(-1*np.cumsum(rates_f['mrate'].iloc[:-1])))
elif tgt0 == 2:
    rates_m['survp'] = np.append(1,np.exp(-1*np.cumsum(rates_m['mrate'].iloc[:-1]-rates_m['mortality'].iloc[:-1]+rates_m['incidence'].iloc[:-1])))
    rates_f['survp'] = np.append(1,np.exp(-1*np.cumsum(rates_f['mrate'].iloc[:-1]-rates_f['mortality'].iloc[:-1]+rates_f['incidence'].iloc[:-1])))

# The baseline date are plotted here.
fig = plt.figure(figsize=(15.0, 12.0))

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)

ax1.scatter(age_base, cmor_m*Const100k, marker="o", facecolor='None', color="deepskyblue", label="Male")
ax1.scatter(age_base, cmor_f*Const100k, marker="o", facecolor='None', color="magenta", label="Female")
ax1.plot(age_grid, rates_m['mortality']*Const100k, color = "deepskyblue", label="Fitted curve")
ax1.plot(age_grid, rates_f['mortality']*Const100k, color = "magenta", label="Fitted curve")
ax1.set_xlabel('Age (Years)')
ax1.set_ylabel('All solid cancer mortality rate per 100,000 PY')
ax1.set_xlim(0,100)
ax1.set_ylim(-0.0,0.03*Const100k)
ax1.grid()
ax1.title.set_text('Mortality rates: All-solid cancer')
ax1.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

ax2.scatter(age_base, cinc_m*Const100k, marker="o", facecolor='None', color="deepskyblue", label="Male")
ax2.scatter(age_base, cinc_f*Const100k, marker="o", facecolor='None', color="magenta", label="Female")
ax2.plot(age_grid, rates_m['incidence']*Const100k, color = "deepskyblue", label="Fitted curve")
ax2.plot(age_grid, rates_f['incidence']*Const100k, color = "magenta", label="Fitted curve")
ax2.set_xlabel('Age (Years)')
ax2.set_ylabel('All solid cancer incidence rate per 100,000 PY')
ax2.set_xlim(0,100)
ax2.set_ylim(-0.0,0.04*Const100k)
ax2.grid()
ax2.title.set_text('Incidence rates: All-solid cancer')
ax2.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

ax3.scatter(age_base, mrates_m*Const100k, marker="o", facecolor='None', color="deepskyblue", label="Male")
ax3.scatter(age_base, mrates_f*Const100k, marker="o", facecolor='None', color="magenta", label="Female")
ax3.plot(age_grid, rates_m['mrate']*Const100k, color = "deepskyblue", label="Fitted curve")
ax3.plot(age_grid, rates_f['mrate']*Const100k, color = "magenta", label="Fitted curve")
ax3.set_xlabel('Age (Years)')
ax3.set_ylabel('All cause mortality rate per 100,000 PY')
ax3.set_xlim(0,100)
ax3.set_ylim(-0.0,0.40*Const100k)
ax3.grid()
ax3.title.set_text('All cause mortality rates')
ax3.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

ax4.plot(age_grid, rates_m['survp'], color = "deepskyblue")
ax4.plot(age_grid, rates_f['survp'], color = "magenta")
ax4.set_xlabel('Age (Years)')
ax4.set_ylabel('Survival probability')
ax4.set_xlim(0,100)
ax4.set_ylim(-0.0,1.0)
ax4.grid()
ax4.title.set_text('Exponential survival function: ' + starget)
ax4.legend(["Male", "Female"], frameon=False, bbox_to_anchor=(0,0), loc='lower left')

#plt.show()
plt.savefig("BaselineData.png")
#sys.exit(0)

#------------------------------------------------------------------------------
# Risk models (ERR and EAR models)
#------------------------------------------------------------------------------
# $$
# f(d,e,a,s) = \beta(0) \cdot d \cdot exp[\beta(1) \cdot \frac{(e-30)}{10} + \beta(2) \cdot log(\frac{a}{70})] \cdot \frac{(1\pm\beta(3))}{1+exp[-(a-e-7.5)]}
# $$
#------------------------------------------------------------------------------

# This is a function to define the risk model.
def func(beta, dose, agex, age, sex):
    if sex == 1:
        sign = -1
    elif sex == 2:
        sign = 1
    else:
        print("ERROR in ERR/EAR function")
        sys.exit(1)

    f = beta[0]*dose * np.exp( beta[1]*(agex-30)/10 + beta[2]*np.log(age/70) ) * (1+sign*beta[3]) / (1+np.exp(-1*(age-agex-lpy0)))
    return f

# This is a function to read the risk-model parameters.
def reading_models(fname):
    with open('data/'+fname, 'r', encoding='UTF-8') as fin:
    
        for line in fin.readlines():
            line = line.replace(' ','')
            line = line.replace('\t','')

            if line.find('n') == 0:
                continue

            line = line.replace('\n','')
            ps = line.find('#')

            if ps == 0:
                continue
            elif ps > 0:
                line = line[0:ps]

            el = len(line)
            if line[el-1] == ',':
                line = line[0:el-1]

            inp.append(line)

    for i in range(0,len(inp)):
        toks = inp[i].split(',')
        for tok in toks:
            try:
                dat = float(tok)
            except ValueError as e:
                print(e, file=sys.stderr)
                sys.exit(1)
            row.append(dat)

        datalist.append(row)

    beta = np.array(row[0:4])
    var  = np.array(row[4:20]).reshape([4,4])
    return beta, var

# Initializing.
beta_err = np.zeros(4); var_err  = np.zeros(16)
beta_ear = np.zeros(4); var_ear  = np.zeros(16)

# Reading the parameters and covariance matrices of the risk models.
# Parameters and covariance matrix for ERR
# Clearing the lists
inp.clear(); datalist.clear(); row.clear(); toks.clear()

if tgt0 == 1:
    beta_err, var_err = reading_models(fname='ERR_mortality.data')
elif tgt0 == 2:
    beta_err, var_err = reading_models(fname='ERR_incidence.data')

# Parameters and covariance matrix for EAR
# Clearing the lists
inp.clear(); datalist.clear(); row.clear(); toks.clear()

if tgt0 == 1:
    beta_ear, var_ear = reading_models(fname='EAR_mortality.data')
elif tgt0 == 2:
    beta_ear, var_ear = reading_models(fname='EAR_incidence.data')

# Plotting.
fig = plt.figure(figsize=(15.0, 6.0))

ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

pred_err = np.zeros((maxage,num_event))
pred_ear = np.zeros((maxage,num_event))

for i in range(0,num_event):
    pred_err[:,i] =func(beta=beta_err, dose=df_input.loc[i,'dose0'], agex=df_input.loc[i,'agex0'], age=ages, sex=sex0)
    pred_ear[:,i] =func(beta=beta_ear, dose=df_input.loc[i,'dose0'], agex=df_input.loc[i,'agex0'], age=ages, sex=sex0)
    ax1.plot(ages, pred_err[:,i], linestyle="solid", label='Event #'+str(i))
    ax2.plot(ages, pred_ear[:,i], linestyle="solid", label='Event #'+str(i))

ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
ax1.set_xlabel('Age (Years)')
ax1.set_ylabel('Excess relative risk')
ax1.set_xlim(0,100)
ax1.set_ylim(-0.0,max(pred_err.ravel()))
ax1.title.set_text('ERR model:' + sgender)
ax1.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper right')

ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
ax2.set_xlabel('Age (Years)')
ax2.set_ylabel('Excess absolute risk')
ax2.set_xlim(0,100)
ax2.set_ylim(-0.0,max(pred_ear.ravel()))
ax2.title.set_text('EAR model:' + sgender)
ax2.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

#plt.show()
plt.savefig("ModelPredictions.png")
#sys.exit(0)

#------------------------------------------------------------------------------
#    Monte-Carlo sampling
#------------------------------------------------------------------------------

# Number of Monte-Carlo sampling
B = 10000

# Setting the initial random seed
np.random.seed(seed=0)

# Mean values of the risk-model parameters
mean_err = beta_err.flatten()
mean_ear = beta_ear.flatten()

# Monte-Carlo sampling assuming normal distributions
mc_paras_err = np.random.multivariate_normal(mean_err, var_err, size=B)
mc_paras_ear = np.random.multivariate_normal(mean_ear, var_ear, size=B)
mc_paras_err = mc_paras_err.T
mc_paras_ear = mc_paras_ear.T

# Creating and initializing the arrays
sim_err = np.zeros((maxage,B,num_event))
sim_ear = np.zeros((maxage,B,num_event))
suv_cond = np.zeros((maxage,num_event))
br = np.zeros((maxage,num_event))
ok_err = np.zeros((maxage,num_event))
ok_ear = np.zeros((maxage,num_event))
ar_err = np.zeros((maxage,B,num_event))
ar_ear = np.zeros((maxage,B,num_event))
tr_err = np.zeros((maxage,B,num_event))
tr_ear = np.zeros((maxage,B,num_event))

# Baseline rate data for male or female
if sex0 == 1:
    rates = rates_m
elif sex0 == 2:
    rates = rates_f

# Extracting the data of the age at exposure
ageatex = df_input.loc[:,'agex0'].copy()

if mode0 >= 0:
    ageatex[:] = minage

#------------------------------------------------------------------------------
# Conditional survival function:
# $$
# S_{i}(a{|}a_{\textrm{min}},s)=\frac{S_{i}(a,s)}{S_{i}(a_{\textrm{min}},s)},
# $$
# where $a_{\textrm{min}}$ is the age at the begining of risk.  
#------------------------------------------------------------------------------

for inum in range(0,num_event):

    # Filter
    ok_err[:,inum] = np.array([True if ageatex[inum]  <= i and i <= maxage else False for i in ages])
    ok_ear[:,inum] = np.array([True if ageatex[inum]  <= i and i <= maxage else False for i in ages])

    # Conditional survival probability
    index_agex = ageatex[inum] - 1 
    index_mage = maxage -1
    suv_cond[:,inum] = rates.loc[:index_mage,"survp"] / rates.loc[index_agex,"survp"] * ok_ear[:,inum]

    # Baseline probability rate
    br[:,inum] =  rates.loc[:index_mage,"crate"] * suv_cond[:,inum]

    # Simulated ERR and EAR
    for imc in range(0,B):
        sim_err[:,imc,inum] = func(beta=mc_paras_err[:,imc], dose=df_input.loc[inum,'dose0'], agex=df_input.loc[inum,'agex0'], age=ages, sex=sex0)
        sim_ear[:,imc,inum] = func(beta=mc_paras_ear[:,imc], dose=df_input.loc[inum,'dose0'], agex=df_input.loc[inum,'agex0'], age=ages, sex=sex0)

    # Attributable probability rate
    for iage in range(0,maxage):
        ar_err[iage,:,inum] = sim_err[iage,:,inum] * br[iage,inum]
        ar_ear[iage,:,inum] = sim_ear[iage,:,inum] * suv_cond[iage,inum]
        # Total probability rate
        for imc in range(0,B):
            tr_err[iage,imc,inum] = ar_err[iage,imc,inum] + br[iage,inum]
            tr_ear[iage,imc,inum] = ar_ear[iage,imc,inum] + br[iage,inum]

### Plotting.
##from matplotlib.ticker import ScalarFormatter
##fig = plt.figure(figsize=(9.5, 4.0))
##
##ax1 = fig.add_subplot(1, 2, 1)
##ax2 = fig.add_subplot(1, 2, 2)
##
##ax1.plot(ages, sim_err[:,1:100,0], linestyle="solid", color = "pink", lw=1)
##ax1.plot(ages, np.mean(sim_err, axis=1), linestyle="solid", color = "tomato", label='MC mean')
##ax1.plot(ages, np.percentile(sim_err, axis=1, q=2.5), linestyle="dashed", color = "tomato", label='2.5th percentile')
##ax1.plot(ages, np.percentile(sim_err, axis=1, q=97.5), linestyle="dashed", color = "tomato", label='97.5th percentile')
##ax1.set_xlabel('Attained age (Years)')
##ax1.set_ylabel('ERR')
##ax1.set_xlim(0,90)
##ax1.set_ylim(-0.00,1.50)
##ax1.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper right')
##
##ax2.plot(ages, sim_ear[:,1:100,0], linestyle="solid", color = "pink", lw=1)
##ax2.plot(ages, np.mean(sim_ear, axis=1), linestyle="solid", color = "tomato", label='MC mean')
##ax2.plot(ages, np.percentile(sim_ear, axis=1, q=2.5), linestyle="dashed", color = "tomato", label='2.5th percentile')
##ax2.plot(ages, np.percentile(sim_ear, axis=1, q=97.5), linestyle="dashed", color = "tomato", label='97.5th percentile')
##ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
##ax2.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
##ax2.set_xlabel('Attained age (Years)')
##ax2.set_ylabel('EAR')
##ax2.set_xlim(0,90)
##ax2.set_ylim(-0.00,3.5e-3)
##ax2.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')
##
###plt.show()
##plt.savefig("MonteCarloSampling.eps",dpi=300)
##sys.exit(1)

# Plotting
fig = plt.figure(figsize=(20.0, 6.0))

ax1 = fig.add_subplot(1, 3, 1)
ax2 = fig.add_subplot(1, 3, 2)
ax3 = fig.add_subplot(1, 3, 3)

ax1.plot(ages, rates.loc[:index_mage,"crate"], linestyle="solid")
ax1.set_xlabel('Age (Years)')
ax1.set_ylabel('m(a,g)')
ax1.set_xlim(0,100)
ax1.set_ylim(1.0e-05,1.0e-01)
ax1.set_yscale('log')
ax1.title.set_text('Baseline '+starget+' rate: ' + sgender)

for inum in range(0,num_event):
    ax2.plot(ages, suv_cond[:,inum], linestyle="solid", label='Event #'+str(inum))
    ax3.plot(ages, br[:,inum], linestyle="solid", label='Event #'+str(inum))

ax2.set_xlabel('Age (Years)')
ax2.set_ylabel('S(a,g)/S(e,g)')
ax2.set_xlim(0,100)
ax2.set_ylim(-0.00,1.1)
ax2.title.set_text('Conditional survival probability: ' + sgender)
ax2.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper right')

ax3.set_xlabel('Age (Years)')
ax3.set_ylabel('m(a,g)*S(a,g)/S(e,g)')
ax3.set_xlim(0,100)
ax3.set_ylim(-0.00,0.020)
ax3.title.set_text('Baseline probability rate for '+starget+': ' + sgender)
ax3.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

#plt.show()
plt.savefig("BaselineRisk.png")
#sys.exit(1)

#------------------------------------------------------------------------------
#    Attributable probability rate
#------------------------------------------------------------------------------
# $$
# AR_{i}(d,e,a|a_{\textrm{min},s}) = M_{i}(d,e,a,s){\times}S_{i}(a|a_{\textrm{min}},s),
# $$
# where $M_{i}(d,e,a,s)$ is the excess cancer mortality/incidence rates.  
#------------------------------------------------------------------------------

# Adding the attributable probability rate per year for each exposure events
ar_err_esum = np.sum(ar_err, axis=2)
ar_ear_esum = np.sum(ar_ear, axis=2)
tr_err_esum = np.sum(tr_err, axis=2)
tr_ear_esum = np.sum(tr_ear, axis=2)

# Plotting.
fig = plt.figure(figsize=(15.0, 6.0))

ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

ax1.plot(ages, ar_err_esum[:,1:100], linestyle="solid", color = "pink", lw=1)
ax1.plot(ages, np.mean(ar_err_esum, axis=1), linestyle="solid", color = "tomato", label='MC mean')
ax1.plot(ages, np.percentile(ar_err_esum, axis=1, q=2.5), linestyle="dashed", color = "tomato", label='2.5th percentile')
ax1.plot(ages, np.percentile(ar_err_esum, axis=1, q=97.5), linestyle="dashed", color = "tomato", label='97.5th percentile')
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
ax1.set_xlabel('Attained age (Years)')
ax1.set_ylabel('M(D,e,a,g)*S(a,g)/S(e,g)')
ax1.set_xlim(0,100)
ax1.set_ylim(-0.00,max(ar_err_esum.ravel()))
ax1.title.set_text(sgender + ', M(D,e,a,g) = ERR(D,e,a,g) * m(a,g)')
ax1.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

ax2.plot(ages, ar_ear_esum[:,1:100], linestyle="solid", color = "pink", lw=1)
ax2.plot(ages, np.mean(ar_ear_esum, axis=1), linestyle="solid", color = "tomato", label='MC mean')
ax2.plot(ages, np.percentile(ar_ear_esum, axis=1, q=2.5), linestyle="dashed", color = "tomato", label='2.5th percentile')
ax2.plot(ages, np.percentile(ar_ear_esum, axis=1, q=97.5), linestyle="dashed", color = "tomato", label='97.5th percentile')

ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
ax2.set_xlabel('Attained age (Years)')
ax2.set_ylabel('M(D,e,a,g)*S(a,g)/S(e,g)')
ax2.set_xlim(0,100)
ax2.set_ylim(-0.00,max(ar_ear_esum.ravel()))
ax2.title.set_text(sgender + ', M(D,e,a,g) = EAR(D,e,a,g)')
ax2.legend(frameon=False, bbox_to_anchor=(0,1), loc='upper left')

#plt.show()
plt.savefig("AttributableRisks.png")
#sys.exit(0)

#------------------------------------------------------------------------------
#    Cumulative excess risk (CER)
#------------------------------------------------------------------------------
# $$
# CER_{i}(d,e,a_{\textrm{min}},s)=\sum_{a=a_{\textrm{min}}}^{a_{\textrm{max}}} AR_{i}(d,e,a|a_{\textrm{min}},s)
# $$
#------------------------------------------------------------------------------

# Integration to the maximum attained age
lar_err = np.sum(ar_err_esum, axis=0)
lar_ear = np.sum(ar_ear_esum, axis=0)

# Cumulative excess risk (CER)
LARs = (1-wgt0)*lar_ear + wgt0*lar_err

#------------------------------------------------------------------------------
#    Output of the calculation results
#------------------------------------------------------------------------------

result = pd.DataFrame(data=[np.mean(LARs, axis=0),np.percentile(LARs, axis=0, q=50),
                            np.percentile(LARs, axis=0, q=2.5),np.percentile(LARs, axis=0, q=97.5)])
result = result.T

result.index=['CER']
result.columns=['Mean','Median','2.5th','97.5th']

print("Results: "+starget)
print("--------------------------------")
print(result)
print("================================")

#------------------------------------------------------------------------------
#    Generating the CSV files for the calculation results
#------------------------------------------------------------------------------

# Input echo
comments = [
           ["# INPUT DATA:"],
           ["# Target = ",starget],
           ["# Sex    = ",sgender],
           ["# Weight = ",wgt0],
           ["# Baseline data = ",sstep],
           ['# Exposure event','Dose','Age at exposure'],
           ['#','Gy','Year']
           ]

for i in range(len(df_input)):
    comments.append(["#"+str(i),df_input.loc[i,"dose0"],df_input.loc[i,'agex0']])
comments.append(["#"])

# Outputing the CERs
with open('results_CER.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(comments)
    writer.writerow(["# OUTPUT:"])
result.to_csv('results_CER.csv',mode='a',float_format='%10.6f')

# Function to output the results.
def f_results(filename,title,age,mean,median,low,hi):
    res = pd.DataFrame(data=[age,mean,median,low,hi])
    res = res.T
    with open(filename, 'w', newline='') as outf:
        writer = csv.writer(outf)
        writer.writerows(comments)
        writer.writerow([title])
    res.columns=['Age','Mean','Median','2.5th','97.5th']
    res.to_csv(filename,mode='a',index=False,float_format='%10.6f')
    return

# Outputing the attributable probability rates (APRs) based on ERR model.
f_results(
    filename="results_APR_ERR.csv",title="# OUTPUT: APR-ERR",age=ages,
    mean=np.mean(ar_err_esum, axis=1),median=np.percentile(ar_err_esum, axis=1, q=50),
    low=np.percentile(ar_err_esum, axis=1, q=2.5),hi=np.percentile(ar_err_esum, axis=1, q=97.5)
    )

# Outputing the APRs based on EAR model.
f_results(
    filename="results_APR_EAR.csv",title="# OUTPUT: APR-EAR",age=ages,
    mean=np.mean(ar_ear_esum, axis=1),median=np.percentile(ar_ear_esum, axis=1, q=50),
    low=np.percentile(ar_ear_esum, axis=1, q=2.5),hi=np.percentile(ar_ear_esum, axis=1, q=97.5)
    )

if num_event > 1:
    print("Finish.")
    sys.exit(0)
elif mode0 == 0 and num_event == 1:
    print("Generating database of CERs for a single acute exposure")
else:
    print("Finish.")
    sys.exit(0)

#------------------------------------------------------------------------------
#    Summary of CERs for male and female
#------------------------------------------------------------------------------

def calculate_LARs_acute(agex, dose, sex, ages, drates, wgt):

    ok_err_acute = np.array([True if agex <= i and i <= max(ages) else False for i in ages])
    ok_ear_acute = np.array([True if agex <= i and i <= max(ages) else False for i in ages])

    index_agex = agex - 1 

    suv_err_acute =  drates["crate"] * drates["survp"] / drates.loc[index_agex,"survp"] * ok_err_acute
    suv_ear_acute =  drates["survp"] / drates.loc[index_agex,"survp"] * ok_ear_acute

    suv_err_value = suv_err_acute.values.reshape(1,-1)
    suv_ear_value = suv_ear_acute.values.reshape(1,-1)

    for imc in range(0,B):
        sim_err_acute[:,imc] = func(beta=mc_paras_err[:,imc], dose=dose, agex=agex, age=ages, sex=sex)
        sim_ear_acute[:,imc] = func(beta=mc_paras_ear[:,imc], dose=dose, agex=agex, age=ages, sex=sex)

    lar_err_acute = np.dot(suv_err_value, sim_err_acute).ravel()
    lar_ear_acute = np.dot(suv_ear_value, sim_ear_acute).ravel()

    LARs_acute = (1-wgt)*lar_ear_acute + wgt*lar_err_acute
    return LARs_acute

# Age at exposure
agexs = np.arange(0,maxage-5) + 1 

# Attained age
ages = np.arange(1,maxage+1)

# Radiation dose
dose = df_input.loc[0,'dose0']

sim_err_acute = np.zeros((len(ages),B))
sim_ear_acute = np.zeros((len(ages),B))
lars_m = np.zeros((B, len(agexs)))
lars_f = np.zeros((B, len(agexs)))
summary_m = np.zeros((4,len(agexs)))
summary_f = np.zeros((4,len(agexs)))

for iagex in range(0,len(agexs)):

    lars_m[:,iagex] = calculate_LARs_acute(agex=agexs[iagex], dose=dose, sex=1, ages=ages, drates=rates_m.iloc[:maxage], wgt=wgt0)
    lars_f[:,iagex] = calculate_LARs_acute(agex=agexs[iagex], dose=dose, sex=2, ages=ages, drates=rates_f.iloc[:maxage], wgt=wgt0)

    summary_m[0,iagex] = np.mean(lars_m[:,iagex], axis=0)
    summary_m[1,iagex] = np.percentile(lars_m[:,iagex], axis=0, q=2.5)
    summary_m[2,iagex] = np.percentile(lars_m[:,iagex], axis=0, q=50)
    summary_m[3,iagex] = np.percentile(lars_m[:,iagex], axis=0, q=97.5)
    summary_f[0,iagex] = np.mean(lars_f[:,iagex], axis=0)
    summary_f[1,iagex] = np.percentile(lars_f[:,iagex], axis=0, q=2.5)
    summary_f[2,iagex] = np.percentile(lars_f[:,iagex], axis=0, q=50)
    summary_f[3,iagex] = np.percentile(lars_f[:,iagex], axis=0, q=97.5)

summary_m = pd.DataFrame({ 'Agex'  : agexs,
                           'Mean'  : summary_m[0,:],
                           '2.5-P' : summary_m[1,:],
                           '50-P'  : summary_m[2,:],
                           '97.5-P': summary_m[3,:] })
summary_f = pd.DataFrame({ 'Agex'  : agexs,
                           'Mean'  : summary_f[0,:],
                           '2.5-P' : summary_f[1,:],
                           '50-P'  : summary_f[2,:],
                           '97.5-P': summary_f[3,:] })

# Generating the CSV files for the calculation results
# Input echo
comments = [
           ["# INPUT DATA:"],
           ["# Target = ",starget],
           ["# Sex    = ",sgender],
           ["# Maximum attained age =",maxage0],
           ["# Weight = ",wgt0],
           ["# Baseline data = ",sstep],
           ['# Exposure event','Dose'],
           ['#','Gy']
           ]

for i in range(len(df_input)):
    comments.append(["#"+str(i),df_input.loc[i,"dose0"]])
comments.append(["#"])

# Outputing the LARs
with open('summary_CER_male.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(comments)
    writer.writerow(["# OUTPUT:"])
summary_m.to_csv('summary_CER_male.csv',mode='a',index=False,float_format='%10.6f')

with open('summary_CER_female.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(comments)
    writer.writerow(["# OUTPUT:"])
summary_f.to_csv('summary_CER_female.csv',mode='a',index=False,float_format='%10.6f')

# Plotting.
fig = plt.figure(figsize=(15.0, 6.0))

ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)

ax1.plot(summary_m["Agex"], summary_m["Mean"]*100, color='blue', label='Mean')
ax1.fill_between(summary_m["Agex"], summary_m["97.5-P"]*100, summary_m["2.5-P"]*100, alpha=.2, color='blue', label='95% confidence interval')

ax1.set_xlabel('Age at exposure (Years)')
ax1.set_ylabel('CER (%)')
ax1.set_xlim(1,max(agexs))
ax1.set_ylim(-0.00,max(lars_f.ravel()*100))
ax1.title.set_text('Acute exposure' + ': Dose =' + str(dose) + ' (Gy)' + '; Male')
ax1.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper right')

ax2.plot(summary_f["Agex"], summary_f["Mean"]*100, color='magenta', label='Mean')
ax2.fill_between(summary_f["Agex"], summary_f["97.5-P"]*100, summary_f["2.5-P"]*100, alpha=.2, color='magenta', label='95% confidence interval')
ax2.set_xlabel('Age at exposure (Years)')
ax2.set_ylabel('CER (%)')
ax2.set_xlim(1,max(agexs))
ax2.set_ylim(-0.00,max(lars_f.ravel()*100))
ax2.title.set_text('Acute exposure' + ': Dose =' + str(dose) + ' (Gy)' + '; Female')
ax2.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper right')

#plt.show()
plt.savefig("CERs.png")