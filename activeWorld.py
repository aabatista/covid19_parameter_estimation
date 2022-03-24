# -*- coding: utf-8 -*-
# Author: Adriano A. Batista, 2020
# activeWorld.py -- Model epidemiological data from any country in the world
#
'''
 command line:
   python3 activeWorld.py configCountry.txt

 Input: 
   csv data files:
      time_series_covid19_confirmed_global_narrow.csv
      time_series_covid19_deaths_global_narrow.csv
      time_series_covid19_recovered_global_narrow.csv
 Esses arquivos podem ser obtidos do site:
    https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases
 configCountry.txt file example:
    country,Brazil
    delay,14
    offset,22
    cutoff,0.43
 Output
 1. confirmed, recovered, death, and active cases plots of comparison 
    between data and corresponding model curve fits.
 2. contagion rate time series, Rt, and lethality probability
'''
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import scipy.integrate
from scipy import stats
from scipy.interpolate import UnivariateSpline
import pandas as pd
import numpy as np
import sys
import time
from datetime import datetime
from collections import OrderedDict
import unidecode
plt.rcParams['axes.grid'] = True
def set_xaxis(ax, major, minor):
  ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
  ax.xaxis.set_major_locator(mdates.DayLocator(interval=major))
  ax.xaxis.set_minor_locator(mdates.DayLocator(interval=minor))


def read_setup(inSetupFile):
  # open SIRD model parameter data 
  # Format of setup file
  # parameter      value
  '''
  country             country name
  delay               tau
  offset              another delay
  cutoff              maximum value of  kappa(t)
  nu                  daily birth rate
  mu                  daily death rate (pre covid pandemic)
  '''
  global country, delay, offset, cutoff
  # Read setup file
  lines = np.genfromtxt(inSetupFile, dtype = 'str', delimiter=',')
  # build a dictionary
  setup_vars = OrderedDict(lines)

  for var in lines:
    if var[0] == 'country':
      setup_vars[var[0]] = var[1]
    else:
      if var[1].isdigit():
        setup_vars[var[0]] = int(var[1])
      else:
        setup_vars[var[0]] = float(var[1])
  # Transform the keys into variables
  globals().update(setup_vars)
read_setup(sys.argv[1])
# create log file
logFile = 'log'+country+'.txt'
flog = open(logFile, 'w')
flog.write('country %s, delay %d, offset %d, cutoff %g\n' %(country, delay, offset,
  cutoff))
##########################################################
# 1. Data processing
csvFiles = ['time_series_covid19_confirmed_global_narrow.csv',
    'time_series_covid19_deaths_global_narrow.csv',
    'time_series_covid19_recovered_global_narrow.csv']
# a. Open confirmed datafile create dataframe
df = pd.read_csv(csvFiles[0], sep =',', skiprows=0, header=0, names=["ProvinceState", "Country/Region", "Lat", "Long", "Date", "Value", "ISO 3166-1 Alpha 3-Codes", "Region Code", "Sub-region Code", "Intermediate Region Code"])
# Filter country data
df = df[df['Country/Region'].str.contains(country)]
df = df.sort_values('Date', ascending=True)
confirmed = df.Value.to_numpy().astype(int)
ind0 = (confirmed!=0).argmax() # first nonzero index of confirmed
ind0 += offset
confirmed = confirmed[ind0:] # cut off the leading zeros
flog.write('ind0 %d, offset %d\n' % (ind0, offset))
dates = pd.to_datetime(df.Date) # covert string dates to datetime format
dates = dates[ind0:]
# b. Open deaths datafile, create dataframe.
df = pd.read_csv(csvFiles[1], sep =',', skiprows=0, header=0, names=["ProvinceState", "Country/Region", "Lat", "Long", "Date", "Value", "ISO 3166-1 Alpha 3-Codes", "Region Code", "Sub-region Code", "Intermediate Region Code"])
df = df[df['Country/Region'].str.contains(country)]
deaths = df.Value.to_numpy().astype(int)
deaths = deaths[::-1]
deaths = deaths[ind0:]
indFirstDeath = (deaths!=0).argmax()
flog.write('indFirstDeath %d, First death on %s, number of deaths %d\n'
  % (indFirstDeath, dates.iloc[indFirstDeath], deaths[indFirstDeath]))
# c. Open recovered datafile, create dataframe.
df = pd.read_csv(csvFiles[2], sep =',', skiprows=0, header=0, names=["ProvinceState", "Country/Region", "Lat", "Long", "Date", "Value", "ISO 3166-1 Alpha 3-Codes", "Region Code", "Sub-region Code", "Intermediate Region Code"])
df = df[df['Country/Region'].str.contains(country)]
recovered = df.Value.to_numpy().astype(int)
recovered = recovered[::-1]
recovered = recovered[ind0:]

# d. generate active cases from data
activeCases = confirmed-deaths-recovered

# e. Obtain country population (World Bank data)
df = pd.read_csv('worldPopulation.csv', sep =',', skiprows=0, header=0,
    names=['Country_Code','Ranking','Country','Population'])
df = df[df['Country'].str.contains(country)]
population = df['Population'].iloc[0]
population = population.replace(' ','')
population = population.replace(',','')
print(df)
P_0 = int(population)*1000
flog.write('População em 2019 %d\n' % (P_0))

# 2. Statistical analysis
# Active cases delay estimate
activeDelay = confirmed[delay:]-confirmed[:-delay]
n_avg=delay
day = 24 # day in hours
dt = 1.0/day  # integration time-step
tau = (1+n_avg*day)*dt # average time duration of infection
# Active cases statistical estimate
q = n_avg/(1+n_avg) # probability to remain sick after 1 day
newConf = np.diff(confirmed) # casos novos confirmados ao dia
newDead = np.diff(deaths) # casos novos confirmados ao dia
N = len(newConf)
activeStat = np.ones(N)
n_pow = np.arange(N)
n_pow_des = n_pow[::-1]
q_to_n = q**n_pow_des
for i in n_pow:
  activeStat[i] = np.sum(newConf[:i]*q_to_n[N-i:])

# 3. Parameter estimation
# A. Contagion rate function
contagionRate = []
N_as = len(activeStat)
for i in np.arange(0, N_as):
  if i+7<=N_as:
    XX = np.arange(i, i+7)
    ZZ = activeStat[i:i+7]
  else:
    XX = np.arange(i-7, i)
    ZZ = activeStat[i-7:i]
  if len(XX)==len(ZZ):
    slope, intercept, r_value, p_value, std_err = stats.linregress(XX, ZZ)
    I_avg = np.sum(ZZ)/7
    S = 1.0-confirmed[i]/P_0
    if I_avg==0:
      ka_t=0
    else:
      ka_t = (slope/I_avg+mu+1.0/tau)/S
    #print(i, ka_t)
    if ka_t<0:
      ka_t = 0
    elif ka_t>cutoff:
      ka_t = cutoff
    contagionRate.append(ka_t)
print(len(dates[1:-7]), len(contagionRate))

##########################################################
def kappa(t):
  ind = int(1.0*t)
  N=len(contagionRate)
  if ind>=N:
    ind = N-1
  ka = contagionRate[ind]
  if ind>0 and ind<N and ka==0:
    ka = (contagionRate[ind-1]+contagionRate[ind+1])/2
  return ka
##########################################################
# vectorize contagion rate
kappa_t_vec = []
n_max = len(confirmed)-1 # in days
tt = np.arange(0.0, n_max, dt)
for ind, t in enumerate(tt[::day]):
  ka = kappa(t)
  kappa_t_vec.append(ka)
kappa_t_vec = np.array(kappa_t_vec)

# C. Lethality probability
P_let = np.zeros(N)
for i in n_pow[:-1]:
    denom = np.sum(newConf[:i]*q_to_n[N-i:])
    if denom>0:
      s = newDead[i]/denom
    else:
      s=0
    P_let[i] = s/(1-q)
print('P_let')
# Lethality and recovery rates
letRate = P_let/tau
rho = 1.0/tau-letRate
P_rho=1-P_let

# 4. Numerical model epidemic evolution (modified SIR model)
##########################################################
def get_l_r(t):
  ind = int(1.0*t)
  N=len(P_let)
  if ind>=N:
    ind = N-1
  la = P_let[ind]/tau
  rho = 1.0/tau - la
  return la, rho
##########################################################
def derivs (x, t): # return derivatives of the array x
   S = x[0] # susceptíveis
   I = x[1] # infectados
   R = x[2] # Recuperados
   M = x[3] # Mortos
   T = S+I+R
   l, r = get_l_r(t)
   dSdt = nu*T-mu*S-kappa(t)*S*I
   dIdt = -mu*I+kappa(t)*S*I-r*I-l*I
   dRdt = r*I-mu*R
   dMdt = l*I
   return [dSdt, dIdt, dRdt, dMdt]
##########################################################
C_0 = confirmed[0]
R_0 = recovered[0]
D_0 = deaths[0]
A_0 = C_0-R_0-D_0
flog.write('C_0 %g, A_0 %g, R_0 %g, D_0 %g\n' % (C_0, A_0, R_0, D_0))
yinit = [1.0-C_0/P_0, A_0/P_0, R_0/P_0, D_0/P_0] # initial values
y = scipy.integrate.odeint(derivs, yinit, tt)
S = y[:, 0]
I = y[:, 1]
R = y[:, 2]
M = y[:, 3]

# 3B. generate Rt
R0_t = tau*kappa_t_vec*S[::day]


# 5. Plot results
# Figure 1
fig1=plt.figure('kappa', figsize=(6, 12))
plt.subplots_adjust(hspace=0.25)
# A. Contagion rate estimate
ax1 = fig1.add_subplot(311)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=28))
ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
plt.gcf().autofmt_xdate()
ax1.set_ylabel('$\kappa(t)$')
ax1.plot(dates[1:], kappa_t_vec, 'k-', ms=2)
ax1.set_title('Contagion rate in %s' % country)
ax1.set_xlim(dates.iloc[0], dates.iloc[-1])
ax1.set_ylim(0., )
ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes, size=12, weight='bold')
plt.gcf().autofmt_xdate()


# B. Dynamic reproductive number R_0(t)
ax2 = fig1.add_subplot(312)
ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes, size=12, weight='bold')
ax2.set_title('Reproduction number')
#ax2.set_xticks(dates[0::7])
ax2.xaxis.set_major_locator(mdates.DayLocator(interval=28))
ax2.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
ax2.yaxis.set_major_locator(MultipleLocator(1))
ax2.yaxis.set_minor_locator(AutoMinorLocator(0.5))
ax2.plot(dates[1:], R0_t, 'k-')
ax2.set_ylabel('$R_0(t)$')
ax2.set_xlim(dates.iloc[0], dates.iloc[-1])
ax2.set_ylim(0., )

# C. Lethality probability, P_\lambda(t), and recovery probability, P_\rho(t)
ax3 = fig1.add_subplot(313)
ax3.set_title('Lethality and recovery probabilities')
fig1.align_ylabels([ax1, ax2, ax3])
ax3.text(-0.1, 1.1, 'C', transform=ax3.transAxes, size=12, weight='bold')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=28))
plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=7))
ax3.plot(dates[1:], P_let, 'k-', linewidth=2.0, label="Lethality probability estimate")
ax3.plot(dates[1:], P_rho, 'b-.', linewidth=1.5,  label="Recovery probability estimate")
ax3.set_xlim(dates.iloc[0], dates.iloc[-1])
ax3.set_ylim(0., 1.)
plt.gcf().autofmt_xdate()
plt.legend(fontsize=10)
plt.gcf().autofmt_xdate()
flog.write('first P_let avg %g\n' % (np.sum(P_let[:21])/21))
flog.write('last P_let avg %g\n' % (np.sum(P_let[-21:])/21))
figure = 'contagionLetRecRates_%s.pdf'% country
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
plt.savefig(figure, bbox_inches='tight')
flog.write('%s\n' % figure)
# Figure 2: confirmed, recovered, deaths, and active
fig2, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

# A. confirmed cases
ax= axs[0,0]
ax.set_title(u'Confirmed cases')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=28))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
legenda = u"Confirmed cases in %s" % (country)
#ax.plot(dates, confirmed, 'k-o', ms=3,  label= legenda)
ax.text(-0.05, 1.02, 'A', transform=ax.transAxes, size=12, weight='bold')
ax.set_xlim(dates.iloc[0], dates.iloc[-1])

# B. recovered cases
ax= axs[0,1]
ax.set_title('Recovered cases')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=28))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
legenda = u"Recovered cases in %s" % (country)
#ax.plot(dates, recovered, 'k-o', ms=3,  label= legenda)
ax.text(-0.05, 1.02, 'B', transform=ax.transAxes, size=12, weight='bold')
ax.set_xlim(dates.iloc[0], dates.iloc[-1])

# C. death cases
ax= axs[1,0]
ax.set_title(u'Death cases')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=28))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
legenda = u"Death cases in %s" % (country)
#ax.plot(dates, deaths, 'k-o', ms=3,  label= legenda)
ax.text(-0.05, 1.02, 'C', transform=ax.transAxes, size=12, weight='bold')
ax.set_xlim(dates.iloc[0], dates.iloc[-1])

# D. active cases
ax= axs[1,1]
ax.set_title(u'Active cases')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=28))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
ax.plot(dates[delay:], activeDelay, 'b+', ms=3,  label= 'estimate, delay %d days' % delay)
ax.plot(dates[1:], activeStat, 'gx', ms=3,  label= 'Statistical estimate')
ax.text(-0.05, 1.02, 'D', transform=ax.transAxes, size=12, weight='bold')
ax.set_xlim(dates.iloc[0], dates.iloc[-1])
plt.gcf().autofmt_xdate()

# Forecasting part
#########################################################################
def getTransMatrix(array):
  dk = np.copy(array) # make a copy
  print('dk', dk)
  dk[dk>0] = 1
  dk[dk<0] = 0
  dk = dk.astype(int)
  print('dk', dk)
  # Create subarrays of 0's and 1's
  subArrays=np.split(dk, np.argwhere(np.diff(dk) != 0)[:,0] + 1)
  print('subArrays', subArrays)
  ones = []
  zeros = []
  for X in subArrays:
    if X[0] == 1:
      ones.append(len(X))
    elif X[0] == 0:
      zeros.append(len(X))
  print('ones', ones)
  print('zeros', zeros)
  ones = np.array(ones)
  zeros = np.array(zeros)
  avgNpos = ones.mean()
  avgNneg = zeros.mean()
  print('avgNpos', avgNpos, 'avgNneg', avgNneg)
  q_pos = avgNpos/(1+avgNpos)
  p_pos = 1-q_pos
  q_neg = avgNneg/(1+avgNneg)
  p_neg = 1-q_neg
  print('q_pos', q_pos, 'p_pos', p_pos, 'q_neg', q_neg, 'p_neg', p_neg)
  return [[q_neg, p_neg], [p_pos, q_pos]]
#######################################################################
def getProbs(dk, Nbins):
  dkPos = dk[dk>=0]
  dkNeg = dk[dk<=0]
  probPos, p_edges = np.histogram(dkPos, bins = Nbins, density=True)
  probNeg, n_edges = np.histogram(dkNeg, bins = Nbins, density=True)
  probPos /= np.sum(probPos)
  probNeg /= np.sum(probNeg)
  # edges at middle of bins
  p_edges = (p_edges[1:]+p_edges[:-1])/2
  n_edges = (n_edges[1:]+n_edges[:-1])/2
  print(len(probPos), len(p_edges))
  print('p_edges', p_edges, 'n_edges', n_edges)
  return probPos, p_edges, probNeg, n_edges
#######################################################################
# Generate forecasts based on Markov Chains for 14 days
Nfcast = 14
nRuns = 1000 # average and confidence interval from these runs
sampleSize = 10
priors = -21
nSamples = int(nRuns/sampleSize)
shift = -14 # it has to be zero or negative, in days
nBinsHist = 40
colors = plt.cm.jet(np.linspace(0, 1, nBinsHist+1))
dia0= dates.iloc[shift]
dtime = pd.date_range(start = dia0, periods=Nfcast)
print("kappa_t", kappa_t_vec[priors+shift:shift])
dk = np.diff(kappa_t_vec[priors+shift:shift])
Nbins = int(len(dk)/2)
print('Nbins', Nbins)
# The state space
states = [0, 1]
# Transition probability matrix
# Markov chain part for contagion rate
W = transitionMatrix = getTransMatrix(dk)
# Get entropy
entropyS = -np.sum(W*np.log(W))
flog.write("Entropia S %g\n" % entropyS)
probPos, p_edges, probNeg, n_edges = getProbs(dk, Nbins)
dl = np.diff(letRate[priors+shift:shift])
letTransMat = getTransMatrix(dl)
letProbPos, lp_edges, letProbNeg, ln_edges = getProbs(dl, Nbins)
# Loop
  # join kappa_t_vec and kappa_for
  # Generate forecasted time series of S(t), I(t), R(t), M(t)
# Obtain most probable forecast and confidence interval

# Plot histograms
fig, Axs = plt.subplots(2, figsize=(8, 8), sharex=True)
fig.align_ylabels(Axs)
fig.tight_layout()
ax= Axs[0]
dkPos = dk[dk>=0]
dkNeg = dk[dk<=0]
ax.hist(dkPos, Nbins, color='blue', label= 'Increments', ec='k')
ax.text(-0.1, 1.00, 'A', transform=ax.transAxes, size=12, weight='bold')
ax.legend()
ax= Axs[1]
ax.hist(dkNeg, Nbins, color='red', label= 'Decrements', ec='k')
ax.text(-0.1, 1.00, 'B', transform=ax.transAxes, size=12, weight='bold')
ax.legend()
figure = 'histograms%s.pdf'% (country)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
plt.savefig(figure, bbox_inches='tight')
flog.write('%s\n' % figure)

# Analytical integration
Sf = np.zeros(Nfcast, dtype=float)
If = np.zeros(Nfcast, dtype=float)
Rf = np.zeros(Nfcast, dtype=float)
Mf = np.zeros(Nfcast, dtype=float)
# allocate space for multiple trajectories
mSfs = np.zeros((nSamples, Nfcast), dtype=float)
mIfs = np.zeros((nSamples, Nfcast), dtype=float)
mRfs = np.zeros((nSamples, Nfcast), dtype=float)
mMfs = np.zeros((nSamples, Nfcast), dtype=float)
S_avg = np.zeros(Nfcast, dtype=float)
I_avg = np.zeros(Nfcast, dtype=float)
R_avg = np.zeros(Nfcast, dtype=float)
M_avg = np.zeros(Nfcast, dtype=float)
ka_forecast = np.zeros(Nfcast)
ka_avg = np.zeros(Nfcast)
#######################################################################
# Contagion rate forecasting
def getForecast(initValue, transMat, p_edges, probPos, n_edges, probNeg):
  up_or_down = np.random.choice(states)
  forecast = np.zeros(Nfcast, dtype = float)
  value = initValue
  for i in range(0, Nfcast):
    # decide to increase or decrease value
    up_or_down = np.random.choice(states, p=transMat[up_or_down])
    if up_or_down == 1:
      # increment value
      value += np.random.choice(p_edges, p=probPos) 
    elif up_or_down == 0:
      # decrement value
      value += np.random.choice(n_edges, p=probNeg) 
      if value<0.0:
        up_or_down == 1
        value = np.random.choice(p_edges, p=probPos) 
    forecast[i] = value
  return forecast
#######################################################################
import colorline
count = 0
for indT in range(nRuns):
  # Contagion rate forecasting
  ka = kappa_t_vec[shift]
  ka_forecast = getForecast(ka, transitionMatrix, p_edges, probPos, n_edges, probNeg)
  ka_avg += ka_forecast
  # Lethality rate forecasting
  la = letRate[shift]
  let_forecast = getForecast(la, letTransMat, lp_edges, letProbPos, ln_edges,
      letProbNeg)
  rec_forecast = 1.0/tau-let_forecast
  #print('kappa forecast', ka_forecast)
  # Obtain forecast of S(t), I(t), R(t), and M(t) based on knowledge of
  # kappa(t)
  # initial conditions
  C_0 = confirmed[shift]
  D_0 = deaths[shift]
  A_0 = activeCases[shift]
  R_0 = recovered[shift]
  # initial values  
  Sf[0] = 1.0-C_0/P_0
  If[0] = A_0/P_0
  Rf[0] = R_0/P_0
  Mf[0] = D_0/P_0
  # Generate new trajectory 
  for i in np.arange(1, Nfcast):
    ka_i = ka_forecast[i-1]
    #if ka_i<0:
    #  ka_i = 0.0
    factor = -1.0/tau+ka_i*Sf[i-1]
    If[i] = If[i-1]*np.exp(factor)
    Delta_I = If[i]-If[i-1]
    Sf[i] = Sf[i-1]*np.exp(-ka_i*Delta_I/factor)
    Rf[i] = Rf[i-1]+rec_forecast[i-1]*Delta_I/factor
    Mf[i] = Mf[i-1]+let_forecast[i-1]*Delta_I/factor
    if ka_i<0:
      print('Error! kappa negativo', ka_i, 'i', i)
      sys.exit(0)
    if Mf[i]<Mf[i-1] or Rf[i]<Rf[i-1]:
      print('Integration error!')
      sys.exit(0)
    if (If[i]+Rf[i]+Mf[i])<(If[i-1]+Rf[i-1]+Mf[i-1]):
      print('indT', indT, 'ka_i', ka_i, 'factor', factor)
      print('Integration error!')
      #sys.exit(0)
  # find average trajectory
  S_avg = S_avg+Sf
  I_avg = I_avg+If
  R_avg = R_avg+Rf
  M_avg = M_avg+Mf
  # pre-average
  if  (indT+1) % sampleSize ==0:
    S_avg /= sampleSize
    I_avg /= sampleSize
    R_avg /= sampleSize
    M_avg /= sampleSize
    mSfs[count,:] = S_avg
    mIfs[count,:] = I_avg
    mRfs[count,:] = R_avg
    mMfs[count,:] = M_avg
    I_avg = S_avg = R_avg= M_avg = 0*I_avg
    count += 1
# Filter out the 5% highest trajectories as given by the final values of mIfs
#zipped = zip(range(nSamples-1), mIfs[:-1, -1])
#sorted_keys = sorted(zipped, key = lambda t: t[1])
#sorted_keys = np.array(sorted_keys)
# 95% confidence interval
#sorted_keys = sorted_keys[:-10]
#print(sorted_keys)
#print(sorted_keys[:, 0])
#mask = np.array(sorted_keys[:, 0])
#mask = mask.astype(int)
#activesFor = P_0*mIfs[mask]
confirmedFor = P_0*(mIfs+mRfs+mMfs)
deathsFor = P_0*mMfs
recoveredFor = P_0*mRfs
activesFor = P_0*mIfs
#for i in range(len(mIfs)):
for i in range(Nfcast):
  distribA = activesFor[:, i]
  distribC = confirmedFor[:, i]
  distribD = deathsFor[:, i]
  distribR = recoveredFor[:, i]
  histA, binEdgesA = np.histogram(distribA, bins=nBinsHist, density=True)
  histC, binEdgesC = np.histogram(distribC, bins=nBinsHist, density=True)
  histD, binEdgesD = np.histogram(distribD, bins=nBinsHist, density=True)
  histR, binEdgesR = np.histogram(distribR, bins=nBinsHist, density=True)
  t = dtime[i].date()
  histAMax = histA.max()
  histCMax = histC.max()
  histDMax = histD.max()
  histRMax = histR.max()
  for j in range(len(histA)-1):
    indCor = int(nBinsHist*histC[j]/histCMax)
    X = [t, t]
    Y = [binEdgesC[j], binEdgesC[j+1]]
    axs[0,0].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
    indCor = int(nBinsHist*histD[j]/histDMax)
    Y = [binEdgesD[j], binEdgesD[j+1]]
    axs[1,0].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
    indCor = int(nBinsHist*histR[j]/histRMax)
    Y = [binEdgesR[j], binEdgesR[j+1]]
    axs[0,1].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
    indCor = int(nBinsHist*histA[j]/histAMax)
    Y = [binEdgesA[j], binEdgesA[j+1]]
    axs[1,1].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
legenda = u"Confirmed cases in %s" % (country)
axs[0,0].plot(dates, confirmed, 'k-o', ms=2,  label= legenda)
legenda = u"Recovered cases in %s" % (country)
axs[0,1].plot(dates, recovered, 'k-o', ms=2,  label= legenda)
legenda = u"Death cases in %s" % (country)
axs[1,0].plot(dates, deaths, 'k-o', ms=2,  label= legenda)
legenda = u"Active cases in %s" % (country)
axs[1,1].plot(dates, activeCases, 'k-o', ms=3,  label= legenda)

# A. confirmed cases
C= P_0*(I+R+M)
ax= axs[0,0]
ax.plot(dates[:-1], C[::day], 'r-', label=u"$P_0(I(t)+R(t)+M(t))$, theoretical model")
ax.legend()

# B. recovered cases
ax= axs[0,1]
ax.plot(dates[:-1], P_0*R[::day], 'r-', label=u"$P_0R(t)$, theoretical model")
ax.legend()

# C. death cases
ax= axs[1,0]
ax.plot(dates[:-1], P_0*M[::day], 'r-', label=u"$P_0M(t)$, theoretical model")
ax.legend()

# D. active cases
ax= axs[1,1]
ax.plot(dates[:-1], P_0*I[::day], 'r-', label=u"$P_0I(t)$, theoretical model")
ax.legend()

# initial conditions
C_0 = confirmed[shift]
D_0 = deaths[shift]
A_0 = activeCases[shift]
R_0 = C_0-A_0-D_0
#A_0 = C_0-R_0-D_0
print('C_0', C_0, 'A_0', A_0, 'R_0', R_0, 'D_0', D_0)
# initial values  
Sf[0] = 1.0-C_0/P_0
If[0] = A_0/P_0
Rf[0] = R_0/P_0
Mf[0] = D_0/P_0
ka_avg /= nRuns
for i in np.arange(1, Nfcast):
  ka_i = ka_avg[i-1]
  factor = -1.0/tau+ka_i*Sf[i-1]
  If[i] = If[i-1]*np.exp(factor)
  Delta_I = If[i]-If[i-1]
  Sf[i] = Sf[i-1]*np.exp(-ka_i*Delta_I/factor)
  Rf[i] = Rf[i-1]+rec_forecast[i-1]*Delta_I/factor
  Mf[i] = Mf[i-1]+let_forecast[i-1]*Delta_I/factor
Cf_avg = P_0*np.sum((mIfs+mRfs+mMfs), axis=0)/nSamples
Mf_avg = P_0*np.sum(mMfs, axis=0)/nSamples
Af_avg = P_0*np.sum(mIfs, axis=0)/nSamples
Rf_avg = P_0*np.sum(mRfs, axis=0)/nSamples
#print(I_avg, I_avg.shape)
axs[0,0].plot(dtime, Cf_avg, 'y-', linewidth=2, label='average forecast')
axs[0,1].plot(dtime, Rf_avg, 'y-', linewidth=2, label='average forecast')
axs[1,0].plot(dtime, Mf_avg, 'y-', linewidth=2, label='average forecast')
axs[1,1].plot(dtime, Af_avg, 'y-', linewidth=2, label='average forecast')
axs[0,0].legend()
axs[0,1].legend()
axs[1,0].legend()
axs[1,1].legend()

figure = 'conRecDeaActive_%s.pdf'% country
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
fig2.tight_layout(pad=2.0)
fig2.savefig(figure, bbox_inches='tight')
flog.write(figure)

# R squared (coefficient of determination)
#y_bar = np.sum(confirmed)/len(confirmed)
#SS_tot = np.sum((confirmed[:-1]-y_bar)**2)
#SS_res = np.sum((confirmed[:-1]-C[::day])**2)
#R2= 1-SS_res/SS_tot
  #print('R2 %g' % R2)
### R squared calculation
# R squared (coefficient of determination)
def Rsquared(y,f):
  y_bar = np.sum(y)/len(y)
  SS_tot = np.sum((y-y_bar)**2)
  #SS_tot = np.sum((confirmed[:-1]-y_bar)**2)
  SS_res = np.sum((y-f)**2)
  R2= 1-SS_res/SS_tot
  return R2
fig3, ax = plt.subplots(2, figsize=(8, 12), sharex=True)
fig3.align_ylabels(axs)
fig3.subplots_adjust(hspace=0.26)

## running window R2 calculation
# confirmed cases
Nwindow = 120
r2Conf = np.zeros(len(confirmed[:-Nwindow]))
for i in np.arange(len(confirmed[:-Nwindow])):
  y = confirmed[i:i+Nwindow]
  f = C[i*day:(i+Nwindow)*day:day]
  r2Conf[i] = Rsquared(y, f)

legenda = u"running %d-day window $R^2$ of confirmed cases in %s" % (Nwindow,
    country)
ax[0].plot(dates[Nwindow:], r2Conf, 'r-', label=legenda)
#ax[0].set_title('Running window $R^2$, %d-day window' % Nwindow)
ax[0].set_ylabel('$R^2$')
ax[0].set_ylim(0, 1.0)
r2D = np.zeros(len(deaths[:-Nwindow]))
# running window R2 calculation
# death cases
for i in np.arange(len(confirmed[:-Nwindow])):
  y = deaths[i:i+Nwindow]
  f = P_0*M[i*day:(i+Nwindow)*day:day]
  r2D[i] = Rsquared(y, f)

legenda = u"running %d-day window $R^2$ of death cases in %s" % (Nwindow,
    country)
ax[0].plot(dates[Nwindow:], r2D, 'k--', label=legenda)
#ax[1].set_title('Running window $R^2$, %d-day window' % Nwindow)
ax[0].set_ylabel('$R^2$')
ax[0].legend(fontsize=10)
ax[0].set_xlim(dates.iloc[Nwindow], dates.iloc[-1])

## Cumulative R2 calculation
Nwindow = 120
r2C = np.zeros(len(confirmed[:-Nwindow]))
for ind, i in enumerate(np.arange(Nwindow, len(confirmed[:-1]))):
  y = confirmed[:i]
  f = C[:i*day:day]
  r2C[ind] = Rsquared(y, f)
legenda = u"Cumulative $R^2$ of confirmed cases in %s" % (country)
ax[1].plot(dates[Nwindow:-1], r2C[:-1], 'r-', label=legenda)
r2D = np.zeros(len(deaths[:-Nwindow]))
for ind, i in enumerate(np.arange(Nwindow, len(confirmed[:-1]))):
  y = deaths[:i]
  f = P_0*M[:i*day:day]
  r2D[ind] = Rsquared(y, f)
legenda = u"Cumulative $R^2$ of death cases in %s" % (country)
ax[1].plot(dates[Nwindow:-1],r2D[:-1], 'k--', label=legenda)
ax[1].set_ylim(0, 1.0)
ax[1].legend()
ax[1].set_xlim(dates.iloc[Nwindow], dates.iloc[-1])
plt.gcf().autofmt_xdate()
set_xaxis(ax[1], 28, 7)

R2 = Rsquared(confirmed[:-1], C[::day])
flog.write('confirmed cases R2 %g' % R2)
R2 = Rsquared(deaths[:-1], P_0*M[::day])
flog.write('death cases R2 %g' % R2)

figure = 'R2confDea%s.pdf'% (country)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
fig3.savefig(figure, bbox_inches='tight')
flog.write('%s\n' % figure)

fig4, ax = plt.subplots(2, figsize=(8, 12), sharex=True)
fig4.align_ylabels(axs)
fig4.subplots_adjust(hspace=0.26)

### NRMSE calculation
# root mean square error
def NRMSE (y, f):
  rmse = np.sqrt(np.sum((y-f)**2)/len(y))
  y_bar = np.sum(y)/len(y)
  nrmse = rmse/y_bar
  return nrmse
## running window NRMSE calculation
# confirmed cases
Nwindow = 90
nrmseC = np.zeros(len(confirmed[:-Nwindow]))
for i in np.arange(len(confirmed[:-Nwindow])):
  y = confirmed[i:i+Nwindow]
  f = C[i*day:(i+Nwindow)*day:day]
  nrmseC[i] = NRMSE(y, f)

legenda = u"running %d-day window $NRMSE$ of confirmed cases in %s" % (Nwindow, country)
ax[0].plot(dates[Nwindow:], nrmseC, 'r-', label=legenda)
#ax[0].set_title('Running window $R^2$, %d-day window' % Nwindow)
ax[0].set_ylabel('$R^2$')
#ax[0].set_ylim(0, 1.0)
## running window NRMSE calculation
# death cases
nrsmeD = np.zeros(len(deaths[:-Nwindow]))
for i in np.arange(len(confirmed[:-Nwindow])):
  y = deaths[i:i+Nwindow]
  f = P_0*M[i*day:(i+Nwindow)*day:day]
  nrsmeD[i] = NRMSE(y, f)

legenda = u"running %d-day window $NRMSE$ of death cases in %s" % (Nwindow, country)
ax[0].plot(dates[Nwindow:], nrsmeD, 'k--', label=legenda)
ax[0].set_ylabel('$NRSME$')
ax[0].legend(fontsize=10)
ax[0].set_xlim(dates.iloc[Nwindow], dates.iloc[-1])
ax[0].set_ylim(0, 0.25)

# Cumulative nrmse calculation
Nwindow = 90
nrmseC = np.zeros(len(confirmed[:-Nwindow]))
for ind, i in enumerate(np.arange(Nwindow, len(confirmed[:-1]))):
  y = confirmed[:i]
  f = C[:i*day:day]
  nrmseC[ind] = NRMSE(y, f)
legenda = u"Cumulative $NRMSE$ of confirmed cases in %s" % (country)
ax[1].plot(dates[Nwindow:-1], nrmseC[:-1], 'r-', label=legenda)
nrmseD = np.zeros(len(deaths[:-Nwindow]))
for ind, i in enumerate(np.arange(Nwindow, len(confirmed[:-1]))):
  y = deaths[:i]
  f = P_0*M[:i*day:day]
  nrmseD[ind] = NRMSE(y, f)
legenda = u"Cumulative $NRMSE$ of death cases in %s" % (country)
ax[1].plot(dates[Nwindow:-1], nrmseD[:-1], 'k--', label=legenda)
ax[0].set_ylabel('$NRMSE$')
ax[1].set_ylabel('$NRMSE$')
ax[0].legend(fontsize=10)
ax[1].legend()
ax[1].set_xlim(dates.iloc[Nwindow], dates.iloc[-1])
ax[1].set_ylim(0, 0.25)
plt.gcf().autofmt_xdate()
set_xaxis(ax[1], 28, 7)

#errConf = np.sqrt(np.sum((confirmed[:-1]-C[::day])**2)/len(confirmed))
#errConf /= P_0
errConf = NRMSE(confirmed[:-1], C[::day])
flog.write('Error in confirmed cases %g' % errConf)
errDea = NRMSE(deaths[:-1], P_0*M[::day])
flog.write('Error in death cases %g' % errDea)

figure = 'RMSEconfDea%s.pdf'% (country)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
fig4.savefig(figure, bbox_inches='tight')
flog.write('%s\n' % figure)
flog.close()
print(logFile)

plt.show()
