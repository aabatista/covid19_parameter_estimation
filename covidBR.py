# -*- coding: utf-8 -*-
# Author: Adriano A. Batista, 2020
# covidBR.py -- dados epidemiológicos e comparação com previsões teóricas
# do Brasil
'''
 command line:
 python3 covidBR.py config.txt

 Input: 
   csv data files:
   caso.csv
 o arquivo caso.csv.gz pode ser obtido do site:
   https://brasil.io/dataset/covid19/files/
 pode ser baixado com o comando:
   wget https://data.brasil.io/dataset/covid19/caso.csv.gz
 config.txt file example:
  state,PB
  city,Campina Grande
  delay,14
  offset,21
  cutoff,0.3653696716410195
 Main Output
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
from datetime import datetime, timedelta
from collections import OrderedDict
import unidecode
plt.rcParams['axes.grid'] = True

def read_setup(inSetupFile):
  # open SIRM model parameter data 
  # Format of setup file
  # parameter      value
  '''
  state               nome do estado
  city                nome da cidade
  delay               atraso
  offset              outro atraso
  cutoff              valor máximo de k(t)
  '''
  global state, city, delay, offset, cutoff
  # Read setup file
  lines = np.genfromtxt(inSetupFile, dtype = 'str', delimiter=',')
  # Crie o dicionario
  setup_vars = OrderedDict(lines)
  #print setup_vars.items()
  for var in lines:
    if var[0] == 'state' or var[0] == 'city':
      setup_vars[var[0]] = var[1]
    else:
      if var[1].isdigit():
        setup_vars[var[0]] = int(var[1])
      else:
        setup_vars[var[0]] = float(var[1])
  # Transform the keys into variables
  globals().update(setup_vars)
read_setup(sys.argv[1])

n_avg=delay 
nu = 3.7844e-05 # taxa de natalidade diária (2018, IBGE) 
mu = 1.6918e-05 # taxa de mortalidade diária (2018, IBGE) 
day = 24 # dia em horas
dt = 1.0/day  # em unidade de dia
tau = (1+n_avg*day)*dt # tempo médio de infecção em dias
print('tau', tau)
estados = {
    'AC':'Acre', 'AL':'Alagoas', 'AM':'Amazônia', 'AP':'Amapá', 
    'BA':'Bahia', 'CE':'Ceará', 'DF':'Distrito Federal', 
    'ES':'Espírito Santo', 'GO':'Goiás', 'MA': 'Maranhão',
    'MG':'Minas Gerais', 'MS':'Mato Grosso do Sul', 'MT':'Mato Grosso',
    'PA':'Pará', 'PB':'Paraíba', 'PR':'Paraná', 'PE':'Pernambuco', 
    'PI':'Piauí', 'RJ':'Rio de Janeiro', 'RN':'Rio Grande do Norte', 
    'RO':'Rondônia', 'RR':'Roraima', 'RS':'Rio Grande do Sul', 
    'SC':'Santa Catarina', 'SE':'Sergipe', 'SP':'São Paulo', 
    'TO':'Tocantins'}
if city == '_':
  city_state = estados[state]
  city = ''
else:
  city_state = city
print (city_state)
# create log file
logFile = 'log'+city+state+'.txt'
# remove white spaces
logFile = logFile.replace(' ','')
f = open(logFile, 'w')
f.write('estado %s, cidade %s, atraso %d, offset %s, cutoff %g\n' 
    % (state,  city, delay, offset, cutoff))
##########################################################
# 1. Data processing
# File caso.csv obtained from the site https://data.brasil.io/dataset/covid19/_meta/list.html
# or download directly with the command:
# wget https://data.brasil.io/dataset/covid19/caso.csv.gz
arquivo_csv = 'caso.csv'
dados = pd.read_csv(arquivo_csv, sep =',', skiprows=0, header=0, parse_dates=['date'])
#dados = pd.read_csv(arquivo_csv, sep =',', skiprows=0, header=0, names=['date','state', 'city', 'place_type', 'confirmed', 'deaths','order_for_place','is_last', 'estimated_population_2019', 'estimated_population', 'city_ibge_code', 'confirmed_per_100k_inhabitants', 'death_rate'], parse_dates=['date'])
dados = dados[dados['state'].str.contains(state)]
if city == '':
  dados = dados[dados['city'].isnull()]
else:
  dados = dados[~dados.isin([np.nan]).any(1)]
  dados = dados[dados['city'].str.contains(city)]

dados = dados.sort_values('date', ascending=True)
confirmed = dados['confirmed'].to_numpy().astype(int)
deaths = dados['deaths'].to_numpy().astype(int)
ind0 = (confirmed!=0).argmax()
dates = pd.to_datetime(dados.date)
dates = dates[ind0:]
indFirstDeath = (deaths!=0).argmax()
f.write('indFirstDeath %d, First death on %s, number of deaths %d\n'
  % (indFirstDeath, dates.iloc[indFirstDeath], deaths[indFirstDeath]))
ind0 += offset
confirmed = confirmed[ind0:]
deaths = deaths[ind0:]
f.write('ind0 %d, offset %d\n' % (ind0, offset))
dates = dates[offset:]
P_0 = dados['estimated_population_2019'].to_numpy(dtype=int)[0]
f.write('População em 2019 %d\n' % (P_0))

n_max = len(confirmed)-1 # em dias
tt = np.arange(0.0, n_max, dt)

# 2. Statistical analysis
# Active cases delay estimate
activeDelay = confirmed[delay:]-confirmed[:-delay]
# Active cases statistical estimate
q = n_avg/(1+n_avg) # probability to remain sick after 1 day
newConf = np.diff(confirmed) # casos novos confirmados ao dia
newDead = np.diff(deaths) # casos novos confirmados ao dia
N = len(newConf)
active_stat = np.ones(N)
n_pow = np.arange(N)
n_pow_des = n_pow[::-1]
q_to_n = q**n_pow_des
for i in n_pow:
  active_stat[i] = np.sum(newConf[:i]*q_to_n[N-i:])

# 3. Parameter estimation
# a. Contagion rate function
contagionRate = []

N_as = len(active_stat)
for i in np.arange(0, N_as):
  if i+7<=N_as:
    XX = np.arange(i, i+7)
    ZZ = active_stat[i:i+7]
  else:
    XX = np.arange(i-7, i)
    ZZ = active_stat[i-7:i]
  if len(XX)==len(ZZ):
    slope, intercept, r_value, p_value, std_err = stats.linregress(XX, ZZ)
    I_avg = np.sum(ZZ)/7
    if I_avg==0:
      ka_t=0
    else:
      ka_t = slope/I_avg+1.0/tau
    if ka_t<0:
      ka_t = 0
    elif ka_t>cutoff:
      ka_t = cutoff
    contagionRate.append(ka_t)
print(len(dates[1:]), len(contagionRate))

##########################################################
# Contagion rate function
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
# Vectorize the contagion rate
kappa_t_vec = []
for ind, t in enumerate(tt[::day]):
  ka = kappa(t)
  kappa_t_vec.append(ka)
kappa_t_vec = np.array(kappa_t_vec)

# Lethality probability
P_let = np.zeros(N)
for i in n_pow[:-1]:
    denom = np.sum(newConf[:i]*q_to_n[N-i:])
    if denom>0:
      s = newDead[i]/denom
    else:
      s=0
    P_let[i] = s/(1-q)
print('P_let')
f.write('first P_let avg %g\n' % (np.sum(P_let[:21])/21))
f.write('last P_let avg %g\n' % (np.sum(P_let[-21:])/21))

# week average
Pl_avg = np.convolve(P_let, np.ones((7,))/7, mode='same')
# Lethality and recovery rates
letRate = Pl_avg/tau
recRate = 1.0/tau-letRate

# Numerical model epidemic evolution (modified SIR model)
##########################################################
# Get lethality and recovery rates
def get_l_r(t):
  ind = int(1.0*t)
  N=len(Pl_avg)
  if ind>=N:
    ind = N-1
  la = Pl_avg[ind]/tau
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
R_0 = 0.0
D_0 = deaths[0]
A_0 = C_0-R_0-D_0
f.write('C_0 %g, A_0 %g, R_0 %g, D_0 %g\n' % (C_0, A_0, R_0, D_0))
yinit = [1.0-C_0/P_0, A_0/P_0, R_0/P_0, D_0/P_0] # initial values
y = scipy.integrate.odeint(derivs, yinit, tt)
S = y[:, 0]
I = y[:, 1]
R = y[:, 2]
M = y[:, 3]

# B. generate Rt
R0_t = tau*kappa_t_vec*S[::day]

# 5. Plot results
# Figure 1

# A. Contagion rate estimate
#fig1 = plt.figure('kappa', figsize=(6, 12))
fig1, axs = plt.subplots(3, figsize=(6, 12), sharex=True)
fig1.align_ylabels(axs)
#fig1.tight_layout()
fig1.subplots_adjust(hspace=0.26)
ax = axs[0]
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.set_ylabel('$\kappa(t)$')
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
ax.plot(dates[1:], kappa_t_vec, 'k-', ms=2)
delta = timedelta(days=14)
ax.set_xlim(dates.iloc[0], dates.iloc[-1]+delta)
ax.set_ylim(0., )
ax.text(-0.1, 1.1, 'A', transform=ax.transAxes, size=12, weight='bold')
ax.set_title('Contagion rate in %s' % city_state)

# B. Dynamic reproductive number R_0(t)
ax = axs[1]
ax.text(-0.1, 1.1, 'B', transform=ax.transAxes, size=12, weight='bold')
ax.set_title('Reproduction number')
#ax.set_xticks(dates[0::7])
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(AutoMinorLocator(0.5))
ax.plot(dates[1:], R0_t, 'k-')
ax.set_ylabel('$R_0(t)$')
ax.set_xlim(dates.iloc[0], dates.iloc[-1]+delta)
ax.set_ylim(0., )

# C. Lethality probability, P_lambda(t), and recovery probability, P_rho(t)
ax = axs[2]
#ax = fig1.add_subplot(313)
ax.set_title('Lethality and recovery rates')
ax.text(-0.1, 1.1, 'C', transform=ax.transAxes, size=12, weight='bold')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=14))
plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=2))
ax.plot(dates[1:], Pl_avg, '-', linewidth=1.5, label="Lethality probability (week average)")
ax.plot(dates[1:], 1-Pl_avg, '-.', linewidth=1.5,  label="Recovery probability (week average)")
ax.set_xlim(dates.iloc[0], dates.iloc[-1]+delta)
ax.set_ylim(0., 1.)
plt.gcf().autofmt_xdate()
plt.legend(fontsize=10)
figure = 'contagionLetRecRates%s_%s.pdf'% (city, state)

# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
fig1.savefig(figure, bbox_inches='tight')
f.write('%s\n' % figure)


# Figure 2: confirmed, deaths, and active
fig2, axs = plt.subplots(3, figsize=(8, 12), sharex=True)
fig2.align_ylabels(axs)
fig2.tight_layout()
#plt.subplots_adjust(hspace=0.15)

# A. confirmed cases
ax= axs[0] 
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
C= P_0*(I+R+M)
ax.plot(dates[1:], C[::day], 'r-', label=u"$P_0[I(t)+R(t)+M(t)]$, theoretical model")
ax.text(-0.1, 1.00, 'A', transform=ax.transAxes, size=12, weight='bold')
ax.set_ylabel(u'number of confirmed cases')
plt.gcf().autofmt_xdate()
ax.set_xlim(dates.iloc[0], dates.iloc[-1]+delta)
ax.legend()

# B. Death cases
ax=axs[1]
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
ax.plot(dates[1:], P_0*M[::day], 'r-', label=u"$P_0M(t)$, theoretical model")
ax.text(-0.1, 1.00, 'B', transform=ax.transAxes, size=12, weight='bold')
ax.legend()
ax.set_ylabel(u'number of death cases')
ax.set_xlim(dates.iloc[0], dates.iloc[-1]+delta)
plt.gcf().autofmt_xdate()

# C. Active cases
ax=axs[2]
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
ax.plot(dates[delay:], activeDelay, 'g+', ms=3,  label= 'estimate, delay %d days' % delay)
print(len(dates), len(I[::day]))
ax.plot(dates[:-1], P_0*I[::day], 'r-', label=u"$P_0I(t)$, theoretical model")
ax.text(-0.1, 1.00, 'C', transform=ax.transAxes, size=12, weight='bold')
ax.set_xlim(dates.iloc[0], dates.iloc[-1]+delta)
ax.set_ylabel(u'number of infected cases')
plt.gcf().autofmt_xdate()

# Analytical integration
Su = np.zeros(len(dates))
In = np.zeros(len(dates))
Re = np.zeros(len(dates))
Mo = np.zeros(len(dates))
# initial conditions
Su[0] = 1.0-C_0/P_0
In[0] = A_0/P_0
Re[0] = R_0/P_0
Mo[0] = D_0/P_0
#kappa = kappa_t_vec
for i in np.arange(1, len(dates)):
  ka_i = kappa(i-1)
  factor = -1.0/tau+ka_i*Su[i-1]
  In[i] = In[i-1]*np.exp(factor)
  Delta_I = In[i]-In[i-1]
  Su[i] = Su[i-1]*np.exp(-ka_i*Delta_I/factor)
  Re[i] = Re[i-1]+recRate[i-1]*Delta_I/factor
  Mo[i] = Mo[i-1]+letRate[i-1]*Delta_I/factor

fig, Axs = plt.subplots(3, figsize=(8, 12), sharex=True)
ax=Axs[0]
ax.plot(dates[:-1], P_0*R[::day], 'r-', label=u"$P_0R(t)$, theoretical model")
ax.plot(dates, P_0*Re, 'k-', label=u"$P_0R(t)$, automatic integration")
ax.legend()
ax=Axs[1]
ax.plot(dates[:-1], P_0*M[::day], 'r-', label=u"$P_0M(t)$, theoretical model")
ax.plot(dates, P_0*Mo, 'k-', label=u"$P_0M(t)$, automatic integration")
ax.legend()
ax=Axs[2]
ax.plot(dates[:-1], P_0*I[::day], 'r-', label=u"$P_0I(t)$, theoretical model")
ax.plot(dates, P_0*In, 'k-', label=u"$P_0I(t)$, automatic integration")
plt.gcf().autofmt_xdate()
ax.legend()
figure = 'metodosIntegComp%s_%s.pdf'% (city, state)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
fig.savefig(figure, bbox_inches='tight')
f.write('%s\n' % figure)

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
nSamples = int(nRuns/sampleSize)
shift = -14 # it has to be zero or negative, in days
nBinsHist = 40
colors = plt.cm.jet(np.linspace(0, 1, nBinsHist+1))
dia0= dates.iloc[shift]
dtime = pd.date_range(start = dia0, periods=Nfcast)
print("kappa_t", kappa_t_vec[-21+shift:shift])
dk = np.diff(kappa_t_vec[-21+shift:shift])
Nbins = int(len(dk)/2)
print('Nbins', Nbins)
# The state space
states = [0, 1]
# Transition probability matrix
# Markov chain part for contagion rate
W = transitionMatrix = getTransMatrix(dk)
print(transitionMatrix)
# Get entropy
entropyS = -np.sum(W*np.log(W))
f.write("Entropia S %g\n" % entropyS)
probPos, p_edges, probNeg, n_edges = getProbs(dk, Nbins)
dl = np.diff(letRate[-21+shift:shift])
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
figure = 'histograms%s_%s.pdf'% (city, state)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
plt.savefig(figure, bbox_inches='tight')
f.write('%s\n' % figure)

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
timeInter = np.arange(len(In)+shift-1, len(In)+shift-1+Nfcast)
ka_avg = np.zeros(Nfcast)
fig3, Axs = plt.subplots(2, figsize=(8, 8))
ax = Axs[0]
ax1 = Axs[1]
ax.plot(active_stat, 'gx', ms=3,  label= 'Statistical estimate')
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
  A_0 = active_stat[shift]
  R_0 = C_0-A_0-D_0
  # initial values  
  Sf[0] = 1.0-C_0/P_0
  If[0] = A_0/P_0
  Rf[0] = R_0/P_0
  Mf[0] = D_0/P_0
  #Sf[0] = Su[shift]
  #If[0] = In[shift]
  #Rf[0] = Re[shift]
  #Mf[0] = Mo[shift]
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
  #axs[1].plot(dtime, P_0*Mf, 'k-', alpha=0.25)
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
    #C_avg = P_0*(If+Rf+Mf)
    #axs[0].plot(dtime, C_avg, 'k-', alpha=0.25)
    #axs[1].plot(dtime, P_0*M_avg, 'k-', alpha=0.25)
    #print(indT, count, P_0*If[1], P_0*If[2], P_0*If[13])
    #print(indT, count, (P_0*mIfs[count, :]).astype(int))
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
activesFor = P_0*mIfs
#for i in range(len(mIfs)):
#  plt.plot(timeInter, P_0*mIfs[i, :], 'r-')
# create forecast csv file
date = str(dates.iloc[-1])
date = date.replace('00:00:00', '')
forecastFile = 'forecast'+city+state+date+'.csv'
agora = datetime.now()
print(agora)
# remove white spaces
forecastFile = forecastFile.replace(' ','')
ff = open(forecastFile, 'w')
ff.write("date,Confirmed_min,Confirmed_max,deaths_min,deaths_max,active_min,active_max\n")
for i in range(Nfcast):
  distribA = activesFor[:, i]
  distribC = confirmedFor[:, i]
  distribD = deathsFor[:, i]
  histA, binEdgesA = np.histogram(distribA, bins=nBinsHist, density=True)
  histC, binEdgesC = np.histogram(distribC, bins=nBinsHist, density=True)
  histD, binEdgesD = np.histogram(distribD, bins=nBinsHist, density=True)
  if i==1 or i==4 or i==8:
    ax1.hist(distribA, bins = binEdgesA, density=True, alpha=0.6, ec='k', label='%d'%i)
    ax1.axvline(distribA.mean(), ls='dashed', linewidth=1)
  X = timeInter[i]*np.ones(len(histC))
  colorline.colorline(ax, X, binEdgesA, histA, cmap=plt.get_cmap('jet'), linewidth=2)
  t = dtime[i].date()
  histAMax = histA.max()
  histCMax = histC.max()
  histDMax = histD.max()
  for j in range(len(histA)-1):
    indCor = int(nBinsHist*histC[j]/histCMax)
    X = [t, t]
    Y = [binEdgesC[j], binEdgesC[j+1]]
    axs[0].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
    indCor = int(nBinsHist*histD[j]/histDMax)
    Y = [binEdgesD[j], binEdgesD[j+1]]
    axs[1].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
    indCor = int(nBinsHist*histA[j]/histAMax)
    Y = [binEdgesA[j], binEdgesA[j+1]]
    axs[2].plot(X, Y, color=colors[indCor], linewidth=6, alpha=0.5)
  ff.write("%s, %d, %d, %d, %d, %d, %d\n" % (t, binEdgesC[0], binEdgesC[-1], binEdgesD[0], binEdgesD[-1], binEdgesA[0], binEdgesA[-1]))
ff.close()
print(forecastFile)
# Save last two weeks of epidemiological data
lastDaysDataFile = 'lastDays'+city+state+date+'.csv'
# remove white spaces
lastDaysDataFile = lastDaysDataFile.replace(' ','')
ldf = pd.DataFrame({"date": dates[-18:], "Confirmed": confirmed[-18:], "Deaths":
  deaths[-18:], "Active": active_stat[-18:].astype(int)})
ldf.to_csv(lastDaysDataFile, index=False)
print(lastDaysDataFile)
##############################################################################
legenda = u"Confirmed cases in %s" % (city_state)
axs[0].plot(dates, confirmed, 'k-o', ms=3,  label= legenda)
legenda = u"Death cases in %s" % (city_state)
axs[1].plot(dates, deaths, 'k-o', ms=3,  label= legenda)
axs[2].plot(dates[1:], active_stat, 'kx', ms=4,  label= 'Statistical estimate')
ax1.legend(loc='upper right')
ax1.set_title('Histograms')
# initial conditions
C_0 = confirmed[shift]
D_0 = deaths[shift]
A_0 = active_stat[shift]
R_0 = C_0-A_0-D_0
#A_0 = C_0-R_0-D_0
print('C_0', C_0, 'A_0', A_0, 'R_0', R_0, 'D_0', D_0)
# initial values  
Sf[0] = 1.0-C_0/P_0
If[0] = A_0/P_0
Rf[0] = R_0/P_0
Mf[0] = D_0/P_0
#Sf[0] = Su[-1]
#If[0] = In[-1]
#Rf[0] = Re[-1]
#Mf[0] = Mo[-1]
ka_avg /= nRuns
for i in np.arange(1, Nfcast):
  ka_i = ka_avg[i-1]
  factor = -1.0/tau+ka_i*Sf[i-1]
  If[i] = If[i-1]*np.exp(factor)
  Delta_I = If[i]-If[i-1]
  Sf[i] = Sf[i-1]*np.exp(-ka_i*Delta_I/factor)
  Rf[i] = Rf[i-1]+rec_forecast[i-1]*Delta_I/factor
  Mf[i] = Mf[i-1]+let_forecast[i-1]*Delta_I/factor
#plt.plot(timeInter, P_0*If, 'k-', label='average $\kappa(t)$ forecast')
Cf_avg = P_0*np.sum((mIfs+mRfs+mMfs), axis=0)/nSamples
Mf_avg = P_0*np.sum(mMfs, axis=0)/nSamples
Af_avg = P_0*np.sum(mIfs, axis=0)/nSamples
#print(I_avg, I_avg.shape)
ax.plot(timeInter, P_0*I_avg, 'k--', label='average forecast')
axs[0].plot(dtime, Cf_avg, 'y-', linewidth=2, label='average forecast')
axs[1].plot(dtime, Mf_avg, 'y-', linewidth=2, label='average forecast')
axs[2].plot(dtime, Af_avg, 'y-', linewidth=2, label='average forecast')
axs[0].legend()
axs[1].legend()
axs[2].legend()
ax.set_ylabel('Active cases')
ax.set_title('Modeling and forecasting active cases in %s %s' % (city, state))
#ax.set_title(dates.iloc[-1].strftime('%Y-%m-%d'))
ax.legend(loc='upper left')
figure = 'forecast%s_%s.pdf'% (city, state)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
plt.savefig(figure, bbox_inches='tight')
f.write('%s\n' % figure)

figure = 'confDeActive%s_%s.pdf'% (city, state)
# remove accents
figure = unidecode.unidecode(figure)
# remove white spaces
figure = figure.replace(' ','')
fig2.savefig(figure, bbox_inches='tight')
f.write('%s\n' % figure)
f.close()
print(logFile)
plt.show()
