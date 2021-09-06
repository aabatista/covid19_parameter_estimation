# -*- coding: utf-8 -*-
#!/usr/bin/python
# activeWorld.py -- dados epidemiológicos de casos ativos de qqr país do mundo
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

def read_setup(inSetupFile):
  # open SIRM model parameter data 
  # Format of setup file
  # parameter      value
  '''
  country             nome do país
  delay               atraso
  offset              outro atraso
  cutoff              valor máximo de k(t)
  '''
  global country, delay, offset, cutoff
  # Read setup file
  lines = np.genfromtxt(inSetupFile, dtype = 'str', delimiter=',')
  # Crie o dicionario
  setup_vars = OrderedDict(lines)
  # Save setup file
  outSetupFile = 'config'+setup_vars['country']+'.txt'
  print(outSetupFile)
  f = open(outSetupFile, 'w')
  #print setup_vars.items()
  for key,val in setup_vars.items():
    line = key+','+val+'\n'
    f.write(line)
  f.close()	
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
print('country', country, 'atraso', delay, 'offset', offset, 'cutoff', cutoff)

n_avg=delay
day = 24 # dia em horas
dt = 1.0/day  # passo de integração
tau = (1+n_avg*day)*dt # tempo médio de infecção em dias
nu = mu= 0.0
##########################################################
csvFile = ['time_series_covid19_confirmed_global_narrow.csv',
    'time_series_covid19_deaths_global_narrow.csv',
    'time_series_covid19_recovered_global_narrow.csv']
# Comparação com dados reais 
dados = pd.read_csv(csvFile[0], sep =',', skiprows=0, header=0, names=["ProvinceState", "Country/Region", "Lat", "Long", "Date", "Value", "ISO 3166-1 Alpha 3-Codes", "Region Code", "Sub-region Code", "Intermediate Region Code"])
dados = dados[dados['Country/Region'].str.contains(country)]
dados = dados.sort_values('Date', ascending=True)
confirmed = dados.Value.to_numpy().astype(int)
ind0 = (confirmed!=0).argmax()
ind0 += offset
confirmed = confirmed[ind0:]
print('ind0', ind0, 'offset', offset)
dates = pd.to_datetime(dados.Date)
dates = dates[ind0:]
dados = pd.read_csv(csvFile[1], sep =',', skiprows=0, header=0, names=["ProvinceState", "Country/Region", "Lat", "Long", "Date", "Value", "ISO 3166-1 Alpha 3-Codes", "Region Code", "Sub-region Code", "Intermediate Region Code"])
dados = dados[dados['Country/Region'].str.contains(country)]
deaths = dados.Value.to_numpy().astype(int)
deaths = deaths[::-1]
deaths = deaths[ind0:]
indFirstDeath = (deaths!=0).argmax()
print(indFirstDeath, 'First death on', dates.iloc[indFirstDeath], deaths[indFirstDeath])
dados = pd.read_csv(csvFile[2], sep =',', skiprows=0, header=0, names=["ProvinceState", "Country/Region", "Lat", "Long", "Date", "Value", "ISO 3166-1 Alpha 3-Codes", "Region Code", "Sub-region Code", "Intermediate Region Code"])
dados = dados[dados['Country/Region'].str.contains(country)]
recovered = dados.Value.to_numpy().astype(int)
recovered = recovered[::-1]
recovered = recovered[ind0:]

dados = pd.read_csv('worldPopulation.csv', sep =',', skiprows=0, header=0,
    names=['Country_Code','Ranking','Country','Population'])
dados = dados[dados['Country'].str.contains(country)]
population = dados['Population'].iloc[0]
population = population.replace(' ','')
population = population.replace(',','')
print(dados)
P_0 = int(population)*1000
print(P_0)
n_max = len(confirmed)-1 # em dias
tt = np.arange(0.0, n_max, dt)
dia0= dates.iloc[0]
dia1= dates.iloc[-1]
dtime = pd.date_range(start = dia0, end = dia1, freq = '360min')

# plot confirmed
fig1, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
ax= axs[0,0]
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
legenda = u"Confirmed cases in %s" % (country)
ax.plot(dates, confirmed, 'b-o', ms=3,  label= legenda)
ax.set_title(u'Confirmed cases')

# plot recovered
ax= axs[0,1]
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
legenda = u"Recovered cases in %s" % (country)
ax.plot(dates, recovered, 'b-o', ms=3,  label= legenda)
ax.set_title('Recovered cases')

# plot deaths
ax= axs[1,0]
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
legenda = u"Death cases in %s" % (country)
ax.plot(dates, deaths, 'b-o', ms=3,  label= legenda)
ax.set_title(u'Death cases')

# plot active
ax= axs[1,1]
activeCases = confirmed-deaths-recovered
#print(activeCases)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=2))
legenda = u"Active cases in %s" % (country)
ax.plot(dates, activeCases, 'b-o', ms=3,  label= legenda)
activeTest = confirmed[delay:]-confirmed[:-delay]
ax.plot(dates[delay:], activeTest, 'k+', ms=3,  label= 'estimate, delay %d days' % delay)
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
plt.plot(dates[1:], active_stat, 'gx', ms=3,  label= 'Statistical estimate')
ax.set_title(u'Active cases')
plt.gcf().autofmt_xdate()

# Contagion rate function
fig=plt.figure('kappa', figsize=(6, 12))
plt.subplots_adjust(hspace=0.25)
ax1 = fig.add_subplot(311)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
plt.gcf().autofmt_xdate()
#plt.xlabel('time (in days)')
ax1.set_ylabel('$\kappa(t)$')
contagionRate = []
#active_stat = (activeCases[1:]+active_stat)/2
#for i in np.arange(0, len(active_stat), 7):
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
    #print(i, ka_t)
    if ka_t<0:
      ka_t = 0
    elif ka_t>cutoff:
      ka_t = cutoff
    contagionRate.append(ka_t)
#print(contagionRate)
ax1.set_title('Contagion rate in %s' % country)
print(len(dates[1:-7]), len(contagionRate))
#ax1.bar(dates[1:], contagionRate, color='red', width=1.0, alpha=0.8)
# spline smoothing
X = np.arange(len(contagionRate))
sp = UnivariateSpline(X, contagionRate)
sp.set_smoothing_factor(0.0)
xs = np.linspace(X[0], X[-1], len(contagionRate))
ys = sp(xs)
print(len(dates[1:len(ys)+1]), len(ys))
#ax1.plot(dates[1:len(ys)+1], ys, '-', label='spline-smoothed contagion rate')
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
#Plot contagion rate
kappa_t_vec = []
for ind, t in enumerate(tt[::day]):
  ka = kappa(t)
  kappa_t_vec.append(ka)
kappa_t_vec = np.array(kappa_t_vec)
ax1.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax1.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax1.plot(dates[:-1], kappa_t_vec, 'k-', ms=2)
ax1.set_xlim(dtime[0], dtime[-1])
ax1.set_ylim(0., )
ax1.text(-0.1, 1.1, 'A', transform=ax1.transAxes, size=12, weight='bold')
plt.gcf().autofmt_xdate()
#plt.legend()

# Plot R_0(t)
ax2 = fig.add_subplot(312)
ax2.text(-0.1, 1.1, 'B', transform=ax2.transAxes, size=12, weight='bold')
ax2.set_title('Reproduction number')
#ax2.set_xticks(dates[0::7])
ax2.xaxis.set_major_locator(mdates.DayLocator(interval=14))
ax2.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax2.yaxis.set_major_locator(MultipleLocator(1))
ax2.yaxis.set_minor_locator(AutoMinorLocator(0.5))
R0_t = tau*kappa_t_vec
ax2.plot(dates[:-1], R0_t, 'k-')
ax2.set_ylabel('$R_0(t)$')
ax2.set_xlim(dtime[0], dtime[-1])
ax2.set_ylim(0., )
#ax2.annotate('(b)', xy=(.025, .975), xycoords='figure fraction',
#              horizontalalignment='left', verticalalignment='top', fontsize=12)
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
#print(len(P_let), len(dates), P_let)

# Plot lambda(t) and rho(t)
ax3 = fig.add_subplot(313)
ax3.set_title('Lethality and recovery rates')
fig.align_ylabels([ax1, ax2, ax3])
ax3.text(-0.1, 1.1, 'C', transform=ax3.transAxes, size=12, weight='bold')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=14))
plt.gca().xaxis.set_minor_locator(mdates.DayLocator(interval=1))
la = P_let/tau
rho = 1.0/tau-la
ax3.plot(dates[1:], P_let, 'k-', linewidth=2.0, label="Lethality probability estimate")
P_rho=1-P_let
ax3.plot(dates[1:], P_rho, 'b-.', linewidth=1.5,  label="Recovery probability estimate")
#ax3.plot(dates[1:], la, 'k-', linewidth=2.0, label="Lethality rate estimate")
#ax3.plot(dates[1:], rho, 'b-.', linewidth=1.5,  label="Recovery rate estimate")
# week average
#Pl_avg = np.convolve(P_let, np.ones((7,))/7, mode='same')
#la = Pl_avg/tau
#ax3.plot(dates[1:], la, '-.', linewidth=1.5,  label="Lethality rate week average")
#rho = 1.0/tau-la
#ax3.plot(dates[1:], rho, '-.', linewidth=1.5,  label="Recovery rate week average")
# spline smoothing
#X = np.arange(len(P_let))
#sp = UnivariateSpline(X, P_let/tau)
#sp.set_smoothing_factor(0.05)
#xs = np.linspace(X[0], X[-1], len(X))
#ys = sp(xs)
#print(len(dates[1:len(ys)+1]), len(ys))
##plt.plot(dates[1:len(ys)+1], ys, '-', label='spline-smoothed lethality rate')
ax3.set_xlim(dtime[0], dtime[-1])
ax3.set_ylim(0., 1.)
#ax3.set_ylim(0., 0.1)
plt.gcf().autofmt_xdate()
plt.legend(fontsize=10)
#plt.title(u'Lethality and recovery rates in %s' % country)
plt.gcf().autofmt_xdate()
#plt.tight_layout()
figura = 'contagionLetRecRates_%s.pdf'% country
# remove accents
figura = unidecode.unidecode(figura)
# remove white spaces
figura = figura.replace(' ','')
plt.savefig(figura, bbox_inches='tight')
print(figura)

# Numerical model epidemic evolution (modified SIR model)
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
print('C_0', C_0, 'A_0', A_0, 'R_0', R_0, 'D_0', D_0)
yinit = [1.0-C_0/P_0, A_0/P_0, R_0/P_0, D_0/P_0] # initial values
y = scipy.integrate.odeint(derivs, yinit, tt)
S = y[:, 0]
I = y[:, 1]
R = y[:, 2]
M = y[:, 3]

# plot confirmed cases
ax= axs[0,0]
C= P_0*(I+R+M)
ax.plot(dates[:-1], C[::day], 'r-', label=u"$P_0(I(t)+R(t)+M(t))$, theoretical model")
ax.text(-0.05, 1.02, 'A', transform=ax.transAxes, size=12, weight='bold')
ax.legend()

#Plot Recovered cases
ax= axs[0,1]
ax.plot(dates[:-1], P_0*R[::day], 'r-', label=u"$P_0R(t)$, theoretical model")
ax.text(-0.05, 1.02, 'B', transform=ax.transAxes, size=12, weight='bold')
ax.legend()

#Plot death cases
ax= axs[1,0]
ax.plot(dates[:-1], P_0*M[::day], 'r-', label=u"$P_0M(t)$, theoretical model")
ax.text(-0.05, 1.02, 'C', transform=ax.transAxes, size=12, weight='bold')
ax.legend()

# plot active cases
ax= axs[1,1]
print(len(dates), len(I[::day]))
ax.plot(dates[:-1], P_0*I[::day], 'r-', label=u"$P_0I(t)$, theoretical model")
ax.text(-0.05, 1.02, 'D', transform=ax.transAxes, size=12, weight='bold')
plt.gcf().autofmt_xdate()
ax.legend()
figura = 'conRecDeaActive_%s.pdf'% country
# remove accents
figura = unidecode.unidecode(figura)
# remove white spaces
figura = figura.replace(' ','')
fig1.tight_layout(pad=2.0)
fig1.savefig(figura, bbox_inches='tight')
print(figura)

plt.show()
