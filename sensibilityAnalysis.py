# -*- coding: utf-8 -*-
#!/usr/bin/python
# sensibilityAnalysis.py -- Check the dependence on the contagion rate and
# lethality rate
import matplotlib.pyplot as plt
import scipy.integrate
import numpy as np

# parameters
P_0 = 2.1105e+8 # Initial population
#P_0 = 409730.0 # população inicial
k = kappa=0.34  # taxa de contágio máxima
mu = 1.6918e-05 # taxa de mortalidade diária (2018, IBGE) 
nu = 3.7844e-05 # taxa de natalidade diária (2018, IBGE) 
day = 24 # day in hours
dt = 1.0/day  # integration time-step
tau = (1+14*day)*dt # tempo médio de infecção em dias
n_med = tau*day # número médio de minutos de infecção
print("tau=", tau)
##########################################################
def derivs (x, t): # return derivatives of the array x
   S = x[0] # susceptible
   I = x[1] # infected
   R = x[2] # Recovered
   M = x[3] # Deaths
   T = S+I+R
   dSdt = nu*T-mu*S-k*S*I
   dIdt = -mu*I+k*S*I-r*I-l*I
   dRdt = r*I-mu*R
   dMdt = l*I
   return [dSdt, dIdt, dRdt, dMdt]
##########################################################
n_max = 365.0 # em dias
tt = np.arange(0.0, n_max, dt)
day = 24 # day in hours
plt.figure()
N = 128
ds = np.linspace(1.0/N, 1, N, endpoint=True)
totalM= 0*ds
lines = np.array(['-', '--', '-.', ':','-'])
colors= np.array(['b', 'r', 'k', 'g', 'orange'])
lws = np.array([1.0, 1.5, 1.5, 1.5, 2.0])
for indP, P_l in enumerate([0.09, 0.05, 0.03, 0.01, 0.001]):
  print(indP, P_l)
  P_r = 1-P_l # recovery probability
  l = P_l/tau # lethality rate
  r =  P_r/tau # recovery rate
  print('lambda=', l)
  print('r=', r)
  for ind, d in enumerate(ds):
    k= d*kappa
    yinit = [1.0, 1.0/P_0, 0.0, 0.0] # initial values
    y = scipy.integrate.odeint(derivs, yinit, tt)
    M = y[:, 3]
    Mpdia = M[::day] # M(T), M(2T), M(3T), ... T=day
    Mpdia = np.diff(Mpdia) # Deaths by day
    # Total number of deaths until 1 year after outbreak of the pandemic
    totalM[ind] = np.sum(Mpdia)*P_0 
  plt.plot(ds*kappa, totalM, linestyle=lines[indP], color=colors[indP], lw =
      lws[indP], label='$P_\lambda=%g$' % P_l)
plt.suptitle(u'Deaths due to COVID-19 after one year of outbreak')
plt.legend()
plt.ylim(1.0, 3.0e+7)
plt.title(r'$P_0=%d$, $\tau=%.2g$' % (P_0, tau), fontsize=10)
plt.yscale('log')
plt.gca().yaxis.set_ticks_position('both')
plt.xlim(0, k)
plt.grid()
plt.xlabel('$\kappa$')
plt.ylabel('Deaths')
figure = 'sensibilityTest.pdf'
plt.savefig(figure, bbox_inches='tight')
print(figure)
plt.show()


