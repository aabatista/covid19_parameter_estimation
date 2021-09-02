# covid19_parameter_estimation
Python code of the arxiv preprint https://arxiv.org/abs/2011.06515
entitled **"An epidemiological compartmental model with automated parameter
estimation and forecasting of the spread of COVID-19 with analysis of data
from Germany and Brazil"**

The code of this project is distributed in following files:
* sensibilityAnalysis.py (figure 1)
* activeWorld.py (figures 2, 3, 4, and 5)
* covidBR.py (the remaining figures)

# **How to run the code**

The command line to run activeWorld.py from a Linux terminal is

$ python3 activeWorld.py config.txt

The file config.txt has the following lines in the case of Germany:

country,Germany

delay,13

offset,41

cutoff,0.41

The command line to run activeWorld.py from a Linux terminal is

> $ python3 activeWorld.py config.txt

The file config.txt has the following lines in the case of Germany:

country,Germany  
delay,13  
offset,41  
cutoff,0.41  

The file config.txt has the following lines in the case of Brazil:

country,Brazil  
delay,14  
offset,22  
cutoff,0.43  

The command line to run covidBR.py from a Linux terminal is

> $ python3 covidBR.py config.txt

In the case of Paraíba, the file config.txt has the following lines:

state,PB  
city,\_  
delay,14  
offset,0  
cutoff,0.38  

The file configCampinaGrandePB.txt has the following lines in the case of Campina Grande:

state,PB  
city,Campina Grande  
delay,14  
offset,21  
cutoff,0.35  
