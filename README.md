# covid19_parameter_estimation
Python code of the arxiv preprint https://arxiv.org/abs/2011.06515
entitled **"An epidemiological compartmental model with automated parameter
estimation and forecasting of the spread of COVID-19 with analysis of data
from Germany and Brazil"**

The code of this project is distributed in the following files:
* sensibilityAnalysis.py (figure 2)
* activeWorld.py (figures 3, 4, 5, and 6)
* covidBR.py (the remaining figures)

![image](flowchartAW.png )

<div align="center">

*Flow chart of the data processing, analysing, plotting in activeWorld.py.
covidBR.py has a very similar flow chart*

</div>

## **How to run the code**

The command line to run sensibilityAnalysis.py from a Linux terminal is simply

> **$ python3 sensibilityAnalysis.py**

The command line to run activeWorld.py from a Linux terminal is

> **$ python3 activeWorld.py configCountry.txt**

where Country is the alpha-3 ISO 3166 country code. 
The file configDEU.txt has the following lines in the case of Germany:

country,Germany  
delay,14  
offset,41  
cutoff,0.43  
nu,2.591e-5  
mu,3.169e-5  

The parameter **delay** is tau in the paper.
The parameter **offset** indicates the initial datapoint to be used by the
simulation from the date of the first recorded case.
The parameter **cutoff** is the upper bound of the contagion rate. 
**nu** is the daily birth rate.
**mu** is the daily death rate.
The file configBRA.txt has the following lines in the case of Brazil:

country,Brazil  
delay,14  
offset,23  
cutoff,0.29  
nu,3.7844e-5  
mu,1.6918e-5 

The command line to run covidBR.py from a Linux terminal is

> **$ python3 covidBR.py configCityState.txt**

In the case of Paraíba, the file configPB.txt has the following lines:

state,PB  
city,\_  
delay,14  
offset,19
cutoff,0.26
nu,4.0139e-5
mu,1.7956e-5

The file configCampinaGrandePB.txt has the following lines in the case of Campina Grande:

state,PB  
city,Campina Grande  
delay,14  
offset,21  
cutoff,0.33  
nu,4.2943e-5  
mu,1.9616e-5  

## DATA

All data in the manuscript are contained in the following databases. The data
from Germany and Brazil, was obtained from the site
https://data.humdata.org/dataset/novel-coronavirus-2019-ncov-cases which
contains the data compiled by the Johns Hopkins University Center for Systems
Science and Engineering (JHU CSSE). The specific CSV files for the confirmed,
recovered, and deceased are:

* time_series_covid19_confirmed_global_narrow.csv
* time_series_covid19_deaths_global_narrow.csv
* time_series_covid19_recovered_global_narrow.csv

All these three csv files have to be in the same directory of activeWorld.py
file. The header line of all these three csv files is:

Province/State,Country/Region,Lat,Long,Date,Value,ISO 3166-1 Alpha 3-Codes,
Region Code,Sub-region Code,Intermediate Region Code

The data from Brazilian states and cities are obtained from the url: 

https://data.brasil.io/dataset/covid19/caso.csv.gz

The file *caso.csv* has to be in the same directory as covidBR.py
The header line of *caso.csv* is  
date,state,city,place_type,confirmed,deaths,order_for_place,is_last,
estimated_population_2019,city_ibge_code,confirmed_per_100k_inhabitants,
death_rate
