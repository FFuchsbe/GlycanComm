# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 10:14:00 2020

@author: Fuchsberger
"""


import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt     
import os
from FlowCytometryTools import FCMeasurement
from scipy.stats.mstats import gmean 
from statsmodels import robust

dir_path = os.path.dirname(os.path.realpath(__file__))   # set the working directory to the location of the Python script file. 

folder_name = "06092019 comp/TNF1"
file_names = os.listdir(folder_name)  #folder including fcs files
elongated_path_name = dir_path +"\\" +folder_name



sample = list(range(0,len(file_names)))

for i in range(0, len(file_names)):
    sample[i] = FCMeasurement(ID='Test Sample', datafile=r''+elongated_path_name + "\\" + file_names[i])



low_percentile = 0.001
high_percentile = 0.999



input_channel = "BL1-A"
output_channel = "VL2-A"  
#adjustable channel
total_io = pd.DataFrame()
for i in range(0, len(file_names)):
    temp_type = i
    sample_temp = sample[i].data
    temp_input = sample_temp[input_channel]
    temp_output = sample_temp[output_channel]
    temp_io = pd.concat([temp_input, temp_output], axis = 1)
    temp_io["Inputs"] = temp_type
    temp_io = temp_io[temp_io[input_channel]>temp_io[input_channel].quantile(low_percentile)]
    temp_io = temp_io[temp_io[input_channel]<temp_io[input_channel].quantile(high_percentile)]
    temp_io = temp_io[temp_io[output_channel]>temp_io[output_channel].quantile(low_percentile)]
    temp_io = temp_io[temp_io[output_channel]<temp_io[output_channel].quantile(high_percentile)]
    temp_io = temp_io[temp_io[output_channel]>1]
    temp_io = temp_io[temp_io[input_channel]>1]
    total_io = total_io.append(temp_io)


#this time with sklearn
from sklearn.linear_model import LinearRegression
# Create linear regression object.
lr_total_io = LinearRegression()
# Fit linear regression.
lr_total_io.fit(total_io[['BL1-A']], total_io['VL2-A'])
# Get the slope and intercept of the line best fit.
#y=k*x+d
print ("intercept is")
print(lr_total_io.intercept_)
print ("coefficient is")
print(lr_total_io.coef_)
#r_sq = lr_total_io.score(total_io['BL1-A'], total_io['VL2-A'])
#print ("R² is")
#print(r_sq)


#now for the transformation of catesian into polar coordinates

rho = np.sqrt(total_io['BL1-A']**2+total_io['VL2-A']**2)
theta = np.arctan2(total_io['VL2-A'],total_io['BL1-A'])
#print(rho)
#print(theta)  #both checks out

# then substract a from phi or actually arctan_value from t
# so we turn all data points by the angle a(alpha) to correct for mAmenrine correaltaion
import math
print ("slope :", lr_total_io.coef_)
arctan_Value = np.arctan(lr_total_io.coef_) 
print ("Inverse Tangent of slope :", arctan_Value)
alpha = ((arctan_Value)*(180/(math.pi)))       
print ("Alpha is:", alpha)
            
thetaMo = (theta-arctan_Value)


GFP_to_slope = rho * np.cos(thetaMo)
mAmetrine_to_slope = rho * np.sin(thetaMo)


#showing the data pre turning 
plt.scatter(total_io['BL1-A'],total_io['VL2-A'],)
plt.show()


#showing the data post turning 
plt.scatter(GFP_to_slope, mAmetrine_to_slope,)
plt.show()

#creating summary plot of both


plt.scatter(total_io['BL1-A'],total_io['VL2-A'], color='g', marker=".", alpha=0.1, )
plt.scatter(GFP_to_slope, mAmetrine_to_slope, color='b', marker=".",alpha=0.1 )
plt.xscale('log')
plt.yscale('log')
plt.xlim([1000, 200000])
plt.ylim([200, 100000])
plt.title('Scatter Plot')
plt.xlabel('GFP')
plt.ylabel('mAmetrine')
plt.legend()
plt.show()

#calculating CV
noise_gfp = (np.std(total_io['BL1-A']))/(np.mean(total_io['BL1-A']))
noise_mAmetrine = (np.std(total_io['VL2-A']))/(np.mean(total_io['VL2-A']))
noise_with_corr = (np.std(GFP_to_slope))/(np.mean(GFP_to_slope))
noise_normal_corr = (np.std(mAmetrine_to_slope))/(np.mean(mAmetrine_to_slope))



#calculating noise (CV²)
print ("noise normal to correlation is")
print((noise_normal_corr)*(noise_normal_corr))
print ("noise_GFP is")
print((noise_gfp)*(noise_gfp))
print ("noise_mAmetrine-A is")
print((noise_mAmetrine)*(noise_mAmetrine))
print ("noise of the Expression is")
print((noise_with_corr)*(noise_with_corr))