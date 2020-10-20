import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt     
import os
import FlowCytometryTools
from FlowCytometryTools import FCMeasurement



dir_path = os.path.dirname(os.path.realpath(__file__))   # set the working directory to the location of the Python script file. 


folder_name = "Channel capacity/Dose TNF-a"  # drop .fcs files into (sub)directory. The script will use all .fcs files found in a folder
file_names = os.listdir(folder_name)  #folder including fcs files
elongated_path_name = dir_path +"\\" +folder_name

sample = list(range(0,len(file_names)))
for i in range(0, len(file_names)):
    sample[i] = FCMeasurement(ID='Test Sample', datafile=r''+elongated_path_name + "\\" + file_names[i])




###############################################################    MI calculation   ###############################################################


output_channel = "BL1-A"
#adjustable channel, check your .fcs files for channel names could be e.g. GFP/FITC channel
#output channel is the response of your cells
input_channel = "Time"  # another channel can be used as input, if unlabelled input is used, put any urelated parameter here.


low_percentile = 0.02
high_percentile = 0.98
total_io = pd.DataFrame()
for i in range(0, len(file_names)):
    temp_type = i
    sample_temp = sample[i].data
    temp_input = sample_temp[output_channel]
    temp_output = sample_temp[input_channel]
    temp_io = pd.concat([temp_input, temp_output], axis = 1)
    temp_io["Inputs"] = temp_type
    temp_io = temp_io[temp_io[output_channel]>temp_io[output_channel].quantile(low_percentile)]
    temp_io = temp_io[temp_io[output_channel]<temp_io[output_channel].quantile(high_percentile)]
    temp_io = temp_io[temp_io[input_channel]>temp_io[input_channel].quantile(low_percentile)]
    temp_io = temp_io[temp_io[input_channel]<temp_io[input_channel].quantile(high_percentile)]
    temp_io = temp_io[temp_io[input_channel]>1]
    temp_io = temp_io[temp_io[output_channel]>1]
    total_io = total_io.append(temp_io)


data = np.array(total_io)
data = np.column_stack([data[:,2], data[:,0]])
i_o = data.astype(int)
p = np.repeat(1/len(file_names), len(file_names))



def MI(p):    
    
    global history
    global i_o
    data = i_o
    bin_num = 25  #binning number for y. note that x binning is deternmined by experiment (i,e., dose)
    #usually something between 16-32 is fine as binning number
    starting_data = np.min(data[:,1])  #Low cut-off from blue laser
    ending_data =np.max(data[:,1])  #High cut-off from blue laser
    binning_number = bin_num #The number of bins
    binning_ratio =      (ending_data/ starting_data) ** (1/binning_number)          

    
    #bin range calulation
    next_value = starting_data
    bins_x = range(1,len(p)+1)     #binning for x is already determined by doses
    bins_y =[starting_data]  #initial value for y binning
    
    #generation of y binning interval. it is log scale
    for i in range(0,binning_number):
        next_value = next_value * binning_ratio
        bins_y.append(next_value) 
    
    
    bins_y = np.array(bins_y)
    bins_y = np.delete(bins_y, [0])
    
    x = data[:,0]
    y = data[:,1]
    x = x.astype(int)
    y = y.astype(int)
    probability_x = []
    
    bincount_x = np.bincount(x)
    
    for i in range(0, len(bincount_x)):
        initial_repeat = np.repeat(p[i],bincount_x[i]).tolist()
        probability_x = probability_x + initial_repeat
    
    probability_x = np.array(probability_x)
    
    data = np.column_stack((data,probability_x.T))
    
    
     
    # calculation of probability distribution px
    Px = np.bincount(np.digitize(x,bins_x), weights = data[:,2])
    Px = Px/sum(Px)
    
    # calculation of entropy x    
    entropy_x = 0
    for i in range(0, len(Px)):
        entropy_Px = -Px[i]*np.nan_to_num(np.log2(Px[i]))
        entropy_x = entropy_x + entropy_Px
        
    # calculation of probability distribution py 
    Py = np.bincount(np.digitize(y,bins_y), weights = data[:,2])
    Py = Py/sum(Py)
       
    # calculation of entropy y    
    entropy_y = 0
    for i in range(0, len(Py)):
        entropy_Py = -Py[i]*np.nan_to_num(np.log2(Py[i] + 0.0000001))
        entropy_y = entropy_y + entropy_Py
    
    
    # calculation of p(x,y)
    Jentropy = 0
    
    for i in range(0,len(Px)):
        a = np.digitize(y,bins_y)[np.digitize(x,bins_x) ==i]  # X that correspond to ith y binned region
        w = data[:,2][np.digitize(x,bins_x) ==i] 
        b = np.bincount(a, weights = w)   #count the number of Xs 
        Pyx = np.nan_to_num(b/float(sum(b))) #Marginal probability *
        c =sum( -Px[i]*Pyx *np.nan_to_num(np.log2(Px[i]*Pyx + 0.0000001)) )     #joint entropy calculation
        Jentropy = Jentropy + c  
        print(Jentropy)
    
    MI = entropy_x + entropy_y - Jentropy
    print("Entropy of input = " + str(entropy_x))
    print("Entropy of output = " + str(entropy_y))
    print("Joint Entropy = " + str(Jentropy))
    print("------------------------------------------------------------------")
    print("maximal Mutual Information (channel capacity) = " + str(MI))
    print("------------------------------------------------------------------")
    
    values = p
    history.append(values)
    print(values)
    return -MI




bin_vs_MI = []
for i in range(2,1):
    value = MI(i)
    bin_vs_MI.append(value)
    
#maximizing mutual information 
history = []


def constraint1(p):
    return sum(p)- 1


    
b = np.array([[0.0, 10.0]])
bnds = np.repeat(b, len(p),axis=0)

#con1 = {"type": "ineq", "fun": constraint1}
con1 = {"type": "eq", "fun": constraint1}

cons = [con1]


    
minimize(MI, p,  bounds = bnds, constraints = cons,options={ 'disp': True})

history = np.array(history)
pd.history = pd.DataFrame(history)
history

fig, ax = plt.subplots()

im = ax.imshow(history[:,0:len(p)].T, aspect = "auto")
plt.colorbar(im,orientation="horizontal")
ax.set_title("Mutual information")
plt.show()



















