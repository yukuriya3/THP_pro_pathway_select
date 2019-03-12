import numpy as np
from numpy.random import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import csv
import pandas as pd

# model import
import test_model
func = [test_model.model1]



number = 0


#df = []

eval = []
p1_pv = []
p2_pv = []

DHPAA_Dopamin_ratio = []


for number in range(1):
 
 if number == 0:
  seed(103)
  # file open
  with open('result_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_product_inhibition.csv','a',newline='') as f:
   writer = csv.writer(f)
   header = ['k0', 'Vmax1', 'KmA', 'Vmax2', 'KmB', 'Kmo2', 'O2', 'k3', 'k4', 'Inh', 'KiB', 'KiC']
   writer.writerow(header)
   
   Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 1.0, 10, 10, 10]
   
   for var in range(0, 10000):
    
    # initial values
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # setting
    dt = 0.1
    T = 100.0
    times = np.arange(0.0, T, dt)
    
    Param = Param_ub * rand(12)
    
    k0 = Param[0]
    Vmax1 = Param[1]
    KmA = Param[2]
    Vmax2 = Param[3]
    KmB = Param[4]
    Kmo2 = Param[5]
    O2 = Param[6]
    k3 = Param[7]
    k4 = Param[8]
    Inh = Param[9]
    KiB = Param[10]
    KiC = Param[11]
    
    
    # numerical solution using scipy.integrate
    args = (Param[0], Param[1], Param[2], Param[3], Param[4], Param[5], Param[6], Param[7], Param[8], Param[9], Param[10], Param[11])
    orbit = odeint(func[0], u0, times, args)
    
    
    # output of parameter sets
    writer.writerow(Param)
    eval.append(orbit[len(orbit)-1, 5])
    
    
  eval2 = np.c_[eval] / 50 * 100
  np.savetxt('eval_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_product_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
  df = eval2
  
  # Histgram using matplotlib
  
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  ax.hist(eval2, bins=50)
  ax.set_title('Histgram of Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_xlabel('Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_ylabel('Frequency (-)', fontsize=14)
  
  plt.tick_params(labelsize=14)
  
  plt.tight_layout()
  
  plt.savefig("output_histgram_eval_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_product_inhibition.png")
  #plt.show()
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['DDC+DHPAAS\nno drain\n+inhibition'])
  
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  #plt.xlabel('DDC+DHPAAS\nno drain\n+inhibition')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  plt.setp(bp['medians'][0], color='green', linewidth=2) #メディアン
  
  plt.savefig("output_boxplot_eval_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_product_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []


