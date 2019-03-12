import numpy as np
from numpy.random import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import csv
import pandas as pd

# model import
import model
func = [model.model1, model.model2, model.model3, model.model4, model.model5, model.model6]



number = 0


#df = []

eval = []
p1_pv = []
p2_pv = []

DHPAA_Dopamin_ratio = []


for number in range(6):
 
 if number == 0:
  seed(103)
  # file open
  with open('result_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_inhibition.csv','a',newline='') as f:
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
  np.savetxt('eval__DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
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
  
  plt.savefig("output_histgram_eval_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_inhibition.png")
  #plt.show()
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['Path2'])
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  plt.xlabel('Pathway')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  
  plt.savefig("output_boxplot_eval_DDC_plus_MAO_pathway_wo_DHPAA_drain_wi_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []
  
 elif number == 1:
  
  seed(103)
  
  # file open
  with open('result_DDC_plus_MAO_pathway_wi_DHPAA_drain_wi_inhibition.csv','a',newline='') as f:
   writer = csv.writer(f)
   header = ['k0', 'Vmax1', 'KmA', 'Vmax2', 'KmB', 'Kmo2', 'O2', 'k3', 'k4', 'Inh', 'KiB', 'KiC', 'k5']
   writer.writerow(header)
   
   Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 1.0, 10, 10, 10, 1.0]
   
   for var in range(0, 10000):
    
    # initial values
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    # setting
    dt = 0.1
    T = 100.0
    times = np.arange(0.0, T, dt)
    
    Param = Param_ub * rand(13)
    
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
    k5 = Param[12]
    
    
    # numerical solution using scipy.integrate
    args = (Param[0], Param[1], Param[2], Param[3], Param[4], Param[5], Param[6], Param[7], Param[8], Param[9], Param[10], Param[11], Param[12])
    orbit = odeint(func[1], u0, times, args)
    
    # output of parameter sets
    writer.writerow(Param)
    eval.append(orbit[len(orbit)-1, 5])
    
  # output of THP yield
  eval2 = np.c_[eval] / 50 * 100
  np.savetxt('eval_DDC_plus_MAO_pathway_wi_DHPAA_drain_wi_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
  df1 = eval2
  
  
  # Histgram using matplotlib
  
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  ax.hist(eval2, bins=50)
  ax.set_title('Histgram of Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_xlabel('Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_ylabel('Frequency (-)', fontsize=14)
  
  plt.tick_params(labelsize=14)
  
  plt.savefig("output_histgram_eval_DDC_plus_MAO_pathway_wi_DHPAA_drain_wi_inhibition.png")
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['Path2'])
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  plt.xlabel('Pathway')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  
  plt.savefig("output_boxplot_eval_DDC_plus_MAO_pathway_wi_DHPAA_drain_wi_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []
  
 elif number == 2:
  
  seed(103)
  
  # file open
  with open('result_DDC_plus_MAO_pathway_wi_DHPAA_drain_wo_inhibition.csv','a',newline='') as f:
   writer = csv.writer(f)
   header = ['k0', 'Vmax1', 'KmA', 'Vmax2', 'KmB', 'Kmo2', 'O2', 'k3', 'k4', 'k5']
   writer.writerow(header)
   
   Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 1.0, 1.0]
   
   for var in range(0, 10000):
    
    # initial values
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    # setting
    dt = 0.1
    T = 100.0
    times = np.arange(0.0, T, dt)
    
    Param = Param_ub * rand(10)
    
    k0 = Param[0]
    Vmax1 = Param[1]
    KmA = Param[2]
    Vmax2 = Param[3]
    KmB = Param[4]
    Kmo2 = Param[5]
    O2 = Param[6]
    k3 = Param[7]
    k4 = Param[8]
    k5 = Param[9]
    
    
    # numerical solution using scipy.integrate
    args = (Param[0], Param[1], Param[2], Param[3], Param[4], Param[5], Param[6], Param[7], Param[8], Param[9])
    orbit = odeint(func[2], u0, times, args)
    
    # output of parameter sets
    writer.writerow(Param)
    eval.append(orbit[len(orbit)-1, 5])
    
  # output of THP yield
  eval2 = np.c_[eval] / 50 * 100
  np.savetxt('eval_DDC_plus_MAO_pathway_wi_DHPAA_drain_wo_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
  df2 = eval2
  
  
  # Histgram using matplotlib
  
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  ax.hist(eval2, bins=50)
  ax.set_title('Histgram of Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_xlabel('Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_ylabel('Frequency (-)', fontsize=14)
  
  plt.tick_params(labelsize=14)
  
  plt.savefig("output_histgram_eval_DDC_plus_MAO_pathway_wi_DHPAA_drain_wo_inhibition.png")
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['Path2'])
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  plt.xlabel('Pathway')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  
  plt.savefig("output_boxplot_eval_DDC_plus_MAO_pathway_wi_DHPAA_drain_wo_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []
  
  
 elif number == 3:
  
  seed(110)
  
  # file open
  with open('result_DDC_plus_DHPAAS_wo_DHPAA_drain_wi_product_inhibition.csv','a',newline='') as f:
   writer = csv.writer(f)
   header = ['k0', 'Vmax1', 'KmA', 'Vmax2', 'KmA2', 'Kmo2', 'O2', 'k3', 'k4', 'KiB', 'KiC']
   writer.writerow(header)
   
   Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 100, 10.0, 10.0]
   
   for var in range(0, 10000):
    #Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 100, 10.0, 10.0, 1.0]
    #Param_ub = [100] * 9  # initialization value: 100, elements: 9 
    #Param_ub[6] = 5.5e-5
    #Param = Param_ub * rand(12)
    #Param = [10.0, 8.0/3.0, 28.0]
    
    
    # initial values
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # setting
    dt = 0.1
    T = 100.0
    times = np.arange(0.0, T, dt)
    
    Param = Param_ub * rand(11)
    
    k0 = Param[0]
    Vmax1 = Param[1]
    KmA = Param[2]
    Vmax2 = Param[3]
    KmA2 = Param[4]
    Kmo2 = Param[5]
    O2 = Param[6]
    k3 = Param[7]
    k4 = Param[8]
    KiB = Param[9]
    KiC = Param[10]
    
    
    # numerical solution using scipy.integrate
    args = (Param[0], Param[1], Param[2], Param[3], Param[4], Param[5], Param[6], Param[7], Param[8], Param[9], Param[10])
    orbit = odeint(func[3], u0, times, args)
    
    # output of parameter sets
    writer.writerow(Param)
    eval.append(orbit[len(orbit)-1, 5])
    
  # output of THP yield
  eval2 = np.c_[eval] / 50 * 100
  np.savetxt('eval_DDC_plus_DHPAAS_wo_DHPAA_drain_wi_product_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
  df3 =eval2
  
  # Histgram using matplotlib
  
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  ax.hist(eval2, bins=50)
  ax.set_title('Histgram of Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_xlabel('Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_ylabel('Frequency (-)', fontsize=14)
  
  plt.tick_params(labelsize=14)
  
  plt.savefig("output_histgram_eval_DDC_plus_DHPAAS_wo_DHPAA_drain_wi_product_inhibition.png")
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['Path2'])
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  plt.xlabel('Pathway')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  
  plt.savefig("output_boxplot_eval_DDC_plus_DHPAAS_wo_DHPAA_drain_wi_product_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []
  
  
 elif number == 4:
  
  seed(110)
  
  with open('result_DDC_plus_DHPAAS_wi_DHPAA_drain_wi_product_inhibition.csv','a',newline='') as f:
   writer = csv.writer(f)
   header = ['k0', 'Vmax1', 'KmA', 'Vmax2', 'KmA2', 'Kmo2', 'O2', 'k3', 'k4', 'KiB', 'KiC', 'k5']
   writer.writerow(header)
   
   Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 100, 10.0, 10.0, 1.0]
   
   
   #eval = []
   #p1_pv = []
   #p2_pv = []
   
   #DHPAA_Dopamin_ratio = []
   
   
   
   for var in range(0, 10000):
    #Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 100, 10.0, 10.0, 1.0]
    #Param_ub = [100] * 9  # initialization value: 100, elements: 9 
    #Param_ub[6] = 5.5e-5
    #Param = Param_ub * rand(12)
    #Param = [10.0, 8.0/3.0, 28.0]
    
    
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
    KmA2 = Param[4]
    Kmo2 = Param[5]
    O2 = Param[6]
    k3 = Param[7]
    k4 = Param[8]
    KiB = Param[9]
    KiC = Param[10]
    k5 = Param[11]
    
    # numerical solution using scipy.integrate
    args = (Param[0], Param[1], Param[2], Param[3], Param[4], Param[5], Param[6], Param[7], Param[8], Param[9], Param[10], Param[11])
    orbit = odeint(func[4], u0, times, args)
    
    
    
    # output of parameter sets
    writer.writerow(Param)
    eval.append(orbit[len(orbit)-1, 5])
  
  
  # output of THP yield
  eval2 = np.c_[eval] / 50 * 100
  np.savetxt('eval_DDC_plus_DHPAAS_wi_DHPAA_drain_wi_product_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
  df4 = eval2
  
  
  # Histgram using matplotlib
  
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  ax.hist(eval2, bins=50)
  ax.set_title('Histgram of Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_xlabel('Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_ylabel('Frequency (-)', fontsize=14)
  
  plt.tick_params(labelsize=14)
  
  plt.savefig("output_histgram_eval_DDC_plus_DHPAAS_wi_DHPAA_drain_wi_product_inhibition.png")
  #plt.show()
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['Path2'])
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  plt.xlabel('Pathway')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  
  plt.savefig("output_boxplot_eval_DDC_plus_DHPAAS_wi_DHPAA_drain_wi_product_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []

 elif number == 5:
  
  seed(213)
  
  with open('result_DDC_plus_DHPAAS_wi_DHPAA_drain_wo_product_inhibition.csv','a',newline='') as f:
   writer = csv.writer(f)
   header = ['k0', 'Vmax1', 'KmA', 'Vmax2', 'KmA2', 'Kmo2', 'O2', 'k3', 'k4','k5']
   writer.writerow(header)
   
   Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 100, 1.0]
   
   
   #eval = []
   #p1_pv = []
   #p2_pv = []
   
   #DHPAA_Dopamin_ratio = []
   
   
   
   for var in range(0, 10000):
    #Param_ub = [10.0, 1.0e+4, 10, 1.0e+4, 10, 1.0, 0.236, 1.0, 100, 10.0, 10.0, 1.0]
    #Param_ub = [100] * 9  # initialization value: 100, elements: 9 
    #Param_ub[6] = 5.5e-5
    #Param = Param_ub * rand(12)
    #Param = [10.0, 8.0/3.0, 28.0]
    
    
    # initial values
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # setting
    dt = 0.1
    T = 100.0
    times = np.arange(0.0, T, dt)
    
    Param = Param_ub * rand(10)
    
    k0 = Param[0]
    Vmax1 = Param[1]
    KmA = Param[2]
    Vmax2 = Param[3]
    KmA2 = Param[4]
    Kmo2 = Param[5]
    O2 = Param[6]
    k3 = Param[7]
    k4 = Param[8]
    k5 = Param[9]
    
    # numerical solution using scipy.integrate
    args = (Param[0], Param[1], Param[2], Param[3], Param[4], Param[5], Param[6], Param[7], Param[8], Param[9])
    orbit = odeint(func[5], u0, times, args)
    
    
    
    # output of parameter sets
    writer.writerow(Param)
    eval.append(orbit[len(orbit)-1, 5])
  
  
  # output of THP yield
  eval2 = np.c_[eval] / 50 * 100
  np.savetxt('eval_DDC_plus_DHPAAS_wi_DHPAA_drain_wo_product_inhibition.csv', eval2, delimiter=',',header='Eval_index', fmt='%.5f')
  
  df5 = eval2
  
  data = np.concatenate((df, df1, df2, df3, df4, df5), 1)
  
  print(data)
  
  # Histgram using matplotlib
  
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  ax.hist(eval2, bins=50)
  ax.set_title('Histgram of Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_xlabel('Conversion efficiency from Substrate to Target (%)', fontsize=14)
  ax.set_ylabel('Frequency (-)', fontsize=14)
  
  plt.tick_params(labelsize=14)
  
  plt.savefig("output_histgram_eval_DDC_plus_DHPAAS_wi_DHPAA_drain_wo_product_inhibition.png")
  #plt.show()
  plt.close()
  
  ## output of statistics
  print(pd.DataFrame(pd.Series(eval2.ravel()).describe()).transpose())
  
  ## Boxplot using matplotlib
  fig = plt.figure()
  ax = fig.add_subplot(1,1,1)
  
  bp = ax.boxplot(eval2)
  ax.set_xticklabels(['Path2'])
  
  plt.title('Box plot of Conversion efficiency from Substrate to Target (%)')
  plt.grid()
  plt.xlabel('Pathway')
  plt.ylabel('Conversion efficiency from Substrate to Target (%)')
  plt.ylim([0,100])
  
  plt.savefig("output_boxplot_eval_DDC_plus_DHPAAS_wi_DHPAA_drain_wo_product_inhibition.png")
  #plt.show()
  plt.close()
  
  eval = []



np.savetxt('For_eval_models.csv', data, delimiter=',', header="DDC+MAO-drain+inhibition','DDC+MAO+drain+inhibition','DDC+MAO+drain-inhibition','DDC+DHPAAS-drain+inhibition','DDC+DHPAAS+drain+inhibition','DDC+DHPAAS+drain-inhibition'", fmt='%.5f', comments = '')

## Boxplot using matplotlib
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#csv_input.plot(ax=ax, kind='box')
#bp=csv_input.plot.box(ax=ax)
#bp=data.plot.box(ax=ax)
#bp=csv_input.plot.box(ax=ax, patch_artist=True)
#ap = pd.DataFrame(data, columns=['A', 'B', 'C', 'D'])
ap = pd.DataFrame(data)
bp = ap.plot.box(ax=ax)


#bp = ax.boxplot(data)
#bp = ax.boxplot(csv_input)
ax.set_xticklabels(['DDC+MAO\n -drain\n+inhibition', 'DDC+MAO\n +drain\n+inhibition', 'DDC+MAO\n +drain\n-inhibition', 'DDC+DHPAAS\n -drain\n+inhibition','DDC+DHPAAS\n +drain\n+inhibition', 'DDC+DHPAAS\n +drain\n-inhibition'])

plt.title('Comparison among models\n wi/wo DHPAA drain and product inhibition', fontsize=14)
#plt.grid()
plt.xlabel('Model', fontsize=14)
plt.ylabel('THP Yield (%)', fontsize=14)
plt.ylim([0,100])
plt.tick_params(labelsize=14)
#plt.setp(bp['medians'][0], color='grean', linewidth=2) #Median

plt.setp(ax.get_xticklabels(), rotation=0)
plt.tight_layout()

plt.savefig("output_boxplot_for_comparison_among_6_models.png")

#plt.show()
plt.close()

