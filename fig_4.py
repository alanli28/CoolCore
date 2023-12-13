#grid_studies
import matplotlib.pyplot as plt
import numpy as np
import time
from solver_CM import *

elapsed_times = []
stepsize = []
steady_temp = []
iterations = []
factors = [1,2,3,4,5,6,7] #to adjust step size

for i in factors:
    print('Current iteration', i)
    print(elapsed_times)
    print(stepsize)
    print(steady_temp)
    print(iterations)
    start_time = time.time()
    mean_temp_M, mean_temp_C, step_size, tot_iter = solve(row_num_M=1*i, 
                                                          forced_convection=False, 
                                                          convergence=0.00005,
                                                          temp_track = True,
                                                          show_plot=False,
                                                          ini_temp=6540)
    end_time = time.time()
    elapsed_time = end_time - start_time  # Calculate the elapsed time

    stepsize.append(step_size)
    print(mean_temp_M)
    steady_temp.append(mean_temp_M[-1])
    iterations.append(tot_iter)
    elapsed_times.append(elapsed_time)


#step_size = np.array([1.0, 0.5, 0.3333333333333333, 0.25, 0.2, 0.16666666666666666, 0.14285714285714288])
#steady_temp = np.array([22.490440175097905, 22.489394551878814, 22.48927570394604, 22.489243973724463, 22.48925516099185, 22.489268009892882, 22.48927449568771])

stepsize = np.array(stepsize)
steady_temp = 293 * np.array(steady_temp)

#plotting
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'

plt.figure(figsize=(5, 3))
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
x = stepsize
y = steady_temp

def quadratic_function(x, a, b):
    return a * x**2 + b

coefficients = np.polyfit(x, y, 2)
a, b, c = coefficients
x_fit = np.linspace(min(x), max(x), 100)
y_fit = a * x_fit**2 + b*x_fit + c
plt.plot(x_fit, y_fit, label='Fit')
plt.plot(x, y, 'o')

plt.xlabel('Step Size (mm)')
plt.ylabel('Mean Temperature (K)')
plt.legend()
plt.tight_layout()
#plt.savefig('path', dpi=500)
plt.show()
