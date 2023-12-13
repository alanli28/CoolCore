from solver_all import *

steady_temp_num = []
num_fins = [int(i) for i in range(10,26)]

steady_temp_height = []
fin_height = [i * 1e-3 for i in range(1, 41)]

for item in num_fins:
   print("Current fin number: ",item)
   mean_temp_M, mean_temp_C, mean_temp_S, step_size, iter_num = solve(num_fins = item, 
                                                         fin_separation = 2e-3, 
                                                         row_num_M = 1, 
                                                         fin_height=10e-3, 
                                                         forced_convection = False,
                                                         convergence=0.1,
                                                         ini_temp=293,
                                                         temp_track = True,
                                                         show_plot=False,
                                                         interval = 100000)
   steady_temp_num.append(mean_temp_M[-1])

for item in fin_height:
   print("Current fin height: ",item)
   mean_temp_M, mean_temp_C, mean_temp_S, step_size, iter_num = solve(num_fins = 12, 
                                                         fin_separation = 2e-3, 
                                                         row_num_M = 1, 
                                                         fin_height=item, 
                                                         forced_convection = False,
                                                         convergence=0.1,
                                                         ini_temp=293,
                                                         temp_track = True,
                                                         show_plot=False,
                                                         interval = 100000)
   steady_temp_height.append(mean_temp_M[-1])

print(num_fins)
print(steady_temp_num)
print(fin_height)
print(steady_temp_height)

#fin_height = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009000000000000001, 0.01, 0.011, 0.012, 0.013000000000000001, 0.014, 0.015, 0.016, 0.017, 0.018000000000000002, 0.019, 0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026000000000000002, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036000000000000004, 0.037, 0.038, 0.039, 0.04]
#steady_temp_height = [11.829475768869234, 10.287736409631965, 9.212515273495123, 8.37125359646327, 7.689460432249462, 7.120185276088421, 6.714145659667463, 6.365704193165994, 5.996986733296116, 5.73756348524891, 5.507810765543837, 5.302820997255043, 5.118721480618765, 4.896282359843774, 4.748352626181175, 4.613273740423309, 4.489412820761751, 4.375405684340991, 4.270103480807239, 4.172531583805297, 4.081857575131836, 3.9973660568841347, 3.918438649372752, 3.844537969031483, 3.77519469075868, 3.709997022366042, 3.648582081340288, 3.5906287837465847, 3.535851944040877, 3.4839973513139872, 3.4348376380450323, 3.3881687960490514, 3.3438072240214205, 3.3015872141276614, 3.2613588030857277, 3.2229859273385624, 3.186344833110648, 3.1513227010541978, 3.11781645232443, 3.085731708665669]
steady_temp_height = 293 * np.array(steady_temp_height)
fin_height = 1000 * np.array(fin_height)

#num_fins = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
#steady_temp_num = [6.4206240624047535, 6.05259691742553, 5.73756348524891, 5.4650471978445925, 5.22645660604538, 5.071187600547651, 4.882504935379643, 4.713740601151902, 4.561556573751186, 4.423916145759276, 4.298536654003371, 4.184128668001464, 4.118190695172832, 4.020978774960813, 3.9310795760207493, 3.8479414190896355]
steady_temp_num = 293 * np.array(steady_temp_num)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5), gridspec_kw={'hspace': 0.4})

ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
ax1.plot(num_fins, steady_temp_num, 'o')
ax1.set_xlabel('Number of Fins')
ax1.set_ylabel('Mean Temperature (K)')

ax2.grid(True, which='both', linestyle='--', linewidth=0.5)
ax2.plot(fin_height, steady_temp_height, 'o', color = 'dimgrey')
ax2.set_xlabel('Fin Height (mm)')
ax2.set_ylabel('Mean Temperature (K)')

ax1.annotate('(a)', (0.47, 1.05), xycoords='axes fraction')
ax2.annotate('(b)', (0.47, 1.05), xycoords='axes fraction')
#plt.savefig('path', dpi = 500)
plt.show()

