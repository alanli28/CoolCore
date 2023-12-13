from solver_all import *

steady_temp = []
windspeed = [i * 5 for i in range(2, 17)]
for item in windspeed:
    print("Current Wind Speed, ",item)
    mean_temp_M, mean_temp_C, mean_temp_S, step_size, iter_num,  = solve( 
                                                          num_fins = 12, 
                                                          fin_separation = 1e-3, 
                                                          row_num_M = 1, 
                                                          fin_height=10e-3, 
                                                          forced_convection = True,
                                                          convergence=0.01,
                                                          ini_temp=293,
                                                          wind_speed=item,
                                                          show_plot=False,
                                                          temp_track=True)
    steady_temp.append(mean_temp_M[-1])
print(windspeed)
print(steady_temp)

# windspeed=[10, 20, 30, 40, 50, 60, 70, 80, 15, 25, 35, 45, 55, 65, 75]
# steady_temp = np.array([2.2121577054419177, 1.6864155818588205, 1.491061723785757, 
#                         1.3872691488632305, 1.32301968536175, 1.2799762181733814, 
#                         1.249597830333194, 1.2255071988605637, 1.8721306477108055, 
#                         1.571281765997275, 1.43175251657887, 1.3519149327228, 
#                         1.3003682602976288, 1.26375195128394, 1.236256945832753])

steady_temp = 293 * np.array(steady_temp)

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'
plt.figure(figsize=(6, 3))
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.plot(windspeed, steady_temp, 'o')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Mean Temperature (K)')
plt.tight_layout()
#plt.savefig('path', dpi=500)
plt.show()

