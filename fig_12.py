from solver_all import *
solve(iteration = 100, num_fins = 35, 
      fin_separation = 1e-3, 
      row_num_M = 1, 
      fin_height=40e-3, 
      forced_convection = True,
      convergence=0.01,
      ini_temp=330,
      temp_track = True,
      show_plot=True,
      interval = 100000,
      wind_speed=40)