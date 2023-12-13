from solver_all import *
#Figure 7
solve(iteration = 100, num_fins = 35, 
      fin_separation = 1e-3, 
      row_num_M = 5, 
      fin_height=25e-3, 
      forced_convection = True,
      convergence=0.01,
      ini_temp=315,
      temp_track = True,
      show_plot=True,
      interval = 100000)

#Figure 8
solve(iteration = 100, num_fins = 20, 
      fin_separation = 1e-3, 
      row_num_M = 10, 
      fin_height=25e-3, 
      forced_convection = True,
      convergence=0.01,
      ini_temp=351,
      temp_track = True,
      show_plot=True,
      interval = 100000)