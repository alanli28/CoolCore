from solver_CM import*
import time
convection = [False, False, False, False, True, True, True, True]
iterations = [10000, 10000, 1000000, 1000000, 10000, 10000, 1000000, 1000000]
stepsize = [1, 10, 1, 10, 1, 10, 1, 10]
time_taken = []
for i in range(2):
    start_time = time.time()
    mean_temp_M, mean_temp_C, step_size, tot_iter = solve(
        iteration=iterations[i], 
        row_num_M=int(stepsize[i]), 
        forced_convection=convection[i], 
        convergence=0.0001, 
        temp_track = True,
        show_plot=False,
        interval=100)
    end_time = time.time()
    elapsed_time = end_time - start_time 
    time_taken.append(elapsed_time)
print("Time taken for each task in seconds:", time_taken)
