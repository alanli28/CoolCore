import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def create_mesh(column_num_M, row_num_M, column_num_C, row_num_C, ini_temp):
    M_width = 14e-3 #in meters
    M_height = 1e-3 #in meters
    C_width = 20e-3 #in meters
    C_height = 2e-3 #in meters
    
    M_dx = M_width / (column_num_M)
    M_dy = M_height / (row_num_M)
    C_dx = C_width / (column_num_C)
    C_dy = C_height / (row_num_C)

    if C_dx != C_dy or M_dx != M_dy or C_dx != M_dy:
      print("matrix dimension wrong!!")
    
    step_size = C_dx #in meters

    microprocessor_mesh = np.pad(np.full((row_num_M, column_num_M), ini_temp, 
                                         dtype=float), ((1, 1), (1, 1)), 
                                         mode='constant', constant_values=0)
    ceramic_mesh = np.pad(np.full((row_num_C, column_num_C), ini_temp, 
                                  dtype=float), ((1, 1), (1, 1)), 
                                  mode='constant', constant_values=0)
    return microprocessor_mesh, ceramic_mesh, step_size


def update_C(matrix, step_size, k, forced, wind_speed = 20, length_scale = 0.001):
    rows, columns = matrix.shape

    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k

    top_diff = matrix[1, 1:columns-1] - 1
    bottom_diff = matrix[-2, 1:columns-1] - 1
    left_diff = matrix[1:rows-1, 1] - 1
    right_diff = matrix[1:rows-1, -2] - 1

    matrix[0, 1:columns-1] = matrix[1, 1:columns-1] - constant * \
        np.sign(top_diff) * np.abs(top_diff) ** power_factor
    matrix[-1, 1:columns-1] = matrix[-2, 1:columns-1] - constant * \
        np.sign(bottom_diff) * np.abs(bottom_diff) ** power_factor
    matrix[1:rows-1, 0] = matrix[1:rows-1, 1] - constant * np.sign(left_diff) \
        * np.abs(left_diff) ** power_factor
    matrix[1:rows-1, -1] = matrix[1:rows-1, -2] - constant * \
        np.sign(right_diff) * np.abs(right_diff) ** power_factor
    return matrix

def update_M(matrix, step_size, k, forced, wind_speed = 20, length_scale = 0.001):
    rows, columns = matrix.shape

    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k

    bottom_diff = matrix[-2, 1:columns-1] - 1
    left_diff = matrix[1:rows-1, 1] - 1
    right_diff = matrix[1:rows-1, -2] - 1

    matrix[-1, 1:columns-1] = matrix[-2, 1:columns-1] - constant * \
         np.sign(bottom_diff) * np.abs(bottom_diff) ** power_factor
    matrix[1:rows-1, 0] = matrix[1:rows-1, 1] - constant * np.sign(left_diff)\
          * np.abs(left_diff) ** power_factor
    matrix[1:rows-1, -1] = matrix[1:rows-1, -2] - constant * \
         np.sign(right_diff) * np.abs(right_diff) ** power_factor
    return matrix

def update_interface(small_matrix, big_matrix, center_start, center_end):
    # replace ghost points at interface
    big_matrix[-1, center_start:center_end+1] = small_matrix[1, 1:-1]
    small_matrix[0, 1:-1] = big_matrix[-2, center_start:center_end+1]

def bottom_edge_stencil(buffer, center_start, center_end, matrix):
    buffer[-2, 1:center_start] = 0.25 * (matrix[-1, 1:center_start] + 
                                         matrix[-3, 1:center_start] + 
                                         matrix[-2, 2:center_start+1] + 
                                         matrix[-2, 0:center_start-1])
    buffer[-2, center_end+1:-1] = 0.25 * (matrix[-1, center_end+1:-1] + 
                                          matrix[-3, center_end+1:-1] + 
                                          matrix[-2, center_end+2:] + 
                                          matrix[-2, center_end:-2])


def progress_update(iter_num, progress_interval, iteration, step_size, k_eff_CM, 
                    matrix_M, matrix_C, center_start_C, center_end_C, 
                    mean_temp_M, forced, length_scale = 0.001, 
                    k_M = 150, wind_speed = 20):
    power_difference_percent = 1000
    temp_difference_percent = 1000
    if iter_num % progress_interval == 0 or iter_num == iteration - 1:
        
        m = matrix_M[0:-1, 1:-1]
        c = matrix_C[1:-1, 1:-1]

        if forced == False:
            power_factor = 4 / 3
            constant = 1.31 * step_size * length_scale * 293**(1/3)/ k_M
        else: 
            power_factor = 1
            constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k_M

        M_bottom = np.sum((m[-1,:] - 1)**(power_factor)*constant)
        M_top = np.sum((m[1,:] - m[0,:]))*k_eff_CM/k_M
        M_left = np.sum((m[1:,0] - 1)**(power_factor)*constant)
        M_right = np.sum((m[1:,-1] - 1)**(power_factor)*constant)

        ceramic_top = np.sum((c[0,:] - 1)**(power_factor)*constant)
        ceramic_left = np.sum((c[:,0] - 1)**(power_factor)*constant)
        ceramic_right = np.sum((c[:,-1] - 1)**(power_factor)*constant)
        ceramic_bottom = np.sum((matrix_C[-2, 1:center_start_C] - 1)**
                                (power_factor)*constant) + \
                                    np.sum((matrix_C[-2, center_end_C+1:-1] - 1)
                                           **(power_factor)*constant)

        microprocessor_power_out = M_bottom + M_left + M_right + M_top
        overall_power_out =  ceramic_top + ceramic_left + ceramic_right + \
            ceramic_bottom + M_bottom + M_left + M_right
        print('===============================================================================')
        print(f"Iteration {iter_num+1}/{iteration} complete ({(iter_num+1)/iteration*100:.0f}%)")
        print('-----------------------------------------------------------------')
        print('Theoretical power in     ', 5e8*0.001*0.014/(293*150))
        print('Microprocessor power out ', f"{microprocessor_power_out:.10g}")
        print('Overall power out        ', f"{overall_power_out:.10g}")
        if iter_num>5:
            power_difference_percent = (overall_power_out - 5e8*0.001*0.014/(293*150))*100/overall_power_out
            temp_difference_percent = (mean_temp_M[-1] - mean_temp_M[-2])*100/mean_temp_M[-1]
            print('-----------------------------------------------------------------')
            print("Microprocessor temperature change from last iteration", f"{(mean_temp_M[-1] - mean_temp_M[-2])*100/mean_temp_M[-1]:.6g}", "%")
            print("Percentage difference from theoretical power out     ", f"{(overall_power_out - 5e8*0.001*0.014/(293*150))*100/overall_power_out:.6g}", "%")
            print("Mean temperature in Microprocessor                   ", f"{mean_temp_M[-1]*293:.6g}", "K")
        print(" ")
    return power_difference_percent/100, temp_difference_percent/100

def create_plot(matrix_M, matrix_C, step_size, grid=False):
    
    M = matrix_M[1:-1, 1:-1]
    C = matrix_C[1:-1, 1:-1]

    M_row_num, M_column_num = M.shape
    C_row_num, C_column_num = C.shape

    total_rows = M_row_num + C_row_num
    total_cols = C_column_num

    big_mesh = np.full((total_rows, total_cols), np.nan)

    def insert_centered(big_mesh, matrix, start_row):
        row_num, col_num = matrix.shape
        col_start = (big_mesh.shape[1] - col_num) // 2
        big_mesh[start_row:start_row + row_num, col_start:col_start + col_num] = matrix

    insert_centered(big_mesh, C, 0)
    insert_centered(big_mesh, M, C_row_num)

    heatmap = big_mesh * 293
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 4))
    im1 = ax1.imshow(heatmap, cmap='turbo', extent = (0,int(len(big_mesh[0])*step_size),0,int(len(big_mesh)*step_size)), interpolation='nearest')
    cbar = plt.colorbar(im1, orientation='horizontal')
    cbar.set_label('Temperature Scale (K)')
    ax1.set_xlabel('x, mm')
    ax1.set_ylabel('y, mm')


    if grid: 
        ax1.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        ax1.set_xticks(np.arange(-0.5, total_cols, 1), minor=True)
        ax1.set_yticks(np.arange(-0.5, total_rows, 1), minor=True)

        ax1.tick_params(which='minor', size=0)
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
    #plt.savefig('path', dpi=500)
    plt.show()



def solve(row_num_M, iteration = 50000000, k_M = 150, k_C = 230, q = 5e8, 
                        ini_temp = 293, length_scale = 0.001, plot_grid = False, 
                        forced_convection = False, convergence = 0.001,
                        temp_track = False, show_plot = True, interval = 10000):
    ini_temp = ini_temp/293
    column_num_M = int(14 * row_num_M)
    row_num_C = int(2 * row_num_M)
    column_num_C = int(column_num_M*20/14)
    matrix_M, matrix_C, step_size = create_mesh(column_num_M, row_num_M, 
                                                column_num_C, row_num_C, 
                                                ini_temp)
    step_size = step_size/length_scale
    heat_term = 0.25 * length_scale ** 2 * step_size ** 2 * q / (k_M * 293)
    k_eff_CM = 2*k_M*k_C/(k_M +k_C)
    interface_heat_term = (length_scale ** 2 * step_size ** 2 * q
                           )/((k_eff_CM + 3 * k_M) * 293)
    center_start_C = (len(matrix_C[0]) - len(matrix_M[0])) // 2 + 1
    center_end_C = (len(matrix_C[0]) + len(matrix_M[0])) // 2 - 2

    mean_temp_M = []
    mean_temp_C = []

    progress_interval = max(iteration // interval, 1)

    for i in range(iteration):
        #update ghost points for all components
        update_C(matrix_C, step_size, k_C, length_scale = length_scale, 
                 forced = forced_convection)
        update_M(matrix_M, step_size, k_M, length_scale = length_scale, 
                 forced = forced_convection)
        
        power_frac_diff, temp_frac_diff = progress_update(iter_num = i, progress_interval = progress_interval, 
                        iteration = iteration, step_size = step_size, 
                        k_eff_CM = k_eff_CM, matrix_M = matrix_M, 
                        matrix_C = matrix_C,
                        center_start_C = center_start_C, 
                        center_end_C = center_end_C, 
                        mean_temp_M = mean_temp_M, 
                        forced = forced_convection)

        # replace ghost points at Ceramic Microprocessor interface
        update_interface(small_matrix = matrix_M, big_matrix = matrix_C, 
                         center_start=center_start_C, center_end=center_end_C)

        buffer_M = np.zeros_like(matrix_M)
        buffer_C = np.zeros_like(matrix_C)

        # Apply the stencil to the internal matrix_M cells, except second row
        buffer_M[2:-1, 1:-1] = 0.25 * (matrix_M[3:, 1:-1] + 
                                       matrix_M[1:-2, 1:-1] + 
                                       matrix_M[2:-1, 2:] + 
                                       matrix_M[2:-1, :-2]) + heat_term

        # Apply the stencil to the internal matrix_C cells, except the second 
        # last row
        buffer_C[1:-2, 1:-1] = 0.25 * (matrix_C[2:-1, 1:-1] + 
                                       matrix_C[:-3, 1:-1] + 
                                       matrix_C[1:-2, 2:] + 
                                       matrix_C[1:-2, :-2])
        
        # Apply the stencil to specific parts of the second last row of matrix_C
        bottom_edge_stencil(buffer = buffer_C, center_start = center_start_C, 
                            center_end = center_end_C, matrix = matrix_C)

        #Apply the stencil to second row of matrix_M with special stencil
        buffer_M[1, 1:-1] = 1/(k_eff_CM + 3 * k_M) * (k_eff_CM * 
                                                      matrix_M[0, 1:-1] + 
                                                      k_M * matrix_M[2, 1:-1] + 
                                                      k_M * matrix_M[1, 2:] + 
                                                      k_M * matrix_M[1, :-2]) +\
                                                        interface_heat_term

        #Apply the stencil to last row of matrix_C with special stencil
        buffer_C[-2, center_start_C:center_end_C+1] = 1/(k_eff_CM + 3 * k_C) *\
              (k_eff_CM * matrix_C[-1, center_start_C:center_end_C+1] + k_C *\
                matrix_C[-3, center_start_C:center_end_C+1] + k_C *\
                      matrix_C[-2, center_start_C+1:center_end_C+2] + k_C *\
                          matrix_C[-2, center_start_C-1:center_end_C])


        matrix_C[1:-1, 1:-1] = buffer_C[1:-1, 1:-1]
        matrix_M[1:-1, 1:-1] = buffer_M[1:-1, 1:-1]
        # Calculate and store the mean temperature
        mean_temp_M.append(np.mean(matrix_M[1:-1, 1:-1]))
        mean_temp_C.append(np.mean(matrix_C[1:-1, 1:-1]))
        if np.abs(power_frac_diff) < convergence and np.abs(temp_frac_diff) < convergence:
            print("Convergence Criteria Met!")
            break
    if show_plot == True:
        create_plot(matrix_M, matrix_C, step_size, grid = plot_grid)

    if temp_track == True:
        return mean_temp_M, mean_temp_C, step_size, i
    else:
        return step_size