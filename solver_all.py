import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def create_mesh(column_num_M, row_num_M, column_num_C, row_num_C, ini_temp):
    M_width = 14e-3 #in meters
    M_height = 1e-3 #in meters
    C_width = 20e-3 #in meters
    C_height = 2e-3 #in meters
    
    # Calculate M step size
    M_dx = M_width / (column_num_M)
    M_dy = M_height / (row_num_M)
    # Calculate C step size
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

def create_heatsink_mesh(step_size, num_fins, fin_separation, 
                         fin_height = 30e-3, ini_temp = 293):

    # Constants
    S_height, S_width = 4e-3, num_fins * 1e-3 + (num_fins - 1) * fin_separation
    F_width = 1e-3

    # Calculate cell dimensions
    row_num_F, column_num_F = int(fin_height/step_size), int(F_width/step_size)
    row_num_S, column_num_S = int(S_height / step_size), \
        int(np.round(S_width / step_size))
    if (column_num_S - 20e-3/step_size)%2 != 0:
      print("Asymmetrical heatsink placement!")

    # Function to create and pad mesh
    def create_and_pad(rows, columns):
        return np.pad(np.full((rows, columns), ini_temp, dtype=float), 
                      ((1, 1), (1, 1)), mode='constant', constant_values=0)

    # Creating meshes
    sink_mesh = create_and_pad(row_num_S, column_num_S)

    # Creating fin meshes
    fin_meshes = {n: create_and_pad(row_num_F, column_num_F) for n in \
                  range(0, num_fins)}

    return sink_mesh, fin_meshes, row_num_F, column_num_F, row_num_S, \
         column_num_S

def update_C(matrix, step_size, k, forced, wind_speed, length_scale = 0.001):
    rows, columns = matrix.shape
    # Precompute common values
    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k

    # Compute the difference and its sign/absolute value once
    top_diff = matrix[1, 1:columns-1] - 1
    bottom_diff = matrix[-2, 1:columns-1] - 1
    left_diff = matrix[1:rows-1, 1] - 1
    right_diff = matrix[1:rows-1, -2] - 1

    # Apply the update rules using vectorized operations
    matrix[0, 1:columns-1] = matrix[1, 1:columns-1] - constant * \
        np.sign(top_diff) * np.abs(top_diff) ** power_factor
    matrix[-1, 1:columns-1] = matrix[-2, 1:columns-1] - constant * \
        np.sign(bottom_diff) * np.abs(bottom_diff) ** power_factor
    matrix[1:rows-1, 0] = matrix[1:rows-1, 1] - constant * np.sign(left_diff) \
        * np.abs(left_diff) ** power_factor
    matrix[1:rows-1, -1] = matrix[1:rows-1, -2] - constant * \
        np.sign(right_diff) * np.abs(right_diff) ** power_factor
    return matrix

def update_M(matrix, step_size, k, forced, wind_speed, length_scale = 0.001):
    rows, columns = matrix.shape

    # Precompute common values
    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k

    # Compute the difference and its sign/absolute value once
    bottom_diff = matrix[-2, 1:columns-1] - 1
    left_diff = matrix[1:rows-1, 1] - 1
    right_diff = matrix[1:rows-1, -2] - 1

    # Apply the update rules using vectorized operations
    matrix[-1, 1:columns-1] = matrix[-2, 1:columns-1] - constant * \
         np.sign(bottom_diff) * np.abs(bottom_diff) ** power_factor
    matrix[1:rows-1, 0] = matrix[1:rows-1, 1] - constant * np.sign(left_diff)\
          * np.abs(left_diff) ** power_factor
    matrix[1:rows-1, -1] = matrix[1:rows-1, -2] - constant * \
         np.sign(right_diff) * np.abs(right_diff) ** power_factor
    return matrix

def update_fin(dictionary, rows, columns, step_size, k, forced, wind_speed,
                length_scale = 0.001):
    # Precompute common values
    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k

    for key, matrix in dictionary.items():
        # Compute the difference and its sign/absolute value once
        top_diff = matrix[1, 1:columns+1] - 1
        left_diff = matrix[1:rows+1, 1] - 1
        right_diff = matrix[1:rows+1, -2] - 1
        # Apply the update rules using vectorized operations
        matrix[0, 1:columns+1] = matrix[1, 1:columns+1] - constant * \
            np.sign(top_diff) * np.abs(top_diff) ** power_factor
        matrix[1:rows+1, 0] = matrix[1:rows+1, 1] - constant * \
            np.sign(left_diff) * np.abs(left_diff) ** power_factor
        matrix[1:rows+1, -1] = matrix[1:rows+1, -2] - constant * \
            np.sign(right_diff) * np.abs(right_diff) ** power_factor
    
    return dictionary

def update_sink(matrix, rows, columns, step_size, k, forced, wind_speed,
                 length_scale = 0.001):
    # Precompute common values
    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k

    # Compute the difference and its sign/absolute value once
    top_diff = matrix[1, 1:columns+1] - 1
    bottom_diff = matrix[-2, 1:columns+1] - 1
    left_diff = matrix[1:rows+1, 1] - 1
    right_diff = matrix[1:rows+1, -2] - 1

    # Apply the update rules using vectorized operations
    matrix[0, 1:columns+1] = matrix[1, 1:columns+1] - constant * \
        np.sign(top_diff) * np.abs(top_diff) ** power_factor
    matrix[-1, 1:columns+1] = matrix[-2, 1:columns+1] - constant * \
        np.sign(bottom_diff) * np.abs(bottom_diff) ** power_factor
    matrix[1:rows+1, 0] = matrix[1:rows+1, 1] - constant * np.sign(left_diff) \
        * np.abs(left_diff) ** power_factor
    matrix[1:rows+1, -1] = matrix[1:rows+1, -2] - constant * \
        np.sign(right_diff) * np.abs(right_diff) ** power_factor
    
    return matrix

def update_interface(small_matrix, big_matrix, center_start, center_end):
    # replace ghost points at interface
    big_matrix[-1, center_start:center_end+1] = small_matrix[1, 1:-1]
    small_matrix[0, 1:-1] = big_matrix[-2, center_start:center_end+1]

def pair_sink_fin(step_size, dictionary, matrix, num_fins, fin_separation, 
                  column_num_F, length_scale = 0.001):
    separation_num = int(fin_separation / (step_size*length_scale))
    for i in range(num_fins):
        matrix[0, 1+i*(column_num_F+separation_num):
               1+i*(column_num_F+separation_num)+column_num_F] \
                =dictionary[i][-2, 1:-1]
        dictionary[i][-1, 1:-1] = matrix[1, 1+i*(column_num_F+separation_num):
                                         1+i*(column_num_F+separation_num)\
                                            +column_num_F]
    return matrix, dictionary

def power_out_fins(fin_meshes, step_size, forced, wind_speed, k_M = 150, 
                   length_scale = 0.001):
    fins = [mesh[1:-1, 1:-1] for mesh in fin_meshes.values()]
    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k_M
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k_M
    
    total = 0
    for fin in fins:
        fin_top = np.sum((fin[0,:] - 1)**(power_factor)*constant)
        fin_left = np.sum((fin[:,0] - 1)**(power_factor)*constant)
        fin_right = np.sum((fin[:,-1] - 1)**(power_factor)*constant)
        total += fin_top + fin_left + fin_right
    return total

def sink_top_power_out(matrix_S, column_num_F, separation_num, step_size, 
                       forced, wind_speed, k_M = 150, length_scale = 0.001):
    power_out = 0
    current_col = column_num_F
    top_row = matrix_S[1,1:-1]
    if forced == False:
        power_factor = 4 / 3
        constant = 1.31 * step_size * length_scale * 293**(1/3)/ k_M
    else: 
        power_factor = 1
        constant = (11.4+5.7*wind_speed) * step_size * length_scale/ k_M

    while current_col < len(top_row):
        # Sum the power over the gap width
        end_col = current_col + separation_num
        power_out += np.sum((top_row[current_col:end_col] - 1)**\
                            (power_factor)*constant)
        # Move past this gap and the subsequent fin
        current_col += separation_num + column_num_F

    return power_out

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
                    matrix_M, matrix_C, matrix_S, fin_meshes, center_start_C, 
                    center_end_C, center_start_S, center_end_S, column_num_F, 
                    separation_num, mean_temp_M, forced, wind_speed, 
                    length_scale = 0.001, k_M = 150):
    power_difference_percent = 1000
    temp_difference_percent = 1000
    if iter_num % progress_interval == 0 or iter_num == iteration - 1:
        
        m = matrix_M[0:-1, 1:-1]
        c = matrix_C[1:-1, 1:-1]
        s = matrix_S[1:-1, 1:-1]
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

        ceramic_left = np.sum((c[:,0] - 1)**(power_factor)*constant)
        ceramic_right = np.sum((c[:,-1] - 1)**(power_factor)*constant)
        ceramic_bottom = np.sum((matrix_C[-2, 1:center_start_C] - 1)**
                                (power_factor)*constant) + \
                                    np.sum((matrix_C[-2, center_end_C+1:-1] - 1)
                                           **(power_factor)*constant)

        sink_bottom = np.sum((matrix_S[-2, 1:center_start_S] - 1)**\
            (power_factor)*constant) + np.sum((matrix_S[-2, center_end_S+1:-1]\
                                                - 1)**(power_factor)*constant)
        sink_left = np.sum((s[:,0] - 1)**(power_factor)*constant)
        sink_right = np.sum((s[:,-1] - 1)**(power_factor)*constant)
        sink_top = sink_top_power_out(matrix_S, column_num_F, separation_num, 
                                      step_size, forced = forced,
                                      wind_speed = wind_speed)

        fins = power_out_fins(fin_meshes, step_size, forced = forced, 
                              wind_speed = wind_speed)
        microprocessor_power_out = M_bottom + M_left + M_right + M_top
        overall_power_out = fins + sink_bottom + sink_left + sink_right + \
            sink_top + ceramic_left + ceramic_right + ceramic_bottom + \
                M_bottom + M_left + M_right
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
            print("Microprocessor temperature change from last iteration", f"{temp_difference_percent:.6g}", "%")
            print("Percentage difference from theoretical power out     ", f"{power_difference_percent:.6g}", "%")
            print("Mean temperature in Microprocessor                   ", f"{mean_temp_M[-1]*293:.6g}", "K")
        print(" ")
    return power_difference_percent/100, temp_difference_percent/100

def create_plot(matrix_M, matrix_C, matrix_S, fin_meshes, fin_separation, 
                step_size, grid, length_scale = 0.001):
    separation_num = int(fin_separation / (step_size*length_scale))
    # Extracting central parts of matrices
    M = matrix_M[1:-1, 1:-1]
    C = matrix_C[1:-1, 1:-1]
    S = matrix_S[1:-1, 1:-1]
    fins = [mesh[1:-1, 1:-1] for mesh in fin_meshes.values()]  

    # Determining dimensions
    M_row_num, M_column_num = M.shape
    C_row_num, C_column_num = C.shape
    S_row_num, S_column_num = S.shape
    F_row_num, F_column_num = fins[0].shape

    # Total dimensions for the big mesh
    total_rows = F_row_num + M_row_num + C_row_num + S_row_num
    total_cols = S_column_num

    # Initialize the big mesh
    big_mesh = np.full((total_rows, total_cols), np.nan)

    # Insert fin meshes horizontally
    current_col = 0
    for fin in fins:
        big_mesh[:F_row_num, current_col:current_col + F_column_num] = fin
        current_col += F_column_num + separation_num

    # Function to center-align and insert a matrix
    def insert_centered(big_mesh, matrix, start_row):
        row_num, col_num = matrix.shape
        col_start = (big_mesh.shape[1] - col_num) // 2
        big_mesh[start_row:start_row + row_num, col_start:col_start + 
                 col_num] = matrix

    # Insert other matrices
    insert_centered(big_mesh, S, F_row_num)
    insert_centered(big_mesh, C, F_row_num + S_row_num)
    insert_centered(big_mesh, M, F_row_num + S_row_num + C_row_num)

    heatmap = big_mesh * 293
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))
    im1 = ax1.imshow(heatmap, cmap='turbo', extent = (0,int(len(big_mesh[0])*step_size),0,int(len(big_mesh)*step_size)), interpolation='nearest')
    cbar = plt.colorbar(im1, shrink=0.8)
    cbar.set_label('Temperature Scale (K)')
    ax1.set_xlabel('x, mm')
    ax1.set_ylabel('y, mm')

    if grid == True: 
    # Set up the grid
        ax1.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        ax1.set_xticks(np.arange(-0.5, big_mesh.shape[1], 1), minor=True)
        ax1.set_yticks(np.arange(-0.5, big_mesh.shape[0], 1), minor=True)

        # Customize the grid appearance
        ax1.tick_params(which='minor', size=0)
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))

    plt.savefig('/Users/alanli/Desktop/Comp Phys Proj/Submission_folder/allcomp1.png', dpi=500)
    plt.show()
    return big_mesh



def solve(num_fins, row_num_M, fin_separation, iteration = 50000000, 
                        k_M = 150, k_C = 230, k_S = 250, q = 5e8, 
                        fin_height = 30e-3, ini_temp = 293, 
                        length_scale = 0.001, plot_grid = False, 
                        forced_convection = False, convergence = 0.001, 
                        temp_track = False, show_plot = True, interval = 10000,
                        wind_speed = 20):
    ini_temp = ini_temp/293
    column_num_M = int(14 * row_num_M)
    row_num_C = int(2 * row_num_M)
    column_num_C = int(column_num_M*20/14)
    matrix_M, matrix_C, step_size = create_mesh(column_num_M, row_num_M, 
                                                column_num_C, row_num_C, 
                                                ini_temp)
    matrix_S, fin_meshes, row_num_F, column_num_F, row_num_S, column_num_S = \
        create_heatsink_mesh(step_size = step_size, num_fins = num_fins, 
                             fin_height = fin_height, 
                             fin_separation = fin_separation, 
                             ini_temp = ini_temp)
    step_size = step_size/length_scale
    heat_term = 0.25 * length_scale ** 2 * step_size ** 2 * q / (k_M * 293)
    k_eff_CM = 2*k_M*k_C/(k_M +k_C)
    k_eff_SC = 2*k_S*k_C/(k_S +k_C)
    interface_heat_term = (length_scale ** 2 * step_size ** 2 * q
                           )/((k_eff_CM + 3 * k_M) * 293)
    center_start_C = (len(matrix_C[0]) - len(matrix_M[0])) // 2 + 1
    center_end_C = (len(matrix_C[0]) + len(matrix_M[0])) // 2 - 2
    center_start_S = (len(matrix_S[0]) - len(matrix_C[0])) // 2 + 1
    center_end_S = (len(matrix_S[0]) + len(matrix_C[0])) // 2 - 2

    mean_temp_M = []
    mean_temp_C = []
    mean_temp_S = []

    progress_interval = max(iteration // interval, 1)

    for i in range(iteration):
        update_C(matrix_C, step_size, k_C, length_scale = length_scale, 
                 forced = forced_convection, wind_speed = wind_speed)
        update_M(matrix_M, step_size, k_M, length_scale = length_scale, 
                 forced = forced_convection, wind_speed = wind_speed)
        update_sink(matrix_S, row_num_S, column_num_S, step_size, k_S, 
                    length_scale = length_scale, forced = forced_convection,
                    wind_speed = wind_speed)
        update_fin(fin_meshes, row_num_F, column_num_F, step_size, k_S, 
                   length_scale = length_scale, forced = forced_convection, 
                   wind_speed = wind_speed)

        power_frac_diff, temp_frac_diff = progress_update\
            (iter_num = i, progress_interval = progress_interval, 
             iteration = iteration, step_size = step_size, 
            k_eff_CM = k_eff_CM, matrix_M = matrix_M, 
            matrix_C = matrix_C, matrix_S = matrix_S, 
            fin_meshes = fin_meshes, 
            center_start_C = center_start_C, 
            center_end_C = center_end_C, 
            center_start_S = center_start_S, 
            center_end_S = center_end_S, 
            column_num_F = column_num_F, 
            separation_num = int(fin_separation/(step_size*length_scale)), 
            mean_temp_M = mean_temp_M, 
            forced = forced_convection,
            wind_speed = wind_speed)

        # replace ghost points at CM and SC interface
        update_interface(small_matrix = matrix_M, big_matrix = matrix_C, 
                         center_start=center_start_C, center_end=center_end_C)
        update_interface(small_matrix = matrix_C, big_matrix = matrix_S, 
                         center_start=center_start_S, center_end=center_end_S)
        pair_sink_fin(step_size = step_size, dictionary = fin_meshes, 
                      matrix = matrix_S, num_fins = num_fins, 
                      fin_separation = fin_separation, 
                      column_num_F = column_num_F, 
                      length_scale = length_scale)

        buffer_M = np.zeros_like(matrix_M)
        buffer_C = np.zeros_like(matrix_C)
        buffer_S = np.zeros_like(matrix_S)

        # Apply the stencil to the internal matrix_M cells, except second row
        buffer_M[2:-1, 1:-1] = 0.25 * (matrix_M[3:, 1:-1] + 
                                       matrix_M[1:-2, 1:-1] + 
                                       matrix_M[2:-1, 2:] + 
                                       matrix_M[2:-1, :-2]) + heat_term

        # Apply the stencil to the internal matrix_C cells, except the second 
        # and second last row
        buffer_C[2:-2, 1:-1] = 0.25 * (matrix_C[3:-1, 1:-1] + 
                                       matrix_C[1:-3, 1:-1] + matrix_C[2:-2, 2:]
                                         + matrix_C[2:-2, :-2])
        
        # Apply the stencil to the heatsink, except the second last row
        buffer_S[1:-2, 1:-1] = 0.25 * (matrix_S[2:-1, 1:-1] + 
                                       matrix_S[:-3, 1:-1] + matrix_S[1:-2, 2:]
                                         + matrix_S[1:-2, :-2])
        for i in range(num_fins):
            fin_meshes[i][1:-1, 1:-1] = 0.25 * (fin_meshes[i][2:, 1:-1] + 
                                                fin_meshes[i][:-2, 1:-1] + 
                                                fin_meshes[i][1:-1, 2:] + 
                                                fin_meshes[i][1:-1, :-2])
            
        # Apply the stencil to specific parts of the second last row of 
        # matrix_C and second last row of matrix_S
        bottom_edge_stencil(buffer = buffer_C, center_start = center_start_C, 
                            center_end = center_end_C, matrix = matrix_C)
        bottom_edge_stencil(buffer = buffer_S, center_start = center_start_S, 
                            center_end = center_end_S, matrix = matrix_S)

        #Apply the stencil to second row of matrix_M with special stencil
        buffer_M[1, 1:-1] = 1/(k_eff_CM + 3 * k_M) * (k_eff_CM * 
                                                      matrix_M[0, 1:-1] + 
                                                      k_M * matrix_M[2, 1:-1] + 
                                                      k_M * matrix_M[1, 2:] + 
                                                      k_M * matrix_M[1, :-2]) +\
                                                        interface_heat_term

        #Apply the stencil to second and second last row of matrix_C with 
        # special stencil
        buffer_C[1, 1:-1] = 1/(k_eff_SC + 3 * k_C) * \
            (k_eff_SC * matrix_C[0, 1:-1] + k_C * matrix_C[2, 1:-1] + k_C * \
             matrix_C[1, 2:] + k_C * matrix_C[1, :-2])
        buffer_C[-2, center_start_C:center_end_C+1] = 1/(k_eff_CM + 3 * k_C) *\
              (k_eff_CM * matrix_C[-1, center_start_C:center_end_C+1] + k_C *\
                matrix_C[-3, center_start_C:center_end_C+1] + k_C *\
                      matrix_C[-2, center_start_C+1:center_end_C+2] + k_C *\
                          matrix_C[-2, center_start_C-1:center_end_C])

        #Apply the stencil to second last row of matrix_S with special stencil
        buffer_S[-2, center_start_S:center_end_S+1] = 1/(k_eff_SC + 3 * k_S) *\
              (k_eff_SC * matrix_S[-1, center_start_S:center_end_S+1] + k_S *\
               matrix_S[-3, center_start_S:center_end_S+1] + k_S *\
                  matrix_S[-2, center_start_S+1:center_end_S+2] + k_S *\
                      matrix_S[-2, center_start_S-1:center_end_S])

        matrix_C[1:-1, 1:-1] = buffer_C[1:-1, 1:-1]
        matrix_M[1:-1, 1:-1] = buffer_M[1:-1, 1:-1]
        matrix_S[1:-1, 1:-1] = buffer_S[1:-1, 1:-1]
        # Calculate and store the mean temperature
        mean_temp_M.append(np.mean(matrix_M[1:-1, 1:-1]))
        mean_temp_C.append(np.mean(matrix_C[1:-1, 1:-1]))
        mean_temp_S.append(np.mean(matrix_S[1:-1, 1:-1]))
        if np.abs(power_frac_diff) < convergence and np.abs(temp_frac_diff) < convergence:
            print("Convergence Criteria Met!")
            break
    if show_plot == True:
        create_plot(matrix_M, matrix_C, matrix_S, fin_meshes, fin_separation, 
                    step_size, grid = plot_grid)
    if temp_track == True:
        return mean_temp_M, mean_temp_C, mean_temp_S, step_size, i
    else:
        return step_size