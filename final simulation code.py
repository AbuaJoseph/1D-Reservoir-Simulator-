# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 15:54:50 2023

@author: ukehj
"""

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import csv


# Number of grid blocks
n = int(input("Enter the number of grid blocks: "))
# Well parameters
well_index = int((input("Enter the well index of the well with the gridblock:")).strip())  # Grid block with the well (0-based index)
well_rate = float((input("Enter the well rate of the gridblock with the well:")).strip())  # Production rate (negative for production)

# Grid block dimensions
dx = float(input("Enter the length of a grid block (dx) in feet: "))
dy = float(input("Enter the width of a grid block (dy) in feet: "))
dz = float(input("Enter the height of a grid block (dz) in feet: "))

# Reservoir and fluid properties
viscosity = float(input("Enter the fluid viscosity (cp): "))
porosity = float(input("Enter the porosity: "))
k = float(input("permeability:"))
FVF = float(input("FVF:"))
comp_total = float(input("Enter the compressibility (psi^-1): "))

# Initialize arrays to store pressure data
pressure_data = []


# Initial reservoir pressure
initial_pressure = np.ones(n) * float(input("Enter the initial reservoir pressure (psi): "))
boundary_condition = input("Choose the type of boundary condition (Dirichlet, Neumann, No flow, Mixed(left_boundary_pressure, right_pressure_gradient), Mixed(left_pressure_gradient, right_boundary_pressure )): ").strip()

# Time parameters
num_time_steps = int(input("Enter the number of time steps: "))
delta_t = float(input("Enter the time step size (days): "))
if boundary_condition == "Dirichlet":
    left_boundary_pressure = float(input("Enter the pressure at the left boundary (psi): "))
    right_boundary_pressure = float(input("Enter the pressure at the right boundary (psi): "))
elif boundary_condition == "Neumann":
    left_pressure_gradient = float(input("Enter the left boundary pressure gradient (psi/ft): "))
    right_pressure_gradient = float(input("Enter the right boundary pressure gradient (psi/ft): "))
elif boundary_condition ==  "Mixed(left_boundary_pressure, right_pressure_gradient)":
    left_boundary_pressure = float(input("Enter the pressure at the left boundary (psi): "))
    right_pressure_gradient = float(input("Enter the right boundary pressure gradient (psi/ft): "))
elif boundary_condition == "Mixed(left_pressure_gradient, right_boundary_pressure )":
    left_pressure_gradient = float(input("Enter the left boundary pressure gradient (psi/ft): "))
    right_pressure_gradient = float(input("Enter the right boundary pressure gradient (psi/ft): "))
    
    

# Create arrays to store the values
lower_diag = np.zeros(n - 1)
main_diag = np.zeros(n)
upper_diag = np.zeros(n - 1)
# Calculate transmissibility and accumulation terms
transmissibility = (0.001127*k * dy*dz) / (viscosity * FVF * dx)
accumulation = (dx * dy * dz * comp_total) * porosity / (5.615 * FVF * delta_t)
# Time-stepping loop
for step in range(num_time_steps):
    # Populate the lower and upper diagonals with transmissibility
    lower_diag.fill(transmissibility)
    upper_diag.fill(transmissibility)
    # populate the interior gridblocks of the main diagonal
    main_diag[1:-1].fill(-accumulation-2*transmissibility)
    

    
    # initialize matrix B
    B = np.zeros(n)
    
    # Handle user-specified boundary conditions
    if boundary_condition == "Dirichlet":
        main_diag[0] = (-3*transmissibility - accumulation)
        main_diag[-1] = (-3*transmissibility - accumulation)
        
    elif boundary_condition == "Neumann":
        main_diag[0] = -(transmissibility + accumulation)
        main_diag[-1] = -(transmissibility + accumulation)
    
    elif boundary_condition ==  "Mixed(left_boundary_pressure, right_pressure_gradient)":
        main_diag[0] =  (-3*transmissibility - accumulation)
        main_diag[-1] = -(transmissibility + accumulation)
        
    elif boundary_condition == "Mixed(left_pressure_gradient, right_boundary_pressure )":
        main_diag[0] = -(transmissibility + accumulation)
        main_diag[-1] = (-3*transmissibility - accumulation)
    
        
        
    else:    # no flow
        main_diag[0] = -(transmissibility + accumulation)
        main_diag[-1] = -(transmissibility + accumulation)
        
    # Handle boundary conditions and well term
    if well_index == 0:
        # Left boundary with well
        if boundary_condition == "Dirichlet":
            B[0] = (-accumulation * initial_pressure[0]) - (2 * transmissibility * left_boundary_pressure) - well_rate
        elif boundary_condition == "Neumann":
            B[0] = (-accumulation * initial_pressure[0]) + (transmissibility * dx * left_pressure_gradient) - well_rate
        elif boundary_condition ==  "Mixed(left_boundary_pressure, right_pressure_gradient)":
            B[0] = (-accumulation * initial_pressure[0]) - (2 * transmissibility * left_boundary_pressure) - well_rate
        elif boundary_condition == "Mixed(left_pressure_gradient, right_boundary_pressure )":
            B[0] = (-accumulation * initial_pressure[0]) + (transmissibility * dx * left_pressure_gradient) - well_rate
            
        elif boundary_condition == "No flow":
             B[0] = (-accumulation * initial_pressure[0])  - well_rate
    else: # left boundary without well
        if boundary_condition == "Dirichlet":
            B[0] =(-accumulation * initial_pressure[0]) - (2 * transmissibility * left_boundary_pressure)
        elif boundary_condition == "Neumann":
            B[0] = (-accumulation * initial_pressure[0]) + (transmissibility * dx * left_pressure_gradient)
        elif boundary_condition ==  "Mixed(left_boundary_pressure, right_pressure_gradient)":
            B[0] = (-accumulation * initial_pressure[0]) - (2 * transmissibility * left_boundary_pressure)
        elif boundary_condition == "Mixed(left_pressure_gradient, right_boundary_pressure )":
            B[0] = (-accumulation * initial_pressure[0]) + (transmissibility * dx * left_pressure_gradient)
        elif boundary_condition == "No flow":
            B[0] = (-accumulation * initial_pressure[0]) 
    if well_index == n - 1:
         # Right boundary with well
        if boundary_condition == "Dirichlet":
            B[-1] = (-accumulation * initial_pressure[-1]) - (2 * transmissibility * right_boundary_pressure) - well_rate
        elif boundary_condition == "Neumann":
            B[-1] = (-accumulation * initial_pressure[-1]) - (transmissibility * dx * right_pressure_gradient) - well_rate
        elif boundary_condition ==  "Mixed(left_boundary_pressure, right_pressure_gradient)":
            B[-1] = (-accumulation * initial_pressure[-1]) - (transmissibility * dx * right_pressure_gradient) - well_rate
        elif boundary_condition == "Mixed(left_pressure_gradient, right_boundary_pressure )":
            B[-1] = (-accumulation * initial_pressure[-1]) - (2 * transmissibility * right_boundary_pressure) - well_rate
        elif boundary_condition == "No flow":
            B[-1] = (-accumulation * initial_pressure[-1])  - well_rate
    else:
        
        # Right boundary without well
        if boundary_condition == "Dirichlet":
            B[-1] = (-accumulation * initial_pressure[-1]) - (2 * transmissibility * right_boundary_pressure)
        elif boundary_condition == "Neumann":
            B[-1] = (-accumulation * initial_pressure[-1]) - (transmissibility * dx * right_pressure_gradient)
        elif boundary_condition ==  "Mixed(left_boundary_pressure, right_pressure_gradient)":
            B[-1] = (-accumulation * initial_pressure[-1]) - (transmissibility * dx * right_pressure_gradient)
        elif boundary_condition == "Mixed(left_pressure_gradient, right_boundary_pressure )":
            B[-1] =  (-accumulation * initial_pressure[-1]) - (2 * transmissibility * right_boundary_pressure)
        elif boundary_condition == "No flow":
            B[-1] = (-accumulation * initial_pressure[-1])
    
    # Interior grid blocks (non-well terms)
    B[1:-1] = -accumulation * initial_pressure[1:-1]

    # Add the well term to the specific grid block
    B[well_index] = -accumulation * initial_pressure[well_index] - well_rate

    
        
        
    
                
    
  
    
    # Create the tridiagonal matrix A
    A = diags([lower_diag, main_diag, upper_diag], [-1, 0, 1], shape=(n, n), format='csr')
    
    
    # update matrix B
    B_0 = np.array([B[0]])
    B_interior = B[1:-1]
    B_n_minus_1 = np.array([B[-1]])
    B_combined = np.concatenate((B_0, B_interior, B_n_minus_1)).reshape(-1, 1)
    
    
    # Solve for pressure using implicit formulation: Ap = B
    pressure = spsolve(A, B_combined)

    # Update initial pressure for the next time step
    initial_pressure = pressure

    # Print the pressure distribution at each time step
    print(f"Time Step {step + 1}: Pressure = {pressure}")

    # Final pressure distribution
    #print("Final Pressure Distribution:")
    #print(pressure)

 # Store pressure data for this time step
    pressure_data.append(pressure)

# Convert pressure_data to a NumPy array for contour plot
pressure_data = np.array(pressure_data)

# Create a contour plot
plt.figure(figsize=(8, 6))
contour = plt.contourf(pressure_data, cmap='viridis', levels=20)
plt.colorbar(contour, label='Pressure (psi)')
plt.xlabel('Grid Block')
plt.ylabel('Time Step')
plt.title('Pressure Distribution in 1D Reservoir')
plt.show()
# Define the output CSV file name
output_file = "pressure_distribution.csv"

# Write pressure data to CSV
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write header row with grid block indices
    header_row = ["Time Step"] + [f"Grid Block {i+1}" for i in range(n)]
    writer.writerow(header_row)

    # Write pressure data for each time step
    for time_step, pressure_values in enumerate(pressure_data):
        data_row = [time_step] + pressure_values.tolist()
        writer.writerow(data_row)

print(f"Pressure distribution data has been saved to {output_file}")