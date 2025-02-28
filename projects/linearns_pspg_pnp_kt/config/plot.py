import matplotlib.pyplot as plt

# Lists to store data
time = []
vx_error = []
vy_error = []
p_error = []
phi_error = []
c1_error = []
c2_error = []

# Read the data from the text file
with open('Error.dat', 'r') as file:
    for line in file:
        # Split the line into values
        values = line.split()
        # Convert values to float and append to respective lists
        time.append(float(values[0]))
        vx_error.append(float(values[1]))
        vy_error.append(float(values[2]))
        p_error.append(float(values[3]))
        phi_error.append(float(values[4]))
        c1_error.append(float(values[5]))
        c2_error.append(float(values[6]))

# Create a single plot with legends
plt.figure(figsize=(12, 8))

# Set x-axis to logarithmic scale
plt.xscale('log')

# Plot vx_error versus time with legend
plt.plot(time, vx_error, label='vx_error')
# Plot vy_error versus time with legend
plt.plot(time, vy_error, label='vy_error')
# Plot p_error versus time with legend
plt.plot(time, p_error, label='p_error')
# Plot phi_error versus time with legend
plt.plot(time, phi_error, label='phi_error')
# Plot c1_error versus time with legend
plt.plot(time, c1_error, label='c1_error')
# Plot c2_error versus time with legend
plt.plot(time, c2_error, label='c2_error')

# Add legends
plt.legend()

# Set plot titles and labels
plt.title('Errors vs. Time (Logarithmic Scale)')
plt.xlabel('Time (log scale)')
plt.ylabel('Error Values')

# Show the plot
plt.show()
