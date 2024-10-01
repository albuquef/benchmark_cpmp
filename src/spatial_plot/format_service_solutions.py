# Define the filenames
sol_file = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton.txt'
sol_file = 'data/solutions_service_cinema/test_paca_cinema_commune_p_192_EXACT_CPMP_cover_commune.txt'
sol_file = 'data/solutions_service_cinema/test_paca_cinema_EPCI_p_192_EXACT_CPMP_cover_EPCI.txt'
sol_file = 'data/solutions_service_cinema/test_paca_cinema_p_192_EXACT_CPMP.txt'

cust_file = '/home/falbuquerque/Documents/projects/GeoAvignon/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
#sol file without the extension
sol_file_name = sol_file.split('.')[0]
output_file = sol_file_name + '_table.txt'

# Step 1: Read the selected locations from sol.txt
with open(sol_file, 'r') as f:
    lines = f.readlines()

# Extract "P LOCATIONS" section
p_locations = []
in_p_locations = False

for line in lines:
    line = line.strip()
    if line == "P LOCATIONS":
        in_p_locations = True
        continue
    if in_p_locations:
        if line.startswith("LOCATION USAGES"):
            break  # Stop when we reach the next section
        if line:  # If line is not empty
            p_locations.append(line)

# Step 2: Read customers from cust.txt and filter based on selected locations
selected_customers = []

with open(cust_file, 'r') as f:
    header = f.readline()  # Skip the header line
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:  # Ensure we have enough parts
            id_cust = parts[0]  # ID of the customer
            coord_x = float(parts[2])  # coord_x
            coord_y = float(parts[3])  # coord_y

            # Check if the customer ID matches one of the selected locations
            if id_cust in p_locations:
                selected_customers.append(line.strip())

# Step 3: Write the filtered customers to a new text file with header
with open(output_file, 'w') as f:
    # Write the header with comma separation
    f.write("customer,weight,coord_x,coord_y,id_cust,fid\n")
    for customer in selected_customers:
        # Write each customer, replacing spaces with commas
        f.write(",".join(customer.split()) + '\n')

print(f'Filtered customers written to {output_file}')