isolated_filename = './data/paca_v2/isolated_customers.txt'
customers = './data/paca_v2/customers.txt'
locations = './data/paca_v2/locations.txt'
distances = './data/paca_v2/distances.txt'
output = './data/paca_v3/dist_matrix.txt'

# Load isolated customers
f_isol = open(isolated_filename, "r")
isolated = set()
while True:
    isol = f_isol.readline()
    if not isol:
        break
    isolated.add(int(isol))

# Create distance matrix
f_cust = open(customers, "r")
f_loc = open(locations, "r")
f_dist = open(distances, "r")
f_out = open(output, "w")

f_out.write("customer location distance\n")

while True:
    loc = f_loc.readline()
    if not loc:
        break
    cust = int(f_cust.readline())
    dist = float(f_dist.readline())
    if not cust in isolated:
        string = str(cust) + " " + loc.rstrip() + " " + str("%.2f" % round(dist, 2)) + "\n"
#        print(string)
        f_out.write(string)

