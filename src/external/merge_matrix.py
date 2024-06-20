# customers = './provided/Matrice_cplt/matdistpaca_600_1500_cplt_origin.txt'
# locations = './provided/Matrice_cplt/matdistpaca_600_1500_cplt_destination.txt'
# distances = './provided/Matrice_cplt/matdistpaca_600_1500_cplt_tps.txt'
# output = './data/paca_v4/dist_matrix.txt'


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
    string = str(cust) + " " + loc.rstrip() + " " + str("%.2f" % round(dist, 2)) + "\n"
    f_out.write(string)
