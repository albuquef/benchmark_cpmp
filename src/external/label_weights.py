# weights = './provided/Matrice_cplt/matdistpaca_600_1500_pts_origines.txt'
# output = './data/paca_v4/cust_weights.txt'

f_weights = open(weights, "r")
f_out = open(output, "w")

f_out.write("customer weight\n")

cust = 1
while True:
    weight = f_weights.readline()
    if not weight:
        break
    weight = int(float(weight.rstrip()))
    string = str(cust) + " " + str(weight) + "\n"
    f_out.write(string)
    cust += 1
