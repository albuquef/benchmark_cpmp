# capacities = './provided/Matrice_600_1500_avec_origines_et_destinations_avec_capacités/Points_destination_matrice_600_1500_avec_capacités_par_services.csv'
# output = './data/paca_v4/capacities_urg.txt'
# column = 8

f_cap = open(capacities, "r")
f_out = open(output, "w")

line = f_cap.readline()
f_out.write("location capacity\n")

while True:
    line = f_cap.readline()
    if not line:
        break
    numbers = line.split(';')
    id = int(numbers[0])
    capacity = int(numbers[column].split(',')[0])
    string = str(id) + " " + str(capacity) + "\n"
    f_out.write(string)
