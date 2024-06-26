import os
import math
import matplotlib.pyplot as plt


def format_instance(PATH, filename, data):
    output_filename = f"loc_capacities_{filename}"
    with open(os.path.join(PATH, output_filename), "w") as output_file:
        output_file.write("location capacity coord_x coord_y\n")
        for i, point in enumerate(data):
            output_file.write(f"{i+1} {point['capacity']} {point['x']} {point['y']}\n")
    print(f"File {output_filename} created successfully.")

    output_filename = f"cust_weights_{filename}"
    with open(os.path.join(PATH, output_filename), "w") as output_file:
        output_file.write("customer weight coord_x coord_y\n")
        for i, point in enumerate(data):
            output_file.write(f"{i+1} {point['demand']} {point['x']} {point['y']}\n")
    print(f"File {output_filename} created successfully.")


def format_distance_matrix(PATH, filename, data):

    skip_filenames = [
        "spain737_74_1.txt",
        "spain737_74_2.txt",
        "spain737_148_1.txt",
        "spain737_148_2.txt"
    ]
    
    # spain is not euclidian distance
    if filename in skip_filenames:
        print(f"Skipping {filename} as it is not euclidian distance.")
        return

    distances = []
    for i in range(len(data)):
        for j in range(i+1, len(data)):
            distance = euclidean_distance(data[i], data[j])
            distances.append((i+1, j+1, distance))

    output_filename = f"dist_matrix_{filename}"
    with open(os.path.join(PATH, output_filename), "w") as output_file:
        output_file.write("loc1 loc2 dist\n")
        for loc1, loc2, dist in distances:
            output_file.write(f"{loc1} {loc2} {dist}\n")
    print(f"File {output_filename} created successfully.")


def format_spain_dist_matrix(path,filename):

    with open(path+filename, 'r') as file:
        lines = file.readlines()

    # Convert the lines into a list of integers
    numbers = [int(line.strip()) for line in lines]

    # Determine the size of the square matrix
    size = int(len(numbers) ** 0.5)

    # Check if the number of elements forms a perfect square
    if size * size != len(numbers):
        raise ValueError("The number of elements is not a perfect square, cannot form a square matrix.")

    # Reshape the list into a square matrix
    matrix = [numbers[i*size:(i+1)*size] for i in range(size)]

    # Create output filename and path with .txt extension
    output_filename = f"dist_matrix_{os.path.basename(filename).replace('.grd', '.txt')}"
    output_filepath = os.path.join(path, output_filename)

    # Write the distances to the output file
    with open(output_filepath, "w") as output_file:
        output_file.write("customer location distance\n")
        for i in range(size):
            for j in range(size):
                loc1, loc2, dist = i + 1, j + 1, matrix[i][j]
                output_file.write(f"{loc1} {loc2} {dist}\n")

    print(f"File {output_filename} created successfully.")







def plot_points(x_coordinates, y_coordinates):
    plt.figure(figsize=(8, 6))
    plt.scatter(x_coordinates, y_coordinates, color='blue', label='Points')
    plt.title('Plot of Points with Coordinates')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.legend()
    plt.grid(True)
    plt.show()

def euclidean_distance(point1, point2):
    return math.sqrt((point1['x'] - point2['x'])*2 + (point1['y'] - point2['y'])*2)




# PATH = "./data/group1/"
# # Generating filenames_g1 using a loop
# filenames = [f"cpmp{i:02}.txt" for i in range(1, 21)]

PATH = "./data/Literature/group4/"
# Generating filenames_g1 using a loop
# filenames = ["spain737_74_1.txt","spain737_74_2.txt","spain737_148_1.txt","spain737_148_2.txt"]
filenames = ["spain737_74_1.txt"]

# PATH = "./data/Literature/GB21/"
# filenames = ["XMC10150_100.txt","XMC10150_1000.txt", "XMC10150_2000.txt","SRA104814_100.txt", "SRA104814_1000.txt", "SRA104814_2000.txt",\
#                 "FNA52057_100.txt", "FNA52057_1000.txt", "FNA52057_2000.txt","LRA498378_100.txt", "LRA498378_1000.txt", "LRA498378_2000.txt"]
# filenames_g6 = ["SRA104814_100.txt", "SRA104814_1000", "SRA104814_2000.txt"]
# filenames_g6 = ["FNA52057_100.txt", "FNA52057_1000.txt", "FNA52057_2000.txt"]
# filenames_g6 = ["LRA498378_100.txt", "LRA498378_1000.txt", "LRA498378_2000.txt"]

# Loop through each filename
for filename in filenames:
    with open(os.path.join(PATH, filename), "r") as file:
        lines = file.readlines()

    # Extracting the number of points and number of medians from the first line
    num_points, num_medians = map(int, lines[0].split())

    # Extracting the coordinates, capacities, and demands for each point
    data = []
    for line in lines[1:]:
        values = list(map(float, line.split()))
        point_data = {
            "x": values[0],
            "y": values[1],
            "capacity": values[2],
            "demand": values[3]
        }
        data.append(point_data)

    # Printing the extracted data
    print("File:", filename)
    print("Number of Points:", num_points)
    print("Number of Medians:", num_medians)
    # format_instance(PATH, filename, data)
    # format_spain_dist_matrix(PATH,filename)
    # format_distance_matrix(PATH, filename, data)

    PATH = "./data/Literature/group4/"  
    filename = "spain.grd"
    format_spain_dist_matrix(PATH,filename)

    # Unpacking the data into separate lists for x and y coordinates
    x_coordinates = [point["x"] for point in data]
    y_coordinates = [point["y"] for point in data]
    plot_points(x_coordinates, y_coordinates)