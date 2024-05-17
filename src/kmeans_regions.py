import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import os
import re



# Function to read data from CSV
def read_location_identifiers(txt_file):
    location_identifiers = {}
    with open(txt_file, 'r') as file:
        next(file)  # Skip the header
        for line in file:
            id_, identif = line.strip().split()
            location_identifiers[identif] = int(id_)
    return location_identifiers

def read_csv(csv_file, location_identifiers):
    locations = []
    coordinates = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            identif = row['identif']
            if identif in location_identifiers:
                x = float(row['x'])
                y = float(row['y'])
                locations.append(location_identifiers[identif])
                coordinates.append((x, y))
    return locations, coordinates

def read_location_identifiers(txt_file):
    location_identifiers = {}
    with open(txt_file, 'r') as file:
        next(file)  # Skip the header
        for line in file:
            id_, identif = line.strip().split()
            location_identifiers[identif] = int(id_)
    return location_identifiers

# Function to read data from a text file
def read_data_txt(input_txt_file):
    locations = []
    data = []
    with open(input_txt_file, 'r') as txtfile:
        next(txtfile)  # Skip the header line
        for line in txtfile:
            parts = line.split()
            location = int(parts[0])
            coord_x = float(parts[2])
            coord_y = float(parts[3])
            locations.append(location)
            data.append([coord_x, coord_y])
    return np.array(locations), np.array(data)

# Function to perform K-means clustering
def perform_kmeans(data, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(data)
    return kmeans

# Function to plot clusters and boundaries
def plot_clusters(data, kmeans):
    centroids = kmeans.cluster_centers_
    labels = kmeans.labels_

    # Define grid parameters
    x_min, x_max = data[:, 0].min() - 1, data[:, 0].max() + 1
    y_min, y_max = data[:, 1].min() - 1, data[:, 1].max() + 1

    # Adjust step size to reduce memory usage
    x_step = (x_max - x_min) / 500
    y_step = (y_max - y_min) / 500

    # Create a grid of points
    xx, yy = np.meshgrid(np.arange(x_min, x_max, x_step),
                         np.arange(y_min, y_max, y_step))

    # Predict the labels for each point in the grid
    Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()])

    # Reshape the predictions to match the grid shape
    Z = Z.reshape(xx.shape)

    # Plot the decision boundary
    plt.contour(xx, yy, Z, levels=np.arange(0, kmeans.n_clusters + 1) - 0.5, colors='k', linestyles='--')
    # Plot the data points
    plt.scatter(data[:, 0], data[:, 1], c=labels, cmap='viridis', s=20)
    # Plot the centroids
    # plt.scatter(centroids[:, 0], centroids[:, 1], marker='x', s=200, color='red')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('K-Means Clustering with Decision Boundaries')
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.show()

# Function to save points with cluster labels to CSV
def save_points_with_clusters(data, labels, output_csv_file):
    with open(output_csv_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["X", "Y", "Cluster_Label"])
        for point, label in zip(data, labels):
            writer.writerow([point[0], point[1], label])

# Function to save locations with cluster labels and coordinates to a text file
def save_locations_with_clusters(locations, data, labels, output_txt_file):
    with open(output_txt_file, 'w') as txtfile:
        writer = csv.writer(txtfile, delimiter=' ')
        writer.writerow(["location", "subarea", "coor_x", "coord_y"])
        for location, coord, label in zip(locations, data, labels):
            writer.writerow([location, (label+1), coord[0], coord[1]])


# Main function to run the entire process
def main():

    # # input_csv_file = "./data/PACA/points_coord.csv"
    # # locations, data = read_data(input_csv_file)
    # group_number = 3
    # directory = f'./data/Literature/group{group_number}/'

    # # List all files in the directory
    # files = os.listdir(directory)

    # # Filter out only files (excluding directories) and remove extensions
    # file_names = [os.path.splitext(file)[0] for file in files if os.path.isfile(os.path.join(directory, file)) and not (file.startswith('loc_') or file.startswith('cust_') or file.startswith('dist_'))]

    # # Output the list of filenames without extensions
    # print(file_names)
    # print(len(file_names))

    # # Pattern to match digits after underscore at the end of the string
    # pattern = r'_(\d+)$'

    # # Loop through the filenames
    # for filename in file_names:
    #     # Extract value from the filename
    #     match = re.search(pattern, filename)
    #     if match:
    #         value_with_zeros = match.group(1)
    #         # Remove leading zeros
    #         value = value_with_zeros.lstrip('0')
    #         # Output the extracted value
    #         print(value)
            
    #         txt_file = filename
    #         input_txt_file = f"{directory}loc_capacities_{txt_file}.txt"
    #         n_clusters = int(value)  # Change this to the desired number of clusters
    #         locations, data = read_data_txt(input_txt_file)
            
    #         # Perform K-means clustering
    #         kmeans = perform_kmeans(data, n_clusters)
    #         # Plot clusters and boundaries
    #         # plot_clusters(data, kmeans)
            
    #         output_txt_file = f"{directory}loc_coverages_kmeans_{txt_file}.txt"
    #         save_locations_with_clusters(locations, data, kmeans.labels_, output_txt_file)


    #         print(f"saved: loc_coverages_kmeans_{txt_file}")
    # exit()
    
    
    group_number = 3
    txt_file='p3038_600'
    # group_number = 5
    # txt_file='rl1304_010'
    # txt_file='pr2392_020'
    # txt_file='fnl4461_0100'
    
    PATH_DATA=f'./data/Literature/group{group_number}/'
    input_txt_file = f"{PATH_DATA}loc_capacities_{txt_file}.txt"
    
    n_clusters = 18  # Change this to the desired number of clusters

    locations, data = read_data_txt(input_txt_file)
    
    PATH_DATA="./data/PACA/"
    input_csv_file = "./data/PACA/points_coord.csv"
    location_identifiers = read_location_identifiers(f"{PATH_DATA}map_id_cust_loc.txt")
    locations, data = read_csv(input_csv_file, location_identifiers)
    txt_file='arrond'
    
    # Perform K-means clustering
    kmeans = perform_kmeans(data, n_clusters)

    # Plot clusters and boundaries
    # plot_clusters(data, kmeans)

    # # Save points with cluster labels to CSV
    # output_csv_file_points = "points_with_clusters.csv"
    # save_points_with_clusters(data, kmeans.labels_, output_csv_file_points)

    output_txt_file = f"{PATH_DATA}loc_coverages_kmeans_{txt_file}.txt"    
    save_locations_with_clusters(locations, data, kmeans.labels_, output_txt_file)


    print(f"saved: loc_coverages_kmeans_{txt_file}")

# Run the main function
if __name__ == "__main__":
    main()