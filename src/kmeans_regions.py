import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import math
from scipy.spatial.distance import cdist
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
    return np.array(locations), np.array(coordinates)

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
def plot_clusters(data, kmeans, title='K-Means Clustering', save_path='clusters_plot'):
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

    plt.xlabel('coord X')
    plt.ylabel('coord Y')
    plt.title(title)  
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Save the figure
    # plt.savefig('./plots/'+save_path+'_HQ.pdf', format='pdf', dpi=1200)
    plt.savefig(save_path+'.pdf', format='pdf')
    plt.savefig(save_path+'_Q.png', format='png', dpi=1200)
    plt.savefig(save_path+'.png', format='png', dpi=300)
    # plt.savefig('./plots/'+save_path+'.svg', format='svg', dpi=1200)
    plt.show()




# Function to create the grid
def create_grid(min_x, max_x, min_y, max_y, k):
    # Calculate the number of rows and columns needed
    cols = math.ceil(math.sqrt(k))
    rows = math.ceil(k / cols)
    
    # Calculate the step size for the grid
    step_x = (max_x - min_x) / cols
    step_y = (max_y - min_y) / rows
    
    # Generate grid points
    x_points = np.arange(min_x, max_x, step_x)
    y_points = np.arange(min_y, max_y, step_y)
    
    grid = []
    
    for i in range(rows):
        for j in range(cols):
            # Create a cell with (x_min, x_max, y_min, y_max)
            x_min = min_x + j * step_x
            x_max = x_min + step_x
            y_min = min_y + i * step_y
            y_max = y_min + step_y
            
            # Ensure the last cell aligns with the max boundaries
            if x_max > max_x:
                x_max = max_x
            if y_max > max_y:
                y_max = max_y
            
            grid.append(((x_min, x_max), (y_min, y_max)))
    
    return grid

def create_varying_grid(min_x, max_x, min_y, max_y, data, k):
    grid = []
    total_area = (max_x - min_x) * (max_y - min_y)
    target_cell_area = total_area / k

    # Initialize the queue with the bounding box
    queue = [((min_x, max_x), (min_y, max_y))]

    while len(queue) > 0:
        # Pop the top bounding box from the queue
        (x_min, x_max), (y_min, y_max) = queue.pop(0)
        
        # Find the point closest to the center of the bounding box
        cell_center = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2])
        distances = cdist([cell_center], data).squeeze()
        closest_point_idx = np.argmin(distances)
        closest_point = data[closest_point_idx]

        # Check if the point is already assigned to a grid cell
        assigned = False
        for cell in grid:
            if x_min <= closest_point[0] < x_max and y_min <= closest_point[1] < y_max:
                assigned = True
                break
        
        if not assigned:
            # Add the grid cell to the grid
            grid.append(((x_min, x_max), (y_min, y_max)))

            # Subdivide the bounding box into four smaller cells
            x_mid = (x_min + x_max) / 2
            y_mid = (y_min + y_max) / 2
            queue.append(((x_min, x_mid), (y_min, y_mid)))  # Top-left
            queue.append(((x_mid, x_max), (y_min, y_mid)))  # Top-right
            queue.append(((x_min, x_mid), (y_mid, y_max)))  # Bottom-left
            queue.append(((x_mid, x_max), (y_mid, y_max)))  # Bottom-right

    return grid
# Function to assign cluster labels to points based on grid cells
def assign_cluster_labels(grid, data):
    labels = np.zeros(len(data), dtype=int)
    
    for idx, (x, y) in enumerate(data):
        for cluster_label, ((x_min, x_max), (y_min, y_max)) in enumerate(grid):
            if x_min <= x < x_max and y_min <= y < y_max:
                labels[idx] = cluster_label
                break
    
    return labels

# Function to plot the grid and points
def plot_grid_and_points(grid, data, labels):
    fig, ax = plt.subplots()
    
    # Generate a colormap with a unique color for each grid cell
    cmap = plt.cm.get_cmap('tab20', len(grid))
    colors = cmap(range(len(grid)))
    
    # Plot the grid with different colors for each cell
    for idx, cell in enumerate(grid):
        (x_min, x_max), (y_min, y_max) = cell
        rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, color=colors[idx], alpha=0.5)
        ax.add_patch(rect)
    
    # Plot the points with their cluster labels
    scatter = ax.scatter(data[:, 0], data[:, 1], c=labels, cmap=cmap, edgecolor='black')
    legend1 = ax.legend(*scatter.legend_elements(), title="Clusters")
    ax.add_artist(legend1)
    
    plt.xlim(grid[0][0][0], grid[-1][0][1])
    plt.ylim(grid[0][1][0], grid[-1][1][1])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

# # Function to plot the grid and points
# def plot_grid_and_points(grid, data, labels, k):
#     plt.figure(figsize=(10, 8))
    
#     # Generate a colormap with a unique color for each grid cell
#     cmap = plt.cm.get_cmap('tab20', k)
#     colors = cmap(range(k))
    
#     # Plot the grid with different colors for each cell
#     for idx, cell in enumerate(grid):
#         (x_min, x_max), (y_min, y_max) = cell
#         rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, color=colors[idx % k], alpha=0.5)
#         plt.gca().add_patch(rect)
    
#     # Plot the points with their cluster labels
#     for idx in range(k):
#         cluster_points = data[labels == idx]
#         plt.scatter(cluster_points[:, 0], cluster_points[:, 1], color=colors[idx % k], label=f'Cluster {idx+1}', edgecolor='black')
    
#     plt.legend()
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.show()


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
    
    # # exit()

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
    
    
    n_clusters = 20  # Change this to the desired number of clusters

    # Literature plot clusters
    # group_number = 3
    # txt_file='p3038_600'
    group_number = 5
    # txt_file='rl1304_010'
    txt_file='pr2392_020'
    # txt_file='fnl4461_0100'
    PATH_DATA=f'./data/Literature/group{group_number}/'
    input_txt_file = f"{PATH_DATA}loc_capacities_{txt_file}.txt"
    locations, data = read_data_txt(input_txt_file)
    
    
    # Perform K-means clustering
    # kmeans = perform_kmeans(data, n_clusters)
    # Plot clusters and boundaries
    # title = f'K-Means Clustering {txt_file}'
    # plot_clusters(data, kmeans, title, f'plots/plots_lit/clusters_plot_{txt_file}')
    # Save points with cluster labels to txt
    # output_txt_file = f"{PATH_DATA}loc_coverages_kmeans_{txt_file}.txt"    
    # save_locations_with_clusters(locations, data, kmeans.labels_, output_txt_file)
    # print(f"saved: loc_coverages_kmeans_{txt_file}")


    
    # Determine the bounding box of the points
    min_x, max_x = np.min(data[:, 0]), np.max(data[:, 0])
    min_y, max_y = np.min(data[:, 1]), np.max(data[:, 1])
    # Create the grid
    # grid = create_grid(min_x, max_x, min_y, max_y, n_clusters)
    # Create a grid with varying cell sizes
    grid = create_varying_grid(min_x, max_x, min_y, max_y, data, n_clusters)
    # Assign cluster labels to points based on grid cells
    labels = assign_cluster_labels(grid, data)
    # Plot the grid and points
    plot_grid_and_points(grid, data, labels)
    output_txt_file = f"{PATH_DATA}loc_coverages_rand_grid_{txt_file}.txt" 
    # Save points with cluster labels to CSV
    save_points_with_clusters(data, labels, output_txt_file)
    # Save locations with cluster labels and coordinates to a text file
    save_locations_with_clusters(locations, data, labels, output_txt_file)

    exit()

    
    PATH_DATA="./data/PACA/"
    input_csv_file = "./data/PACA/points_coord.csv"
    location_identifiers = read_location_identifiers(f"{PATH_DATA}map_id_cust_loc.txt")
    locations, data = read_csv(input_csv_file, location_identifiers)
    txt_file='commune'
    # Perform K-means clustering
    kmeans = perform_kmeans(data, n_clusters)

    # Plot clusters and boundaries
    title = f'K-Means Clustering {txt_file} PACA'
    plot_clusters(data, kmeans, title, f'clusters_plot_{txt_file}')

    # # Save points with cluster labels to CSV
    # output_csv_file_points = "points_with_clusters.csv"
    # save_points_with_clusters(data, kmeans.labels_, output_csv_file_points)

    output_txt_file = f"{PATH_DATA}loc_coverages_kmeans_{txt_file}.txt"    
    save_locations_with_clusters(locations, data, kmeans.labels_, output_txt_file)


    print(f"saved: loc_coverages_kmeans_{txt_file}")

# Run the main function
if __name__ == "__main__":
    main()