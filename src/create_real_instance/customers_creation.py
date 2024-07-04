import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import re
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors


# points population 1km paca with demand weight
df_points_pop_paca = pd.read_csv('data/data_qgis/data_instance_paca/points_population_1km_paca_table.csv')
# grid population 5km paca without demand weight
df_points_grid_5km = pd.read_csv('data/data_qgis/data_instance_paca/points_grid_5km_paca_table.csv')

# Function to extract coordinates from the 'idcar_1km' column
def extract_coordinates(idcar):
    match = re.search(r'N(\d+)E(\d+)', idcar)
    if match:
        northing = int(match.group(1))
        easting = int(match.group(2))
        return easting, northing
    else:
        return None, None

def load_population_data(filepath):
    """Loads population data and extracts coordinates."""
    df = pd.read_csv(filepath)
    df['easting'], df['northing'] = zip(*df['idcar_1km'].apply(extract_coordinates))
    df = df.dropna(subset=['easting', 'northing'])
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['easting'], df['northing']))
    gdf.set_crs(epsg=3035, inplace=True)
    return gdf

def load_grid_points_data(filepath):
    """Loads grid data and converts CRS."""
    df = pd.read_csv(filepath, low_memory=False)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']), crs='ESRI:102100') # coordinates in WGS 84
    gdf = gdf.to_crs(epsg=3035) # Convert to EPSG:3035
    return gdf

def calculate_weights(gdf_pop, gdf_grid):
    """Performs spatial join to find nearest neighbors and calculate weights."""
    print('-' * 50)
    print("[INFO] Calculating weights...")

    gdf_grid['weight'] = 0.0  # Ensure weight column is of type float
    gdf_pop_coords = list(zip(gdf_pop.geometry.x, gdf_pop.geometry.y))
    gdf_grid_coords = list(zip(gdf_grid.geometry.x, gdf_grid.geometry.y))
    
    nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(gdf_grid_coords)
    distances, indices = nn.kneighbors(gdf_pop_coords)
    gdf_pop['closest_idx'] = indices.flatten()
    
    for idx in range(len(gdf_grid)):
        gdf_grid.at[idx, 'weight'] = gdf_pop[gdf_pop['closest_idx'] == idx]['ind'].sum()
    
    # compare the sum of ind column with the sum of weight column
    sum_ind = gdf_pop['ind'].sum()
    sum_weight = gdf_grid['weight'].sum()
    
    if sum_ind == sum_weight:
        print(f"The sums of the weights are equal: {sum_ind}.")
    else:
        print("[WARN] The sums of the weights are not equal.")
        print(f"Sum of 'ind' pop: {sum_ind}")
        print(f"Sum of 'weight' grid: {sum_weight}")
        
    print('-' * 50)
    
    return gdf_grid

def filter_weights(gdf_grid, k):
    """Eliminates points with zero weight and transfers weights smaller than k to the nearest non-zero point."""
    print('-' * 50)
    print(f"[INFO] Filtering points with zero weight and transferring weights smaller than {k}...")
    print(f"Number of elements before remove the zeros weights: {len(gdf_grid)}")
    
    # Initialize color column
    gdf_grid['color'] = 'green'
    # Separate points with zero weight for later plotting
    zero_weight_points = gdf_grid[gdf_grid['weight'] == 0].copy()
    zero_weight_points['color'] = 'red'
    print(f'Number of elements with zero weight: {len(gdf_grid[gdf_grid["weight"] == 0])}')
    
    # Filter out points with zero weight
    gdf_grid = gdf_grid[gdf_grid['weight'] != 0].copy()
    
    
    # Identify points with weights smaller than k
    small_weight_mask = gdf_grid['weight'] < k
    large_weight_mask = gdf_grid['weight'] >= k
    
    small_weight_points = gdf_grid[small_weight_mask]
    large_weight_points = gdf_grid[large_weight_mask]
    
    if not small_weight_points.empty and not large_weight_points.empty:
        # Get coordinates for small and large weight points
        small_weight_coords = list(zip(small_weight_points.geometry.x, small_weight_points.geometry.y))
        large_weight_coords = list(zip(large_weight_points.geometry.x, large_weight_points.geometry.y))
        
        # Find the nearest large weight point for each small weight point
        nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(large_weight_coords)
        distances, indices = nn.kneighbors(small_weight_coords)
        
        # Transfer small weights to the nearest large weight point
        for i, idx in enumerate(small_weight_points.index):
            nearest_large_idx = large_weight_points.index[indices[i][0]]
            gdf_grid.at[nearest_large_idx, 'weight'] += gdf_grid.at[idx, 'weight']
            gdf_grid.at[idx, 'weight'] = 0  # Set small weight to zero after transfer
        
            # Identify and separate points that now have zero weight after transfer
        transferred_zero_weight_points = gdf_grid[gdf_grid['weight'] == 0].copy()
        transferred_zero_weight_points['color'] = 'red'
        
        # Combine both sets of zero-weight points
        zero_weight_points = zero_weight_points.append(transferred_zero_weight_points)
        
        
        # Remove points that now have zero weight after transfer
        gdf_grid = gdf_grid[gdf_grid['weight'] != 0]
    
    
    # Print how many ellement before and now
    print(f"Number of elements before: {len(gdf_grid) + len(small_weight_points)}")
    print(f"Number of elements after filtering less than {k}: {len(gdf_grid)}")
    if k != 0:
        print(f'Number of elements removed in total: {len(small_weight_points) + len(transferred_zero_weight_points)}')   
    print('-' * 50) 
    
    create_plot = True
    # Plotting
    if create_plot:
        fig, ax = plt.subplots()
        gdf_grid[gdf_grid['color'] == 'green'].plot(ax=ax, color='green', label=f'Remaining (size: {len(gdf_grid)})')
        if not zero_weight_points.empty:
            zero_weight_points.plot(ax=ax, color='red', label=f'Removed (size: {len(zero_weight_points)})')
        plt.legend()
        plt.title(f"Grid points with the filter to reallocate weights smaller than {k}") 
        plt.savefig('plots/creation_instance/customers/points_grid_5km_paca_table_filter.png')
        plt.show()
    
    return gdf_grid

def compare_grid_with_old_dataset(gdf_grid, gdf_old_dataset):
    
    print('-' * 50)
    print("[INFO] Comparing the new dataset with the old one...")
    print(f"Number of elements in the new dataset: {len(gdf_grid)}")
    print(f"Number of elements in the old dataset: {len(gdf_old_dataset)}")
    # compare dist max min avg and std of the two datasets
    print('New dataset:')
    print(f"Sum of 'weight' column in the new dataset: {gdf_grid['weight'].sum()}")
    print(f"Max value in the 'weight' column of the new dataset: {gdf_grid['weight'].max()}")
    print(f"Min value in the 'weight' column of the new dataset: {gdf_grid['weight'].min()}")
    print(f"Mean value in the 'weight' column of the new dataset: {gdf_grid['weight'].mean()}")
    print(f"Standard deviation of the 'weight' column of the new dataset: {gdf_grid['weight'].std()}")
    print('Old dataset:')
    print(f"Sum of 'weight' column in the old dataset: {gdf_old_dataset['weight'].sum()}")
    print(f"Max value in the 'weight' column of the old dataset: {gdf_old_dataset['weight'].max()}")
    print(f"Min value in the 'weight' column of the old dataset: {gdf_old_dataset['weight'].min()}")
    print(f"Mean value in the 'weight' column of the old dataset: {gdf_old_dataset['weight'].mean()}")
    print(f"Standard deviation of the 'weight' column of the old dataset: {gdf_old_dataset['weight'].std()}")
    
    create_plot = True  
    if create_plot:
        # compare the histogram and heatmap of the two datasets
        fig, ax = plt.subplots(1, 1, figsize=(20, 10))
        gdf_grid['weight'].plot(kind='hist', bins=5, alpha=0.5, label=f'new dataset ({len(gdf_grid)})')
        gdf_old_dataset['weight'].plot(kind='hist', bins=5, alpha=0.5, label=f'old dataset ({len(gdf_old_dataset)})')
        plt.title("Distribution of Weights")
        plt.xlabel("Weight")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()
        
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        gdf_grid.plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
        gdf_old_dataset.plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
        ax[0].set_title(f"Heat Map of Grid Population (new dataset)\nnum elements= {len(gdf_grid)}, sum = {gdf_grid['weight'].sum()}")
        ax[0].set_xlabel("Easting")
        ax[0].set_ylabel("Northing")
        ax[1].set_title(f"Heat Map of Grid Population (old dataset)\nnum elements={len(gdf_old_dataset)}  , sum wi= {gdf_old_dataset['weight'].sum()}")
        ax[1].set_xlabel("Easting")
        ax[1].set_ylabel("Northing")
        plt.savefig('plots/creation_instance/customers/heatmaps_comparison.png')
        plt.show()


    print('-' * 50)


def plot_pop_and_grid_points(gdf_pop, gdf_grid):
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    gdf_pop.plot(ax=ax, marker='x', color='blue', markersize=5, label=f'Population (size: {gdf_pop.shape[0]})')
    gdf_grid.plot(ax=ax, marker='o', color='green', markersize=10, label=f'Grid points (size: {gdf_grid.shape[0]})')
    plt.title("Plot of Points from Both Datasets")
    plt.xlabel("Easting")
    plt.ylabel("Northing")
    plt.legend()
    plt.savefig('plots/creation_instance/customers/points_pop_grid_5km_paca_table.png')
    plt.show()
    
    # """Plots datasets side by side."""
    # fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    # gdf_grid.plot(ax=ax[0], marker='x', color='blue', markersize=5, label='Grid')
    # gdf_pop.plot(ax=ax[1], marker='o', color='red', markersize=5, label='Population')
    # ax[0].set_title("Grid Points")
    # ax[0].set_xlabel("Easting")
    # ax[0].set_ylabel("Northing")
    # ax[1].set_title("Population Points")
    # ax[1].set_xlabel("Easting")
    # ax[1].set_ylabel("Northing")
    # plt.show()

def plot_heatmap_grid_points(gdf_grid):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    gdf_grid.plot(ax=ax, column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
    plt.title("Heat Map of Grid Points with Weights")
    plt.xlabel("Easting")
    plt.ylabel("Northing")
    plt.savefig('plots/creation_instance/customers/points_grid_5km_paca_table_heatmap.png')
    plt.show()

def create_final_table_instance(gdf_grid):
    # create one file txt with 4 columns:id, weight, coord_x, coord_y, fid
    filename = './outputs/PACA_Jun2024/cust_weights_Jun2024.txt'

    print('-' * 50)
    print(f"[INFO] Creating final table instance in {filename}...")
    
    cont = 1
    # print columns names  of df_grid
    with open(filename, 'w') as f:
        # first line with the columns names
        f.write("customer weight coord_x coord_y id_grid5km fid\n")
        for idx, row in gdf_grid.iterrows():
            f.write(f"{cont} {row['weight']} {row['geometry'].x} {row['geometry'].y} {idx} {row['fid']}\n")
            cont += 1

    
    print('-' * 50)
    
    

gdf_points_pop_paca = load_population_data('data/data_qgis/data_instance_paca/points_population_1km_paca_table.csv')
gdf_points_grid_5km = load_grid_points_data('data/data_qgis/data_instance_paca/points_grid_5km_paca_table.csv')
gdf_points_grid_5km = calculate_weights(gdf_points_pop_paca, gdf_points_grid_5km)
plot_pop_and_grid_points(gdf_points_pop_paca, gdf_points_grid_5km)

gdf_points_grid_5km = filter_weights(gdf_points_grid_5km, 0) # Filter weights smaller than k and remove points with zero weight
plot_heatmap_grid_points(gdf_points_grid_5km)


create_final_table_instance(gdf_points_grid_5km)



#  -----    compare with old dataset ------
df_cust_weights = pd.read_csv('data/PACA/cust_weights.txt', delim_whitespace=True, names=['customer', 'weight'])
df_map_id_cust_loc = pd.read_csv('data/PACA/map_id_cust_loc.txt', delim_whitespace=True, names=['id', 'identif'])
df_locations = pd.read_csv('data/PACA/locations_paca_2017_coord.csv')

# Ensure that 'identif' columns are of the same type (string)
df_map_id_cust_loc['identif'] = df_map_id_cust_loc['identif'].astype(str)
df_locations['identif'] = df_locations['identif'].astype(str)

# Merge the DataFrames to associate weights with coordinates
df_merged = pd.merge(df_map_id_cust_loc, df_cust_weights, left_on='id', right_on='customer')
df_old_dataset = pd.merge(df_merged, df_locations, left_on='identif', right_on='identif')
# convert gdf_final['weight'] to float
df_old_dataset['weight'] = df_old_dataset['weight'].astype(float)

# Create a GeoDataFrame with the coordinates and weights in WGS 84 (EPSG:4326)
df_old_dataset = gpd.GeoDataFrame(df_old_dataset, geometry=gpd.points_from_xy(df_old_dataset.x, df_old_dataset.y), crs='EPSG:3035')
compare_grid_with_old_dataset(gdf_points_grid_5km, df_old_dataset)




exit()

# Read the tables into DataFrames
df_cust_weights = pd.read_csv('data/PACA/cust_weights.txt', delim_whitespace=True, names=['customer', 'weight'])
df_map_id_cust_loc = pd.read_csv('data/PACA/map_id_cust_loc.txt', delim_whitespace=True, names=['id', 'identif'])
df_locations = pd.read_csv('data/PACA/locations_paca_2017_coord.csv')

# Ensure that 'identif' columns are of the same type (string)
df_map_id_cust_loc['identif'] = df_map_id_cust_loc['identif'].astype(str)
df_locations['identif'] = df_locations['identif'].astype(str)

# Merge the DataFrames to associate weights with coordinates
df_merged = pd.merge(df_map_id_cust_loc, df_cust_weights, left_on='id', right_on='customer')
df_final = pd.merge(df_merged, df_locations, left_on='identif', right_on='identif')

# Create a GeoDataFrame with the coordinates and weights in WGS 84 (EPSG:4326)
gdf_final = gpd.GeoDataFrame(df_final, geometry=gpd.points_from_xy(df_final.x, df_final.y), crs='EPSG:3035')

# Plot the points
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_grid_5km.plot(ax=ax, marker='x', color=gdf_points_grid_5km['color'], markersize=5, label='Dataset 2')
gdf_final.plot(ax=ax, marker='o', color='red', markersize=5, label='old dataset')
plt.title("Plot of Points with Weights (ESRI:102100)")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

# convert gdf_final['weight'] to float
gdf_final['weight'] = gdf_final['weight'].astype(float)

# compare the weights of the two datasets
sum_weight_gdf_points_grid_5km = gdf_points_grid_5km['weight'].sum()
sum_weight_gdf_final = gdf_final['weight'].sum()
print(f"Sum of 'weight' column in gdf_points_grid_5km: {sum_weight_gdf_points_grid_5km}")
print(f"Sum of 'weight' column in gdf_final: {sum_weight_gdf_final}")
# graphic comparation of weights of the two datasets
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_grid_5km['weight'].plot(kind='hist', bins=5, alpha=0.5, label='Dataset 2')
gdf_final['weight'].plot(kind='hist', bins=5, alpha=0.5, label='old dataset')
plt.title("Distribution of Weights")
plt.xlabel("Weight")
plt.ylabel("Frequency")
plt.legend()
plt.show()

# now graphically compare the two datasets plot side by side the heat map of the two datasets
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
gdf_points_grid_5km.plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
gdf_final.plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
ax[0].set_title("Heat Map of Points with Weights (Dataset 2)")
ax[0].set_xlabel("Easting")
ax[0].set_ylabel("Northing")
ax[1].set_title("Heat Map of Points with Weights (old dataset)")
ax[1].set_xlabel("Easting")
ax[1].set_ylabel("Northing")
plt.show()

#eliminate columns that are not necessary gdf_points_grid_5km
gdf_points_grid_5km = gdf_points_grid_5km.drop(columns=['color'])

# Save the updated df_points_grid_5km with the 'weight' column to a new CSV file
gdf_points_grid_5km.drop(columns='geometry').to_csv('test.csv', index=False)



# print the slowest last five weights of the gdf_points_grid_5km
# the slowest five weights of the gdf_points_grid_5km with unique values
unique_weights = gdf_points_grid_5km['weight'].unique()
# order
unique_weights.sort()
# print the slowest five weights
print(unique_weights[:10])

# same for the gdf_final
unique_weights = gdf_final['weight'].unique()
unique_weights.sort()
print(unique_weights[:10])


# compare the distribution of the smallest hundred unique weights of the two datasets
# get the smallest hundred weights of the gdf_points_grid_5km
num_of_smallest_weights = 1
unique_weights = gdf_points_grid_5km['weight'].unique()
unique_weights.sort()
smallest_weights_gdf_points_grid_5km = unique_weights[:num_of_smallest_weights]
# get the smallest hundred weights of the gdf_final
unique_weights = gdf_final['weight'].unique()
unique_weights.sort()
smallest_weights_gdf_final = unique_weights[:num_of_smallest_weights]
# plot the distribution of the smallest hundred weights of the two datasets
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
gdf_points_grid_5km[gdf_points_grid_5km['weight'].isin(smallest_weights_gdf_points_grid_5km)].plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
gdf_final[gdf_final['weight'].isin(smallest_weights_gdf_final)].plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})

ax[0].set_title(f"Heat Map of Points with Smallest Hundred Weights (Dataset 2) [size: {len(smallest_weights_gdf_points_grid_5km)}]")
ax[0].set_xlabel("Easting")
ax[0].set_ylabel("Northing")
ax[1].set_title(f"Heat Map of Points with Smallest Hundred Weights (old dataset) [size: {len(smallest_weights_gdf_final)}]")
ax[1].set_xlabel("Easting")
ax[1].set_ylabel("Northing")
plt.show()







# import pandas as pd
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import re
# from sklearn.neighbors import NearestNeighbors

# def extract_coordinates(idcar):
#     """Extracts coordinates from the 'idcar' column."""
#     match = re.search(r'N(\d+)E(\d+)', idcar)
#     if match:
#         northing = int(match.group(1))
#         easting = int(match.group(2))
#         return easting, northing
#     else:
#         return None, None

# def load_population_data(filepath):
#     """Loads population data and extracts coordinates."""
#     df = pd.read_csv(filepath)
#     df['easting'], df['northing'] = zip(*df['idcar_1km'].apply(extract_coordinates))
#     df = df.dropna(subset=['easting', 'northing'])
#     gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['easting'], df['northing']))
#     gdf.set_crs(epsg=3035, inplace=True)
#     return gdf

# def load_grid_data(filepath):
#     """Loads grid data and converts CRS."""
#     df = pd.read_csv(filepath)
#     gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']), crs='EPSG:4326')
#     gdf = gdf.to_crs(epsg=3035)
#     return gdf

# def find_nearest_weights(gdf_pop, gdf_grid):
#     """Performs spatial join to find nearest neighbors and calculate weights."""
#     gdf_grid['weight'] = 0
#     gdf_pop_coords = list(zip(gdf_pop.geometry.x, gdf_pop.geometry.y))
#     gdf_grid_coords = list(zip(gdf_grid.geometry.x, gdf_grid.geometry.y))
#     nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(gdf_grid_coords)
#     distances, indices = nn.kneighbors(gdf_pop_coords)
#     gdf_pop['closest_idx'] = indices.flatten()
#     for idx in range(len(gdf_grid)):
#         gdf_grid.at[idx, 'weight'] = gdf_pop[gdf_pop['closest_idx'] == idx]['ind'].sum()
#     return gdf_grid

# def plot_datasets(gdf_pop, gdf_grid, title):
#     """Plots datasets side by side."""
#     fig, ax = plt.subplots(1, 2, figsize=(20, 10))
#     gdf_grid.plot(ax=ax[0], marker='x', color=gdf_grid['color'], markersize=5, label='Grid')
#     gdf_pop.plot(ax=ax[1], marker='o', color='red', markersize=5, label='Population')
#     ax[0].set_title(f"{title} - Grid")
#     ax[0].set_xlabel("Easting")
#     ax[0].set_ylabel("Northing")
#     ax[1].set_title(f"{title} - Population")
#     ax[1].set_xlabel("Easting")
#     ax[1].set_ylabel("Northing")
#     plt.show()

# def analyze_weights(gdf_grid):
#     """Analyzes and visualizes weights."""
#     num_zeros = (gdf_grid['weight'] == 0).sum()
#     num_non_zeros = (gdf_grid['weight'] != 0).sum()
#     sum_weight = gdf_grid['weight'].sum()
#     print(f"Number of zeros in the 'weight' column: {num_zeros}")
#     print(f"Number of non-zero values in the 'weight' column: {num_non_zeros}")
#     print(f"Sum of 'weight' column: {sum_weight}")
#     print(f"Max value in the 'weight' column: {gdf_grid['weight'].max()}")
#     print(f"Min value in the 'weight' column: {gdf_grid['weight'].min()}")
#     print(f"Mean value in the 'weight' column: {gdf_grid['weight'].mean()}")
#     print(f"Standard deviation of the 'weight' column: {gdf_grid['weight'].std()}")
#     gdf_grid['weight'].plot(kind='hist', bins=5)
#     plt.title("Distribution of Weights")
#     plt.xlabel("Weight")
#     plt.ylabel("Frequency")
#     plt.show()
#     fig, ax = plt.subplots(1, 1, figsize=(10, 10))
#     gdf_grid.plot(ax=ax, column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
#     plt.title("Heat Map of Points with Weights")
#     plt.xlabel("Easting")
#     plt.ylabel("Northing")
#     plt.show()

# def main():
#     # Load and process population data
#     gdf_pop = load_population_data('data/data_qgis/data_instance_paca/points_population_1km_paca_table.csv')
    
#     # Load and process grid data
#     gdf_grid = load_grid_data('data/data_qgis/data_instance_paca/points_grid_5km_paca_table.csv')
    
#     # Perform spatial join and calculate weights
#     gdf_grid = find_nearest_weights(gdf_pop, gdf_grid)
    
#     # Plot datasets side by side
#     plot_datasets(gdf_pop, gdf_grid, "Population vs Grid")
    
#     # Analyze weights
#     analyze_weights(gdf_grid)

#     # Additional operations if needed
#     # ...

# if __name__ == "__main__":
#     main()