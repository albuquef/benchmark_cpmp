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

# Apply the function to extract coordinates and create new columns in the DataFrame
df_points_pop_paca['easting'], df_points_pop_paca['northing'] = zip(*df_points_pop_paca['idcar_1km'].apply(extract_coordinates))

# Drop rows where coordinates couldn't be extracted
df_points_pop_paca = df_points_pop_paca.dropna(subset=['easting', 'northing'])

# Create a GeoDataFrame for the first dataset
gdf_points_pop_paca = gpd.GeoDataFrame(df_points_pop_paca, geometry=gpd.points_from_xy(df_points_pop_paca.easting, df_points_pop_paca.northing))

# Set the coordinate reference system (CRS) to the appropriate CRS (EPSG:3035)
gdf_points_pop_paca.set_crs(epsg=3035, inplace=True)


# Create a GeoDataFrame for the second dataset with initial CRS as EPSG:4326 (WGS 84)
gdf_points_grid_5km = gpd.GeoDataFrame(df_points_grid_5km, geometry=gpd.points_from_xy(df_points_grid_5km.longitude, df_points_grid_5km.latitude), crs='ESRI:102100')

# Convert the CRS to EPSG:3035
gdf_points_grid_5km = gdf_points_grid_5km.to_crs(epsg=3035)

# Perform a spatial join to associate points in df_points_pop_paca with the closest point in df_points_grid_5km
gdf_points_grid_5km['weight'] = 0  # Initialize the 'weight' column

# Extract x and y coordinates
gdf_points_pop_paca_coords = list(zip(gdf_points_pop_paca.geometry.x, gdf_points_pop_paca.geometry.y))
gdf_points_grid_5km_coords = list(zip(gdf_points_grid_5km.geometry.x, gdf_points_grid_5km.geometry.y))

# Use NearestNeighbors to find the nearest point in gdf_points_grid_5km for each point in gdf_points_pop_paca
nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(gdf_points_grid_5km_coords)
distances, indices = nn.kneighbors(gdf_points_pop_paca_coords)

# Add the closest point indices to gdf_points_pop_paca
gdf_points_pop_paca['closest_idx'] = indices.flatten()

# Sum the 'ind' values for each point in gdf_points_grid_5km
for idx in range(len(gdf_points_grid_5km)):
    gdf_points_grid_5km.at[idx, 'weight'] = gdf_points_pop_paca[gdf_points_pop_paca['closest_idx'] == idx]['ind'].sum()

# Plot both sets of points on the same map
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_pop_paca.plot(ax=ax, marker='o', color='red', markersize=5, label='Dataset 1')
gdf_points_grid_5km.plot(ax=ax, marker='x', color='blue', markersize=5, label='Dataset 2')
plt.title("Plot of Points from Both Datasets")
plt.xlabel("Easting")
plt.ylabel("Northing")
plt.legend()
plt.show()


# count the number of zeros in the 'weight' column
num_zeros = (gdf_points_grid_5km['weight'] == 0).sum()
print(f"Number of zeros in the 'weight' column: {num_zeros}")
# number of non-zero values in the 'weight' column
num_non_zeros = (gdf_points_grid_5km['weight'] != 0).sum()
print(f"Number of non-zero values in the 'weight' column: {num_non_zeros}")


# compare the sum of ind column with the sum of weight column
sum_ind = gdf_points_pop_paca['ind'].sum()
sum_weight = gdf_points_grid_5km['weight'].sum()
print(f"Sum of 'ind' column: {sum_ind}")
print(f"Sum of 'weight' column: {sum_weight}")
if sum_ind == sum_weight:
    print("The sums are equal.")
else:
    print("The sums are not equal.")


# change the color the points with zero weight
gdf_points_grid_5km['color'] = 'green'
gdf_points_grid_5km.loc[gdf_points_grid_5km['weight'] == 0, 'color'] = 'red'
# plot the points with different colors
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_pop_paca.plot(ax=ax, marker='o', color='red', markersize=5, label='Dataset 1')
gdf_points_grid_5km.plot(ax=ax, marker='x', color=gdf_points_grid_5km['color'], markersize=5, label='Dataset 2')
plt.title("Plot of Points from Both Datasets")
plt.xlabel("Easting")
plt.ylabel("Northing")
plt.legend()
plt.show()

# eliminate the points with zero weight
gdf_points_grid_5km = gdf_points_grid_5km[gdf_points_grid_5km['weight'] != 0]


# print max min mean stddev of the 'weight' column
print(f"Max value in the 'weight' column: {gdf_points_grid_5km['weight'].max()}")
print(f"Min value in the 'weight' column: {gdf_points_grid_5km['weight'].min()}")
print(f"Mean value in the 'weight' column: {gdf_points_grid_5km['weight'].mean()}")
print(f"Standard deviation of the 'weight' column: {gdf_points_grid_5km['weight'].std()}")

# plot graphic dist of weight with better visualization
gdf_points_grid_5km['weight'].plot(kind='hist', bins=5)
plt.title("Distribution of Weights")
plt.xlabel("Weight")
plt.ylabel("Frequency")
plt.show()

# plot the heat map of the df_points_grid_5km using the 'weight' column
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_grid_5km.plot(ax=ax, column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
plt.title("Heat Map of Points with Weights")
plt.xlabel("Easting")
plt.ylabel("Northing")
plt.show()


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