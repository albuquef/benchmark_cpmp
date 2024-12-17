import pandas as pd
import numpy as np
import geopandas as gpd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from shapely.geometry import Point

# Function to parse coord_5km from file
def load_coord_5km(file_path):
    return pd.read_csv(
        file_path, 
        sep=r'\s+', 
        names=['location', 'capacity', 'coord_x', 'coord_y', 'id_loc', 'fid'],
        skiprows=1
    )

# Function to parse coord_2km from file
def load_coord_2km(file_path):
    return pd.read_csv(
        file_path, 
        sep=r'\s+', 
        names=['location', 'weight', 'coord_x', 'coord_y', 'id_grid2km', 'fid'],
        skiprows=1
    )

# Function to extract p locations from solution file
def load_p_locations(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_idx = lines.index('P LOCATIONS\n') + 1
    locations = []
    for line in lines[start_idx:]:
        if line.strip() == '':
            break
        locations.append(int(line.strip()))
    return locations

# Load the files
coord_5km_file = 'outputs/PACA_5km/loc_capacities_cinema.txt'
coord_2km_file = 'outputs/PACA_2km/loc_capacities_cinema.txt'
solution_5km_file = 'outputs/solutions_files/solution_cinema_5km_p_135.txt'

coord_5km = load_coord_5km(coord_5km_file)
coord_2km = load_coord_2km(coord_2km_file)
p_locations = load_p_locations(solution_5km_file)

# Ensure numeric coordinates
coord_2km['coord_x'] = pd.to_numeric(coord_2km['coord_x'], errors='coerce')
coord_2km['coord_y'] = pd.to_numeric(coord_2km['coord_y'], errors='coerce')
coord_2km = coord_2km.dropna(subset=['coord_x', 'coord_y'])

# Check for CRS transformation
print(f"Number of 5km locations: {len(coord_5km)}")
print(f"Number of 2km locations: {len(coord_2km)}")
print(f"Number of p locations: {len(p_locations)}")

# Create GeoDataFrame for 5km and 2km locations
gdf_5km = gpd.GeoDataFrame(
    coord_5km, 
    geometry=gpd.points_from_xy(coord_5km['coord_x'], coord_5km['coord_y']),
    crs="EPSG:900913"  # Specify the CRS of your 5km points, e.g., EPSG:2154 or EPSG:4326
)

# create a gdf for the p_locations
gdf_p_locations = gpd.GeoDataFrame(
    {'location': p_locations},
    geometry=[Point(gdf_5km.loc[gdf_5km['location'] == loc][['coord_x', 'coord_y']].values[0]) for loc in p_locations],
    crs="EPSG:900913"
)


gdf_2km = gpd.GeoDataFrame(
    coord_2km, 
    geometry=gpd.points_from_xy(coord_2km['coord_x'], coord_2km['coord_y']),
    crs="EPSG:3035"  # Specify the CRS of your 2km points, e.g., EPSG:2154 or EPSG:4326
)

# Transform both GeoDataFrames to the same CRS (e.g., EPSG:4326 or any common CRS)
gdf_5km = gdf_5km.to_crs("EPSG:3035")  # Replace with the desired CRS
gdf_2km = gdf_2km.to_crs("EPSG:3035")  # Replace with the desired CRS
gdf_p_locations =  gdf_p_locations.to_crs("EPSG:3035")  # Replace with the desired CRS

# plot the 5km locations and the p_locations
fig, ax = plt.subplots(figsize=(10, 8))
gdf_5km.plot(ax=ax, marker='o', color='blue', label='5km Locations')
gdf_p_locations.plot(ax=ax, marker='o', color='red', label='P Locations')
ax.set_aspect('auto')  # Adjust to 'auto' or specific ratio
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('5km Locations and P Locations')
ax.legend()
plt.show()




# debug gdfs
print(gdf_5km.head())
print(gdf_2km.head())


# Initialize list for selected 2km locations
selected_locations = []
used_indices = set()

# For each p_location, find the closest 2km location
for loc in p_locations:
    loc_5km = gdf_5km.loc[gdf_5km['location'] == loc]
    if loc_5km.empty:
        print(f"Location {loc} not found in coord_5km.")
        continue
    coord = loc_5km[['coord_x', 'coord_y']].values

    # Compute distances between the current 5km location and all 2km points
    distances = cdist(coord, gdf_2km[['coord_x', 'coord_y']].values).flatten()
    
    # Find the nearest unused 2km location
    for idx in np.argsort(distances):
        if idx not in used_indices:
            used_indices.add(idx)
            selected_locations.append(gdf_2km.iloc[idx])
            break



# Convert selected locations to GeoDataFrame
selected_2km = gpd.GeoDataFrame(
    selected_locations, 
    geometry=gpd.points_from_xy([loc['coord_x'] for loc in selected_locations], [loc['coord_y'] for loc in selected_locations]),
    crs="EPSG:3035"  # Use the same CRS
).to_crs("EPSG:3035")  # Replace with the desired CRS

selected_2km = selected_2km.to_crs("EPSG:3035")  # Replace with the desired CRS

# Ensure geometries are valid
selected_2km = selected_2km[selected_2km.is_valid]

print(gdf_5km.crs)
print(gdf_2km.crs)
print(selected_2km.crs)
print(gdf_5km[['coord_x', 'coord_y']].describe())
print(gdf_2km[['coord_x', 'coord_y']].describe())
# debug selected_2km
print(selected_2km.head())

print(f"Number of selected 2km locations: {len(selected_2km)}")
print(selected_2km[['coord_x', 'coord_y']].describe())

# Extract coordinates from selected_2km GeoDataFrame
x_coords = selected_2km['coord_x']
y_coords = selected_2km['coord_y']


# Create a plot with a fixed aspect ratio
fig, ax = plt.subplots(figsize=(10, 8))

# Plot the 5km locations
gdf_5km.plot(ax=ax, marker='o', color='blue', label='Original 5km Locations')

# Plot selected 2km locations
selected_2km.plot(ax=ax, marker='o', color='red', label='Selected 2km Locations')

# Adjust aspect ratio manually if needed
ax.set_aspect('auto')  # Adjust to 'auto' or specific ratio
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Original 5km Locations and Selected 2km Locations')
ax.legend()

plt.show()

# # Plot side-by-side the original 5km points and the selected 2km locations
# fig, axs = plt.subplots(1, 2, figsize=(15, 8))

# # Plot the original 5km locations
# axs[0].scatter(coord_5km['coord_x'], coord_5km['coord_y'], c='blue', label='Original 5km Locations')
# axs[0].set_xlabel('X Coordinate')
# axs[0].set_ylabel('Y Coordinate')
# axs[0].set_title('Original 5km Locations')
# axs[0].legend()

# # Plot the selected 2km locations
# axs[1].scatter(coord_2km['coord_x'], coord_2km['coord_y'], c='lightgrey', label='All 2km Points')
# axs[1].scatter(selected_2km['coord_x'], selected_2km['coord_y'], c='red', label='Selected Locations')
# axs[1].set_xlabel('X Coordinate')
# axs[1].set_ylabel('Y Coordinate')
# axs[1].set_title('Selected 2km Locations')
# axs[1].legend()


# # Display the plots
# plt.tight_layout()
# plt.show()

