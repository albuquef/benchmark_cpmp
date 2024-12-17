import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt


def extract_coordinates(coord_string):
    try:
        # Ensure coord_string is a string
        if not isinstance(coord_string, str):
            raise ValueError(f"Invalid format for coordinates: {coord_string}")
        
        # Remove the prefix "CRS3035RES1000m"
        coord_string = coord_string.replace("CRS3035RES1000m", "")
        
        # Ensure the string starts with "N" and contains 'E'
        if not coord_string.startswith("N") or 'E' not in coord_string:
            raise ValueError(f"Invalid format for coordinates: {coord_string}")
        
        # Find positions of 'N' and 'E'
        pos_n = coord_string.index('N') + 1
        pos_e = coord_string.index('E')
        
        # Extract northing and easting parts
        northing_str = coord_string[pos_n:pos_e]
        easting_str = coord_string[pos_e + 1:]
        
        # Convert to integers
        northing = int(northing_str)
        easting = int(easting_str)
        
        return pd.Series([northing, easting])
    
    except ValueError as ve:
        raise ve  # Let ValueError propagate to handle specific errors
    except Exception as e:
        raise ValueError(f"Error processing {coord_string}: {e}")


# Read the CSV file
df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/data_Felipe26_06/carreaux_1km_PACA.csv')


# Apply the function to the DataFrame
df[['northing', 'easting']] = df['idcar_1km'].apply(extract_coordinates)

# Define the CRS
crs = "EPSG:3035"

# Create a geometry column
geometry = [Point(xy) for xy in zip(df['easting'], df['northing'])]

# Create a GeoDataFrame
gdf = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)

# # Plot the points
gdf.plot(marker='o', color='red', markersize=5)
plt.title(f'Points from CSV {gdf.shape}')
plt.xlabel('Easting')
plt.ylabel('Northing')
plt.show()

exit()

#################################################################################################


# Define boundaries for the grid (example values, adjust as needed)
min_easting = gdf['easting'].min()
max_easting = gdf['easting'].max()
min_northing = gdf['northing'].min()
max_northing = gdf['northing'].max()

# Define grid cell dimensions (width and height in meters)
cell_width = 5000
cell_height = 5000

# Create grid cells as polygons
grid_cells = []
for x in range(int(min_easting), int(max_easting), cell_width):
    for y in range(int(min_northing), int(max_northing), cell_height):
        polygon = Polygon([
            (x, y),
            (x + cell_width, y),
            (x + cell_width, y + cell_height),
            (x, y + cell_height)
        ])
        grid_cells.append(polygon)

# Create GeoDataFrame for grid cells
grid_gdf = gpd.GeoDataFrame({'geometry': grid_cells}, crs="EPSG:3035")

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))
gdf.plot(ax=ax, marker='o', color='red', markersize=50, label='Points')
grid_gdf.plot(ax=ax, edgecolor='blue', facecolor='none', linewidth=0.5, label='Grid')
plt.title('Points with Custom Grid')
plt.xlabel('Easting')
plt.ylabel('Northing')
plt.legend()
plt.show()



exit()











#######################################################################################################

# Define boundaries for the grid (example values, adjust as needed)
min_easting = gdf['easting'].min()
max_easting = gdf['easting'].max()
min_northing = gdf['northing'].min()
max_northing = gdf['northing'].max()

# Define grid cell dimensions (width and height in meters)
cell_width = 5000
cell_height = 5000

# Create grid cells as polygons
# Create grid cells as polygons and calculate points inside each grid
grid_cells = []
for x in range(int(min_easting), int(max_easting), cell_width):
    for y in range(int(min_northing), int(max_northing), cell_height):
        polygon = Polygon([
            (x, y),
            (x + cell_width, y),
            (x + cell_width, y + cell_height),
            (x, y + cell_height)
        ])
        points_inside = gdf[gdf.within(polygon)]
        if not points_inside.empty:
            grid_cells.append((polygon, points_inside.shape[0], polygon.centroid))

# Create GeoDataFrame for grid cells
grid_gdf = gpd.GeoDataFrame(geometry=[gc[0] for gc in grid_cells], crs="EPSG:3035")
grid_gdf['points_inside'] = [gc[1] for gc in grid_cells]
grid_gdf['centroid'] = [gc[2] for gc in grid_cells]

# Calculate central point only if there are points inside any grid cell
if grid_gdf['points_inside'].sum() > 0:
    central_point_x = grid_gdf.centroid.x.mean()
    central_point_y = grid_gdf.centroid.y.mean()
    central_point = Point(central_point_x, central_point_y)
else:
    central_point = None

# Count of grid cells
num_grids = len(grid_gdf)

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))
gdf.plot(ax=ax, marker='o', color='red', markersize=50, label='Points')
grid_gdf.plot(ax=ax, column='points_inside', cmap='Blues', edgecolor='blue', linewidth=0.5, legend=True)
# grid_gdf.plot(ax=ax, edgecolor='blue', facecolor='none', linewidth=0.5, label='Grid')
for idx, row in grid_gdf.iterrows():
    ax.scatter(row['centroid'].x, row['centroid'].y, color='green', s=50)
# if central_point:
#     ax.scatter(central_point_x, central_point_y, color='green', s=100, label='Central Point')
plt.title(f'Points with Custom Grid - {num_grids} grids')
plt.xlabel('Easting')
plt.ylabel('Northing')
plt.legend()
plt.show()


# # Example usage:
# coord_string = "CRS3035RES1000mN2216000E4010000"
# northing, easting = extract_coordinates(coord_string)
# print(f"Northing: {northing}, Easting: {easting}")


