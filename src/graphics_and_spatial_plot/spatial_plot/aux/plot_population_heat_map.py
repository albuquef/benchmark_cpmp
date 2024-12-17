import pandas as pd
import geopandas as gpd
import plotly.express as px
import plotly.graph_objects as go  # Importing the graph_objects module
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry import box  # Import box from shapely
import folium
import numpy as np
# https://python-charts.com/spatial/spatial-heatmap-plotly/


# # Example 1
# Data with latitude/longitude and values
# df = pd.read_csv('https://raw.githubusercontent.com/R-CoderDotCom/data/main/sample_datasets/population_galicia.csv')

# fig = px.density_mapbox(df, lat = 'latitude', lon = 'longitude', z = 'tot_pob',
#                         radius = 7,
#                         center = dict(lat = 42.83, lon = -8.35),
#                         zoom = 6,
#                         mapbox_style = 'open-street-map',
#                         color_continuous_scale = 'rainbow',
#                         opacity = 0.5)
# fig.show()


### Example 2
# # Load the data
# df = pd.read_csv('https://raw.githubusercontent.com/R-CoderDotCom/data/main/sample_datasets/population_galicia.csv')

# # Set up the figure and the projection (Limbert Cylindrical)
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertCylindrical())

# # Scatter plot on the map
# ax.scatter(df['longitude'], df['latitude'], transform=ccrs.PlateCarree(), s=df['tot_pob']/100, alpha=0.5)

# # Add coastlines for reference
# ax.coastlines()

# plt.show()

# def add_point_to_map(point, crs_name, color='red'):
#     """Create a Folium map with a point."""
#     # Create a GeoDataFrame for the point
#     gdf_point = gpd.GeoDataFrame(geometry=[Point(point)], crs="EPSG:2154")  # Start with known CRS

#     # Reproject to WGS 84 (latitude/longitude)
#     gdf_reprojected = gdf_point.to_crs(epsg=4326)
    
#     # Create a Folium map centered at the point
#     folium_map = folium.Map(location=[gdf_reprojected.geometry.y[0], gdf_reprojected.geometry.x[0]], zoom_start=14)

#     # Add the point to the map
#     folium.Marker(
#         location=[gdf_reprojected.geometry.y[0], gdf_reprojected.geometry.x[0]],
#         popup=f'CRS: {crs_name}',
#         icon=folium.Icon(color=color)
#     ).add_to(folium_map)

#     return folium_map

# # Example point (replace with your coordinates)
# test_point = (478861.72591274, 5390061.79565943)

# # List of different CRS to test
# crs_list = ["EPSG:4326", "EPSG:2154", "EPSG:27572"]  # Add more CRS as needed

# # Create a Folium map and add points for each CRS
# folium_map = None
# for crs in crs_list:
#     folium_map = add_point_to_map(test_point, crs, color='red' if crs == "EPSG:2154" else 'blue')

# # Display the map
# folium_map.save('map_with_points.html')
# folium_map  # If using Jupyter Notebook, this will display the map inline




# Read the data
# df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt', sep=' ')

# # Create a GeoDataFrame from the Lambert 94 coordinates
# gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs="EPSG:2154")

# # Reproject to WGS 84 (latitude/longitude)
# gdf = gdf.to_crs(epsg=4326)

# # Extract the transformed coordinates
# df['latitude'] = gdf.geometry.x
# df['longitude'] = gdf.geometry.y

# # Plot using Plotly
# fig = px.density_mapbox(df, lat='latitude', lon='longitude', z='weight',
#                         radius=7, zoom=6,
#                         mapbox_style='open-street-map',
#                         color_continuous_scale='rainbow',
#                         opacity=0.5)
# fig.show()




# Load your data
# df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt', sep=' ')

# # Create a GeoDataFrame assuming you suspect Lambert 93
# gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs="EPSG:2154")

# # Plot
# gdf.plot()
# plt.show()

def categorize_weights(df, weight_col, num_bins, labels=None):
    """
    Categorizes weights in a DataFrame into specified bins.

    Parameters:
    - df: The DataFrame containing the weight data.
    - weight_col: The column name containing weight values.
    - num_bins: The number of bins to create.
    - labels: Optional list of labels for the bins. If None, generates labels.

    Returns:
    - The DataFrame with a new column for weight categories.
    """
    # If labels are not provided, create default labels
    if labels is None:
        labels = [f'Category {i + 1}' for i in range(num_bins)]

    # Use pd.qcut to categorize weights into quantiles
    df['weight_category'] = pd.qcut(df[weight_col], q=num_bins, labels=labels, duplicates='drop')

    return df


########################### WORKING EXAMPLE ###########################

# # Read the data
# # df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037_shuffle.txt', sep=' ')
# df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt', sep=' ')


# # Create a GeoDataFrame with Lambert 93 coordinates (EPSG:2154)
# gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs="EPSG:900913")

# # Reproject to WGS84 (EPSG:4326)
# gdf = gdf.to_crs(epsg=4326)

# # Add transformed coordinates back to the original DataFrame
# df['longitude'] = gdf.geometry.x
# df['latitude'] = gdf.geometry.y

# # # Weight categorization
# # quantiles = df['weight'].quantile([0, 0.25, 0.5, 0.75, 1.0]).values
# # labels = ['Very Low', 'Low', 'Medium', 'High']
# # df['weight_category'] = pd.cut(df['weight'], bins=quantiles, labels=labels, include_lowest=True)

# # Categorize weights into 10 bins
# # df = categorize_weights(df, weight_col='weight', num_bins=50)

# # Categorize into quantile-based bins
# # df['weight_category'] = pd.qcut(df['weight'], q=10, labels=[f'Bin {i+1}' for i in range(10)], duplicates='drop')

# # Example data
# # df['log_weight'] = np.log1p(df['weight'])  # log1p for log(x + 1) to avoid log(0)
# # # Create bins based on the log-transformed weights
# # log_bins = np.linspace(df['log_weight'].min(), df['log_weight'].max(), 6)  # 10 bins
# # df['weight_category'] = pd.cut(df['log_weight'], bins=log_bins, labels=[f'Bin {i+1}' for i in range(5)], include_lowest=True)

# # Define bins that cover the entire range of weights
# bins = [0, 1000, 5000, 10000, 50000, 100000, 200000]  # Adjust as needed based on your data
# labels = ['Very Low', 'Low', 'Medium', 'High', 'Very High', 'Extreme']

# # Ensure that every weight has a corresponding category
# df['weight_category'] = pd.cut(df['weight'], bins=bins, labels=labels, include_lowest=True)

# # Count the number of occurrences in each category
# category_counts = df['weight_category'].value_counts()

# # Check which categories are missing
# missing_categories = set(labels) - set(category_counts.index)

# # Create a DataFrame for missing categories and add them to the original DataFrame
# for category in missing_categories:
#     # Create a new row with NaN values for weight and the missing category
#     new_row = {'weight': np.nan, 'weight_category': category}
#     df = df.append(new_row, ignore_index=True)

# # Display the DataFrame with all categories
# print(df[['weight', 'weight_category']].dropna().head())


# # Calculate the center point for the map
# center_lat = df['latitude'].mean()
# center_lon = df['longitude'].mean()

# # Load the shapefile
# shapefile_path = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
# region_gdf = gpd.read_file(shapefile_path, crs="EPSG:2154")

# # Reproject the shapefile to WGS84 (EPSG:4326)
# region_gdf = region_gdf.to_crs(epsg=4326)

# # Plot using Plotly with customized center
# fig = px.density_mapbox(df, lat='latitude', lon='longitude', z='weight',
#                         #  radius=16,  # Adjusted radius for 5 km
#                          radius=22,  # Adjusted radius for 5 km
#                          zoom=7,
#                          center=dict(lat=center_lat, lon=center_lon),
#                          mapbox_style='open-street-map',
#                         #  mapbox_style = 'carto-darkmatter',
#                         #  mapbox_style = 'stamen-terrain',
#                         #  color_continuous_scale='rainbow',
#                          color_continuous_scale=px.colors.sequential.Rainbow,
#                          opacity=0.5)

# # # Add the shapefile boundaries to the plot
# # for _, row in region_gdf.iterrows():
# #     if row.geometry.geom_type == 'Polygon':
# #         x, y = row.geometry.exterior.xy
# #         fig.add_scattermapbox(lon=x, lat=y, mode='lines', line=dict(width=2, color='black'))

# # Show the plot
# fig.show()



#################################################################################################################################################


# Read the data
df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037_shuffle.txt', sep=' ')
# df = pd.read_csv('/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt', sep=' ')

# Create a GeoDataFrame with Lambert 93 coordinates (EPSG:2154)
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs="EPSG:900913")
# Reproject to WGS84 (EPSG:4326)
gdf = gdf.to_crs(epsg=4326)

# Add transformed coordinates back to the original DataFrame
df['longitude'] = gdf.geometry.x
df['latitude'] = gdf.geometry.y

# Define bins that cover the entire range of weights
bins = [0, 1000, 5000, 10000, 50000, 100000, 200000]  # Adjust as needed based on your data
labels = ['Very Low', 'Low', 'Medium', 'High', 'Very High', 'Extreme']
# Ensure that every weight has a corresponding category
df['weight_category'] = pd.cut(df['weight'], bins=bins, labels=labels, include_lowest=True)

# Count the number of occurrences in each category
category_counts = df['weight_category'].value_counts()

# Check which categories are missing
missing_categories = set(labels) - set(category_counts.index)
# Create a DataFrame for missing categories and add them to the original DataFrame
for category in missing_categories:
    # Create a new row with NaN values for weight and the missing category
    new_row = {'weight': np.nan, 'weight_category': category}
    df = df.append(new_row, ignore_index=True)
# Display the DataFrame with all categories
print(df[['weight', 'weight_category']].dropna().head())

# Calculate the center point for the map
center_lat = df['latitude'].mean()
center_lon = df['longitude'].mean()

# Load the shapefile
shapefile_path = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
region_gdf = gpd.read_file(shapefile_path, crs="EPSG:2154")
# Reproject the shapefile to WGS84 (EPSG:4326)
region_gdf = region_gdf.to_crs(epsg=4326)


# Assuming the shapefile contains a column 'region_name' for filtering
highlighted_region = region_gdf[region_gdf['layer'] == 'Region_to_Highlight']

# Plot the density map for the data points
fig = px.density_mapbox(df, lat='latitude', lon='longitude', z='weight',
                        radius=22,
                        zoom=7,
                        center=dict(lat=center_lat, lon=center_lon),
                        # mapbox_style='open-street-map',
                        # mapbox_style='carto-darkmatter',
                        mapbox_style='carto-positron',
                        color_continuous_scale=px.colors.sequential.Rainbow,
                        opacity=0.5)

# Add the shapefile boundaries in color
for _, row in region_gdf.iterrows():
    geom = row.geometry
    if geom.geom_type == 'Polygon':
        x, y = geom.exterior.xy
        fig.add_scattermapbox(lon=list(x), lat=list(y), mode='lines', line=dict(width=2, color='black'), showlegend=False)
    elif geom.geom_type == 'MultiPolygon':
        for part in geom.geoms:
            x, y = part.exterior.xy
            fig.add_scattermapbox(lon=list(x), lat=list(y), mode='lines', line=dict(width=2, color='black'), showlegend=False)

# Show the plot
fig.show()
