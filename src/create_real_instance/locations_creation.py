import pandas as pd
import geopandas as gpd
import re
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, cKDTree
from geovoronoi import voronoi_regions_from_coords, points_to_region
import osmnx as ox



def check_crs_consistency(gdf1, gdf2):
    if gdf1.crs != gdf2.crs:
        raise ValueError("CRS mismatch between GeoDataFrames. Ensure both GeoDataFrames use the same CRS.")

# Function to extract coordinates from the 'idcar_1km' column
def extract_coordinates(idcar):
    match = re.search(r'N(\d+)E(\d+)', idcar)
    if match:
        northing = float(match.group(1))
        easting = float(match.group(2))
        return easting, northing
    else:
        return None, None

def load_population_data(filepath):
    """Loads population data and extracts coordinates."""
    df = pd.read_csv(filepath)
    df['easting'], df['northing'] = zip(*df['idcar_1km'].apply(extract_coordinates))
    df = df.dropna(subset=['easting', 'northing'])
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['easting'], df['northing']))
    gdf.set_crs(epsg=3395, allow_override=True)
    return gdf

def filter_by_type_equipment(df, type_equipment):
    return df[df['TYPEQU'] == type_equipment]

def load_service_points_data(filepath, type_equipment):
    """Loads grid data and converts CRS."""
    df = pd.read_csv(filepath, low_memory=False)
    df = filter_by_type_equipment(df, type_equipment)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['LAMBERT_X'], df['LAMBERT_Y']), crs='IGNF:RGF93LAMB93') # coordinates in LAMB93
    gdf = gdf.to_crs(epsg=3035) # Convert to EPSG:3035
    return gdf

# not covering all services points
def create_voronoi(gdf_population, gdf_service, delta=10000):
    
    # area of interest
    # extent of the area of interest with a buffer of delta
    xmin, ymin, xmax, ymax = gdf_population.total_bounds
    lat_point_list = [ymin-delta, ymax+delta, ymax+delta, ymin-delta, ymin-delta]
    lon_point_list = [xmin-delta, xmin-delta, xmax+delta, xmax+delta, xmin-delta]
    polygon_geom = Polygon(zip(lon_point_list, lat_point_list))
    
    # Create a GeoDataFrame with the boundary polygon
    boundary_gdf = gpd.GeoDataFrame(geometry=[polygon_geom], crs='EPSG:3395')
    
    # Set the CRS for gdf_service
    # gdf_service.set_crs(epsg=3395, allow_override=True)
    
    # Create Voronoi diagram
    coords = np.array([[p.x, p.y] for p in gdf_service.geometry])
    print(f'Number of service points: {len(coords)}')
    vor = Voronoi(coords)
  
    polygons = []
    for region in vor.regions:
        if not -1 in region and len(region) > 0:  # Ensure the region is valid
            polygon = Polygon([vor.vertices[i] for i in region])
            polygons.append(polygon)        
            
    voronoi_gdf = gpd.GeoDataFrame(geometry=polygons, crs=gdf_service.crs)
    
    # Clip Voronoi polygons to the boundary
    clipped_gdf = gpd.clip(voronoi_gdf, boundary_gdf)
    
    # Plotting
    fig, ax = plt.subplots(figsize=(10, 10))
    clipped_gdf.plot(ax=ax, edgecolor='black', facecolor='red', alpha=0.5)
    boundary_gdf.boundary.plot(ax=ax, color='red')
    voronoi_gdf.plot(ax=ax, edgecolor='black', facecolor='green', alpha=0.5)
    gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
    gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
    plt.title(f'Voronoi Diagram ({len(voronoi_gdf)}) with Service Points')
    plt.xlabel('Easting')
    plt.ylabel('Northing')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    
    return voronoi_gdf
   
def create_buffer_regions(gdf_service, buffer_radius):
    buffer_regions = []

    for idx, service_point in gdf_service.iterrows():
        buffer_region = service_point.geometry.buffer(buffer_radius)
        buffer_regions.append(buffer_region)

    # Convert to GeoDataFrame
    buffer_gdf = gpd.GeoDataFrame(geometry=buffer_regions, crs=gdf_service.crs)

    return buffer_gdf

def plot_population_and_service_points(gdf_population, gdf_service):
    fig, ax = plt.subplots()
    gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
    gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
    plt.legend()
    plt.show()

def plot_voronoi_polygons(gdf_population, voronoi_gdf, gdf_service, delta=10000):
    fig, ax = plt.subplots()
    gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
    voronoi_gdf.plot(ax=ax, edgecolor='black', facecolor='green', alpha=0.5, label='Voronoi Polygons')
    gdf_service.plot(ax=ax, color='orange', markersize=20, marker='^', label='Service Points')
    # cut the size of the plot using the population points
    plt.xlim([gdf_population.total_bounds[0]-delta, gdf_population.total_bounds[2]+delta])
    plt.ylim([gdf_population.total_bounds[1]-delta, gdf_population.total_bounds[3]+delta])
    plt.legend()
    plt.title('Voronoi Polygons for Service Points')
    plt.show()

def plot_buffer_regions(gdf_population, buffer_gdf, gdf_service):
    fig, ax = plt.subplots()
    gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
    buffer_gdf.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1.5, label='Buffer Regions')
    gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
    plt.legend()
    plt.title('Buffer Regions around Service Points')
    plt.xlabel('Easting')
    plt.ylabel('Northing')
    plt.show()
    
def create_locations_final_table(cap_polygons):
    # read a txt file with the locations of the service points
    # customer weight coord_x coord_y fid
    # 0 11.0 3859533.886436741 2284448.02565863 51
    # get the coordinates of the service points
    df = pd.read_csv('outputs/PACA_Jun2024/cust_weights_Jun2024', sep=' ')   
    
    # create intervals for the capacities using cap_polygons
    # cap_intervals = []
    # cap_intervals.append(0)
    # for cap in cap_polygons:
    #     cap_intervals.append(cap)
    
    
    
# Load population data
gdf_population = load_population_data('data/data_qgis/data_instance_paca/points_population_1km_paca_table.csv')
# # Load service points data
gdf_service = load_service_points_data('data/data_qgis/data_instance_paca/bpe21_sport_loisir_xy_paca_table.csv', 'F303')
# plot_population_and_service_points(gdf_population, gdf_service)

# create_voronoi(gdf_population, gdf_service, delta=10000) # not covering all services points
# read a polygon shapefile as a GeoDataFrame
gdf_voronoi = gpd.read_file('/home/felipe/Documents/Projects/GeoAvigon/qgis/PACA_instance_qgis/voronoi_paca_bpe_F303.shp')
gdf_voronoi = gdf_voronoi.to_crs(epsg=3035)

# # plotting
fig, ax = plt.subplots()
gdf_voronoi.plot(ax=ax, edgecolor='black', facecolor='green', alpha=0.5)
gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
# plt.title('Voronoi Polygons for Service Points')
plt.show()

# compare sizes
print(len(gdf_voronoi))
print(len(gdf_service))

sum_points = 0
sum_weight = 0
cap_polygons = []
for idx, voronoi_polygon in gdf_voronoi.iterrows():
    cont = 0
    sum_weight_polygon = 0
    for idx_pop, point_pop in gdf_population.iterrows() :
        if voronoi_polygon.geometry.contains(point_pop.geometry):
            cont += 1
            # access the point information "ind" in the gdf_population
            sum_weight_polygon += float(point_pop.ind)
    print(f"Number of points in Voronoi Polygon {idx}: {cont}")
    sum_points += cont
    sum_weight += sum_weight_polygon
    cap_polygons.append(sum_weight_polygon)
    
print(f"Total number of points in Voronoi Polygons: {sum_points}")
print(f"Total weight of points in Voronoi Polygons: {sum_weight}")
print(f"Total number of points in Population Points: {len(gdf_population)}")
print(f"Total weight of points in Population Points: {gdf_population['ind'].sum()}")

# plot cap_polygon distribution
fig, ax = plt.subplots()
ax.hist(cap_polygons, bins=30)
plt.title('Distribution Capacities in Voronoi Polygons')
plt.xlabel('Capacities')
plt.ylabel('Frequency')
plt.show()

# plot each polygon voronoi with the label of the cap_polygon
fig, ax = plt.subplots()
gdf_voronoi.plot(ax=ax, edgecolor='black', facecolor='green', alpha=0.5)
# gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
# gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
# add cap_polygon label
gdf_voronoi['cap_polygon'] = cap_polygons
for idx, voronoi_polygon in gdf_voronoi.iterrows():
    ax.text(voronoi_polygon.geometry.centroid.x, voronoi_polygon.geometry.centroid.y, voronoi_polygon['cap_polygon'], fontsize=8)
plt.show()

