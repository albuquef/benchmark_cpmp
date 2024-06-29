import pandas as pd
import geopandas as gpd
import re
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi
import osmnx as ox

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

def filter_by_type_equipment(df, type_equipment):
    return df[df['TYPEQU'] == type_equipment]

def load_service_points_data(filepath, type_equipment):
    """Loads grid data and converts CRS."""
    df = pd.read_csv(filepath, low_memory=False)
    df = filter_by_type_equipment(df, type_equipment)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['LAMBERT_X'], df['LAMBERT_Y']), crs='IGNF:RGF93LAMB93') # coordinates in LAMB93
    gdf = gdf.to_crs(epsg=3035) # Convert to EPSG:3035
    return gdf

def plot_population_and_service_points(gdf_population, gdf_service):
    fig, ax = plt.subplots()
    gdf_population.plot(ax=ax, color='blue', markersize=5, marker='o', label='Population Points')
    gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
    plt.legend()
    plt.show()

def create_voronoi_polygons(gdf_service):
    coords = np.array([[point.x, point.y] for point in gdf_service.geometry])
    vor = Voronoi(coords)
    
    def voronoi_regions(vor, shape):
        regions = []
        for region_index in vor.point_region:
            vertices = vor.regions[region_index]
            if all(v >= 0 for v in vertices):
                region = [vor.vertices[i] for i in vertices]
                poly = Polygon(region)
                if shape.contains(poly):
                    regions.append(poly)
        return regions
    
    min_x, min_y, max_x, max_y = coords[:, 0].min(), coords[:, 1].min(), coords[:, 0].max(), coords[:, 1].max()
    bounding_box = Polygon([(min_x, min_y), (min_x, max_y), (max_x, max_y), (max_x, min_y)])
    
    polygons = voronoi_regions(vor, bounding_box)
    voronoi_gdf = gpd.GeoDataFrame(geometry=polygons)
    return voronoi_gdf

def create_isochrones(gdf_service, travel_times):
    isochrones = []
    for idx, service_point in gdf_service.iterrows():
        center_point = (service_point.geometry.y, service_point.geometry.x)
        
        G = ox.graph_from_point(center_point, dist=max(travel_times) * 100, network_type='drive')
        G = ox.add_edge_speeds(G)
        G = ox.add_edge_travel_times(G)

        for tt in travel_times:
            subgraph = ox.truncate.truncate_graph_bbox(G, north=center_point[0]+0.1, south=center_point[0]-0.1, east=center_point[1]+0.1, west=center_point[1]-0.1)
            nodes = ox.distance.nearest_nodes(G, center_point[1], center_point[0])
            subgraph = ox.truncate.truncate_graph_bbox(G, north=center_point[0], south=center_point[0], east=center_point[1], west=center_point[1])
            isochrone = ox.plot_figure_ground(subgraph, bgcolor='k', edge_color='w', node_alpha=0.1, edge_linewidth=0.2, save=False, show=False)
            isochrones.append(isochrone)

    return isochrones

# Load population data
gdf_population = load_population_data('data/data_qgis/data_instance_paca/points_population_1km_paca_table.csv')

# Load service points data
gdf_service = load_service_points_data('data/data_qgis/data_instance_paca/bpe21_sport_loisir_xy_paca_table.csv', 'F303')

# Plot population and service points
plot_population_and_service_points(gdf_population, gdf_service)

# Create Voronoi polygons
voronoi_gdf = create_voronoi_polygons(gdf_service)

# Plot Voronoi polygons
fig, ax = plt.subplots()
voronoi_gdf.plot(ax=ax, edgecolor='black', facecolor='none', label='Voronoi Polygons')
gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
plt.legend()
plt.show()

# Create isochrones
travel_times = [300, 600, 900]  # in seconds (5, 10, 15 minutes)
isochrones = create_isochrones(gdf_service, travel_times)

# Plot isochrones (Note: This plotting might require adaptation based on how isochrones are returned)
fig, ax = plt.subplots()
for isochrone in isochrones:
    isochrone.plot(ax=ax)
gdf_service.plot(ax=ax, color='orange', markersize=30, marker='^', label='Service Points')
plt.legend()
plt.show()
