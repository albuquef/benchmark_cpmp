import geopandas as gpd
import numpy as np
from shapely.geometry import Point, Polygon, box
from scipy.spatial import Voronoi
from shapely.ops import unary_union
import matplotlib.pyplot as plt
import logging

# Function to generate K subdivisions inside a region
def divide_region(shapefile_path, K, gap=0.1, seed=None, max_iterations=10):
    # Set the seed for reproducibility
    if seed is not None:
        np.random.seed(seed)
    
    # Load the region shapefile
    region_gdf = gpd.read_file(shapefile_path)
    
    # Ensure the region is in a projected coordinate system for accurate distances
    region_gdf = region_gdf.to_crs(epsg=3395)  # EPSG:3395 for meter-based coordinates

    # Get the bounding box of the region for random point generation
    bounds = region_gdf.geometry.total_bounds
    minx, miny, maxx, maxy = bounds
    
    # Add a gap to the bounding box
    minx -= gap
    miny -= gap
    maxx += gap
    maxy += gap
    
    # Generate initial K random points inside the bounding box
    points = []
    while len(points) < K:
        x = np.random.uniform(minx, maxx)
        y = np.random.uniform(miny, maxy)
        point = Point(x, y)
        if region_gdf.contains(point).any():  # Ensure point is inside the region
            points.append(point)
    
    # Add additional points around the bounding box to ensure all Voronoi regions are bounded
    extra_points = [
        Point(minx, miny), Point(minx, maxy), Point(maxx, miny), Point(maxx, maxy),
        Point(minx, (miny + maxy) / 2), Point(maxx, (miny + maxy) / 2),
        Point((minx + maxx) / 2, miny), Point((minx + maxx) / 2, maxy)
    ]
    points.extend(extra_points)
    
    # Iteratively add points to cover the entire region
    for _ in range(max_iterations):
        # Convert points into a numpy array for Voronoi computation
        points_array = np.array([[p.x, p.y] for p in points])
        
        # Create Voronoi polygons
        vor = Voronoi(points_array)
        
        # Create a bounding box polygon with the added gap
        bbox_polygon = Polygon([(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)])
        
        # List to hold valid polygons
        polygons = []
        
        # Create polygons from Voronoi regions
        for region in vor.regions:
            if not -1 in region and len(region) > 0:
                polygon = Polygon([vor.vertices[i] for i in region])
                clipped_polygon = polygon.intersection(bbox_polygon)
                if not clipped_polygon.is_empty:
                    polygons.append(clipped_polygon)
        
        # Clip the polygons to the region
        clipped_polygons = [polygon.intersection(region_gdf.unary_union) for polygon in polygons]
        
        # Check if the entire region is covered
        uncovered_area = region_gdf.difference(unary_union(clipped_polygons))
        if uncovered_area.is_empty.any():
            break
        
        # Add points in the uncovered areas
        for geom in uncovered_area.geometry:
            if not geom.is_empty:
                x, y = geom.centroid.x, geom.centroid.y
                points.append(Point(x, y))
    
    # Ensure no intersections between polygons
    non_intersecting_polygons = []
    for polygon in clipped_polygons:
        non_intersecting_polygon = polygon.difference(unary_union(non_intersecting_polygons))
        if not non_intersecting_polygon.is_empty:
            non_intersecting_polygons.append(non_intersecting_polygon)
    
    # Create a GeoDataFrame from the non-intersecting polygons
    subdivisions_gdf = gpd.GeoDataFrame(geometry=non_intersecting_polygons, crs=region_gdf.crs)
    
    return region_gdf, subdivisions_gdf, points

# Function to generate a grid of subdivisions
def divide_region_grid(shapefile_path, cell_size, gap=0.1):
    region_gdf = gpd.read_file(shapefile_path).to_crs(epsg=3395)
    
    # Get bounding box with gap
    minx, miny, maxx, maxy = region_gdf.geometry.total_bounds
    minx -= gap
    miny -= gap
    maxx += gap
    maxy += gap

    # Create grid cells
    x_coords = np.arange(minx, maxx, cell_size)
    y_coords = np.arange(miny, maxy, cell_size)
    grid_cells = [
        box(x, y, x + cell_size, y + cell_size) for x in x_coords for y in y_coords
    ]

    # Clip grid to region
    grid_gdf = gpd.GeoDataFrame(geometry=grid_cells, crs="EPSG:3395")
    clipped_grid = gpd.overlay(grid_gdf, region_gdf, how='intersection')

    return clipped_grid


def count_regions_with_points(subdivisions_gdf, points_gdf):
    count = subdivisions_gdf.intersects(points_gdf.unary_union).sum()
    return count



# save log file
save_log = True
if save_log:
    log_filename = "subdivision_log_grid.txt"
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_filename),
                            logging.StreamHandler()  # Prints to console
                        ])

seed = 42  # Seed for reproducibility
shapefile_path = '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/PACA_region.shp'  # Path to your shapefile
# Load points and region shapefile
POINTS_2KM = gpd.read_file("/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/grid_2km/grid_paca_points_2km.shp").to_crs(epsg=3395)
POINTS_5KM = gpd.read_file("/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/grid_5km/points_5km_2037_paca.shp").to_crs(epsg=3395)


seed = 42  # For Voronoi reproducibility

# for K in range(1040, 1045):  
#     # Voronoi-based subdivision
#     region_gdf, subdivisions_gdf_voronoi, points = divide_region(shapefile_path, K, seed=seed)

#     # Save and count points
#     subdivisions_gdf_voronoi.to_file(f"outputs/decoupages_voronoi/PACA_subdivisions_voronoi_{K}.shp")
#     count_2km_voronoi = count_regions_with_points(subdivisions_gdf_voronoi, POINTS_2KM)
#     count_5km_voronoi = count_regions_with_points(subdivisions_gdf_voronoi, POINTS_5KM)

#     print('-------------------------------------------------')
#     print(f"Voronoi (K={K}): {len(subdivisions_gdf_voronoi)} subdivisions")
#     print(f"Points 2KM: {count_2km_voronoi},\n Points 5KM: {count_5km_voronoi}")
#     print('-------------------------------------------------')

#     if save_log:
#         # Logging information
#         logging.info('-------------------------------------------------')
#         logging.info(f"Voronoi (K={K}): {len(subdivisions_gdf_voronoi)} subdivisions")
#         logging.info(f"Points 2KM: {count_2km_voronoi}")
#         logging.info(f"Points 5KM: {count_5km_voronoi}")
#         logging.info('-------------------------------------------------')


grid_cells = [41.4, 19.9, 7.958, 8.2712]
# multiply by 1000 to get km
grid_cells = [x*1000 for x in grid_cells]

# start = 8.271 # 40 km
# finish = 8.37  # 90 km
# interval = 0.001 #10km
# for cell_size in range(int(start*1000), int(finish*1000),int(interval*1000)):
for cell_size in grid_cells:
    # Grid-based subdivision
    # subdivisions_gdf_grid = divide_region_grid(shapefile_path, cell_size+0.28123)
    subdivisions_gdf_grid = divide_region_grid(shapefile_path, cell_size)

    # Save and count points
    subdivisions_gdf_grid.to_file(f"outputs/decoupages_grid/PACA_subdivisions_grid_{int(cell_size/1000)}km.shp")
    count_2km_grid = count_regions_with_points(subdivisions_gdf_grid, POINTS_2KM)
    count_5km_grid = count_regions_with_points(subdivisions_gdf_grid, POINTS_5KM)

    print(f"Grid (Cell Size={cell_size}m): {len(subdivisions_gdf_grid)} subdivisions")
    print(f"Points 2KM: {count_2km_grid},\n Points 5KM: {count_5km_grid}")

    # Logging information
    logging.info('-------------------------------------------------')
    logging.info(f"Grid (Cell Size={cell_size}m): {len(subdivisions_gdf_grid)} subdivisions")
    logging.info(f"Points 2KM: {count_2km_grid}")
    logging.info(f"Points 5KM: {count_5km_grid}")
    logging.info('-------------------------------------------------')

