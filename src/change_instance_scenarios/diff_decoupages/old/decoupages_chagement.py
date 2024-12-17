import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.validation import explain_validity


# load the grids files
def load_shapefiles(shapefiles):
    gdfs = []
    for shapefile in shapefiles:
        gdf = gpd.read_file(shapefile)
        gdf = gdf.to_crs(epsg=3035)
        gdf['index_polygon'] = range(1, len(gdf) + 1)
        gdfs.append(gdf)
    return gdfs

def load_locs_points(shapefile):
    gdf = gpd.read_file(shapefile)
    gdf['index_point'] = range(1, len(gdf) + 1)
    gdf = gdf.to_crs(epsg=3035)
    return gdf

def plot_polygons_and_points(gdf_polygons, gdf_points, representation_name):
    fig, ax = plt.subplots()
    gdf_polygons.plot(ax=ax, color='white', edgecolor='black')
    gdf_points.plot(ax=ax, color='orange', markersize=10, marker='o')   
    plt.title(f"Grid representing {representation_name} ({len(gdf_polygons)})")  
    plt.savefig(f"plots/diff_decoupages/grid/{representation_name}.png")
    plt.show()

def check_and_fix_invalid_geometries(gdf):
    valid_geometries = gdf.geometry.apply(lambda geom: geom.is_valid)
    if not valid_geometries.all():
        invalid_geoms = gdf[~valid_geometries]
        for idx, row in invalid_geoms.iterrows():
            print(f"Invalid geometry at index {idx}: {explain_validity(row.geometry)}")
        gdf.geometry = gdf.geometry.buffer(0)  # Fix invalid geometries
    return gdf

def create_table_vet_sum_loc_coverages(gdf_polygons, gdf_points,typepolygon):
    
    
    # Ensure the same CRS
    if gdf_points.crs != gdf_polygons.crs:
        gdf_polygons = gdf_polygons.to_crs(gdf_points.crs)
    
    # Check and fix geometries
    gdf_polygons = check_and_fix_invalid_geometries(gdf_polygons)
    
    print(f"gdf_polygons columns: {gdf_polygons.columns}")
    print(f"head gdf_polygons:\n {gdf_polygons.head()}")
    
    vet_sum_loc_coverages = []
    for i in range(len(gdf_polygons)):
        vet_sum_loc_coverages.append(gdf_points.within(gdf_polygons.loc[i, 'geometry']).sum())
    
    # compare the cum loc coverages points and the number of points
    sum_vet_sum_loc_coverages = sum(vet_sum_loc_coverages)
    num_points = len(gdf_points)
    if sum_vet_sum_loc_coverages != num_points:
        print(f"Error: the sum of loc coverages ({sum_vet_sum_loc_coverages}) is different from the number of points ({num_points})")

    print(f"gdf_points columns: {gdf_points.columns}")

    filename = f"data/PACA/loc_coverages_{typepolygon}.txt"
    with open(filename, 'w') as f:
        
        # Write the header
        f.write("location subarea coord_x coord_y\n")
        # Iterate over each point in gdf_points
        for idx_point, point in gdf_points.iterrows():
            point_geometry = point.geometry

            # Iterate over each polygon in gdf_polygons
            for idx_polygon, polygon in gdf_polygons.iterrows():
                polygon_geometry = polygon.geometry
        
                # Check if the point is inside the polygon
                if polygon_geometry.contains(point_geometry):
                    # Write the point information to the file
                    f.write(f"{int(point['index_point'])} {int(polygon['index_polygon'])} {point['geometry'].x} {point['geometry'].y}\n")
                    # Optionally, you can break the loop if you only need to know if it is inside any polygon
                    # break
    


shapefile_grid_arrond = "/home/felipe/Documents/Projects/GeoAvigon/qgis/PACA_instance_qgis/arrondissement_grid60.shp"
shapefile_grid_epci = "/home/felipe/Documents/Projects/GeoAvigon/qgis/PACA_instance_qgis/EPCI_grid30.shp"
shapefile_grid_canton = "/home/felipe/Documents/Projects/GeoAvigon/qgis/PACA_instance_qgis/canton_grid15.shp"
shapefile_grid_commune = "/home/felipe/Documents/Projects/GeoAvigon/qgis/PACA_instance_qgis/commune_grid6_3.shp"
shapefiles_grids = [shapefile_grid_arrond, shapefile_grid_epci, shapefile_grid_canton, shapefile_grid_commune]
names_grids = ["arrond", "EPCI", "canton", "commune"]
# shapefiles_grids = [shapefile_grid_commune]
# names_grids = ["commune"]

gdf_grid_locs_points = load_locs_points("/home/felipe/Documents/Projects/GeoAvigon/qgis/PACA_instance_qgis/location.shp")


gdfs_grid = load_shapefiles(shapefiles_grids)
cont = 0
# plot each grid
for gdf_grid in gdfs_grid:
    
    # create_table_vet_sum_loc_coverages(gdf_grid, gdf_grid_locs_points,f"grid_{names_grids[cont]}")
    plot_polygons_and_points(gdf_grid, gdf_grid_locs_points, names_grids[cont])
    cont += 1