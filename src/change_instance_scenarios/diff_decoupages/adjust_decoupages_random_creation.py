import geopandas as gpd
import os

def count_subareas_per_points_file(shapefile_paths, points_gdfs):
    """
    Print the number of subareas in each shapefile covered by each points GeoDataFrame.

    Parameters:
    shapefile_paths (list of str): List of paths to the shapefiles (subregions).
    points_gdfs (list of GeoDataFrames): List of GeoDataFrames representing different sets of points.
    """
    for shapefile_path in shapefile_paths:
        subregions_gdf = gpd.read_file(shapefile_path).to_crs(epsg=3395)
        print(f"\nProcessing shapefile: {shapefile_path}")

        for i, points_gdf in enumerate(points_gdfs):
            joined = gpd.sjoin(points_gdf, subregions_gdf, predicate='within')
            # joined = gpd.sjoin(subregions_gdf, points_gdf, how="inner", predicate='contains')
            subareas_with_points = set(joined['index_right'].unique())
            print(f"  Points file {i+1} covers {len(subareas_with_points)} subareas.")


def filter_and_save_subareas(shapefile_paths, points_gdfs):
    """
    Filter subareas with no points for each points file and save to a new shapefile.

    Parameters:
    shapefile_paths (list of str): List of paths to the shapefiles (subregions).
    points_gdfs (list of GeoDataFrames): List of GeoDataFrames representing different sets of points.
    """
    point_labels = ["2km", "5km"]  # Labels for the points files
    
    for i, points_gdf in enumerate(points_gdfs):
        label = point_labels[i]
        
        for shapefile_path in shapefile_paths:
            subregions_gdf = gpd.read_file(shapefile_path).to_crs(epsg=3395)

            # Reset index to avoid issues with non-matching indices
            subregions_gdf = subregions_gdf.reset_index(drop=True)

            # Perform spatial join and filter subareas with points
            # joined = gpd.sjoin(subregions_gdf, points_gdf, how="inner", predicate='contains')
            # joined = gpd.sjoin(points_gdf, subregions_gdf, predicate='within')

            # Validate geometries
            points_gdf = points_gdf[points_gdf.is_valid]
            subregions_gdf = subregions_gdf[subregions_gdf.is_valid]

            # # Print CRS and geometry validity
            # print(f"Points CRS: {points_gdf.crs}")
            # print(f"Subregions CRS: {subregions_gdf.crs}")
            # print(f"Points valid: {points_gdf.is_valid.all()}")
            # print(f"Subregions valid: {subregions_gdf.is_valid.all()}")

            # Perform spatial join
            joined = gpd.sjoin(points_gdf, subregions_gdf, predicate='within')

            # # Print the result to debug
            # print(joined.head())
            # print(joined.columns)


            if 'index_right' in joined.columns:
                # Get unique indices
                indices = joined['index_right'].unique()

                # Filter out indices that are out of bounds
                valid_indices = [idx for idx in indices if 0 <= idx < len(subregions_gdf)]

                if valid_indices:
                    # Use iloc to filter valid positions
                    filtered_subregions = subregions_gdf.iloc[valid_indices]

                    # Determine the output filename
                    num_areas_covered = len(filtered_subregions)
                    output_filename = f"voronoi_region_{label}_{num_areas_covered}.shp"
                    output_path = os.path.join("outputs/filtered_voronoi", output_filename)

                    # Save filtered subareas to a new file
                    filtered_subregions.to_file(output_path)
                    print(f"Saved {output_filename} with {num_areas_covered} subareas.")
                else:
                    print(f"Warning: No valid subareas covered for {shapefile_path} with points file {label}.")
            else:
                print(f"Warning: No subareas covered for {shapefile_path} with points file {label}.")

# Example usage
shapefiles = [
    '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_voronoi/voronoi_Arrondisement.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_50.shp',
    '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_voronoi/voronoi_EPCI.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_188.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_192.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_voronoi_978.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_voronoi_1041.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_voronoi_1058.shp',
    'outputs/decoupages_voronoi/PACA_subdivisions_voronoi_1133.shp'
]

points = [
    gpd.read_file("/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/grid_2km/grid_paca_points_2km.shp").to_crs(epsg=3395),
    gpd.read_file("/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/grid_5km/points_5km_2037_paca.shp").to_crs(epsg=3395)
]

count_subareas_per_points_file(shapefiles, points)

filter_and_save_subareas(shapefiles, points)