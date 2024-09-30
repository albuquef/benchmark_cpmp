import pandas as pd
import geopandas as gpd
import plotly.express as px
import numpy as np

# Global variable for Mapbox style
MAPBOX_STYLE = 'carto-positron'

# Read the data
def load_data(filepath):
    return pd.read_csv(filepath, sep=' ')

# Create GeoDataFrame and reproject to WGS84
def create_geodataframe(df):
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs="EPSG:900913")
    return gdf.to_crs(epsg=4326)

# Add transformed longitude and latitude columns to the DataFrame
def add_transformed_coordinates(df, gdf):
    df['longitude'] = gdf.geometry.x
    df['latitude'] = gdf.geometry.y

# Add shapefile boundaries with color highlighting to the plot
# Add shapefile boundaries without filling polygons, just highlighting the boundaries
def add_shapefile_boundaries(fig, region_gdf, line_color='black'):
    for _, row in region_gdf.iterrows():
        geom = row.geometry
        
        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            fig.add_scattermapbox(
                lon=list(x), 
                lat=list(y), 
                mode='lines', 
                line=dict(width=2, color=line_color), 
                showlegend=False
            )
        elif geom.geom_type == 'MultiPolygon':
            for part in geom.geoms:
                x, y = part.exterior.xy
                fig.add_scattermapbox(
                    lon=list(x), 
                    lat=list(y), 
                    mode='lines', 
                    line=dict(width=2, color=line_color), 
                    showlegend=False
                )

# Load shapefile and reproject to WGS84
def load_shapefile(filepath):
    region_gdf = gpd.read_file(filepath, crs="EPSG:2154")
    return region_gdf.to_crs(epsg=4326)


# Spatial join to calculate the sum of weights inside each region
# Spatial join to calculate the sum of weights inside each region
def sum_weights_in_regions(df, points_gdf, regions_gdf):
    # Ensure the CRS of both GeoDataFrames match
    points_gdf = points_gdf.set_crs("EPSG:4326")
    regions_gdf = regions_gdf.set_crs("EPSG:4326")
    
    # Perform spatial join using `predicate` instead of deprecated `op`
    joined = gpd.sjoin(points_gdf, regions_gdf, how="inner", predicate="within")
    
    # Sum weights inside each region (assuming 'weight' is a column in your points data)
    region_weight_sum = joined.groupby('index_right').agg({'weight': 'sum'}).reset_index()
    
    # Merge the weight sum back to the regions_gdf
    regions_gdf['weight_sum'] = region_weight_sum.set_index('index_right')['weight']
    
    # Fill NaN values with 0 in case some regions have no points inside
    regions_gdf['weight_sum'].fillna(0, inplace=True)
    
    return regions_gdf

# Plot the choropleth map
def plot_choropleth_mapbox(regions_gdf, center_lat, center_lon):
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,  # Convert GeoDataFrame geometries to GeoJSON
        locations=regions_gdf.index,  # Use the index as the location identifier
        color='weight_sum',  # Use the summed weight for coloring
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        mapbox_style=MAPBOX_STYLE,
        color_continuous_scale="Viridis",  # Choose color scale
        opacity=0.6
    )
    
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
    fig.show()

# Main function to run the entire process
def main():
    # File paths
    data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    # data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037_shuffle.txt'
    shapefile_paca = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
    
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/Arrondissement.shp'  # Replace with actual path to your regions shapefile
    shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/EPCI.shp'
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/canton.shp'
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/commune.shp'
    
    
        # Load and prepare data
    df = load_data(data_path)
    gdf = create_geodataframe(df)
    add_transformed_coordinates(df, gdf)

    # Load shapefiles
    region_gdf_paca = load_shapefile(shapefile_paca)
    region_gdf_regions = load_shapefile(shapefile_regions)

    # Calculate the center point for the map
    center_lat = df['latitude'].mean()
    center_lon = df['longitude'].mean()

    plot_points = False
    # Plot points if `plot_points` is True
    if plot_points:
        fig = px.scatter_mapbox(
            df,
            lat='latitude',
            lon='longitude',
            hover_name='customer',  # Replace 'customer' with the actual column name for hover info
            zoom=7,  # Adjust the zoom level as needed
            center={"lat": center_lat, "lon": center_lon},
            mapbox_style=MAPBOX_STYLE
        )

        # Add shapefile boundaries for PACA
        add_shapefile_boundaries(fig, region_gdf_paca, line_color='black')

        # Add shapefile boundaries for the additional 'shp_regions'
        add_shapefile_boundaries(fig, region_gdf_regions, line_color='red')  # Different line color for distinction

        # Show the plot
        fig.show()

    # Calculate the sum of weights in each region and update the region GeoDataFrame
    region_gdf_regions = sum_weights_in_regions(df, gdf, region_gdf_regions)

    # Plot the choropleth map with the summed weights
    plot_choropleth_mapbox(region_gdf_regions, center_lat, center_lon)
    


# Execute the main function
if __name__ == "__main__":
    main()