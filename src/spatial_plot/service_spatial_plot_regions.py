import pandas as pd
import geopandas as gpd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np

# Global variable for Mapbox style
MAPBOX_STYLE = 'carto-positron'

# Read the data
def load_data(filepath):
    return pd.read_csv(filepath, sep=' ')

# Create GeoDataFrame and reproject to WGS84
def create_geodataframe(df, CRS="EPSG:900913"):
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs=CRS)
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
    # Define quantile intervals for 'weight_sum'
    # q = np.quantile(regions_gdf['weight_sum'], [0.25, 0.5, 0.75])
    
    # # Add a new column with the quantile classification
    # regions_gdf['weight_quantile'] = pd.cut(
    #     regions_gdf['weight_sum'],
    #     bins=[regions_gdf['weight_sum'].min(), q[0], q[1], q[2], regions_gdf['weight_sum'].max()],
    #     labels=["Low", "Medium-Low", "Medium-High", "High"]
    # )
    
    # num_bins = 8
    # # Define equal interval breaks for 'weight_sum'
    # min_val = regions_gdf['weight_sum'].min()
    # max_val = regions_gdf['weight_sum'].max()
    # bins = np.linspace(min_val, max_val, num_bins + 1)
    
    # # Add a new column with the equal interval classification
    # regions_gdf['weight_interval'] = pd.cut(
    #     regions_gdf['weight_sum'],
    #     bins=bins,
    #     labels=[f"Interval {i+1}" for i in range(num_bins)],
    #     include_lowest=True
    # )
    
    # bin_size = 50000
    # # Define intervals starting from min value and increasing by bin_size (100k)
    # min_val = regions_gdf['weight_sum'].min()
    # max_val = regions_gdf['weight_sum'].max()
    
    # # Create bins from min value to max value with 100k intervals
    # bins = np.arange(min_val, max_val + bin_size, bin_size)
    
    # # Add a new column with the custom interval classification
    # regions_gdf['weight_interval'] = pd.cut(
    #     regions_gdf['weight_sum'],
    #     bins=bins,
    #     include_lowest=True
    # )
    
    # Plot the choropleth map
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,  # Convert GeoDataFrame geometries to GeoJSON
        locations=regions_gdf.index,  # Use the index as the location identifier
        # color='weight_quantile',  # Use the quantile classification for coloring
        # color='weight_interval',  # Use the quantile classification for coloring
        color='weight_sum',
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        color_continuous_scale="Viridis",  # Choose a color scale
        # color_discrete_map={
        #     "Low": "lightblue",
        #     "Medium-Low": "lightgreen",
        #     "Medium-High": "orange",
        #     "High": "red"
        # },  # Define color for each quantile
        mapbox_style=MAPBOX_STYLE,
        opacity=0.6
    )
    
    
    # Update layout to include colorbar
    # fig.update_layout(
    #     margin={"r": 0, "t": 0, "l": 0, "b": 0},
    #     coloraxis_colorbar={
    #         'title': 'Weight Sum',
    #         'tickvals': bins,  # Show the bin edges as ticks on the colorbar
    #         'ticktext': [f'{int(val/1000)}k' for val in bins]  # Format tick labels in thousands (k)
    #     }
    # )
    
    
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
    fig.show()
    
# Read the service data (new file)
def load_service_data(filepath):
    df = pd.read_csv(filepath, sep=',', header=0)
    return df

# Plot the number of points inside each region
def count_points_in_regions(df, points_gdf, regions_gdf):
    # Ensure the CRS of both GeoDataFrames match
    points_gdf = points_gdf.set_crs("EPSG:4326")
    regions_gdf = regions_gdf.set_crs("EPSG:4326")
    
    
    points_gdf = points_gdf.to_crs(epsg=4326)
    regions_gdf = regions_gdf.to_crs(epsg=4326)

    
    # Perform spatial join to count points within regions
    # joined = gpd.sjoin(points_gdf, regions_gdf, how="inner", predicate="within")
    
    joined = gpd.sjoin(points_gdf, regions_gdf, how="inner", predicate="intersects")

    
    # Count the number of points inside each region
    region_point_count = joined.groupby('index_right').size().reset_index(name='point_count')
    
    # Merge the point count back to the regions_gdf
    regions_gdf['point_count'] = region_point_count.set_index('index_right')['point_count']
    
    # Fill NaN values with 0 in case some regions have no points inside
    regions_gdf['point_count'].fillna(0, inplace=True)
    
    return regions_gdf

def plot_point_count_mapbox(regions_gdf, center_lat, center_lon):
    # Create custom bins for point_count
    max_points = regions_gdf['point_count'].max()
    bins = [0, 1, 5] + list(range(10, int(max_points) + 5, 5))  # Create intervals (0-0, 1-4, 5-9, etc.)
    
    # Create bins and assign labels
    regions_gdf['point_count_bins'] = pd.cut(
        regions_gdf['point_count'], 
        bins=bins, 
        include_lowest=True, 
        labels=[f'{bins[i]}-{bins[i+1]-1}' for i in range(len(bins)-1)]
    )

    # Map colors to bins, ensuring '0-0' is gray
    color_discrete_map = {'0-0': 'gray'}
    
    # Use a subset of the Viridis color scale for the other bins
    viridis_colors = px.colors.sequential.Viridis[1:]  # Skip the first color
    for i, label in enumerate(regions_gdf['point_count_bins'].unique()[1:], 1):
        color_discrete_map[label] = viridis_colors[i % len(viridis_colors)]

    # Fill NaNs with '0-0' for regions without points
    regions_gdf['point_count_bins'].fillna('0-0', inplace=True)

    # Convert 'point_count_bins' to a categorical type
    regions_gdf['point_count_bins'] = regions_gdf['point_count_bins'].astype('category')

    # Plot with discrete colors
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,
        locations=regions_gdf.index,
        color='point_count_bins',  # Use the binned point count for coloring
        color_discrete_map=color_discrete_map,  # Apply custom color mapping
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        mapbox_style=MAPBOX_STYLE,
        opacity=0.6
    )
    
    # Update layout and fix the legend order
    fig.update_layout(
        margin={"r": 0, "t": 0, "l": 0, "b": 0},
        legend_title_text='Number of Points',
        coloraxis_colorbar=dict(
            title="Number of Points",
            tickvals=[i for i in range(len(bins))],  # Tick values correspond to bins
            ticktext=[f'{bins[i]}-{bins[i+1]-1}' for i in range(len(bins)-1)]  # Show intervals in the legend
        )
    )

    fig.show()


def plot_point_percentage_mapbox(regions_gdf, center_lat, center_lon, show_percent_labels=True):
    # Calculate the percentage of points inside each region
    total_points = regions_gdf['point_count'].sum()
    regions_gdf['point_percentage'] = (regions_gdf['point_count'] / total_points) * 100  # Calculate percentage
    
    # Set up a continuous color scale from light blue to dark blue
    color_scale = px.colors.sequential.Blues
    
    # Plot with continuous colors based on point percentage
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,
        locations=regions_gdf.index,
        color='point_percentage',  # Use the percentage of points for coloring
        color_continuous_scale=color_scale,  # Apply continuous color scale
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        mapbox_style=MAPBOX_STYLE,
        opacity=0.6
    )
    
    # Reproject geometries to a projected CRS (e.g., EPSG:2154 for Lambert 94) before calculating centroids
    regions_projected = regions_gdf.to_crs(epsg=2154)
    regions_projected['centroid'] = regions_projected.geometry.centroid  # Calculate centroids in the projected CRS

    # Reproject centroids back to the original CRS (EPSG:4326) for plotting
    centroids_geo = regions_projected['centroid'].to_crs(epsg=4326)
    
    if show_percent_labels:
        # Create a scatter layer for text labels using go.Scattermapbox
        scatter_data = go.Scattermapbox(
            lat=centroids_geo.y,  # Latitude from reprojected centroids
            lon=centroids_geo.x,  # Longitude from reprojected centroids
            mode='text',
            text=regions_gdf['point_percentage'].apply(lambda x: f"{x:.2f}%"),  # Format percentages
            textfont=dict(size=300000000, color='black'),  # Customize font size and color
            textposition="middle center",  # Place text at the center of each region
            marker=dict(size=100000000, color='black')
        )
        fig.add_trace(scatter_data)
    
    # Update layout and color bar settings
    fig.update_layout(
        margin={"r": 0, "t": 0, "l": 0, "b": 0},
        legend_title_text='Percentage of Points',
        coloraxis_colorbar=dict(
            title="Percentage of Service Points",
            ticks="outside",
            tickformat=".2f"  # Show two decimal places for percentages
        )
    )

    fig.show()

# Plot the choropleth map along with points
def plot_choropleth_mapbox_with_points(regions_gdf, points_gdf, center_lat, center_lon):
    # Plot the choropleth map
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,  # Convert GeoDataFrame geometries to GeoJSON
        locations=regions_gdf.index,  # Use the index as the location identifier
        color='point_count',  # Use the point count for coloring
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        color_continuous_scale="Viridis",  # Choose a color scale
        mapbox_style=MAPBOX_STYLE,
        opacity=0.6
    )
    
    # Add points on top of the choropleth
    fig.add_scattermapbox(
        lat=points_gdf.geometry.x,
        lon=points_gdf.geometry.y,
        mode='markers',
        # marker=px.scatter_mapbox.Marker(size=7, color='red', opacity=0.7),
        hoverinfo='text',
        text=points_gdf['fid'],  # Change this to whatever column contains point info
        showlegend=False
    )

    # Show the plot
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})
    fig.show()

# Main function to run the entire process
def main():
    # File paths
    # data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    data_path = '/home/falbuquerque/Documents/projects/GeoAvignon/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    # Load and prepare data
    df = load_data(data_path)
    gdf = create_geodataframe(df)
    add_transformed_coordinates(df, gdf)
    # Calculate the center point for the map
    center_lat = df['latitude'].mean()
    center_lon = df['longitude'].mean()
    
    # data_service_path = 'data/solutions_service_cinema/points_BPE23_F303_paca_table.csv'
    # df_service = load_service_data(data_service_path)
    # print(df_service.head())    
    # gdf_service = create_geodataframe(df_service,"EPSG:2154")  # Create a GeoDataFrame from service data
    
    data_service_path = 'data/solutions_service_cinema/test_paca_cinema_p_192_EXACT_CPMP_table.txt'
    # data_service_path = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton_table.txt'
    # data_service_path = 'data/solutions_service_cinema/test_paca_cinema_EPCI_p_192_EXACT_CPMP_cover_EPCI_table.txt'
    # data_service_path = 'data/solutions_service_cinema/test_paca_cinema_commune_p_192_EXACT_CPMP_cover_commune_table.txt'
    df_service = load_service_data(data_service_path)
    print(df_service.head())    
    gdf_service = create_geodataframe(df_service)  # Create a GeoDataFrame from service data
    
    
    
    # shapefile_paca = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
    shapefile_paca = '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/PACA_region.shp'

    shapefile_regions = '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/canton.shp'
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/Arrondissement.shp'  # Replace with actual path to your regions shapefile
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PAeCA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/EPCI.shp'
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/canton.shp'
    # shapefile_regions = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/commune.shp'

    # Load shapefiles
    region_gdf_paca = load_shapefile(shapefile_paca)
    region_gdf_regions = load_shapefile(shapefile_regions)


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


    # # Calculate the sum of weights in each region and update the region GeoDataFrame
    # region_gdf_regions = sum_weights_in_regions(df, gdf, region_gdf_regions)
    # # Plot the choropleth map with the summed weights
    # plot_choropleth_mapbox(region_gdf_regions, center_lat, center_lon)
    
    
    # Calculate the number of points in each region and update the region GeoDataFrame
    region_gdf_regions = count_points_in_regions(df_service, gdf_service, region_gdf_regions)
    # print(region_gdf_regions)
    # Plot the choropleth map with the point counts
    # plot_point_count_mapbox(region_gdf_regions, center_lat, center_lon)
    # plot_choropleth_mapbox_with_points(region_gdf_regions, gdf_service, center_lat, center_lon)
    plot_point_percentage_mapbox(region_gdf_regions, center_lat, center_lon)


# Execute the main function
if __name__ == "__main__":
    main()
