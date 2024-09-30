import pandas as pd
import geopandas as gpd
import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt

# Global variable for Mapbox style
MAPBOX_STYLE = 'carto-positron'

# Read the data
def load_data(filepath):
    return pd.read_csv(filepath, sep=' ')

# Create GeoDataFrame and reproject to WGS84
def create_geodataframe(df):
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs="EPSG:900913")
    return gdf.to_crs(epsg=4326)

# Add transformed coordinates to DataFrame
def add_transformed_coordinates(df, gdf):
    df['longitude'] = gdf.geometry.x
    df['latitude'] = gdf.geometry.y

# Define weight categories
def define_weight_categories(df):
    bins = [0, 1000, 5000, 10000, 50000, 100000, 200000]
    labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]
    df['weight_category'] = pd.cut(df['weight'], bins=bins, labels=labels, include_lowest=True)
    df['weight_category'] = pd.Categorical(df['weight_category'], categories=labels, ordered=True)

    # Add missing categories
    for category in set(labels) - set(df['weight_category'].dropna().unique()):
        df = df.append({'weight': np.nan, 'weight_category': category}, ignore_index=True)
    return df

# Load shapefile and reproject to WGS84
def load_shapefile(filepath):
    region_gdf = gpd.read_file(filepath, crs="EPSG:2154")
    return region_gdf.to_crs(epsg=4326)

# Create heatmap
def create_heatmap(df, center_lat, center_lon, color_mode='continuous', color_scale=None):
    if color_scale is None:
        color_scale = ['violet', 'blue', 'green', 'yellow', 'orange', 'red']

    if color_mode == 'continuous':
        fig = px.density_mapbox(df, lat='latitude', lon='longitude', z='weight',
                                radius=22, zoom=7, center=dict(lat=center_lat, lon=center_lon),
                                mapbox_style=MAPBOX_STYLE, color_continuous_scale=color_scale, 
                                opacity=0.5)
    else:  # 'interval'
        fig = px.density_mapbox(df, lat='latitude', lon='longitude', z='weight_category',
                                radius=22, zoom=7, center=dict(lat=center_lat, lon=center_lon),
                                mapbox_style=MAPBOX_STYLE, color_discrete_sequence=color_scale, 
                                opacity=0.5)

    # Update the legend title
    fig.update_layout(coloraxis_colorbar=dict(title='Weight Range (Numeric Intervals)'))
    return fig

# Create bubble map
def create_bubble_map(df, center_lat, center_lon, color_mode='continuous', color_scale=None):
    if color_scale is None:
        color_scale = ['violet', 'blue', 'green', 'yellow', 'orange', 'red']

    if color_mode == 'continuous':
        fig = px.scatter_mapbox(df, lat='latitude', lon='longitude', size='weight',
                                zoom=7, center=dict(lat=center_lat, lon=center_lon),
                                mapbox_style=MAPBOX_STYLE, color='weight', 
                                color_continuous_scale=color_scale, opacity=0.5)
    else:  # 'interval'
        fig = px.scatter_mapbox(df, lat='latitude', lon='longitude', size='weight',
                                zoom=7, center=dict(lat=center_lat, lon=center_lon),
                                mapbox_style=MAPBOX_STYLE, color='weight_category',
                                color_discrete_sequence=color_scale, opacity=0.5)

    # Update the legend title
    fig.update_layout(coloraxis_colorbar=dict(title='Weight Range (Numeric Intervals)'))
    return fig

# Add shapefile boundaries to the plot
def add_shapefile_boundaries(fig, region_gdf):
    for _, row in region_gdf.iterrows():
        geom = row.geometry
        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            fig.add_scattermapbox(lon=list(x), lat=list(y), mode='lines', line=dict(width=2, color='black'), showlegend=False)
        elif geom.geom_type == 'MultiPolygon':
            for part in geom.geoms:
                x, y = part.exterior.xy
                fig.add_scattermapbox(lon=list(x), lat=list(y), mode='lines', line=dict(width=2, color='black'), showlegend=False)

# Function to calculate total weight per region
def calculate_weight_per_region(df, region_gdf):
    # Perform a spatial join to associate points with regions
    points_gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']), crs="EPSG:4326")
    merged = gpd.sjoin(points_gdf, region_gdf, predicate='within')  # Use 'predicate'

    # Group by region and calculate sum of weights
    weight_per_region = merged.groupby('index_right')['weight'].sum().reset_index()

    # Merge with the region GeoDataFrame to preserve geometry
    result = region_gdf.merge(weight_per_region, left_index=True, right_on='index_right', how='left')

    # Fill NaN weights with 0 for regions with no points
    result['weight'] = result['weight'].fillna(0)
    
    print(result)
    
    return result

def plot_weight_by_region(center_lat, center_lon, weighted_region_gdf, color_scale=None, binning_method='quantile', num_bins=5):
    """
    Plots the regions with weights and automatically sets intervals using either quantiles or evenly spaced intervals.
    Deals with NaN values in the weights by assigning them to 'No Data' category.

    Parameters:
    - weighted_region_gdf: GeoDataFrame containing region geometries and weight data.
    - color_scale: List of colors for the categories (default is a cold-to-hot scale).
    - binning_method: 'even' for evenly spaced intervals, 'quantile' for quantile-based intervals.
    - num_bins: Number of bins/intervals for coloring the regions.

    Returns:
    - fig: Plotly choropleth mapbox figure.
    """
    # Set the color scale if not provided
    
    color_scale = ['violet', 'blue', 'green', 'yellow', 'orange', 'red']


    # Remove rows where 'weight' is NaN or 0
    weighted_region_gdf = weighted_region_gdf.dropna(subset=['weight'])  # Remove NaN values in 'weight'
    weighted_region_gdf = weighted_region_gdf[weighted_region_gdf['weight'] != 0]  # Remove rows where 'weight' is 0

    # Remove duplicate geometries
    weighted_region_gdf = weighted_region_gdf[~weighted_region_gdf.geometry.duplicated(keep='first')]


    # Now you can proceed with your existing logic
    if weighted_region_gdf.empty:
        print("Warning: The DataFrame is empty after filtering out NaN and zero weights.")
        return None  # or handle it accordingly

    # Proceed with the rest of your code, like applying pd.qcut()
    if len(weighted_region_gdf) < num_bins:
        print("Warning: Not enough data to create the requested number of bins.")
        return None  # or adjust num_bins accordingly

    # Apply qcut
    weighted_region_gdf['weight_category'] = pd.qcut(weighted_region_gdf['weight'], q=num_bins, labels=False, duplicates='drop')



    # Create bins based on the chosen binning method
    if binning_method == 'quantile':
        # Create quantile bins
        weighted_region_gdf['weight_category'] = pd.qcut(weighted_region_gdf['weight'], q=num_bins, labels=False, duplicates='drop')
    elif binning_method == 'even':
        # Create evenly spaced bins
        bins = pd.cut(weighted_region_gdf['weight'], bins=num_bins, include_lowest=True)
        weighted_region_gdf['weight_category'] = bins
    else:
        raise ValueError("Invalid binning method. Use 'quantile' or 'even'.")

    # Convert to categorical type
    weighted_region_gdf['weight_category'] = weighted_region_gdf['weight_category'].astype('category')

    # Add 'No Data' category
    weighted_region_gdf['weight_category'] = weighted_region_gdf['weight_category'].cat.add_categories('No Data')
    weighted_region_gdf['weight_category'] = weighted_region_gdf['weight_category'].fillna('No Data')

    # Convert the geometries to GeoJSON format
    geojson_data = weighted_region_gdf.__geo_interface__
    
    
    
    print(weighted_region_gdf)
    print(weighted_region_gdf.crs)
    
    invalid_geometries = weighted_region_gdf[~weighted_region_gdf.is_valid]

    print("Invalid geometries:")
    print(invalid_geometries)
    
    
    for column in weighted_region_gdf.select_dtypes(['category']).columns:
        weighted_region_gdf[column] = weighted_region_gdf[column].astype(str)
    weighted_region_gdf.to_file("output.shp")

    # Create the plot using the correct geometry and coloring
    fig = px.choropleth_mapbox(weighted_region_gdf,
                               geojson=geojson_data,  # Use converted GeoJSON data
                               locations=weighted_region_gdf.index,
                               color='weight_category',
                               color_continuous_scale=color_scale,
                               category_orders={'weight_category': sorted(weighted_region_gdf['weight_category'].unique(), key=lambda x: str(x))},
                               mapbox_style=MAPBOX_STYLE,
                               center=dict(lat=center_lat, lon=center_lon),  # Pass center here
                               zoom=7,
                               opacity=0.5)

    # Apply layout changes to use the given map style
    fig.update_layout(mapbox=dict(center=dict(lat=center_lat, lon=center_lon)))

    return fig


def plot_lines_from_geodataframe(gdf, center_lat, center_lon):
    """
    Plots only the lines from the given GeoDataFrame.
    
    Parameters:
    gdf (GeoDataFrame): A GeoDataFrame containing LineString geometries.
    center_lat (float): The latitude to center the plot.
    center_lon (float): The longitude to center the plot.
    """
    # Create a new figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Set the limits of the plot around the specified center
    ax.set_xlim(center_lon - 0.1, center_lon + 0.1)  # Adjust as needed
    ax.set_ylim(center_lat - 0.1, center_lat + 0.1)  # Adjust as needed
    
    # Plot only LineString geometries
    gdf[gdf.geometry.type == 'LineString'].plot(ax=ax, color='blue', linewidth=2)

    # Set title and labels
    ax.set_title("Line Plot from GeoDataFrame")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    
    # Show grid
    ax.grid(True)

    # Show the plot
    plt.show()


# Main function to run the entire process
def main():
    # File paths
    data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037_shuffle.txt'
    shapefile_paca = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
    
    # arrondissement.shp
    # shapefile_region = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/Arrondissement.shp'
    # epci.shp
    shapefile_region = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/EPCI.shp'
    # canton.shp
    # shapefile_region = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/canton.shp'
    # commune.shp
    # shapefile_region = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/commune.shp'
    
    
    # Load and prepare data
    df = load_data(data_path)
    gdf = create_geodataframe(df)
    add_transformed_coordinates(df, gdf)
    df = define_weight_categories(df)
    
    # Load shapefile
    region_gdf = load_shapefile(shapefile_paca)
    region_within_paca = load_shapefile(shapefile_region)
    
    # Calculate the center point for the map
    center_lat = df['latitude'].mean()
    center_lon = df['longitude'].mean()

    # Toggle plot type and color mode
    plot_type = 'bubble'  # Change to 'heatmap' or 'bubble'
    color_mode = 'interval'  # Change to 'continuous' or 'interval'

    plot_type = 'region' 


    # Create the plot based on the selected type
    if plot_type == 'heatmap':
        fig = create_heatmap(df, center_lat, center_lon, color_mode=color_mode)
    elif plot_type == 'bubble':
        fig = create_bubble_map(df, center_lat, center_lon, color_mode=color_mode)
    elif plot_type == 'region':
        weighted_region_gdf = calculate_weight_per_region(df, region_within_paca)
        plot_lines_from_geodataframe(weighted_region_gdf, center_lat, center_lon)
        fig = plot_weight_by_region(center_lat, center_lon,weighted_region_gdf, binning_method='quantile', num_bins=5)
    else:
        raise ValueError('Invalid plot type. Choose between "heatmap", "bubble", and "region".')

    # Add shapefile boundaries to the plot
    if plot_type != 'region':
        add_shapefile_boundaries(fig, region_gdf)

    # Show the plot
    fig.show()

# Execute the main function
if __name__ == "__main__":
    main()
