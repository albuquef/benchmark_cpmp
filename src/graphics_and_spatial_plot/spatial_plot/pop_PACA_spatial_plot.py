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

# Add transformed coordinates to DataFrame
def add_transformed_coordinates(df, gdf):
    df['longitude'] = gdf.geometry.x
    df['latitude'] = gdf.geometry.y

# Define weight categories
# def define_weight_categories(df):
#     bins = [0, 1000, 5000, 10000, 50000, 100000, 200000]
#     labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]
#     df['weight_category'] = pd.cut(df['weight'], bins=bins, labels=labels, include_lowest=True)
#     df['weight_category'] = pd.Categorical(df['weight_category'], categories=labels, ordered=True)

#     # Add missing categories
#     for category in set(labels) - set(df['weight_category'].dropna().unique()):
#         df = df.append({'weight': np.nan, 'weight_category': category}, ignore_index=True)
#     return df

def define_weight_categories(df):
    
    bins = [0, 1000, 5000, 10000, 50000, 100000, 200000]
    labels = [f'{bins[i]}-{bins[i + 1]}' for i in range(len(bins) - 1)]
    df['weight_category'] = pd.cut(df['weight'], bins=bins, labels=labels, include_lowest=True)
    df['weight_category'] = pd.Categorical(df['weight_category'], categories=labels, ordered=True)
    
    # Add missing categories
    for category in set(labels) - set(df['weight_category'].dropna().unique()):
        new_row = pd.DataFrame({'weight': [np.nan], 'weight_category': [category]})
        df = pd.concat([df, new_row], ignore_index=True)
        
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
    fig.update_layout(coloraxis_colorbar=dict(title='Number of inhabitants'))
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
    fig.update_layout(coloraxis_colorbar=dict(title='Number of inhabitants'))
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

# Main function to run the entire process
def main():
    # File paths
    # data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037_shuffle.txt'
    # data_path = '/home/falbuquerque/Documents/projects/GeoAvignon/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    data_path = '/home/falbuquerque/Documents/projects/GeoAvignon/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    # data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    # shapefile_paca = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
    shapefile_paca = '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/PACA_region.shp'

    # Load and prepare data
    df = load_data(data_path)
    gdf = create_geodataframe(df)
    add_transformed_coordinates(df, gdf)
    df = define_weight_categories(df)
    
    # Load shapefile
    region_gdf = load_shapefile(shapefile_paca)

    # Calculate the center point for the map
    center_lat = df['latitude'].mean()
    center_lon = df['longitude'].mean()

    # Toggle plot type and color mode
    plot_type = 'bubble'  # Change to 'heatmap' or 'bubble'
    color_mode = 'continuous'  # Change to 'continuous' or 'interval'

    # Create the plot based on the selected type
    if plot_type == 'heatmap':
        fig = create_heatmap(df, center_lat, center_lon, color_mode=color_mode)
    elif plot_type == 'bubble':
        fig = create_bubble_map(df, center_lat, center_lon, color_mode=color_mode)
    else:
        raise ValueError('Invalid plot type. Choose between "heatmap", "bubble", and "region".')

    # Add shapefile boundaries to the plot
    add_shapefile_boundaries(fig, region_gdf)

    # Show the plot
    fig.show()

# Execute the main function
if __name__ == "__main__":
    main()
