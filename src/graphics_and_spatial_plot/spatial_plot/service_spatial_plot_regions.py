import pandas as pd
import geopandas as gpd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Global variable for Mapbox style
MAPBOX_STYLE = 'carto-positron'

# Read the data
def load_data(filepath):
    return pd.read_csv(filepath, sep=' ')

# Create GeoDataFrame and reproject to WGS84
def create_geodataframe(df, CRS="EPSG:900913"):
    # print(df.head())
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['coord_x'], df['coord_y']), crs=CRS)
    return gdf.to_crs(epsg=4326)

# Add transformed longitude and latitude columns to the DataFrame
def add_transformed_coordinates(df, gdf):
    df['longitude'] = gdf.geometry.x
    df['latitude'] = gdf.geometry.y
    


def filter_shapefile_fields(shapefile_gdf, shapefile_path):

    if shapefile_path.endswith("Arrondissement.shp"):
        shapefile_gdf = shapefile_gdf.rename(columns={"codearondo": "region_name"})
    elif shapefile_path.endswith("EPCI.shp"):
        shapefile_gdf = shapefile_gdf.rename(columns={"NOM_EPCI": "region_name"})
    elif shapefile_path.endswith("canton.shp"):
        shapefile_gdf = shapefile_gdf.rename(columns={"numcant": "region_name"})
    elif shapefile_path.endswith("commune.shp"):
        shapefile_gdf = shapefile_gdf.rename(columns={"NOM_COMM": "region_name"})

    return shapefile_gdf

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

def plot_percentage_choropleth_mapbox(regions_gdf, center_lat, center_lon, show_zero_percentage=False):
    
    # Calculate the percentage of points inside each region
    regions_gdf['point_count'] = regions_gdf['point_count'].fillna(0)  # Avoid inplace=True

    total_points = regions_gdf['point_count'].sum()
    regions_gdf['point_percentage'] = (regions_gdf['point_count'] / total_points) * 100  # Calculate percentage
    
    # Define custom colorscale
    colorscale = ["#f7fbff", "#ebf3fb", "#deebf7", "#d2e3f3", "#c6dbef", "#b3d2e9", "#9ecae1",
                  "#85bcdb", "#6baed6", "#57a0ce", "#4292c6", "#3082be", "#2171b5", "#1361a9",
                  "#08519c", "#0b4083", "#08306b"]

    # Create the choropleth map
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,
        locations=regions_gdf.index,  # Ensure this is a unique identifier
        color='point_percentage',
        mapbox_style='carto-positron',  # Use a valid Mapbox style
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        opacity=0.6,
        color_continuous_scale=colorscale,
        labels={'point_percentage': '% of Points'},
        hover_name='region_name',  # Add the region names here
        hover_data={"point_percentage": True, "point_count": True}  # Add point count to hover data
    )

    # Add static annotations for point count
    for _, row in regions_gdf.iterrows():
        centroid = row.geometry.centroid  # Get the centroid of each region
        
        # Format the text based on show_zero_percentage
        if row['point_percentage'] == 0 and not show_zero_percentage:
            text = ""  # Do not show text for zero percentage
        else:
            text = f"{row['point_percentage']:.1f}%"  # Format the percentage with one decimal place
        
        fig.add_trace(
            go.Scattermapbox(
                lon=[centroid.x], lat=[centroid.y],  # Position of the centroid
                mode='text',
                text=text,  # Use the formatted text
                textfont=dict(size=14, color="black"),  # Customize the font size and color
                showlegend=False
            )
        )

    # Update layout
    fig.update_layout(
        title_text='Region by Percentage of Points (with Static Point Count Labels)',
        legend_title='% of Points',
        mapbox=dict(
            style="carto-positron",
            zoom=7,
            center={"lat": center_lat, "lon": center_lon},
        ),
        margin={"r":0,"t":0,"l":0,"b":0}
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



def plot_percentage_choropleth_mapbox_with_points(regions_gdf, points_gdf, center_lat, center_lon, show_zero_percentage=True):
    
    # Calculate the percentage of points inside each region
    regions_gdf['point_count'] = regions_gdf['point_count'].fillna(0)  # Avoid inplace=True
    # Calculate bounding box of the regions
    bounds = regions_gdf.total_bounds  # [minx, miny, maxx, maxy]
    min_lon, min_lat, max_lon, max_lat = bounds


    total_points = regions_gdf['point_count'].sum()
    regions_gdf['point_percentage'] = (regions_gdf['point_count'] / total_points) * 100  # Calculate percentage
    
    # Define custom colorscale
    # colorscale = ["#f7fbff", "#ebf3fb", "#deebf7", "#d2e3f3", "#c6dbef", "#b3d2e9", "#9ecae1",
    #               "#85bcdb", "#6baed6", "#57a0ce", "#4292c6", "#3082be", "#2171b5", "#1361a9",
    #               "#08519c", "#0b4083", "#08306b"]
    colorscale = px.colors.sequential.OrRd  # A perceptually distinct color scale


    # Create the choropleth map
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,
        locations=regions_gdf.index,  # Ensure this is a unique identifier
        color='point_percentage',
        # mapbox_style='carto-positron',  # Use a valid Mapbox style
        center={"lat": center_lat, "lon": center_lon},
        zoom=5,
        opacity=0.6,
        color_continuous_scale=colorscale,
        labels={'point_percentage': '% of Service Points'},
        hover_name='region_name',  # Add the region names here
        hover_data={"point_percentage": True, "point_count": True}  # Add point count to hover data
    )

        # Add points to the map
    fig.add_trace(
        go.Scattermapbox(
            lon=points_gdf.geometry.x,  # Longitude of points
            lat=points_gdf.geometry.y,  # Latitude of points
            mode='markers',
            marker=dict(size=8, color='blue', opacity=0.3),  # Customize marker size and color
            name='Points',  # Legend entry for points
            hoverinfo='text',
            showlegend=False,
            text=points_gdf['point_name'] if 'point_name' in points_gdf.columns else None  # Optional point labels
        )
    )

    # Add static annotations for point count
    for _, row in regions_gdf.iterrows():
        centroid = row.geometry.centroid  # Get the centroid of each region
        
        # Format the text based on show_zero_percentage
        if row['point_percentage'] == 0 and not show_zero_percentage:
            text = ""  # Do not show text for zero percentage
        else:
            text = f"{row['point_percentage']:.1f}%"  # Format the percentage with one decimal place
        
        fig.add_trace(
            go.Scattermapbox(
                lon=[centroid.x], lat=[centroid.y],  # Position of the centroid
                mode='text',
                text=text,  # Use the formatted text
                textfont=dict(size=15, color="black"),  # Customize the font size and color
                showlegend=False
            )
        )
    
    fig.update_layout(
        # title_text='Percentage of Service Points in Arrondissement divisions (PACA) - Model Cover Arrond',
        title_text=f'Model {model_type} - {decoupage_type} divisions (PACA)',
        title_font=dict(size=20),  # Increase title font size
        legend_title='% of Points',
        mapbox=dict(
            style="carto-positron",
            zoom=7,
            center={"lat": center_lat, "lon": center_lon},
            # bounds={"west": min_lon, "east": max_lon, "south": min_lat, "north": max_lat}
        ),
        margin={"r":0,"t":50,"l":0,"b":0},  # Adjust top margin for title
        annotations=[
            dict(
                text='Data Source: [Your Data Source Here]',  # Add data source if needed
                xref='paper', yref='paper',
                x=0.5, y=-0.1,  # Position below the map
                showarrow=False,
                font=dict(size=12, color='black')
            )
        ]
    )

    add_shapefile_boundaries(fig, region_gdf_paca, line_color='black')

    fig.show()


def plot_percentage_choropleth_mapbox_with_regions_only(regions_gdf, center_lat, center_lon, show_zero_percentage=True):
    
    # Calculate the percentage of points inside each region
    regions_gdf['point_count'] = regions_gdf['point_count'].fillna(0)  # Avoid inplace=True

    total_points = regions_gdf['point_count'].sum()
    regions_gdf['point_percentage'] = (regions_gdf['point_count'] / total_points) * 100  # Calculate percentage
    
    # Define custom colorscale
    colorscale = px.colors.sequential.OrRd  # A perceptually distinct color scale

    # Create the choropleth map for regions only
    fig = px.choropleth_mapbox(
        regions_gdf,
        geojson=regions_gdf.geometry.__geo_interface__,
        locations=regions_gdf.index,  # Ensure this is a unique identifier
        color='point_percentage',
        mapbox_style='carto-positron',  # Use a valid Mapbox style
        center={"lat": center_lat, "lon": center_lon},
        zoom=7,
        opacity=0.6,
        color_continuous_scale=colorscale,
        labels={'point_percentage': '% of Service Points'},
        hover_name='region_name',  # Add the region names here
        hover_data={"point_percentage": True, "point_count": True}  # Add point count to hover data
    )

    # Add static annotations for point count
    for _, row in regions_gdf.iterrows():
        centroid = row.geometry.centroid  # Get the centroid of each region
        
        # Format the text based on show_zero_percentage
        if row['point_percentage'] == 0 and not show_zero_percentage:
            text = ""  # Do not show text for zero percentage
        else:
            text = f"{row['point_percentage']:.1f}%"  # Format the percentage with one decimal place
        
        fig.add_trace(
            go.Scattermapbox(
                lon=[centroid.x], lat=[centroid.y],  # Position of the centroid
                mode='text',
                text=text,  # Use the formatted text
                textfont=dict(size=15, color="black"),  # Customize the font size and color
                showlegend=False
            )
        )
    
    fig.update_layout(
        title_text=f'Model {model_type} - {decoupage_type} divisions (PACA)',
        title_font=dict(size=20),  # Increase title font size
        legend_title='% of Points',
        mapbox=dict(
            style="carto-positron",
            zoom=7,
            center={"lat": center_lat, "lon": center_lon},
        ),
        margin={"r":0,"t":50,"l":0,"b":0},  # Adjust top margin for title
        annotations=[
            dict(
                text='Data Source: [Your Data Source Here]',  # Add data source if needed
                xref='paper', yref='paper',
                x=0.5, y=-0.1,  # Position below the map
                showarrow=False,
                font=dict(size=12, color='black')
            )
        ]
    )

    add_shapefile_boundaries(fig, region_gdf_paca, line_color='black')

    fig.show()

def plot_points_with_boundaries(points_gdf, region_gdf_paca, region_gdf_regions, center_lat, center_lon, 
                                point_color='blue', point_size=8, point_opacity=1.2, output_file=None):
    """
    Plot selected points on a map and add shapefile boundaries.

    Parameters:
        points_gdf (GeoDataFrame): GeoDataFrame containing points to plot.
        region_gdf_paca (GeoDataFrame): GeoDataFrame containing shapefile boundaries to overlay.
        center_lat (float): Latitude to center the map.
        center_lon (float): Longitude to center the map.
        point_color (str): Color of the points.
        point_size (int): Size of the points.
        point_opacity (float): Opacity of the points.

    Returns:
        None
    """
    # Initialize the map figure
    fig = go.Figure()

    # Add points to the map
    fig.add_trace(
        go.Scattermapbox(
            lon=points_gdf.geometry.x,  # Longitude of points
            lat=points_gdf.geometry.y,  # Latitude of points
            mode='markers',
            marker=dict(size=point_size, color=point_color, opacity=point_opacity),
            name='Points',  # Legend entry for points
            hoverinfo='text',
            text=points_gdf['point_name'] if 'point_name' in points_gdf.columns else None  # Optional point labels
        )
    )

    
    # Add shapefile boundaries for the additional 'shp_regions'
    add_shapefile_boundaries(fig, region_gdf_regions, line_color='red')  # Different line color for distinction


    # Add shapefile boundaries
    add_shapefile_boundaries(fig, region_gdf_paca, line_color='black')

    # Configure map layout
    fig.update_layout(
        mapbox=dict(
            style="carto-positron",
            zoom=7,
            center={"lat": center_lat, "lon": center_lon},
        ),
        margin={"r": 0, "t": 0, "l": 0, "b": 0}  # Remove margins
    )

    # Save the figure to a file if output_file is provided
    if output_file:
        fig.write_image(output_file, format='pdf')

    # Display the map
    fig.show()


def plot_percentage_choropleth_static(
    regions_gdf,
    points_gdf,
    show_zero_percentage=True,
    save_path="choropleth_map.pdf",
    colorbar_intervals=None  # New parameter for controlling colorbar intervals
):
    # Calculate the percentage of points inside each region
    regions_gdf['point_count'] = regions_gdf['point_count'].fillna(0)  # Fix chained assignment warning
    total_points = regions_gdf['point_count'].sum()
    regions_gdf['point_percentage'] = (regions_gdf['point_count'] / total_points) * 100

    # Set color scale with continuous colormap
    colorscale = ListedColormap(plt.cm.OrRd(np.linspace(0, 1, 256)))
    # norm = plt.Normalize(vmin=regions_gdf['point_percentage'].min(), vmax=regions_gdf['point_percentage'].max())
    norm = plt.Normalize(vmin=0, vmax=8)  # Use the full range for 0-100%



    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    
    # Plot regions with choropleth based on point percentage
    regions_gdf.plot(column='point_percentage', cmap=colorscale, linewidth=0.8, ax=ax, edgecolor='black', norm=norm)
    
    # Create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=colorscale, norm=norm)
    sm._A = []  # Empty array for ScalarMappable

    # Add colorbar with custom tick intervals if provided
    cbar = fig.colorbar(sm, ax=ax)
    if colorbar_intervals is not None:
        cbar.set_ticks(colorbar_intervals)
    cbar.set_label('% of Service Points')

    # Plot points on the map
    points_gdf.plot(ax=ax, color='blue', markersize=8, alpha=0.3, label='Points')

    add_percent_labels = False

    # Add percentage labels at region centroids
    for _, row in regions_gdf.iterrows():
        centroid = row.geometry.centroid
        percentage = row['point_percentage']
        if add_percent_labels:
            if percentage > 0 or show_zero_percentage:
                ax.annotate(f"{percentage:.1f}%", xy=(centroid.x, centroid.y), ha="center", fontsize=8, color="black")

    # Add title and legend
    ax.set_title("Percentage of Service Points by Region")
    ax.set_axis_off()  # Optional: Remove axis for a cleaner look
    
    # Save as PDF
    plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close()



# Main function to run the entire process
def main():
    # File paths
    # data_path = '/home/felipe/Documents/Projects/GeoAvigon/pmp_code/PMPSolver/data/PACA_jul24/cust_weights_PACA_2037.txt'
    data_path = '/home/falbuquerque/Documents/projects/GeoAvignon/jupyter_analysis/data/PACA_5km/cust_weights_5km.txt'
    # Load and prepare data
    df = load_data(data_path)
    gdf = create_geodataframe(df)
    add_transformed_coordinates(df, gdf)
    # Calculate the center point for the map
    center_lat = df['latitude'].mean()
    center_lon = df['longitude'].mean()
    
    global model_type 
    # model_type = 'Real BPE23'
    model_type = 'Cover Canton'  # 'Cover Arrond', 'Cover EPCI', 'Cover Canton', 'Cover Commune'
    # model_type = 'Cover Commune'  # 'Cover Arrond', 'Cover EPCI', 'Cover Canton', 'Cover Commune'
    global decoupage_type
    decoupage_type = 'canton' # 'Arrondissement', 'EPCI', 'canton', 'commune'


    # if model_type == 'Real BPE23':
    #     data_service_path = 'data/solutions_service_cinema/points_BPE23_F303_paca_table.csv'
    #     df_service = load_service_data(data_service_path)
    #     # print(df_service.head())    
    #     gdf_service = create_geodataframe(df_service,"EPSG:2154")  # Create a GeoDataFrame from service data

    # if model_type == 'Cover Arrond' or model_type == 'No Cover':
    #     data_service_path = 'data/solutions_service_cinema/test_paca_cinema_p_192_EXACT_CPMP_table.txt'
    # elif model_type == 'Cover EPCI':
    #     data_service_path = 'data/solutions_service_cinema/test_paca_cinema_EPCI_p_192_EXACT_CPMP_cover_EPCI_table.txt'
    # elif model_type == 'Cover Canton':
    #     data_service_path = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton_table.txt'
    # elif model_type == 'Cover Commune':
    #     data_service_path = 'data/solutions_service_cinema/test_paca_cinema_commune_p_192_EXACT_CPMP_cover_commune_table.txt'    
    
    # data_service_path = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton_table.txt'
    # data_service_path = 'data/solutions_service_cinema/output_paca_grid_cinema_EXACT_CPMP_canton_p_192_EXACT_CPMP_cover_canton_table.txt'



    # data_service_path = 'data/solutions_service_cinema/test_paca_cinema_canton_p_192_EXACT_CPMP_cover_canton_table.txt'
    # shapefile_regions = f'/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/{decoupage_type}.shp'


    # data_service_path = 'data/solutions_service_cinema/output_paca_grid_cinema_EXACT_CPMP_canton_p_192_EXACT_CPMP_cover_canton_table.txt'
    # shapefile_regions = '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_grid/grid_canton.shp'


    data_service_path = 'data/solutions_service_cinema/output_paca_voronoi_cinema_EXACT_CPMP_canton_p_192_EXACT_CPMP_cover_canton_table.txt'
    shapefile_regions = f'/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_voronoi/voronoi_canton_5km.shp'



    df_service = load_service_data(data_service_path)
    gdf_service = create_geodataframe(df_service)

    if model_type == 'Real BPE23':
        gdf_service = create_geodataframe(df_service,"EPSG:2154")  # Create a GeoDataFrame from service data
    
    # shapefile_paca = '/home/felipe/Documents/Projects/GeoAvigon/create_instance_PACA/Create_data_PACA/Create_data_PACA/Creation_Real_Instance/Decoupages_GIS/PACA_region_polygon.shp'
    # /home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_grid/grid_canton.shp
    shapefile_paca = '/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/PACA_region.shp'


    # shapefile_regions = f'/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/{decoupage_type}.shp'

    # /home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_voronoi/voronoi_canton_5km.shp
     # /home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_grid/grid_canton.shp
    
    # shapefile_regions = f'/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_voronoi/voronoi_canton_5km.shp'

    # Load shapefiles
    global region_gdf_paca
    region_gdf_paca = load_shapefile(shapefile_paca)
    region_gdf_regions = load_shapefile(shapefile_regions)   
    region_gdf_regions = filter_shapefile_fields(region_gdf_regions, shapefile_regions)



    plot_points = False
    # Plot points if `plot_points` is True
    if plot_points:
        fig = px.scatter_mapbox(
            df,
            lat='latitude',
            lon='longitude',
            hover_name='customer',  # Replace 'customer' with the actual column name for hover info
            zoom=6,  # Adjust the zoom level as needed
            center={"lat": center_lat, "lon": center_lon},
            mapbox_style=MAPBOX_STYLE
        )

        # Add shapefile boundaries for PACA
        add_shapefile_boundaries(fig, region_gdf_paca, line_color='black')

        # Add shapefile boundaries for the additional 'shp_regions'
        add_shapefile_boundaries(fig, region_gdf_regions, line_color='red')  # Different line color for distinction

        # Show the plot
        fig.show()


    # plot only the selected points and the region boundaries
    plot_points_with_boundaries(gdf_service, region_gdf_paca,region_gdf_regions, center_lat, center_lon, point_color='blue', point_size=8, point_opacity=0.3,output_file="voronoi_canton_p_192.pdf")

    # # Calculate the sum of weights in each region and update the region GeoDataFrame
    # region_gdf_regions = sum_weights_in_regions(df, gdf, region_gdf_regions)
    # # Plot the choropleth map with the summed weights
    # plot_choropleth_mapbox(region_gdf_regions, center_lat, center_lon)
    
    
    # Calculate the number of points in each region and update the region GeoDataFrame
    # region_gdf_regions = count_points_in_regions(df_service, gdf_service, region_gdf_regions)
    # print(region_gdf_regions)
    # Plot the choropleth map with the point counts
    # plot_point_count_mapbox(region_gdf_regions, center_lat, center_lon)

    # plot_choropleth_mapbox_with_points(region_gdf_regions, gdf_service, center_lat, center_lon)
    # plot_point_percentage_mapbox(region_gdf_regions, gdf_service, center_lat, center_lon, show_percent_labels=True)

    # plot_percentage_choropleth_mapbox_with_points(region_gdf_regions, gdf_service, center_lat, center_lon, show_zero_percentage=False)
    # plot_percentage_choropleth_mapbox_with_regions_only(region_gdf_regions, center_lat, center_lon, show_zero_percentage=False)\

    # plot_percentage_choropleth_static(region_gdf_regions, gdf_service, show_zero_percentage=False, save_path="choropleth_map.pdf")

    # colorbar_intervals = [0, 1, 2, 3, 4, 5, 6, 7, 8]  
    # plot_percentage_choropleth_static(
    #     region_gdf_regions,
    #     gdf_service,
    #     show_zero_percentage=False,
    #     save_path="choropleth_map.pdf",
    #     colorbar_intervals=colorbar_intervals  # Set intervals for colorbar display
    # )


# Execute the main function
if __name__ == "__main__":
    main()
