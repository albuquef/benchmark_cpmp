import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import re
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
import shapely.geometry as geom
import random

machine = 'lia'
# machine = 'mine'

#path 1km INSEE population data
path_INSEE_1km='/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/pop_1km_2019/points_population2019_paca_1km_wihtout_islands_table.csv'
if machine == 'mine': path_INSEE_1km='/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/pop_1km_2019/points_population2019_paca_1km_wihtout_islands_table.csv'
# path points grids
GRID_TYPE='2km'
path_points_grid = f'/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/grid_{GRID_TYPE}/points_grid_{GRID_TYPE}_paca_table.csv'
if machine == 'mine': path_points_grid = f'/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Population/grid_{GRID_TYPE}/points_grid_{GRID_TYPE}_paca_table.csv'

DECOUPAGE_FILE = 'commune'
path_to_regions_shapefile =  f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/{DECOUPAGE_FILE}.shp"
if machine == 'mine': path_to_regions_shapefile =  f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/{DECOUPAGE_FILE}.shp"

#define random seed 42
random.seed(42)


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

def load_grid_points_data(filepath):
    """Loads grid data and converts CRS."""
    df = pd.read_csv(filepath, low_memory=False)
    gdf = gpd.GeoDataFrame() # coordinates in WGS 84
    if(GRID_TYPE=='5km'): gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']), crs='ESRI:102100') # coordinates in WGS 84
    if(GRID_TYPE=='2km'): gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x_LAMB93'], df['y_LAMB93']), crs='EPSG:2154') # coordinates in lambert 93
    
    #if GRID_TYPE=='5km add x_LAMB93 and y_LAMB93
    if(GRID_TYPE=='5km'):
        gdf_lamb93 = gdf.to_crs('EPSG:2154')
        gdf['x_LAMB93'] = gdf_lamb93['geometry'].x
        gdf['y_LAMB93'] = gdf_lamb93['geometry'].y
    
    gdf = gdf.to_crs(epsg=3035) # Convert to EPSG:3035
    return gdf

def calculate_weights(gdf_pop, gdf_grid):
    """Performs spatial join to find nearest neighbors and calculate weights."""
    print('-' * 50)
    print("[INFO] Calculating weights...")

    gdf_grid['weight'] = 0.0  # Ensure weight column is of type float
    gdf_pop_coords = list(zip(gdf_pop.geometry.x, gdf_pop.geometry.y))
    gdf_grid_coords = list(zip(gdf_grid.geometry.x, gdf_grid.geometry.y))
    
    nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(gdf_grid_coords)
    distances, indices = nn.kneighbors(gdf_pop_coords)
    gdf_pop['closest_idx'] = indices.flatten()
    
    for idx in range(len(gdf_grid)):
        gdf_grid.at[idx, 'weight'] = gdf_pop[gdf_pop['closest_idx'] == idx]['ind'].sum()
    
    # compare the sum of ind column with the sum of weight column
    sum_ind = gdf_pop['ind'].sum()
    sum_weight = gdf_grid['weight'].sum()
    
    if sum_ind == sum_weight:
        print(f"The sums of the weights are equal: {sum_ind}.")
    else:
        print("[WARN] The sums of the weights are not equal.")
        print(f"Sum of 'ind' pop: {sum_ind}")
        print(f"Sum of 'weight' grid: {sum_weight}")
        
    print('-' * 50)
    
    

    # Step 3: Set up side-by-side histograms
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)  # 1 row, 2 columns

    # Plot the histogram for gdf_pop weights on the first axis
    axes[0].hist(gdf_pop['ind'], bins=30, alpha=0.7, color='blue', edgecolor='black')
    axes[0].set_title('grid 1km')
    axes[0].set_xlabel('Weight')
    axes[0].set_ylabel('Frequency (log)')
    axes[0].set_yscale('log')  # Use logarithmic scale for y-axis

    # Plot the histogram for gdf_grid weights on the second axis
    axes[1].hist(gdf_grid['weight'], bins=30, alpha=0.7, color='green', edgecolor='black')
    axes[1].set_title(f'grid {GRID_TYPE}')
    axes[1].set_xlabel('Weight')
    axes[1].set_yscale('log')  # Use logarithmic scale for y-axis

    # Adding a main title for both histograms
    fig.suptitle(f'Side-by-Side Comparison of Weight Frequency: grid 1km vs grid {GRID_TYPE}', fontsize=16)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f"plot_histogram_grid_1km_vs_grid_{GRID_TYPE}.png")
    plt.show()
    
    return gdf_grid

def filter_weights(gdf_grid, k):
    """Eliminates points with zero weight and transfers weights smaller than k to the nearest non-zero point."""
    print('-' * 50)
    print(f"[INFO] Filtering points with zero weight and transferring weights smaller than {k}...")
    print(f"Number of elements before remove the zeros weights: {len(gdf_grid)}")
    
    # Initialize color column
    gdf_grid['color'] = 'green'
    # Separate points with zero weight for later plotting
    zero_weight_points = gdf_grid[gdf_grid['weight'] == 0].copy()
    zero_weight_points['color'] = 'red'
    print(f'Number of elements with zero weight: {len(gdf_grid[gdf_grid["weight"] == 0])}')
    
    # Filter out points with zero weight
    gdf_grid = gdf_grid[gdf_grid['weight'] != 0].copy()
    
    
    # Identify points with weights smaller than k
    small_weight_mask = gdf_grid['weight'] < k
    large_weight_mask = gdf_grid['weight'] >= k
    
    small_weight_points = gdf_grid[small_weight_mask]
    large_weight_points = gdf_grid[large_weight_mask]
    
    if not small_weight_points.empty and not large_weight_points.empty:
        # Get coordinates for small and large weight points
        small_weight_coords = list(zip(small_weight_points.geometry.x, small_weight_points.geometry.y))
        large_weight_coords = list(zip(large_weight_points.geometry.x, large_weight_points.geometry.y))
        
        # Find the nearest large weight point for each small weight point
        nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(large_weight_coords)
        distances, indices = nn.kneighbors(small_weight_coords)
        
        # Transfer small weights to the nearest large weight point
        for i, idx in enumerate(small_weight_points.index):
            nearest_large_idx = large_weight_points.index[indices[i][0]]
            gdf_grid.at[nearest_large_idx, 'weight'] += gdf_grid.at[idx, 'weight']
            gdf_grid.at[idx, 'weight'] = 0  # Set small weight to zero after transfer
        
            # Identify and separate points that now have zero weight after transfer
        transferred_zero_weight_points = gdf_grid[gdf_grid['weight'] == 0].copy()
        transferred_zero_weight_points['color'] = 'red'
        
        # Combine both sets of zero-weight points
        zero_weight_points = zero_weight_points.append(transferred_zero_weight_points)
        
        
        # Remove points that now have zero weight after transfer
        gdf_grid = gdf_grid[gdf_grid['weight'] != 0]
    
    
    # Print how many ellement before and now
    print(f"Number of elements before: {len(gdf_grid) + len(small_weight_points)}")
    print(f"Number of elements after filtering less than {k}: {len(gdf_grid)}")
    if k != 0:
        print(f'Number of elements removed in total: {len(small_weight_points) + len(transferred_zero_weight_points)}')   
    print('-' * 50) 
    
    create_plot = True
    # Plotting
    if create_plot:
        fig, ax = plt.subplots()
        gdf_grid[gdf_grid['color'] == 'green'].plot(ax=ax, color='green', label=f'Remaining (size: {len(gdf_grid)})')
        if not zero_weight_points.empty:
            zero_weight_points.plot(ax=ax, color='red', label=f'Removed (size: {len(zero_weight_points)})')
        plt.legend()
        plt.title(f"Grid points with the filter to reallocate weights smaller than {k}") 
        plt.savefig(f'plots/creation_instance/customers/points_grid_{GRID_TYPE}_paca_table_filter.png')
        plt.show()
    
    return gdf_grid

def compare_grid_with_old_dataset(gdf_grid, gdf_old_dataset):
    
    print('-' * 50)
    print("[INFO] Comparing the new dataset with the old one...")
    print(f"Number of elements in the new dataset: {len(gdf_grid)}")
    print(f"Number of elements in the old dataset: {len(gdf_old_dataset)}")
    # compare dist max min avg and std of the two datasets
    print('New dataset:')
    print(f"Sum of 'weight' column in the new dataset: {gdf_grid['weight'].sum()}")
    print(f"Max value in the 'weight' column of the new dataset: {gdf_grid['weight'].max()}")
    print(f"Min value in the 'weight' column of the new dataset: {gdf_grid['weight'].min()}")
    print(f"Mean value in the 'weight' column of the new dataset: {gdf_grid['weight'].mean()}")
    print(f"Standard deviation of the 'weight' column of the new dataset: {gdf_grid['weight'].std()}")
    print('Old dataset:')
    print(f"Sum of 'weight' column in the old dataset: {gdf_old_dataset['weight'].sum()}")
    print(f"Max value in the 'weight' column of the old dataset: {gdf_old_dataset['weight'].max()}")
    print(f"Min value in the 'weight' column of the old dataset: {gdf_old_dataset['weight'].min()}")
    print(f"Mean value in the 'weight' column of the old dataset: {gdf_old_dataset['weight'].mean()}")
    print(f"Standard deviation of the 'weight' column of the old dataset: {gdf_old_dataset['weight'].std()}")
    
    create_plot = True  
    if create_plot:
        # compare the histogram and heatmap of the two datasets
        fig, ax = plt.subplots(1, 1, figsize=(20, 10))
        gdf_grid['weight'].plot(kind='hist', bins=5, alpha=0.5, label=f'new dataset ({len(gdf_grid)})')
        gdf_old_dataset['weight'].plot(kind='hist', bins=5, alpha=0.5, label=f'old dataset ({len(gdf_old_dataset)})')
        plt.title("Distribution of Weights")
        plt.xlabel("Weight")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()
        
        fig, ax = plt.subplots(1, 2, figsize=(20, 10))
        gdf_grid.plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
        gdf_old_dataset.plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
        ax[0].set_title(f"Heat Map of Grid Population (new dataset)\nnum elements= {len(gdf_grid)}, sum = {gdf_grid['weight'].sum()}")
        ax[0].set_xlabel("Easting")
        ax[0].set_ylabel("Northing")
        ax[1].set_title(f"Heat Map of Grid Population (old dataset)\nnum elements={len(gdf_old_dataset)}  , sum wi= {gdf_old_dataset['weight'].sum()}")
        ax[1].set_xlabel("Easting")
        ax[1].set_ylabel("Northing")
        plt.savefig('plots/creation_instance/customers/heatmaps_comparison.png')
        plt.show()


    print('-' * 50)


def plot_pop_and_grid_points(gdf_pop, gdf_grid):
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    gdf_pop.plot(ax=ax, marker='x', color='blue', markersize=5, label=f'Population (size: {gdf_pop.shape[0]})')
    gdf_grid.plot(ax=ax, marker='o', color='green', markersize=10, label=f'Grid points (size: {gdf_grid.shape[0]})')
    plt.title("Plot of Points from Both Datasets")
    plt.xlabel("Easting")
    plt.ylabel("Northing")
    plt.legend()
    plt.savefig(f'plots/creation_instance/customers/points_pop_grid_{GRID_TYPE}_paca_table.png')
    plt.show()
    
    # """Plots datasets side by side."""
    # fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    # gdf_grid.plot(ax=ax[0], marker='x', color='blue', markersize=5, label='Grid')
    # gdf_pop.plot(ax=ax[1], marker='o', color='red', markersize=5, label='Population')
    # ax[0].set_title("Grid Points")
    # ax[0].set_xlabel("Easting")
    # ax[0].set_ylabel("Northing")
    # ax[1].set_title("Population Points")
    # ax[1].set_xlabel("Easting")
    # ax[1].set_ylabel("Northing")
    # plt.show()

def plot_heatmap_grid_points(gdf_grid):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    gdf_grid.plot(ax=ax, column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
    plt.title("Heat Map of Grid Points with Weights")
    plt.xlabel("Easting")
    plt.ylabel("Northing")
    plt.savefig(f'plots/creation_instance/customers/points_grid_{GRID_TYPE}_paca_table_heatmap.png')
    plt.show()


def ensure_coverage_for_regions(shapefile_regions, gdf_points_grid, gdf_points_all_grid):
    """
    Ensure all regions in the shapefile are covered by points in gdf_points_grid. 
    If a region is not covered, either add a random point from gdf_points_all_grid inside the uncovered region
    or generate a random point inside the region.
    
    Additionally, plot the existing points, points added from the 2km grid, and random points added, 
    with different colors for each type of point.
    
    Parameters:
    - shapefile_regions (GeoDataFrame): GeoDataFrame of regions.
    - gdf_points_grid (GeoDataFrame): GeoDataFrame of currently selected points.
    - gdf_points_all_grid (GeoDataFrame): GeoDataFrame of all available points.

    Returns:
    - gdf_points_grid (GeoDataFrame): Updated GeoDataFrame with all regions covered.
    """
    print('Ensure coverage for all regions...') 
    print(f'Region shapefile: {shapefile_regions.shape[0]}')

    # Spatial join to find which regions are covered by points in gdf_points_grid
    joined = gpd.sjoin(shapefile_regions, gdf_points_grid, how='left', predicate='intersects')
    
    # Identify uncovered regions
    uncovered_regions = joined[joined['index_right'].isna()]
    
    if uncovered_regions.empty:
        print("All regions are already covered.")
        return gdf_points_grid  # No update needed

    print(f"Found {len(uncovered_regions)} uncovered regions.")
    
    points_added_from_2km = []  # Track points added from the 2km grid
    random_points_added = []    # Track random points added

    for _, region in uncovered_regions.iterrows():
        # Select points from gdf_points_all_grid inside the uncovered region
        points_in_region = gdf_points_all_grid[gdf_points_all_grid.intersects(region.geometry)]
        
        if not points_in_region.empty:
            # Randomly select one point from the 2km grid
            selected_point = points_in_region.sample(n=1, random_state=42)
            # weight is zero
            selected_point['weight'] = 0
            points_added_from_2km.append(selected_point)
            print(f"Added point {selected_point.iloc[0]['id']} from all grid to cover region {region['subarea_id']}.")
        else:
            print(f"No points available inside region {region['subarea_id']}. Adding a random point.")
            # Generate a random point inside the region
            random_point = generate_random_point_within_polygon(region.geometry)
            random_point_gdf = gpd.GeoDataFrame(
                [{'geometry': random_point, 'id': f'random_{region["subarea_id"]}'}],
                crs=gdf_points_grid.crs
            )
            # weight is zero
            random_point_gdf['weight'] = 0
            # fid ==  id 
            random_point_gdf['fid'] = random_point_gdf['id']
            # Reproject to LAMB93 (EPSG:2154) to extract x_LAMB93 and y_LAMB93
            random_point_gdf_lamb93 = random_point_gdf.to_crs('EPSG:2154')
            random_point_gdf['x_LAMB93'] = random_point_gdf_lamb93.geometry.x
            random_point_gdf['y_LAMB93'] = random_point_gdf_lamb93.geometry.y

            random_points_added.append(random_point_gdf)

    # Combine all points added from 2km grid and random points
    if points_added_from_2km:
        new_points_from_2km_gdf = gpd.GeoDataFrame(pd.concat(points_added_from_2km, ignore_index=True))
        gdf_points_grid = pd.concat([gdf_points_grid, new_points_from_2km_gdf], ignore_index=True)
    
    if random_points_added:
        new_random_points_gdf = gpd.GeoDataFrame(pd.concat(random_points_added, ignore_index=True))
        gdf_points_grid = pd.concat([gdf_points_grid, new_random_points_gdf], ignore_index=True)

    print(f"Total number of points after adding: {len(gdf_points_grid)}")

    # Plot the results
    plot_points_with_regions(
        shapefile_regions,
        gdf_points_grid,
        points_added_from_2km,
        random_points_added
    )
    
    return gdf_points_grid


def generate_random_point_within_polygon(polygon):
    """
    Generate a random point within a given polygon.
    """
    minx, miny, maxx, maxy = polygon.bounds
    while True:
        random_point = geom.Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if polygon.contains(random_point):
            return random_point


def plot_points_with_regions(shapefile_regions, gdf_points_grid, points_added_from_2km, random_points_added):
    """
    Plot the regions with the existing points, points added from 2km grid, and random points added.
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    shapefile_regions.plot(ax=ax, color='lightgrey', edgecolor='black', alpha=0.5, label='Regions')
    gdf_points_grid.plot(ax=ax, color='blue', markersize=5, label='Existing Points')
    
    if points_added_from_2km:
        new_points_from_2km = gpd.GeoDataFrame(pd.concat(points_added_from_2km, ignore_index=True))
        new_points_from_2km.plot(ax=ax, color='green', markersize=20, label='Added from 2km Grid')
    
    if random_points_added:
        new_random_points = gpd.GeoDataFrame(pd.concat(random_points_added, ignore_index=True))
        new_random_points.plot(ax=ax, color='red', markersize=20, label='Random Points')
    
    num_added_points = len(points_added_from_2km) + len(random_points_added)

    plt.legend()
    plt.title(f"Region Coverage with({num_added_points} points added)")
    plt.show()

def create_final_table_instance(gdf_grid):
    # create one file txt with 4 columns:id, weight, coord_x, coord_y, fid
    filename = f'./outputs/PACA_{GRID_TYPE}/cust_weights_{GRID_TYPE}.txt'

    print('-' * 50)
    print(f"[INFO] Creating final table instance in {filename}...")
    
    if (GRID_TYPE=='2km'):
        #      customer weight coord_x coord_y id_grid2km fid x_LAMB93 y_LAMB93
        # line 5085 6552.0 3955570.085641593 2252179.211702804 7489 14508 903606.8 6245458
        # put the weight of the element fid == 14508 put the weight associate to the closest point and remove the point 5085
        # get the point with the x_LAMB93:1027606.8   y_LAMB93: 6277458
        print(f'Reallocation of the point with the x_LAMB93:1027606.8   y_LAMB93: 6277458')
        # just put the point in the last position in table gdg_grid
        gdf_grid = pd.concat([gdf_grid, gdf_grid[(gdf_grid['x_LAMB93'] == 1027606.8) & (gdf_grid['y_LAMB93'] == 6277458)]])
        # Remove the duplicated point
        gdf_grid = gdf_grid.drop_duplicates(subset=['x_LAMB93', 'y_LAMB93'], keep='last')
        # Check the sum of the weight
        print(f"total weight (after): {gdf_grid['weight'].sum()}")


        # point_to_remove = gdf_grid[(gdf_grid['x_LAMB93'] == 1027606.8) & (gdf_grid['y_LAMB93'] == 6277458)]
        # # get closest point to the point to remove
        # point_to_remove_coords = list(zip(point_to_remove.geometry.x, point_to_remove.geometry.y))
        # gdf_grid_coords = list(zip(gdf_grid.geometry.x, gdf_grid.geometry.y))
        # nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(gdf_grid_coords)
        # distances, indices = nn.kneighbors(point_to_remove_coords)
        # idx = indices[0][0]
        # print(f"total weight (before): {gdf_grid['weight'].sum()}")
        # print(f"weight of the point to remove: {point_to_remove['weight'].values[0]}")
        # print(f"weight of the closest point (before): {gdf_grid.at[idx, 'weight']}")
        # gdf_grid.at[idx, 'weight'] += point_to_remove['weight'].values[0]
        # print(f"weight of the closest point (after): {gdf_grid.at[idx, 'weight']}")
        # # remove the point to remove the 
        # gdf_grid = gdf_grid.drop(point_to_remove.index)
        # print(f"total weight (after): {gdf_grid['weight'].sum()}")

    cont = 1
    # print columns names  of df_grid
    with open(filename, 'w') as f:
        # first line with the columns names
        f.write(f"customer weight coord_x coord_y id_grid{GRID_TYPE} fid x_LAMB93 y_LAMB93\n")
        for idx, row in gdf_grid.iterrows():
            f.write(f"{cont} {row['weight']} {row['geometry'].x} {row['geometry'].y} {idx} {row['fid']} {row['x_LAMB93']} {row['y_LAMB93']}\n")
            cont += 1

    
    print(f"Final table instance created in {filename}.")
    print(f'Total number of elements: {cont - 1}')

    print('-' * 50)
    
    
df_points_pop_paca = pd.read_csv(path_INSEE_1km)
gdf_points_pop_paca = load_population_data(path_INSEE_1km)
gdf_points_all_grid = load_grid_points_data(path_points_grid)

gdf_points_all_grid = calculate_weights(gdf_points_pop_paca, gdf_points_all_grid)
plot_pop_and_grid_points(gdf_points_pop_paca, gdf_points_all_grid)

gdf_points_grid = filter_weights(gdf_points_all_grid, 0) # Filter weights smaller than k and remove points with zero weight
# plot_heatmap_grid_points(gdf_points_grid)

# print sum of weights
sum_weight_gdf_points_grid = gdf_points_grid['weight'].sum()
print(f"Sum of 'weight' column in gdf_points_grid: {sum_weight_gdf_points_grid}")


shapefile_regions = gpd.read_file(path_to_regions_shapefile)  # Load regions shapefile
# transform shp to the same crs of gdf_points_grid
shapefile_regions = shapefile_regions.to_crs(gdf_points_grid.crs)

gdf_points_grid = ensure_coverage_for_regions(shapefile_regions, gdf_points_grid, gdf_points_all_grid)

# number points in the grid
print(f"Number of elements in the grid: {len(gdf_points_grid)}")
# sum weight of the grid
print(f"Sum of 'weight' column in the grid: {gdf_points_grid['weight'].sum()}")
# check if any nan value in th gdf_points_grid
print(f"Number of nan values in the grid: {gdf_points_grid.isnull().sum().sum()}")


# print tail
print(gdf_points_grid.tail())

# create_final_table_instance(gdf_points_grid)


exit()
#  -----    compare with old dataset ------
df_cust_weights = pd.read_csv('data/PACA/cust_weights.txt', delim_whitespace=True, names=['customer', 'weight'])
df_map_id_cust_loc = pd.read_csv('data/PACA/map_id_cust_loc.txt', delim_whitespace=True, names=['id', 'identif'])
df_locations = pd.read_csv('data/PACA/locations_paca_2017_coord.csv')

# Ensure that 'identif' columns are of the same type (string)
df_map_id_cust_loc['identif'] = df_map_id_cust_loc['identif'].astype(str)
df_locations['identif'] = df_locations['identif'].astype(str)

# Merge the DataFrames to associate weights with coordinates
df_merged = pd.merge(df_map_id_cust_loc, df_cust_weights, left_on='id', right_on='customer')
df_old_dataset = pd.merge(df_merged, df_locations, left_on='identif', right_on='identif')
# convert gdf_final['weight'] to float
df_old_dataset['weight'] = df_old_dataset['weight'].astype(float)

# Create a GeoDataFrame with the coordinates and weights in WGS 84 (EPSG:4326)
df_old_dataset = gpd.GeoDataFrame(df_old_dataset, geometry=gpd.points_from_xy(df_old_dataset.x, df_old_dataset.y), crs='EPSG:3035')
compare_grid_with_old_dataset(gdf_points_grid, df_old_dataset)


exit()

# Read the tables into DataFrames
df_cust_weights = pd.read_csv('data/PACA/cust_weights.txt', delim_whitespace=True, names=['customer', 'weight'])
df_map_id_cust_loc = pd.read_csv('data/PACA/map_id_cust_loc.txt', delim_whitespace=True, names=['id', 'identif'])
df_locations = pd.read_csv('data/PACA/locations_paca_2017_coord.csv')

# Ensure that 'identif' columns are of the same type (string)
df_map_id_cust_loc['identif'] = df_map_id_cust_loc['identif'].astype(str)
df_locations['identif'] = df_locations['identif'].astype(str)

# Merge the DataFrames to associate weights with coordinates
df_merged = pd.merge(df_map_id_cust_loc, df_cust_weights, left_on='id', right_on='customer')
df_final = pd.merge(df_merged, df_locations, left_on='identif', right_on='identif')

# Create a GeoDataFrame with the coordinates and weights in WGS 84 (EPSG:4326)
gdf_final = gpd.GeoDataFrame(df_final, geometry=gpd.points_from_xy(df_final.x, df_final.y), crs='EPSG:3035')

# Plot the points
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_grid.plot(ax=ax, marker='x', color=gdf_points_grid['color'], markersize=5, label='Dataset 2')
gdf_final.plot(ax=ax, marker='o', color='red', markersize=5, label='old dataset')
plt.title("Plot of Points with Weights (ESRI:102100)")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

# convert gdf_final['weight'] to float
gdf_final['weight'] = gdf_final['weight'].astype(float)

# compare the weights of the two datasets
sum_weight_gdf_points_grid = gdf_points_grid['weight'].sum()
sum_weight_gdf_final = gdf_final['weight'].sum()
print(f"Sum of 'weight' column in gdf_points_grid: {sum_weight_gdf_points_grid}")
print(f"Sum of 'weight' column in gdf_final: {sum_weight_gdf_final}")
# graphic comparation of weights of the two datasets
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
gdf_points_grid['weight'].plot(kind='hist', bins=5, alpha=0.5, label='Dataset 2')
gdf_final['weight'].plot(kind='hist', bins=5, alpha=0.5, label='old dataset')
plt.title("Distribution of Weights")
plt.xlabel("Weight")
plt.ylabel("Frequency")
plt.legend()
plt.show()

# now graphically compare the two datasets plot side by side the heat map of the two datasets
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
gdf_points_grid.plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
gdf_final.plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
ax[0].set_title("Heat Map of Points with Weights (Dataset 2)")
ax[0].set_xlabel("Easting")
ax[0].set_ylabel("Northing")
ax[1].set_title("Heat Map of Points with Weights (old dataset)")
ax[1].set_xlabel("Easting")
ax[1].set_ylabel("Northing")
plt.show()

#eliminate columns that are not necessary gdf_points_grid
gdf_points_grid = gdf_points_grid.drop(columns=['color'])

# Save the updated df_points_grid with the 'weight' column to a new CSV file
gdf_points_grid.drop(columns='geometry').to_csv('test.csv', index=False)



# print the slowest last five weights of the gdf_points_grid
# the slowest five weights of the gdf_points_grid with unique values
unique_weights = gdf_points_grid['weight'].unique()
# order
unique_weights.sort()
# print the slowest five weights
print(unique_weights[:10])

# same for the gdf_final
unique_weights = gdf_final['weight'].unique()
unique_weights.sort()
print(unique_weights[:10])


# compare the distribution of the smallest hundred unique weights of the two datasets
# get the smallest hundred weights of the gdf_points_grid
num_of_smallest_weights = 1
unique_weights = gdf_points_grid['weight'].unique()
unique_weights.sort()
smallest_weights_gdf_points_grid = unique_weights[:num_of_smallest_weights]
# get the smallest hundred weights of the gdf_final
unique_weights = gdf_final['weight'].unique()
unique_weights.sort()
smallest_weights_gdf_final = unique_weights[:num_of_smallest_weights]
# plot the distribution of the smallest hundred weights of the two datasets
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
gdf_points_grid[gdf_points_grid['weight'].isin(smallest_weights_gdf_points_grid)].plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
gdf_final[gdf_final['weight'].isin(smallest_weights_gdf_final)].plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})

ax[0].set_title(f"Heat Map of Points with Smallest Hundred Weights (Dataset 2) [size: {len(smallest_weights_gdf_points_grid)}]")
ax[0].set_xlabel("Easting")
ax[0].set_ylabel("Northing")
ax[1].set_title(f"Heat Map of Points with Smallest Hundred Weights (old dataset) [size: {len(smallest_weights_gdf_final)}]")
ax[1].set_xlabel("Easting")
ax[1].set_ylabel("Northing")
plt.show()







# import pandas as pd
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import re
# from sklearn.neighbors import NearestNeighbors

# def extract_coordinates(idcar):
#     """Extracts coordinates from the 'idcar' column."""
#     match = re.search(r'N(\d+)E(\d+)', idcar)
#     if match:
#         northing = int(match.group(1))
#         easting = int(match.group(2))
#         return easting, northing
#     else:
#         return None, None

# def load_population_data(filepath):
#     """Loads population data and extracts coordinates."""
#     df = pd.read_csv(filepath)
#     df['easting'], df['northing'] = zip(*df['idcar_1km'].apply(extract_coordinates))
#     df = df.dropna(subset=['easting', 'northing'])
#     gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['easting'], df['northing']))
#     gdf.set_crs(epsg=3035, inplace=True)
#     return gdf

# def load_grid_data(filepath):
#     """Loads grid data and converts CRS."""
#     df = pd.read_csv(filepath)
#     gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['longitude'], df['latitude']), crs='EPSG:4326')
#     gdf = gdf.to_crs(epsg=3035)
#     return gdf

# def find_nearest_weights(gdf_pop, gdf_grid):
#     """Performs spatial join to find nearest neighbors and calculate weights."""
#     gdf_grid['weight'] = 0
#     gdf_pop_coords = list(zip(gdf_pop.geometry.x, gdf_pop.geometry.y))
#     gdf_grid_coords = list(zip(gdf_grid.geometry.x, gdf_grid.geometry.y))
#     nn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(gdf_grid_coords)
#     distances, indices = nn.kneighbors(gdf_pop_coords)
#     gdf_pop['closest_idx'] = indices.flatten()
#     for idx in range(len(gdf_grid)):
#         gdf_grid.at[idx, 'weight'] = gdf_pop[gdf_pop['closest_idx'] == idx]['ind'].sum()
#     return gdf_grid

# def plot_datasets(gdf_pop, gdf_grid, title):
#     """Plots datasets side by side."""
#     fig, ax = plt.subplots(1, 2, figsize=(20, 10))
#     gdf_grid.plot(ax=ax[0], marker='x', color=gdf_grid['color'], markersize=5, label='Grid')
#     gdf_pop.plot(ax=ax[1], marker='o', color='red', markersize=5, label='Population')
#     ax[0].set_title(f"{title} - Grid")
#     ax[0].set_xlabel("Easting")
#     ax[0].set_ylabel("Northing")
#     ax[1].set_title(f"{title} - Population")
#     ax[1].set_xlabel("Easting")
#     ax[1].set_ylabel("Northing")
#     plt.show()

# def analyze_weights(gdf_grid):
#     """Analyzes and visualizes weights."""
#     num_zeros = (gdf_grid['weight'] == 0).sum()
#     num_non_zeros = (gdf_grid['weight'] != 0).sum()
#     sum_weight = gdf_grid['weight'].sum()
#     print(f"Number of zeros in the 'weight' column: {num_zeros}")
#     print(f"Number of non-zero values in the 'weight' column: {num_non_zeros}")
#     print(f"Sum of 'weight' column: {sum_weight}")
#     print(f"Max value in the 'weight' column: {gdf_grid['weight'].max()}")
#     print(f"Min value in the 'weight' column: {gdf_grid['weight'].min()}")
#     print(f"Mean value in the 'weight' column: {gdf_grid['weight'].mean()}")
#     print(f"Standard deviation of the 'weight' column: {gdf_grid['weight'].std()}")
#     gdf_grid['weight'].plot(kind='hist', bins=5)
#     plt.title("Distribution of Weights")
#     plt.xlabel("Weight")
#     plt.ylabel("Frequency")
#     plt.show()
#     fig, ax = plt.subplots(1, 1, figsize=(10, 10))
#     gdf_grid.plot(ax=ax, column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
#     plt.title("Heat Map of Points with Weights")
#     plt.xlabel("Easting")
#     plt.ylabel("Northing")
#     plt.show()

# def main():
#     # Load and process population data
#     gdf_pop = load_population_data('data/data_qgis/data_instance_paca/points_population_1km_paca_table.csv')
    
#     # Load and process grid data
#     gdf_grid = load_grid_data('data/data_qgis/data_instance_paca/points_grid_5km_paca_table.csv')
    
#     # Perform spatial join and calculate weights
#     gdf_grid = find_nearest_weights(gdf_pop, gdf_grid)
    
#     # Plot datasets side by side
#     plot_datasets(gdf_pop, gdf_grid, "Population vs Grid")
    
#     # Analyze weights
#     analyze_weights(gdf_grid)

#     # Additional operations if needed
#     # ...

# if __name__ == "__main__":
#     main()