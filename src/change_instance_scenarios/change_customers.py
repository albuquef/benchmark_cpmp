import csv 
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import geopandas as gpd



def create_shuffle_customers(df_customers_original):
    # set a seed 
    seed = 42
    
    df_customers = df_customers_original.copy()
    
    # Shuflle only the weights  and set a seed
    df_customers['weight'] = df_customers['weight'].sample(frac=1, random_state=seed).reset_index(drop=True)
    
    # create a txt with the shuffled weights and columns: customer weight coord_x coord_y indetif
    filename = 'data/PACA/cust_weights_shuffled.txt'
    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(['customer', 'weight', 'coord_x', 'coord_y', 'identif'])
        for index, row in df_customers.iterrows():
            writer.writerow([row['customer'], row['weight'], row['x'], row['y'], row['identif']])   
    
    
    return df_customers

def create_split_customers(df_customers_original):
    # set a seed 
    seed = 42
    
    df_customers = df_customers_original.copy()
    total_weight = df_customers['weight'].sum()
    num_customers = len(df_customers)
    
    # Generate random weights that sum up to total_weight
    random_weights = np.random.rand(num_customers)
    random_weights /= random_weights.sum()  # Normalize to sum up to 1
    random_weights *= total_weight  # Scale to total_weight
    
    
    # change the weights of the customers for the random weights
    df_customers['weight'] = random_weights
    
    # create a txt with the shuffled weights and columns: customer weight coord_x coord_y indetif
    filename = 'data/PACA/cust_weights_split.txt'
    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(['customer', 'weight', 'coord_x', 'coord_y', 'identif'])
        for index, row in df_customers.iterrows():
            writer.writerow([row['customer'], row['weight'], row['x'], row['y'], row['identif']])
    
    return df_customers
    


df_cust_weights = pd.read_csv('data/PACA/cust_weights.txt', delim_whitespace=True, names=['customer', 'weight'])
df_map_id_cust_loc = pd.read_csv('data/PACA/map_id_cust_loc.txt', delim_whitespace=True, names=['id', 'identif'])
df_locations = pd.read_csv('data/PACA/locations_paca_2017_coord.csv')

# Ensure that 'identif' columns are of the same type (string)
df_map_id_cust_loc['identif'] = df_map_id_cust_loc['identif'].astype(str)
df_locations['identif'] = df_locations['identif'].astype(str)

# Merge the DataFrames to associate weights with coordinates
df_merged = pd.merge(df_map_id_cust_loc, df_cust_weights, left_on='id', right_on='customer')
df_dataset = pd.merge(df_merged, df_locations, left_on='identif', right_on='identif')
df_dataset['weight'] = df_dataset['weight'].astype(float)

print(df_dataset.head(5))   

# send a copy of the dataset to the function that shuffles the weights
df_customers_shuffle = create_shuffle_customers(df_dataset)
# df_customers_shuffle = create_shuffle_customers(copy(df_dataset))

print(df_customers_shuffle.head(5))


df_customers_split = create_split_customers(df_dataset)

print(df_customers_split.head(5))


# create heatmaps using geopandas and seaborn comparing the original and shuffled weights side by side
gdf = gpd.GeoDataFrame(df_dataset, geometry=gpd.points_from_xy(df_dataset['x'], df_dataset['y']))
gdf_shuffle = gpd.GeoDataFrame(df_customers_shuffle, geometry=gpd.points_from_xy(df_customers_shuffle['x'], df_customers_shuffle['y']))
gdf_split = gpd.GeoDataFrame(df_customers_split, geometry=gpd.points_from_xy(df_customers_split['x'], df_customers_split['y']))
# set the crs
gdf = gdf.set_crs("EPSG:3035")
gdf_shuffle = gdf_shuffle.set_crs("EPSG:3035")
gdf_split = gdf_split.set_crs("EPSG:3035")
# heatmaps
fig, ax = plt.subplots(1, 3, figsize=(20, 10))
gdf.plot(ax=ax[0], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
gdf_shuffle.plot(ax=ax[1], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
gdf_split.plot(ax=ax[2], column='weight', legend=True, legend_kwds={'label': "Weight", 'orientation': "horizontal"})
ax[0].set_title(f"Heat Map of Grid Population (Original dataset)\nnum elements= {len(gdf)}, sum = {gdf['weight'].sum()}")
ax[0].set_xlabel("Easting")
ax[0].set_ylabel("Northing")
ax[1].set_title(f"Heat Map of Grid Population (Shuffle dataset)\nnum elements={len(gdf_shuffle)}  , sum wi= {gdf_shuffle['weight'].sum()}")
ax[1].set_xlabel("Easting")
ax[1].set_ylabel("Northing")
ax[2].set_title(f"Heat Map of Grid Population (Split dataset)\nnum elements={len(gdf_split)}  , sum wi= {gdf_split['weight'].sum()}")
ax[2].set_xlabel("Easting")
ax[2].set_ylabel("Northing")
plt.show()



