import pandas as pd
import pyproj

# read csv file (with no header)
df_distances = pd.read_csv('/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Reseau/2km_grid/Matrice/mat_full.csv', header=None, sep=';')  

# create a header fid_origin, fid_dest, time 
df_distances.columns = ['fid_origin', 'fid_dest', 'weight_orig', 'weight_dest', 'dist_metres' ,'time_sec']

print(df_distances.head()) # print first 5 rows
# count the number of rows
print(len(df_distances))

# show 100 row to 105
print(df_distances[10410:10415])

# Read the new CSV file
df_customers = pd.read_csv('outputs/PACA_2km/cust_weights_2km.txt', sep='\s+')

# Check if all values in the "fid" column of the new DataFrame are present in "fid_origin" and "fid_dest" of the existing DataFrame
all_fids = set(df_distances['fid_origin']).union(set(df_distances['fid_dest']))
new_fids = set(df_customers['fid'])

# Check if all new_fids are in all_fids
all_present = new_fids.issubset(all_fids)

print("\n\nAll fids customers are present in the distance table:", all_present)

# Define the projection transformations using Transformer
transformer = pyproj.Transformer.from_crs("epsg:2154", "epsg:4326", always_xy=True)

# Convert coordinates
df_customers['x_WGS84'], df_customers['y_WGS84'] = transformer.transform(df_customers['x_LAMB93'].values, df_customers['y_LAMB93'].values)

# Save the updated DataFrame to a new CSV file (sep=\s+)
df_customers.to_csv('outputs/PACA_2km/cust_weights_2km.txt', index=False, sep=' ')

print("\n")

# Print head
print(df_customers.head())
print(df_customers.columns)

# add column "customer_origin" and "customer_dest" to the distance DataFrame
df_distances['customer_origin'] = df_distances['fid_origin'].map(df_customers.set_index('fid')['customer'])
df_distances['customer_dest'] = df_distances['fid_dest'].map(df_customers.set_index('fid')['customer'])



########### Save dist matrix to txt file

# change columns names fid_origin -> customer and fid_dest -> location  dist_metres -> distance
df_distances_metres = pd.DataFrame()
df_distances_metres.loc[:, 'customer'] = df_distances['customer_origin']
df_distances_metres.loc[:, 'location'] = df_distances['customer_dest']
df_distances_metres.loc[:, 'distance'] = df_distances['dist_metres']
df_distances_metres.to_csv('outputs/PACA_2km/dist_matrix_metres_2km.txt', index=False, sep=' ')

# another file fid_origin -> customer and fid_dest -> location  dist_sec -> distance
# generic dataframe
df_distances_sec = pd.DataFrame()
df_distances_sec.loc[:, 'customer'] = df_distances['customer_origin']
df_distances_sec.loc[:, 'location'] = df_distances['customer_dest']
df_distances_sec.loc[:, 'distance'] = df_distances['time_sec']

df_distances_sec.to_csv('outputs/PACA_2km/dist_matrix_seconds_2km.txt', index=False, sep=' ')