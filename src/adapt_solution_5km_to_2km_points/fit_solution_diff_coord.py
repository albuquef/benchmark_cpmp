import pandas as pd
import numpy as np
import geopandas as gpd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from shapely.geometry import Point

# Function to parse coord_5km from file
def load_coord_5km(file_path):
    return pd.read_csv(
        file_path, 
        sep=r'\s+', 
        names=['location', 'capacity', 'coord_x', 'coord_y', 'id_loc', 'fid'],
        skiprows=1
    )

# Function to parse coord_2km from file
def load_coord_2km(file_path):
    return pd.read_csv(
        file_path, 
        sep=r'\s+', 
        names=['location', 'capacity', 'coord_x', 'coord_y', 'id_grid2km', 'fid'],
        skiprows=1
    )

# Function to extract p locations from solution file
def load_p_locations(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    start_idx = lines.index('P LOCATIONS\n') + 1
    locations = []
    for line in lines[start_idx:]:
        if line.strip() == '':
            break
        locations.append(int(line.strip()))
    return locations


for p in [135, 173, 192, 212, 250]:
    # Load the files
    coord_5km_file = 'outputs/PACA_5km/loc_capacities_cinema.txt'
    coord_2km_file = 'outputs/PACA_2km/loc_capacities_cinema.txt'
    solution_5km_file = f'outputs/solutions_files/2024-12-05_PACA_5km_cinema_FORMULATION/Assignments/output_paca_cinema_EXACT_CPMP_p_{p}_EXACT_CPMP.txt'

    coord_5km = load_coord_5km(coord_5km_file)
    coord_2km = load_coord_2km(coord_2km_file)
    p_locations = load_p_locations(solution_5km_file)

    # Ensure numeric coordinates
    coord_2km['coord_x'] = pd.to_numeric(coord_2km['coord_x'], errors='coerce')
    coord_2km['coord_y'] = pd.to_numeric(coord_2km['coord_y'], errors='coerce')
    coord_2km = coord_2km.dropna(subset=['coord_x', 'coord_y'])

    # Create GeoDataFrame for 5km and 2km locations
    gdf_5km = gpd.GeoDataFrame(
        coord_5km, 
        geometry=gpd.points_from_xy(coord_5km['coord_x'], coord_5km['coord_y']),
        crs="EPSG:900913"
    )

    # create a gdf for the p_locations
    gdf_p_locations = gpd.GeoDataFrame(
        {'location': p_locations},
        geometry=[Point(gdf_5km.loc[gdf_5km['location'] == loc][['coord_x', 'coord_y']].values[0]) for loc in p_locations],
        crs="EPSG:900913"
    )

    gdf_2km = gpd.GeoDataFrame(
        coord_2km, 
        geometry=gpd.points_from_xy(coord_2km['coord_x'], coord_2km['coord_y']),
        crs="EPSG:3035"
    )

    # Transform both GeoDataFrames to the same CRS (e.g., EPSG:3035)
    gdf_5km = gdf_5km.to_crs("EPSG:3035")
    gdf_2km = gdf_2km.to_crs("EPSG:3035")
    gdf_p_locations = gdf_p_locations.to_crs("EPSG:3035")

    print(gdf_2km['location'].isna().sum())

    # Initialize a set to track used locations in 2km
    used_2km_locations = set()

    p_locations_2km = []
    total_cap = 0

    for idx, p_location in gdf_p_locations.iterrows():
        # Filter 2km points by capacity
        valid_2km_points = gdf_2km[gdf_2km['capacity'] >= coord_5km.loc[coord_5km['location'] == p_location['location'], 'capacity'].values[0]]
        
        # If no valid points, skip
        if valid_2km_points.empty:
            p_locations_2km.append(None)
            continue

        # Reset index to avoid misalignment
        valid_2km_points = valid_2km_points.reset_index()

        # Calculate distances
        distances = valid_2km_points.geometry.distance(p_location.geometry)

        # Sort valid points by distance and iterate to find an unused one
        sorted_valid_points = valid_2km_points.loc[distances.sort_values().index]
        assigned_location = None
        for _, candidate_point in sorted_valid_points.iterrows():
            if candidate_point['location'] not in used_2km_locations:
                assigned_location = candidate_point
                break

        # If no unused location is found, skip
        if assigned_location is None:
            p_locations_2km.append(None)
            continue

        # Mark the location as used
        used_2km_locations.add(assigned_location['location'])

        # Append the location and update total capacity
        p_locations_2km.append(assigned_location['location'])
        total_cap += assigned_location['capacity']




    print('Total capacity:', total_cap)

    # Add the equivalence as a new column in the p_locations GeoDataFrame
    gdf_p_locations['2km_location'] = p_locations_2km

    print(gdf_p_locations.head())
    print(gdf_2km.loc[gdf_2km['location'].isin(p_locations_2km)].head())  




    # # number point ins gdf_2km
    # print(f'Points in gdf_2km: {len(gdf_2km)}')
    # # number of p_locations_2km
    # print(f'Points in p_locations_2km: {len(p_locations_2km)}')

    # # number point in gdf_2km that are in p_locations_2km
    print(f'Points in gdf_2km that are in p_locations_2km: {len(gdf_2km.loc[gdf_2km["location"].isin(p_locations_2km)])}')
    # # max min location in gdf_2km
    # print(f'Max location in gdf_2km: {gdf_2km["location"].max()}')
    # print(f'Min location in gdf_2km: {gdf_2km["location"].min()}')
    # print(f"Max location in p_locations_2km: {max(p_locations_2km)}")
    # print(f"Min location in p_locations_2km: {min(p_locations_2km)}")



    # print the locations in p_locations_2km that are not in gdf_2km    
    # print(set(p_locations_2km) - set(gdf_2km['location']))


    print(gdf_2km['capacity'].sum())    
    print(gdf_2km.loc[gdf_2km['location'].isin(p_locations_2km)]['capacity'].sum())
    print('Total capacity:', total_cap)
  


    DIR_TO_SAVE = 'outputs/p_location_5km_to_2km'
    NAME_OUTPUT_FILE = f'p_locations_5km_to_2km_p_{p}.txt'
    # save 1 line: number of p_locations 2:capacity sum of 2km locations 3: p_locations 
    with open(f'{DIR_TO_SAVE}/{NAME_OUTPUT_FILE}', 'w') as file:
        file.write(f'{len(p_locations)}\n')
        file.write(f'{total_cap}\n')
        file.write(f'P LOCATIONS\n')
        for p_location_2km in p_locations_2km:
            file.write(f'{p_location_2km}\n')   
    


    # # Plot the 5km locations, 2km locations, and the equivalence of p_locations to 2km locations
    # fig, ax = plt.subplots(figsize=(10, 8))
    # gdf_5km.plot(ax=ax, marker='o', color='blue', label='5km Locations')
    # gdf_2km.plot(ax=ax, marker='x', color='green', label='2km Locations')
    # gdf_p_locations.plot(ax=ax, marker='o', color='red', label='P Locations')
    # ax.set_aspect('auto')
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    # ax.set_title('5km Locations, 2km Locations, and P Locations')
    # ax.legend()
    # plt.show()

    # plot side by side the p_locations_5km with 5km and the other p_locations_2km with 2km points
    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    gdf_5km.plot(ax=ax[0], marker='o', color='blue', label='5km Locations')
    gdf_p_locations.loc[gdf_p_locations['location'].isin(p_locations)].plot(ax=ax[0], marker='o', color='red', label='P Locations')
    ax[0].set_aspect('auto')
    ax[0].set_xlabel('Longitude')
    ax[0].set_ylabel('Latitude')
    ax[0].set_title('5km Locations and P Locations')
    ax[0].legend()

    gdf_2km.plot(ax=ax[1], marker='x', color='green', label='2km Locations')
    gdf_p_locations.loc[gdf_p_locations['2km_location'].isin(p_locations_2km)].plot(ax=ax[1], marker='o', color='red', label='P Locations')
    ax[1].set_aspect('auto')
    ax[1].set_xlabel('Longitude')
    ax[1].set_ylabel('Latitude')
    ax[1].set_title('2km Locations and P Locations')
    ax[1].legend()
    plt.show()




