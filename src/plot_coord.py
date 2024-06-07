import pandas as pd
import matplotlib.pyplot as plt
# import geopandas as gpd
# from shapely.geometry import Point

def plot_points_Lit(file_path, instance_name):
    # with open(file_txt, 'r') as f:
    #     lines = f.readlines()
    #     lines = [line.strip() for line in lines]
    #     lines = [line.split() for line in lines]
    #     df = pd.DataFrame(lines, columns=['location', 'capacity', 'coord_x', 'coord_y'])
    #     df['location'] = df['location'].astype(int)
    #     df['capacity'] = df['capacity'].astype(float)
    #     df['coord_x'] = df['coord_x'].astype(float)
    #     df['coord_y'] = df['coord_y'].astype(float)
    # Specify the path to your text file


    # Read the data into a DataFrame
    df = pd.read_csv(file_path, delim_whitespace=True)
    # plot points df
    plt.figure(figsize=(10, 6))
    plt.scatter(df['coord_x'], df['coord_y'], c='green', marker='o')
    plt.xlabel('coord_x')
    plt.ylabel('coord_y')
    plt.title(f'Points from {instance_name}')
    plt.savefig(f'plots/plots_lit/points_plot_lit_{instance_name}.png')
    plt.show()
    
def plot_points(csv_file, method='matplotlib', title='Points from CSV'):
    df = pd.read_csv(csv_file)
    
    if method == 'matplotlib':
        # Extract coordinates
        x = df['x']
        y = df['y']
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.scatter(x, y, c='green', marker='o')
        
        # # Add labels and title
        # for i, row in df.iterrows():
        #     plt.text(row['x'], row['y'], row['identif'], fontsize=9, ha='right')
        
        plt.xlabel('coord_x')
        plt.ylabel('coord_y')
        plt.title(title)
        # plt.grid(True)
        plt.show()
    
    # elif method == 'geopandas':
    #     # Create a GeoDataFrame
    #     geometry = [Point(xy) for xy in zip(df['x'], df['y'])]
    #     gdf = gpd.GeoDataFrame(df, geometry=geometry)
        
    #     # Plotting
    #     fig, ax = plt.subplots(figsize=(10, 6))
    #     gdf.plot(ax=ax, color='blue', markersize=10)
        
    #     # Add labels
    #     for x, y, label in zip(gdf.geometry.x, gdf.geometry.y, gdf['identif']):
    #         ax.text(x, y, label, fontsize=9, ha='right')
        
    #     plt.xlabel('coord_x')
    #     plt.ylabel('coord_y')
    #     plt.title(title)
    #     # plt.grid(True)
    #     plt.show()

# plot_points('./data/PACA/points_coord.csv', method='matplotlib', title='PACA region')
# instance_name = 'p3038_600'
# plot_points_Lit(f'./data/Literature/group3/loc_capacities_{instance_name}.txt', instance_name)
# instance_name = 'rl1304_010'
# instance_name = 'pr2392_020'
# instance_name = 'fnl4461_0100'
# plot_points_Lit(f'./data/Literature/group5/loc_capacities_{instance_name}.txt', instance_name)
# plot_points_Lit('./data/Literature/group5/loc_capacities_rl1304_010.txt')
# plot_points_Lit('./data/Literature/group5/loc_capacities_pr2392_020.txt')
# plot_points_Lit('./data/Literature/group5/loc_capacities_fnl4461_0100.txt')

instance_name = 'spain737_74_1'
plot_points_Lit(f'./data/Literature/group4/loc_capacities_{instance_name}.txt', instance_name)