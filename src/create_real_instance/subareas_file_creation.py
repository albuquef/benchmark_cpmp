import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import matplotlib.pyplot as plt
# voronoi ="" 
grid = ""
style_decoupage = "grid_" # "voronoi_" "grid_" or ""
# style_decoupage = ""
type_points = "5km"
type_decoupage = "commune_"  # "Arrondissement" "EPCI" "canton" "commune" "departement" "region"
DECOUPAGE_FILE = f"{style_decoupage}{type_decoupage}{type_points}" # "Arrondissement" "EPCI" "canton" "commune" "departement" "region"
# Load the territorial division shapefile

shp_path = f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/{DECOUPAGE_FILE}.shp"
if (style_decoupage == "voronoi_"):
    shp_path = f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_voronoi/{DECOUPAGE_FILE}.shp"
    if(type_decoupage == "Arrondissement_" or type_decoupage == "EPCI_"):
        filename = f"{style_decoupage}{type_decoupage}" # remove the last _ from the filename
        filename = filename[:-1]
        shp_path = f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/voronoi_check/{filename}.shp"
if (style_decoupage == "grid_"):
    shp_path = f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_grid/{DECOUPAGE_FILE}.shp"
    if(type_decoupage == "Arrondissement_" or type_decoupage == "EPCI_" or type_decoupage == "canton_"):
        filename = f"{style_decoupage}{type_decoupage}" # remove the last _ from the filename
        filename = filename[:-1]
        shp_path = f"/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/decoupages_grid/{filename}.shp"



subareas_gdf = gpd.read_file(shp_path)
# subareas_gdf = subareas_gdf.to_crs(epsg=3035)


# Inspect the columns to find the correct subarea name column
print(subareas_gdf.columns)

# Assuming the column that identifies subareas is named 'subarea_name_column' (adjust based on inspection)
# Create sequential subarea names like 'subarea_1', 'subarea_2', etc.
subarea_names = [(int)(i+1) for i in range(len(subareas_gdf))]

# Add the new 'subarea' column with sequential names
subareas_gdf['subarea'] = subarea_names

# print(subareas_gdf.head())  
# print("Number of subareas: ", len(subareas_gdf))

point_df_path = f"outputs/PACA_{type_points}/cust_weights_{type_points}.txt"
points_df = pd.read_csv(point_df_path, sep=" ")
if (type_points == "2km"):  
    points_gdf = gpd.GeoDataFrame(points_df, geometry=gpd.points_from_xy(points_df.x_LAMB93, points_df.y_LAMB93), crs="EPSG:2154")
    points_gdf = points_gdf.to_crs(epsg=4326)
    # points_df = points_gdf.to_crs(epsg=3035)
if (type_points == "5km"):
    # points_df = pd.read_csv("outputs/PACA_5km/cust_weights_5km.txt", sep=" ")
    points_gdf = gpd.GeoDataFrame(points_df, geometry=gpd.points_from_xy(points_df.coord_x, points_df.coord_y), crs="EPSG:900913")
    points_gdf = points_gdf.to_crs(epsg=4326)


# Reproject subareas_gdf to match the CRS of points_gdf
subareas_gdf = subareas_gdf.to_crs(epsg=4326)
print(subareas_gdf.head())  
print("Number of subareas: ", len(subareas_gdf))

# Perform the spatial join (left join)
joined_gdf = gpd.sjoin(points_gdf, subareas_gdf, how="left", predicate="intersects")

subarea_code = 'FID'
if (style_decoupage == ""):
    if type_decoupage == "Arrondissement_":
        subarea_code = 'codearondo'
    elif type_decoupage == "EPCI_":
        subarea_code = 'NOM_EPCI'
    elif type_decoupage == "canton_":
        subarea_code = 'numcant'
    elif type_decoupage == "commune_":
        subarea_code = 'NOM_COMM'

# Select the subarea and customer columns and save the result
output_df = joined_gdf[["subarea", "customer", subarea_code, "geometry"]]

# check if every customer is in a subarea
print("Number of customers: ", len(output_df))
print("Number of subareas: ", len(subareas_gdf))
# print if there are any customers not in a subarea
print("Number of customers not in a subarea: ", len(output_df[output_df["subarea"].isna()]))
print("Number of customers in a subarea: ", len(output_df[output_df["subarea"].notna()]))

print(output_df.head())
# Ensure geometries are valid
output_df = output_df[output_df.is_valid]

# # Ensure subareas_gdf contains polygon geometries
# if all(geom in ['Polygon', 'MultiPolygon'] for geom in subareas_gdf.geom_type.unique()):
#     fig, ax = plt.subplots()
    
#     # Plot the subareas as filled polygons
#     subareas_gdf.plot(ax=ax, color='blue', edgecolor='black')
    
#     # Plot the points that do not have a subarea
#     output_df[output_df["subarea"].isna()].plot(ax=ax, color='red', marker='o', markersize=5)

#     # Set aspect ratio manually if needed
#     ax.set_aspect('equal')
    
#     plt.show()
# else:
#     print("subareas_gdf does not contain polygon geometries.")

# add the points not within a subarea to the closest subarea
if len(output_df[output_df["subarea"].isna()]) > 0:
    print("Number of customers not in a subarea: ", len(output_df[output_df["subarea"].isna()]))
    print('Fitting the points not in a subarea to the closest subarea')
    # get the points not in a subarea
    points_not_in_subarea = output_df[output_df["subarea"].isna()]
    # calculate the distance between each point and each subarea
    for point_index, row in points_not_in_subarea.iterrows():
        point = row["geometry"]
        min_distance = float("inf")
        closest_subarea = None
        subarea_identifier = None
        for subarea_index, subarea in subareas_gdf.iterrows():
            distance = point.distance(subarea["geometry"])
            if distance < min_distance:
                min_distance = distance
                closest_subarea = subarea["subarea"]
                subarea_identifier = subarea[subarea_code]
        output_df.at[point_index, "subarea"] = closest_subarea
        output_df.at[point_index, subarea_code] = subarea_identifier
        print(f"Point {point_index} fitted to subarea {closest_subarea}")

    # check if every customer is in a subarea
    print("Number of customers: ", len(output_df))
    print("Number of subareas: ", len(subareas_gdf))
    # print if there are any customers not in a subarea
    print("Number of customers not in a subarea: ", len(output_df[output_df["subarea"].isna()]))


if (len(output_df[output_df["subarea"].isna()]) == 0):
    # remove geometry column
    output_df = output_df.drop(columns=["geometry"])

    # subare collumn to int
    output_df["subarea"] = output_df["subarea"].astype(int)

    #change name collumn "customer" to "location"  
    output_df = output_df.rename(columns={"customer": "location"})
    
    # change the order of subarea and location
    output_df = output_df[["location", "subarea", subarea_code]]

    # count number of different subareas
    num_subareas = len(output_df["subarea"].unique())
    print(f"Number of subareas final (with point inside): {num_subareas}")

    if DECOUPAGE_FILE == f"{style_decoupage}Arrondissement_{type_points}":
        output_df.to_csv(f"outputs/PACA_{type_points}/loc_coverages_{style_decoupage}arrond_{type_points}.txt", sep=" ", index=False)
    elif DECOUPAGE_FILE == f"{style_decoupage}EPCI_{type_points}":
        output_df.to_csv(f"outputs/PACA_{type_points}/loc_coverages_{style_decoupage}EPCI_{type_points}.txt", sep=" ", index=False)
    elif DECOUPAGE_FILE == f"{style_decoupage}canton_{type_points}":
        output_df.to_csv(f"outputs/PACA_{type_points}/loc_coverages_{style_decoupage}canton_{type_points}.txt", sep=" ", index=False)
    elif DECOUPAGE_FILE == f"{style_decoupage}commune_{type_points}":
        output_df.to_csv(f"outputs/PACA_{type_points}/loc_coverages_{style_decoupage}commune_{type_points}.txt", sep=" ", index=False)


# create a file txt with the subarea and customer
# output_df.to_csv("outputs/PACA_2km/subarea_customer_mapping.txt", sep=" ", index=False)




# output_df.to_csv("subarea_customer_mapping.csv", index=False)




