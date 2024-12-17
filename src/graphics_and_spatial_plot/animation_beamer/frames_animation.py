# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap
# import os

# # Create directory to save frames
# if not os.path.exists('frames'):
#     os.makedirs('frames')

# # Number of frames
# n_frames = 20

# # Initialize the Basemap
# map = Basemap(projection='merc', llcrnrlat=-60, urcrnrlat=80, llcrnrlon=-180, urcrnrlon=180, resolution='c')

# # Create random locations
# np.random.seed(42)
# lats = np.random.uniform(-60, 80, size=50)
# lons = np.random.uniform(-180, 180, size=50)

# for i in range(1, n_frames+1):
#     plt.figure(figsize=(8, 6))
#     map.drawcoastlines()
#     map.drawcountries()

#     # Plot a subset of the random points
#     map.scatter(lons[:i], lats[:i], latlon=True, c='red', s=50, marker='o')

#     # Save frame as PNG
#     plt.savefig(f'frames/frame_{i:03d}.png')
#     plt.close()


import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from shapely.geometry import Point

# Load shapefile
shapefile = gpd.read_file("/home/falbuquerque/Documents/projects/GeoAvignon/Creation_Real_Instance/Decoupages_GIS/PACA_region.shp")

# Create random points for "factories" within the region bounds
import numpy as np
np.random.seed(42)
num_points = 10
points_inside = []

# Generate random points inside the region
while len(points_inside) < num_points:
    lat = np.random.uniform(shapefile.bounds.miny.min(), shapefile.bounds.maxy.max())
    lon = np.random.uniform(shapefile.bounds.minx.min(), shapefile.bounds.maxx.max())
    point = Point(lon, lat)
    
    # Check if the point is within the region
    if shapefile.contains(point).any():
        points_inside.append((lat, lon))

# Create random points for "services" outside the region
points_outside = []
while len(points_outside) < num_points:
    lat = np.random.uniform(shapefile.bounds.miny.min(), shapefile.bounds.maxy.max())
    lon = np.random.uniform(shapefile.bounds.minx.min(), shapefile.bounds.maxx.max())
    point = Point(lon, lat)
    
    # Check if the point is outside the region
    if not shapefile.contains(point).any():
        points_outside.append((lat, lon))

# Function to add emojis to the plot
def add_emoji(ax, lat, lon, emoji_path=None, size=0.4):
    """Add an emoji image to the plot."""
    image = plt.imread(emoji_path)
    im = OffsetImage(image, zoom=size)
    ab = AnnotationBbox(im, (lon, lat), frameon=False)
    ax.add_artist(ab)

# Create plot with a labeled box representing "Service" outside the PACA region
fig, ax = plt.subplots(figsize=(10, 10))
shapefile.plot(ax=ax, color='lightblue')

# Define the position of the "Service" label box
service_box_x = shapefile.bounds.maxx.max() + 0.5
service_box_y = shapefile.bounds.maxy.max() - 0.5

# Add the emoji to represent services
for lat, lon in points_outside:
    add_emoji(ax, lat, lon, emoji_path='/home/falbuquerque/Documents/projects/GeoAvignon/icons8-hospital-48.png')

# Create a box for "Service"
box_width = 1.0
box_height = 0.5
plt.gca().add_patch(plt.Rectangle((service_box_x, service_box_y - box_height), box_width, box_height, fill=True, color='white', edgecolor='black'))
plt.text(service_box_x + box_width / 2, service_box_y, 'Service', ha='center', va='bottom', fontsize=15)

# Remove axes
ax.set_axis_off()
plt.title("Locations Outside the PACA Region")
plt.show()

# Create second plot with emojis inside the region
fig, ax = plt.subplots(figsize=(10, 10))
shapefile.plot(ax=ax, color='lightblue')
for lat, lon in points_inside:
    add_emoji(ax, lat, lon, emoji_path='/home/falbuquerque/Documents/projects/GeoAvignon/icons8-hospital-48.png')

# Remove axes
ax.set_axis_off()
plt.title("Locations Inside the PACA Region")
plt.show()