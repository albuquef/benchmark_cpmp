import sys
import os

# Set the QGIS paths
os.environ['PATH'] += r":/usr/local/qgis/bin"
os.environ['PYTHONPATH'] = r"/usr/local/qgis/share/qgis/python"
os.environ['LD_LIBRARY_PATH'] = r"/usr/local/qgis/lib"

# Add QGIS Python to system path
sys.path.append(r"/usr/local/qgis/share/qgis/python")
sys.path.append(r"/usr/local/qgis/share/qgis/python/plugins")

# Import QGIS modules
from qgis.core import QgsApplication, QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry, QgsPointXY, QgsField
from qgis.PyQt.QtCore import QVariant
from qgis.analysis import QgsVoronoiPolygons

import geopandas as gpd
from shapely.wkb import loads

# Initialize QGIS Application
QgsApplication.setPrefixPath("C:/OSGeo4W64/apps/qgis", True)
qgs = QgsApplication([], False)
qgs.initQgis()

# Function to create a memory layer with points
def create_point_layer(points):
    layer = QgsVectorLayer('Point?crs=EPSG:4326', 'points', 'memory')
    provider = layer.dataProvider()
    
    provider.addAttributes([QgsField('id', QVariant.Int)])
    layer.updateFields()
    
    for i, point in enumerate(points):
        feature = QgsFeature()
        feature.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(point[0], point[1])))
        feature.setAttributes([i])
        provider.addFeature(feature)
    
    layer.updateExtents()
    return layer

# Example points
points = [(10, 10), (20, 20), (30, 10), (40, 30)]

# Create point layer
point_layer = create_point_layer(points)

# Run Voronoi Polygons Algorithm
voronoi_layer = QgsVectorLayer('Polygon?crs=EPSG:4326', 'voronoi', 'memory')
QgsVoronoiPolygons.createVoronoiPolygons(point_layer, voronoi_layer)

# Extract features from Voronoi layer and convert to GeoDataFrame
features = voronoi_layer.getFeatures()
geometries = []
for feature in features:
    geom = feature.geometry()
    wkb = geom.asWkb()
    shapely_geom = loads(wkb)
    geometries.append(shapely_geom)

gdf = gpd.GeoDataFrame(geometry=geometries, crs="EPSG:4326")

# Clean up QGIS application
qgs.exitQgis()

# Display the GeoDataFrame
print(gdf)