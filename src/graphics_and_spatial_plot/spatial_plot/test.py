### test 1

# import geopandas
# import plotly.express as px
# # Shapefile data
# shp = geopandas.read_file("https://github.com/R-CoderDotCom/data/blob/main/shapefile_spain/spain.zip?raw=true")
# shp = shp.to_crs("WGS84")

# fig = px.choropleth_mapbox(
#     data_frame = shp.set_index("ccaa_id"), # Using the id as index of the data
#     geojson = shp.geometry,                # The geometry
#     locations = shp.index,                 # The index of the data
#     color = 'unemp_rate',
#     mapbox_style = 'open-street-map',
#     center = dict(lat = 40.0, lon = -3.72),
#     zoom = 4)

# fig.show()




### test 2


import plotly.graph_objects as go

fig = go.Figure(go.Scattergeo())
fig.update_geos(
    visible=True, resolution=50, scope="usa",
    showcountries=True, countrycolor="Black",
    showsubunits=True, subunitcolor="Blue"
)
fig.update_layout(height=300, margin={"r":0,"t":0,"l":0,"b":0})
fig.show()