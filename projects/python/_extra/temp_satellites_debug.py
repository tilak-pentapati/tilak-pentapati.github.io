# temp_satellites_debug.py
import subprocess
import sys
import importlib

# List of required packages
required_packages = ["requests", "skyfield", "numpy", "pandas", "geopandas", "shapely", "matplotlib", "pyproj", "contextily", "mapclassify", "folium"]

# defining Function to install packages
def install_package(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Checking and installing missing packages
for package in required_packages:
    try:
        importlib.import_module(package)
    except ImportError:
        print(f"{package} not found, installing...")
        install_package(package)

# Importing Required Packages with aliases and functions
import requests
from skyfield.api import load, wgs84, EarthSatellite
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, shape, Polygon, MultiPolygon
import json
from shapely.ops import transform
import pyproj
import matplotlib.pyplot as plt
import contextily as ctx
import folium
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mapclassify import Quantiles
import matplotlib.font_manager as fm

# Using Celes Trak, to get TLE Data for active satellites
try:
    Sat_TLE = requests.get("https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle")
    if Sat_TLE.status_code == 200:
     print("Data Retrieved successfully.")
except requests.exceptions.RequestException as e:
    print(f"An error occurred: {e}")

# first we split after every line, then we strip '\r' at end of every line
tle_lines = [line.strip('\r') for line in Sat_TLE.text.strip().split('\n')]

if len(tle_lines) % 3 != 0:
    print("Warning: The TLE data might be incomplete or corrupted.")

satellites = [tle_lines[i:i+3] for i in range(0, len(tle_lines), 3)]

satellites_data = []
for sat in satellites:
    satellite_info = {
        'name': sat[0],
        'line1': sat[1],
        'line2': sat[2]
    }
    satellites_data.append(satellite_info)

satellites_df = pd.DataFrame(satellites_data)
print("satellites_df head:")
print(satellites_df.head())

ts = load.timescale()

start_datetime = datetime(2024, 3, 4, 0, 0, 0)
times = [(start_datetime + timedelta(minutes=240*i)).strftime('%Y-%m-%dT%H:%M:%SZ') for i in range(6)]

def calculate_positions(name, line1, line2):
    satellite = EarthSatellite(line1, line2, name, ts)
    year, month, day = 2024, 3, 4
    times = ts.utc(year, month, day, 0, range(0, 1440, 240))
    geocentric = satellite.at(times)
    subpoint = wgs84.subpoint(geocentric)
    lons, lats = subpoint.longitude.degrees, subpoint.latitude.degrees
    return [f"{lon},{lat}" for lon, lat in zip(lons, lats)]

positions_dict = {}

for index, row in satellites_df.iterrows():
    positions = calculate_positions(row['name'], row['line1'], row['line2'])
    positions_dict[row['name']] = positions

positions_df = pd.DataFrame.from_dict(positions_dict, orient='index', columns=times)
positions_df.reset_index(inplace=True)
positions_df.rename(columns={'index': 'Satellite'}, inplace=True)

for timestamp in times:
    positions_df[timestamp] = positions_df[timestamp].apply(lambda x: Point([float(coord) for coord in x.split(',')]))

print("positions_df head:")
print(positions_df.head())

base_url = 'https://public.opendatasoft.com/api/explore/v2.1/catalog/datasets/world-administrative-boundaries/records'
total_count = 256
limit = 100
offset = 0
all_countries = []

while offset < total_count:
    response = requests.get(base_url, params={'limit': limit, 'offset': offset})
    if response.status_code == 200:
        data = response.json()
        all_countries.extend(data['results'])
        offset += limit
    else:
        print(f"Failed to fetch data: HTTP {response.status_code}")
        break

results = all_countries
geometries = []
properties = []

for item in results:
    if 'geo_shape' in item and 'geometry' in item['geo_shape'] and 'name' in item:
        geom = shape(item['geo_shape']['geometry'])
        geometries.append(geom)
        properties.append(item['name'])
    else:
        geometries.append(None)
        properties.append(None)

countries_df = pd.DataFrame({'Countries' : properties})
countries_gdf = gpd.GeoDataFrame(countries_df, geometry=geometries)
countries_gdf.set_crs(epsg=4326, inplace=True)

print("countries_gdf head:")
print(countries_gdf.head())

def ensure_multipolygon(geometry):
    if geometry is None:
        return None
    elif isinstance(geometry, Polygon):
        return MultiPolygon([geometry])
    elif isinstance(geometry, MultiPolygon):
        return geometry
    else:
        raise ValueError("Geometry must be a Polygon or MultiPolygon.")

countries_gdf['geometry'] = countries_gdf['geometry'].apply(ensure_multipolygon)

def shift_geom(geometry, offset):
    if geometry is None:
        return None
    def shift_special(x, y, z=None):
        if x < 0:
            new_x = x + 2 * offset
        else:
            new_x = x
        return new_x, y
    return transform(shift_special, geometry)

russia_gdf_bs= countries_gdf[countries_gdf['Countries'] == 'Russian Federation']
russia_gdf = countries_gdf[countries_gdf['Countries'] == 'Russian Federation'].copy()
russia_gdf['geometry'] = russia_gdf['geometry'].apply(lambda geom: shift_geom(geom, 180))

fig, axs = plt.subplots(1, 2, figsize=(15, 14))
russia_gdf_bs.plot(ax=axs[0], color='blue')
axs[0].set_title('Geography of Russia on 2D Map - Before using Shift Function') 
ymin, ymax = axs[0].get_ylim()
axs[0].set_ylim(ymin*0.6, ymax * 1.4)
russia_gdf.plot(ax=axs[1], color='green')
axs[1].set_title('Geography of Russia on 2D Map - After using Shift Function')  
plt.tight_layout() 
plt.savefig("debug_output_1.png")

counts_dict = {}
countries_gdf_copy = countries_gdf.copy()
for timestamp in positions_df.columns[1:]:
    temp_df = positions_df[['Satellite', timestamp]].copy()
    temp_df.rename(columns={timestamp: 'position'}, inplace=True)
    temp_gdf = gpd.GeoDataFrame(temp_df, geometry='position')
    temp_gdf.set_crs(epsg=4326, inplace=True)
    countries_gdf_copy = countries_gdf_copy.dropna(subset=['Countries'])
    joined_gdf = gpd.sjoin(temp_gdf, countries_gdf_copy, how="inner", predicate='intersects')
    counts_series = joined_gdf.groupby('Countries').size()
    counts_series = counts_series.reindex(countries_gdf_copy['Countries'].unique(), fill_value=0)
    counts_dict[timestamp] = counts_series

counts_df = pd.DataFrame(counts_dict)
counts_df_reset = counts_df.reset_index()
counts_df_reset.rename(columns={'index': 'Countries'}, inplace=True)

temporal_satellite_cm = pd.merge(countries_gdf_copy, counts_df_reset, how='left', on='Countries')
print("temporal_satellite_cm head:")
print(temporal_satellite_cm.head())

def add_scalebar(ax, length, location=(0.05, 0.05), linewidth=3,
                 units='km', fontsize=10, color='black'):
    x = location[0]
    y = location[1]
    xend = x + length / 44448
    ax.plot([x, xend], [y, y], color=color, linewidth=linewidth, transform=ax.transAxes)
    ax.text((x + xend) / 1.35, y - 0.01, f'{length} {units}',
            horizontalalignment='center', verticalalignment='top',
            transform=ax.transAxes, fontsize=fontsize, color=color)

fig, axs = plt.subplots(2, 1, figsize=(20, 15))
temporal_satellite_cm.plot(column='2024-03-04T00:00:00Z', ax=axs[0], legend=True,
                           cmap='Greens', edgecolor='black', lw=0.6)
axs[0].set_title(f"No. of Satellites Overhead per country - {times[0]}")
ymin, ymax = axs[0].get_ylim()
axs[0].set_ylim(ymin, ymax * 1.3)
xmin, xmax = axs[0].get_xlim()
axs[0].set_xlim(xmin*0.975, xmax * 1)
add_scalebar(axs[0], length=5556, location=(0.01, 0.05), units='km')

temporal_satellite_cm.plot(column='2024-03-04T00:00:00Z', ax=axs[1], legend=True,
                           cmap='YlGnBu', scheme='naturalbreaks', k=5, 
                           edgecolor='black', lw=0.6, legend_kwds={'loc': 'center right'})
axs[1].set_title(f"No. of Satellites Overhead per country (Natural Breaks) - {times[0]}")
add_scalebar(axs[1], length=5556, location=(0.01, 0.05), units='km')

plt.tight_layout()
plt.savefig("debug_output_2.png")
print("Plots saved and execution finished successfully.")
