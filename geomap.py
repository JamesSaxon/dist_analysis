#!/usr/bin/env python

from fiona.crs import from_epsg
import geopandas as gpd
from shapely.geometry import Point
from glob import glob
import matplotlib.pyplot as plt

def map_format(ax, on = False):

    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.margins(0,0)
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    if not on:
        ax.set_axis_off()
        ax.set_axis_on()
        for a in ["bottom", "top", "right", "left"]:
            ax.spines[a].set_linewidth(0)

    return ax


def make_pa_pittsburgh(f, e = 3364):

    gdf = gpd.read_file(f).to_crs(epsg = e)
    pittsburgh = gpd.GeoDataFrame(crs = from_epsg(4326), geometry = [Point(-79.9959, 40.4406)]).to_crs(epsg = 3364)

    ax = gdf.plot(column = "id", cmap = "nipy_spectral", alpha = 0.4, figsize = (4, 4))
    gdf.plot(color = "none", alpha = 1, edgecolor = "k", linewidth = 0.8, ax = ax)

    pittsburgh.plot(facecolor = "none", edgecolor = "k", markersize = 410, ax = ax)
    pittsburgh.plot(facecolor = "none", edgecolor = "k", markersize = 490, ax = ax)
    pittsburgh.plot(facecolor = "none", edgecolor = "w", markersize = 450, ax = ax)

    ax.set_xlim(290000, 900000)
    map_format(ax)

    # ax.figure.savefig("paper_figs/pa_ex/" + f.split("/")[-1].replace("geojson", "pdf"))
    ax.figure.savefig("paper_figs/pa_ex/" + f.split("/")[-1].replace("geojson", "png"))

    plt.close("all")


for fm in glob("data/pa_ex/*.geojson"):
    make_pa_pittsburgh(fm)

