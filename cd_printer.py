#!/usr/bin/env python 

import pandas as pd
import geopandas as gpd

from fiona.crs import from_epsg

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")

import matplotlib.pyplot as plt

from geomap import map_format

def get_state_cd_map(usps, cd, session = 114, color = "k"):

    g = gpd.GeoDataFrame.from_postgis("""SELECT cd, cd.geom, epsg 
                                         FROM cd JOIN states ON cd.state = fips
                                         WHERE states.usps = UPPER('{}') AND cd = {} AND sessn = {};""".format(usps, cd, session),
                                      con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                             host = "saxon.harris.uchicago.edu", port = 5432),
                                      geom_col = "geom", crs = from_epsg(2163), index_col = "cd")

    g = g.to_crs(epsg = g.loc[cd].epsg)

    ax = g.plot(color = color, alpha = 1, linewidth = 0, figsize = (5, 5))

    ax.set_axis_off()

    color_label = ""
    if color != "k":
        color_label = "_" + color.replace("#" "")

    map_format(ax)
    ax.figure.savefig("paper_figs/bad_districts/{}_{}{}.pdf".format(usps, cd, color_label), bbox_inches = 'tight', pad_inches = 0.05)

    plt.close("all")


for st, cd in [["nc", 12], ["oh", 9], ["fl", 5], ["ny", 10], ["tx", 35]]:
  get_state_cd_map(st, cd)

# get_state_cd_map("nc", 12, color = "#2c72ff")
# get_state_cd_map("oh",  9, color = "#5a15b5")
# get_state_cd_map("fl",  5, color = "#18b240")
# get_state_cd_map("ny", 10, color = "#cc0d0d")
# get_state_cd_map("tx", 35, color = "#baae15")
  
# get_state_cd_map("pa",  7, color = "k")

# get_state_cd_map("tx", 35, color = "b")
# get_state_cd_map("il",  4, color = "k")

