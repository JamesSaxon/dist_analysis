from IPython.display import display

from fiona.crs import from_epsg
import geopandas as gpd
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(rc={"figure.figsize": (6, 3)})
sns.set_style("white")

import psycopg2
from netrc import netrc
user, acct, passwd = netrc().authenticators("harris")

import os, glob

import matplotlib.pyplot as plt

import psycopg2

import jellyfish as jf

def jf_jw_match(x, list_strings):
    return sorted(list_strings, key=lambda y: jf.jaro_winkler(x, y), reverse=True)


def map_seats(name, paths, tract_votes, years, final_table):
    
    if not len(paths): return
    
    outcomes = pd.DataFrame(columns = [str(y) for y in years])
    
    for xi, xpath in enumerate(paths):
        
        map_df = pd.read_csv(xpath, names = ["rn", "region"], header = None)
        mapped_votes = tract_votes.join(map_df)
        cd_seats = mapped_votes.groupby("region").sum()
        cd_seats["map_id"] = xi
        
        for y in years:
            var = y if type(y) is str else ("%02d" % (y%100))
            cd_seats["%s D Fr" % y] = cd_seats["D" + var] / (cd_seats["D" + var] + cd_seats["R" + var])
            cd_seats[str(y)] = cd_seats["D" + var] > cd_seats["R" + var]
        
        outcomes = outcomes.append(cd_seats)
    
    final_table.loc[name] = (outcomes > 0.5).mean()

    outcomes["weights"] = 1 / len(paths)

    outcols = ["map_id", "weights"]
    for y in years: outcols += [str(y), "%s D Fr" % y]
    
    return outcomes[outcols]


def cdmap_seats(name, fips, epsg, sessn, tract_votes, years, final_table):
    
    trcd_df = gpd.GeoDataFrame.from_postgis("""SELECT rn.rn, cd, 
                                                      ST_Centroid(ST_Transform(tr.geom, states.epsg)) geometry
                                               FROM census_tracts_2015 as tr
                                               JOIN (SELECT state, county, tract,
                                                         row_number() over (PARTITION BY state ORDER BY county, tract NULLS LAST) - 1 as rn
                                                     FROM census_tracts_2015) rn ON
                                                  tr.state  = rn.state  AND
                                                  tr.county = rn.county AND
                                                  tr.tract  = rn.tract
                                               JOIN states ON tr.state = states.fips
                                               JOIN cd ON cd.sessn = {} AND
                                                  cd.state = tr.state AND
                                                  ST_Covers(cd.geom, ST_Centroid(tr.geom))
                                               WHERE tr.state = {} ORDER BY rn;
                                               """.format(sessn, fips),
                                               con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                                      host = "saxon.harris.uchicago.edu", port = 5432), 
                                               geom_col = "geometry", crs = from_epsg(epsg), index_col = "rn")

    mapped_votes = tract_votes.join(trcd_df)
    cd_seats = mapped_votes.groupby("cd").sum()
    cd_seats["map_id"] = 0

    for y in years:
        var = y if type(y) is str else ("%02d" % (y%100))
        cd_seats["%s D Fr" % y] = cd_seats["D" + var] / (cd_seats["D" + var] + cd_seats["R" + var])
        cd_seats[str(y)] = cd_seats["D" + var] > cd_seats["R" + var]
    
    final_table.loc[name] = cd_seats.mean()

    outcols = ["map_id"]
    for y in years: outcols += [str(y), "%s D Fr" % y]
    
    return cd_seats[outcols]


def seat_table(usps, years):

    votes = pd.read_csv(usps + "_votes.csv", index_col = "rn")

    st = pd.read_sql("select seats, epsg, lower(usps) usps, fips from states where usps = upper('{}');".format(usps),
                     con = psycopg2.connect(database = "census", user = user, password = passwd,
                                            host = "saxon.harris.uchicago.edu", port = 5432)).ix[0].to_dict()
    
    epsg, fips = st["epsg"], st["fips"]
    
    final_table = pd.DataFrame(columns = [str(y) for y in years],
                               index   = ["Statewide Votes", "107th Congress", "111th Congress", "114th Congress",
                                          "Power Diagram", "Split-Line", "Isoperimeter Quotient", "Rohrbach", 
                                          "Exchange", "Population Hull", "Dynamic Radius", "Inscribed Circles",
                                          "Circumscribing Circles", "Path Fraction"])

    final_table.columns.name = usps.upper()
    
    final_table.loc["Statewide Votes"] = {str(y) :  votes["D{:02d}".format(y % 100)].sum() / \
                                                   (votes["D{:02d}".format(y % 100)].sum() + \
                                                    votes["R{:02d}".format(y % 100)].sum())
                                          for y in years}
    
    stdir = "../cluscious/res/{}/{}/*/*.csv"

    metrics = {"Dynamic Radius"        : "dyn_radius",  "Inscribed Circles"      : "ehrenburg",
               "Exchange"              : "exchange",    "Population Hull"        : "hull_p",
               "Isoperimeter Quotient" : "polsby",      "Circumscribing Circles" : "reock", 
               "Rohrbach"              : "rohrbach",    "Power Diagram"          : "power",
               "Split-Line"            : "split"}


    seat_res = {}
    for ti, m in metrics.items():
        # print(ti, sorted(glob.glob(stdir.format(usps, m))))
        seat_res[m] = map_seats(ti, glob.glob(stdir.format(usps, m)), votes, years, final_table)

    seat_res["107"] = cdmap_seats("107th Congress", fips, epsg, 107, votes, years, final_table)
    seat_res["111"] = cdmap_seats("111th Congress", fips, epsg, 111, votes, years, final_table)
    seat_res["114"] = cdmap_seats("114th Congress", fips, epsg, 114, votes, years, final_table)

    final_table.columns = pd.MultiIndex.from_tuples([(usps.upper(), yr) for yr in years])

    return seat_res, final_table


def plot_share(usps, method, year, data):

    if usps not in data: 
        print(usps, "not found.")
    if method not in data[usps]: 
        print(method, "not found.")

    df = data[usps][method]

    
    var = "{} D Fr".format(year)

    outR = df[df[var] < 0.5]
    outD = df[df[var] > 0.5]

    sns.set(rc={"figure.figsize": (6, 3)})
    sns.set_style("white")
    
    ax = sns.distplot(1-outD[var], hist_kws={'weights':outD.weights.values, "alpha" : 0.7, "color" : "#0F83F4"},
                      bins=np.arange(0.05, 0.96, 0.05), kde=False) 


    sns.distplot(1-outR[var], hist_kws={'weights':outR.weights.values, "alpha" : 0.7, "color" : "#E41214"},
                 bins=np.arange(0.05, 0.96, 0.05), kde=False, ax = ax) 

    sns.despine()
    ax.set_xlabel("Republican Vote Share", size = 12)
    ax.set_ylabel("Average Seats / Plan", size = 12)
    ax.set_xlim([0.05, 0.95])
                      
    ax.get_figure().savefig("{}_{}_{}.pdf".format(usps, year, method), bbox_inches = 'tight', pad_inches = 0.1)

