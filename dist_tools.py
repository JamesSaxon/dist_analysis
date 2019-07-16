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
from matplotlib.patches import Rectangle

import psycopg2

import jellyfish as jf

from sklearn.decomposition.pca import PCA

import scipy
from scipy.stats import norm
import statsmodels.api as sm
import statsmodels.discrete.discrete_model as sm_dm

import shapely
from shapely import wkt
from shapely.geometry import Point

from scipy.spatial import Voronoi, voronoi_plot_2d
from scc import *  ## Nayuki minimum bounding circle, instead of miniball.

import re
import time


us_states = ["al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl", 
             "ga", "hi", "id", "il", "in", "ia", "ks", "ky", "la", "me", 
             "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", 
             "nj", "nm", "ny", "nc", "nd", "oh", "ok", "or", "pa", "ri", 
             "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy"]


def jf_jw_match(x, list_strings):
    return sorted(list_strings, key=lambda y: jf.jaro_winkler(x, y), reverse=True)


def get_tr_rn(usps, epsg):

    query = """SELECT 
                  rn, ST_Transform(tr.geom, epsg) geometry
               FROM census_tracts_2015 AS tr
               JOIN (SELECT state, county, tract,
                            row_number() over 
                              (PARTITION BY state ORDER BY county, tract NULLS LAST) - 1 as rn
                     FROM census_tracts_2015) rn ON
                       tr.state  = rn.state  AND
                       tr.county = rn.county AND
                       tr.tract  = rn.tract
               JOIN states AS st ON st.fips = tr.state
               WHERE st.usps = '{}' ORDER BY rn;
               """
    
    con = psycopg2.connect(database = "census", user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)
    
    tr_rn = gpd.GeoDataFrame.from_postgis(query.format(usps), con, 
                                          geom_col = "geometry", crs = from_epsg(epsg))

    return tr_rn


def get_bg_rn(usps, epsg):

    query = """SELECT 
                  rn, ST_Transform(bg.geom, epsg) geometry
               FROM census_bg_2010 AS bg
               JOIN (SELECT state, county, tract, bgroup,
                            row_number() over 
                              (PARTITION BY state ORDER BY county, tract, bgroup NULLS LAST) - 1 as rn
                     FROM census_bg_2010) rn ON
                       bg.state  = rn.state  AND
                       bg.county = rn.county AND
                       bg.tract  = rn.tract AND 
                       bg.bgroup = rn.bgroup
               JOIN states AS st ON st.fips = bg.state
               WHERE st.usps = UPPER('{}') ORDER BY rn;
               """
    
    con = psycopg2.connect(database = "census", user = user, password = passwd,
                           host = "saxon.harris.uchicago.edu", port = 5432)
    
    bg_rn = gpd.GeoDataFrame.from_postgis(query.format(usps), con, 
                                          geom_col = "geometry", crs = from_epsg(epsg))

    return bg_rn



def merge_tract_number(rndf, vdf, var = "rn"):

    votes = gpd.sjoin(vdf.set_geometry(vdf.centroid), rndf, op = "within", how = "left")

    for vri, row in votes[votes.rn.isnull()].iterrows():
        ctr = row.geometry.centroid
        distances = [(xi, pt.distance(ctr)) for xi, pt in enumerate(rndf.centroid)]
        match = min(distances, key=lambda item:item[1])
        votes.loc[vri, var] = match[0]
      
    votes.rn = votes.rn.astype(int)
    agg_votes = votes.groupby(var).sum().filter(regex = '[RD][901][02468]').astype(float)

    return rndf.set_index("rn", drop = True)[[]].join(agg_votes, how = "left").fillna(0)


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



def cdmap_seats(sessn, usps, votes = None):

    if votes is None: votes = usps

    cd_gdf = gpd.GeoDataFrame.from_postgis("""SELECT cd, cd.geom geometry, epsg
                                              FROM   cd 
                                              JOIN   states ON cd.state = states.fips
                                              WHERE  cd.sessn = {} AND states.usps = UPPER('{}')
                                              ORDER  BY cd;
                                              """.format(sessn, usps),
                                              con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                                     host = "saxon.harris.uchicago.edu", port = 5432),
                                              geom_col = "geometry", crs = from_epsg(2163))
    epsg = cd_gdf.iloc[0]["epsg"]

    cd_gdf = cd_gdf.to_crs(epsg = epsg)

    years = {}
    for y in range(1992, 2020, 2):

        votes_file = os.path.expanduser('~/proj/dist_analysis/voting/mapped/{}_{}.geojson'.format(votes.lower(), y))
        if not os.path.isfile(votes_file): continue

        v = gpd.read_file(votes_file).to_crs(epsg = epsg)
        v.set_geometry(v.centroid, inplace = True)

        v = gpd.sjoin(v, cd_gdf, op = "within", how = "left")

        for vri, vrow in v[v.cd.isnull()].iterrows():
            ctr = vrow.geometry.centroid
            distances = [(cd_row.cd, cd_row.geometry.distance(ctr)) for cdi, cd_row in cd_gdf.iterrows()]
            match = min(distances, key=lambda item:item[1])
            v.loc[vri, "cd"] = match[0]

        v.cd = v.cd.astype(int)

        y2 = "{:02d}".format(y % 100)
        pv = ["D" + y2, "R" + y2]

        cd_votes = v.groupby('cd')[pv].sum()

        years[y] = (cd_votes["D" + y2] > cd_votes["R" + y2]).mean()

    return years



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
    
    stdir = "../chalk/s3/res/{}/{}/s2[78][0-9]/*/*.csv"

    metrics = {"Dynamic Radius"        : "dyn_radius",  "Inscribed Circles"      : "ehrenburg",
               "Exchange"              : "exchange",    "Population Hull"        : "hull_p",
               "Isoperimeter Quotient" : "polsby",      "Circumscribing Circles" : "reock", 
               "Rohrbach"              : "rohrbach",    "Power Diagram"          : "power",
               "Split-Line"            : "split"}


    seat_res = {}
    for ti, m in metrics.items():
        # print(ti, sorted(glob.glob(stdir.format(usps, m))))
        seat_res[m] = map_seats(ti, glob.glob(stdir.format(usps, m)), votes, years, final_table)

    seat_res["107"] = cdmap_seats("107th Congress", 107, epsg, fips, votes, years, final_table)
    seat_res["111"] = cdmap_seats("111th Congress", 111, epsg, fips, votes, years, final_table)
    seat_res["114"] = cdmap_seats("114th Congress", 114, epsg, fips, votes, years, final_table)

    final_table.columns = pd.MultiIndex.from_tuples([(usps.upper(), yr) for yr in years])

    return seat_res, final_table


def plot_share(data, state, method, seats, mark_competitive = True, for_table = False):

    weight = seats / len(data[state][method])

    rwins = [v for v in data[state][method] if v >= 0.5]
    rW    = [weight for v in rwins]

    dwins = [v for v in data[state][method] if v <  0.5]
    dW    = [weight for v in dwins]

    size = (7, 3)
    if for_table:
        size = (5, 1)
        if method == "split":
            size = (5, 1.2)

    sns.set_style("white", rc={"figure.figsize": size, "axes.linewidth" : 3})
    sns.set_context("notebook", font_scale = 2.2)
    # sns.set(rc={"figure.figsize": size, "axes.linewidth" : 4})

    f, ax = plt.subplots(1, 1, figsize = size)

    sns.distplot(dwins, kde=False, ax = ax, bins=np.arange(0, 1.01, 0.025),
                 hist_kws={"linewidth" : 0.2, "rwidth" : 0.85, 'weights' : dW,
                                "alpha" : 1, "color" : "#56A8F7"})


    sns.distplot(rwins, kde=False, ax = ax, bins=np.arange(0, 1.01, 0.025),
                 hist_kws={"linewidth" : 2, "rwidth" : 0.85, 'weights' : rW,
                           "alpha" : 1, "color" : "#EC595A"})

    ax.set_xlabel("Republican Vote Share")
    ax.set_ylabel("Average Seats / Plan")
    ax.set_xlim([0, 1])
    ax.set_xticks([0.25, 0.5, 0.75])
    ax.set_xticklabels([0.25, 0.5, 0.75])

    sns.despine(left = for_table)
    if for_table:
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.set_xlim([0, 1])
        ax.set_xlabel("")
        if method != "split": ax.set_xticks([])

    if mark_competitive:
        ax.add_patch(Rectangle((0.475, 0), 0.05, ax.get_ylim()[1],
                               linewidth = 0, hatch = "///" * (1 + for_table),
                               fill = False, color = "grey", alpha = 1, zorder = -10))

    f.savefig("figs/{}_{}.pdf".format(state.lower(), method), bbox_inches = 'tight', pad_inches = 0.02)

    plt.close("all")




def output_geojson(df, fname, var = None):

    df["fill-opacity"] = 0.25
    df["stroke"] = "black"
    df["stroke-width"] = 0.35
    df["stroke-opacity"] = 0.5

    if var is None: colors = plt.get_cmap("nipy_spectral")
    else: colors = plt.get_cmap("RdBu")

    fill = {}
    for ri, row in df.iterrows():

        if var is None:
            if df.shape[0] == 1: c = 1
            else: c = ri/(df.shape[0]-1)
        else:
            if np.isnan(row[var]):
                df.loc[ri, "fill-opacity"] = 0.0
                fill[ri] = "#000000"
                continue

            c = row[var]

        color = [int(v*255) for v in colors(c)][:3]
        fill[ri] = "#{0:02X}{1:02X}{2:02X}".format(*color)


    df["fill"] = pd.Series(fill)

    with open(fname, "w") as out: out.write(df.to_crs(epsg = 4269).to_json())




def ctr(group, w):
    
    p = group[w]
    xctr = (group["x"] * p).sum() / p.sum()
    yctr = (group["y"] * p).sum() / p.sum()
    
    return Point([xctr, yctr])


def avg_interperson_distance(group):

    g = group[group["pop"] > 0]

    prod = scipy.spatial.distance.cdist(g[["x", "y"]].as_matrix(),
                                        g[["x", "y"]].as_matrix(),
                                        metric = 'euclidean')
    prod = np.dot(g["pop"], prod)
    prod = np.dot(prod, g["pop"])
    
    return prod / np.dot(g["pop"][:,np.newaxis], g["pop"][np.newaxis,:]).sum()



def normed_inertia(group, d, w):

    I = (group[d]**2 * group[w]).sum()
    In = group["area"].sum() * group[w].sum() / 2 / math.pi

    return In / I 


def R_circumscribing(group): return make_circle([xy for xy in zip(group["x"], group["y"])])[2]

def R_mean(group): return (group["area"] * group["dctr"]).sum() / group["area"].sum()
def R_dyn(group):  return np.sqrt((group["area"] * group["dctr"]**2).sum() / group["area"].sum())
def R_harm(group): return group["area"].sum() / (group["area"] / group["dctr"]).sum()

def get_points(x):
    
    l = x if type(x) is shapely.geometry.multipolygon.MultiPolygon else [x]
    
    pts = []
    for poly in l:
        for pt in poly.exterior.coords:
            pts.append(pt)
        for ir in poly.interiors:
            for pt in ir.coords:
                pts.append(pt)

    return pts


def get_lic(poly):
    
    # simp = poly.simplify(100)
    pts = get_points(poly)

    # print(pts)
    vor = Voronoi(pts)

    max_d2, lic_ctr, lic_r = 0, 0, 0
    for vtxi, vtx in enumerate(vor.vertices):

        if not poly.contains(Point(vtx)):
            continue

        region = [ri for ri, r in enumerate(vor.regions) if vtxi in r][0]

        source = None
        for pti, pt_reg in enumerate(vor.point_region):
            if region == pt_reg: 
                source = vor.points[pti]
                break

        d2 = (vtx[0] - source[0])**2 + (vtx[1] - source[1])**2

        if d2 > max_d2:

            max_d2 = d2
            lic_ctr, lic_r = tuple(vtx), math.sqrt(d2)

    return lic_ctr, lic_r


def get_lic_circ(poly):

    ctr, R = get_lic(poly)
    return Point(ctr).buffer(R)
    
def R_inscribed(poly): return get_lic(poly)[1]

def get_exchange_area(poly):
    
    R = math.sqrt(poly.area/math.pi)
    return poly.intersection(Point(poly.centroid).buffer(R))

def exchange(poly):
    
    R = math.sqrt(poly.area/math.pi)
    return poly.intersection(Point(poly.centroid).buffer(R)).area / poly.area

def axis_ratio(group):

    ax_pca = PCA()
    ax_pca.fit(group[["x", "y"]])
    
    return ax_pca.explained_variance_[1]/ax_pca.explained_variance_[0]


def rohrbach(group):
    
    obj = (group["dperim"] * group["area"]).sum()
    R3 = math.pow(group["area"].sum()/math.pi, 3/2)
                  
    return obj / (math.pi * R3/3)

def get_state_info(usps):
    
    st = pd.read_sql("select seats, epsg, lower(usps) usps, fips from states where usps = upper('{}');".format(usps),
                     con = psycopg2.connect(database = "census", user = user, password = passwd,
                                            host = "saxon.harris.uchicago.edu", port = 5432)).loc[0].to_dict()

    return st["epsg"], st["fips"], st["seats"]


def get_state_cd_map(usps, session, year, epsg):

    cd_gdf = gpd.GeoDataFrame.from_postgis("""SELECT cd, ST_Transform(cd.geom, epsg) geometry
                                              FROM cd JOIN states ON cd.state = fips
                                              WHERE states.usps = UPPER('{}') and sessn = {};""".format(usps, session),
                                           con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                                  host = "saxon.harris.uchicago.edu", port = 5432),
                                           geom_col = "geometry", crs = from_epsg(epsg), index_col = "cd")
    
    cd_gdf["usps"] = usps.upper()
    cd_gdf["congress"] = session
    cd_gdf["year"] = year
    
    return cd_gdf


def get_state_ld_map(usps, house, year, epsg):

    cd_gdf = gpd.GeoDataFrame.from_postgis("""SELECT id, ST_Transform(sld.geomland, epsg) geometry
                                              FROM sld JOIN states ON sld.state = fips
                                              WHERE states.usps = UPPER('{}') AND house = UPPER('{}') AND year = {};""".format(usps, house, year),
                                           con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                                  host = "saxon.harris.uchicago.edu", port = 5432),
                                           geom_col = "geometry", crs = from_epsg(epsg), index_col = "id")
    
    cd_gdf["usps"] = usps.upper()
    cd_gdf["congress"] = house
    cd_gdf["year"] = year
    
    return cd_gdf


def get_state_cells(usps, year = 2015, epsg = 2163, cell = "tract"):

    if cell not in ["tract", "tracts", "bg", "bgroup", "bgroups", "blocks"]:
        print(cell, "not a valid cell geometry!")
    
    rn = "county, tract"
    rn_join = ""
    if "tr" in cell: cell = "tracts"
    if "bg" in cell:
        cell = "bg"
        rn += ", bgroup"
        rn_join = "AND rn.bgroup = cell.bgroup"
    if "bl" in cell:
        cell = "blocks"
        rn += ", block"
        rn_join = "AND rn.block = cell.block"

    query = """SELECT 
                 row_number() over (PARTITION BY state ORDER BY {} NULLS LAST) - 1 AS rn,
                 pop, area, vap, bvap, hvap,
                 ST_Centroid(ST_Transform(cell.geom, epsg)) geometry,
                 ST_X(ST_Centroid(ST_Transform(cell.geom, epsg))) x, 
                 ST_Y(ST_Centroid(ST_Transform(cell.geom, epsg))) y
               FROM census_{}_{} cell, states
               WHERE 
                 cell.state = fips AND states.usps = UPPER('{}');"""

    # print(query.format(rn, cell, year, usps))

    pt_gdf = gpd.GeoDataFrame.from_postgis(query.format(rn, cell, year, usps),
                                           con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                                  host = "saxon.harris.uchicago.edu", port = 5432),
                                           geom_col = "geometry", index_col = "rn", crs = from_epsg(epsg))

    return pt_gdf




def evaluate_spatial(cd_gdf, tr_gdf):

    tr_gdf = gpd.tools.sjoin(tr_gdf, cd_gdf.copy().reset_index(), op = "within")
    tr_gdf = tr_gdf.reset_index()[["cd", "pop", "area", "x", "y", "geometry"]] 

    tr_gdf = tr_gdf.merge(pd.DataFrame({"ctr"  : tr_gdf.groupby("cd").apply(ctr, w = "area"),
                                        "pctr" : tr_gdf.groupby("cd").apply(ctr, w = "pop")}).reset_index(), on = "cd")

    tr_gdf["dctr"] = tr_gdf.distance(tr_gdf.set_geometry("ctr"))
    tr_gdf["dpctr"] = tr_gdf.distance(tr_gdf.set_geometry("pctr"))

    boundaries = shapely.ops.unary_union(cd_gdf.boundary)
    tr_gdf["dperim"] = tr_gdf.distance(boundaries)

    tr_cd_group = tr_gdf.groupby("cd")

    cd_gdf["pop_sp"] = tr_cd_group["pop"].sum()

    cd_gdf["ctr"]    = gpd.GeoSeries(tr_cd_group.apply(ctr, w = "area"))
    cd_gdf["pctr"]   = gpd.GeoSeries(tr_cd_group.apply(ctr, w = "pop"))

    cd_gdf["R_isoa"] = np.sqrt(cd_gdf.area / math.pi)
    cd_gdf["R_isop"] = cd_gdf.length / (2 * math.pi)
    cd_gdf["R_LIC"]  = cd_gdf.geometry.apply(R_inscribed)
    cd_gdf["R_SCC"]  = tr_cd_group.apply(R_circumscribing)
    cd_gdf["R_mean"] = tr_cd_group.apply(R_mean)
    cd_gdf["R_harm"] = tr_cd_group.apply(R_harm)
    cd_gdf["R_dyn"]  = tr_cd_group.apply(R_dyn)

    # print(cd_gdf[["geometry"]].set_geometry(cd_gdf.convex_hull).reset_index().head())
    hull_pop = gpd.tools.sjoin(cd_gdf[["geometry"]].copy().set_geometry(cd_gdf.convex_hull).reset_index(), 
                               tr_gdf[["geometry", "pop"]], op = "contains")[["cd", "pop"]]
    
    cd_gdf["pop_hull"]   = hull_pop.groupby("cd").sum()

    state_mp = shapely.ops.unary_union(list(cd_gdf.geometry))
    cd_gdf["A_hull"]     = cd_gdf.convex_hull.intersection(state_mp).area

    cd_gdf["IP_d"]       = tr_cd_group[["x", "y", "pop"]].apply(avg_interperson_distance)

    cd_gdf["obj_ip_dist"]     = (128 * cd_gdf.R_isoa / (45 * math.pi)) / cd_gdf.IP_d
    cd_gdf["obj_polsby"]      = 4 * math.pi * cd_gdf.area / cd_gdf.length**2
    cd_gdf["obj_lic"]         = (math.pi * cd_gdf["R_LIC"] ** 2) / cd_gdf.area
    cd_gdf["obj_scc"]         = cd_gdf.area / (math.pi * cd_gdf["R_SCC"] ** 2)
    cd_gdf["obj_inertia_a"]   = tr_cd_group.apply(normed_inertia, w = "area", d = "dctr")
    cd_gdf["obj_inertia_p"]   = tr_cd_group.apply(normed_inertia, w = "pop", d = "dpctr")
    cd_gdf["obj_mean_radius"] = (2 * cd_gdf.R_isoa / 3) / cd_gdf.R_mean
    cd_gdf["obj_harm_radius"] = (cd_gdf.R_isoa / 2) / cd_gdf.R_harm
    cd_gdf["obj_dyn_radius"]  = (cd_gdf.R_isoa / math.sqrt(2)) / cd_gdf.R_dyn
    cd_gdf["obj_axis"]        = tr_gdf.groupby("cd").apply(axis_ratio)
    cd_gdf["obj_exchange"]    = cd_gdf.geometry.apply(exchange)
    cd_gdf["obj_hull_pop"]    = cd_gdf["pop_sp"] / cd_gdf["pop_hull"]
    cd_gdf["obj_hull_area"]   = cd_gdf.area / cd_gdf["A_hull"]
    cd_gdf["obj_rohrbach"]    = tr_cd_group.apply(rohrbach)

    cd_gdf['obj_pca1']        = pca.transform(cd_gdf[pca_cols])[:,0]
    cd_gdf['obj_pca2']        = pca.transform(cd_gdf[pca_cols])[:,1]


    cd_gdf["a_sq_mi"]         = cd_gdf.area * 3.8610216e-7


def evaluate_demographics(cd_gdf, pt_gdf, usps, epsg):

    race = gpd.sjoin(pt_gdf, cd_gdf.reset_index(), op = "within").groupby("cd").sum()

    cd_gdf["pop"] = race["pop"]
    cd_gdf["pop_ratio"] = cd_gdf["pop"] * cd_gdf.shape[0] / cd_gdf["pop"].sum()

    cd_gdf["Black VAP Frac"]    = race["bvap"] / race["vap"]
    cd_gdf["Hispanic VAP Frac"] = race["hvap"] / race["vap"]

    for y in range(1992, 2020, 4):

        votes_file = os.path.expanduser('~/proj/dist_analysis/voting/mapped/{}_{}.geojson'.format(usps.lower(), y))
        if not os.path.isfile(votes_file): continue
            
        v = gpd.read_file(votes_file).to_crs(epsg = epsg)
        v.set_geometry(v.centroid, inplace = True)

        v = gpd.sjoin(v, cd_gdf.reset_index(), op = "within", how = "left")

        for vri, vrow in v[v.cd.isnull()].iterrows():
            ctr = vrow.geometry.centroid
            distances = [(cd_row.cd, cd_row.geometry.distance(ctr))
                         for cdi, cd_row in cd_gdf.reset_index().iterrows()]
            match = min(distances, key=lambda item:item[1])
            v.loc[vri, "cd"] = match[0]

        v.cd = v.cd.astype(int)

        y2 = "{:02d}".format(y % 100)
        pv = ["D" + y2, "R" + y2]

        cd_votes = v.groupby('cd')[pv].sum()
        
        cd_gdf["D" + y2 + " Votes"] = cd_votes["D" + y2] 
        cd_gdf["R" + y2 + " Votes"] = cd_votes["R" + y2] 
        cd_gdf["D" + y2 + " Share"] = cd_votes["D" + y2] / (cd_votes["D" + y2] + cd_votes["R" + y2])
        cd_gdf["R" + y2 + " Share"] = cd_votes["R" + y2] / (cd_votes["D" + y2] + cd_votes["R" + y2])


def efficiency_gap_rep(cd_gdf):
    
    years = [x[1:] for x in list(cd_gdf.filter(regex = "^D[901][02468] Votes").columns)]
    if not len(years): return 0
    
    EG = {}
    for vy in years:

        cd_gdf["D" + vy + " wasted"] = np.where(cd_gdf["D" + vy] < cd_gdf["R" + vy], cd_gdf["D" + vy],
                                                cd_gdf["D" + vy] - 0.5 * (cd_gdf["R" + vy] + cd_gdf["D" + vy]))
        cd_gdf["R" + vy + " wasted"] = np.where(cd_gdf["R" + vy] < cd_gdf["D" + vy], cd_gdf["D" + vy],
                                                cd_gdf["R" + vy] - 0.5 * (cd_gdf["D" + vy] + cd_gdf["D" + vy]))

        egy = cd_gdf.filter(regex = vy + " wasted").sum()/(cd_gdf["D" + vy].sum() + cd_gdf["R" + vy].sum())
        EG[vy] = egy["R" + vy + " wasted"] - egy["D" + vy + " wasted"]
    
    return sum(EG.values()) / len(years)


def party_seats(cd_gdf, seats):
    
    years = [x[1:] for x in list(cd_gdf.filter(regex = "^D[901][02468] Votes$").columns)]
    if not len(years): return 0, 0

    if seats: seat_scale = seats / cd_gdf.shape[0]
    else: seat_scale = 1
    
    exp_seats = {vy : seat_scale * (cd_gdf["D" + vy] > cd_gdf["R" + vy]).sum() for vy in years}
        
    dem_seats = sum(exp_seats.values()) / len(exp_seats)

    return dem_seats, seats - dem_seats



min_rep = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + "/min_rep.csv")
min_rep = min_rep[min_rep.Session == 115].copy(deep = True)
min_rep.reset_index(drop = True, inplace = True)

probit_b = sm_dm.Probit(min_rep["BRep"], min_rep[["BFrac", "const"]]).fit(disp=0)
probit_h = sm_dm.Probit(min_rep["HRep"], min_rep[["HFrac", "const"]]).fit(disp=0)

def minority_seats(cd_gdf, seats = 0):
    
    cd_gdf["BSeats"] = norm.sf(- (cd_gdf["Black VAP Frac"]    * probit_b.params["BFrac"] + probit_b.params.const))
    cd_gdf["HSeats"] = norm.sf(- (cd_gdf["Hispanic VAP Frac"] * probit_h.params["HFrac"] + probit_h.params.const))

    if seats: seat_scale = seats / cd_gdf.shape[0]
    else: seat_scale = 1

    return cd_gdf["BSeats"].sum() * seat_scale, cd_gdf["HSeats"].sum() * seat_scale
    

dec = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + "/decennial_compactness.csv")

pca = PCA(n_components = 2)
pca.fit(dec.filter(regex = "^obj"))
pca_cols = list(dec.filter(regex = "^obj").columns)


