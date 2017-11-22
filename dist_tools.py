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


us_states = ["al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl", 
             "ga", "hi", "id", "il", "in", "ia", "ks", "ky", "la", "me", 
             "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh", 
             "nj", "nm", "ny", "nc", "nd", "oh", "ok", "or", "pa", "ri", 
             "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy"]


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


def cdmap_seats(name, sessn, epsg, fips, tract_votes, years, final_table = None):

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
    
    if final_table is not None: final_table.loc[name] = cd_seats.mean()

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
    
    stdir = "../chalk/s3/res/{}/{}/*/*/*.csv"

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
                      bins=np.arange(0.00, 1.01, 0.025), kde=False) 


    sns.distplot(1-outR[var], hist_kws={'weights':outR.weights.values, "alpha" : 0.7, "color" : "#E41214"},
                 bins=np.arange(0.00, 1.01, 0.025), kde=False, ax = ax) 

    sns.despine()
    ax.set_xlabel("Republican Vote Share", size = 12)
    ax.set_ylabel("Average Seats / Plan", size = 12)
    ax.set_xlim([0.05, 0.95])
                      
    ax.get_figure().savefig("{}_{}_{}.pdf".format(usps, year, method), bbox_inches = 'tight', pad_inches = 0.1)



def output_geojson(df, fname):

    spectral = plt.get_cmap("nipy_spectral")

    fill = {}
    for ri, row in df.iterrows():
        c = ri/(df.shape[0]-1)
        color = [int(v*255) for v in spectral(c)][:3]
        fill[ri] = "#{0:02X}{1:02X}{2:02X}".format(*color)

    df["fill"] = pd.Series(fill)
    df["fill-opacity"] = 0.25
    df["stroke"] = "black"
    df["stroke-width"] = 0.35
    df["stroke-opacity"] = 0.5

    with open(fname, "w") as out: out.write(df.to_crs(epsg = 4269).to_json())




def ctr(group, w):
    
    p = group[w]
    xctr = (group["x"] * p).sum() / p.sum()
    yctr = (group["y"] * p).sum() / p.sum()
    
    return Point([xctr, yctr])


def avg_interperson_distance(group):

    dist2 =  scipy.spatial.distance.cdist(group[["x", "y"]].as_matrix(),
                                          group[["x", "y"]].as_matrix(),
                                          metric = 'euclidean')

    interp_dist = np.dot(np.dot(group["pop"], dist2), group["pop"]) / \
                  np.dot(group["pop"][:,np.newaxis], group["pop"][np.newaxis,:]).sum()

    return interp_dist # * 0.000621371 # meters to miles


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

    pca = PCA()
    pca.fit(group[["x", "y"]])
    
    return pca.explained_variance_[1]/pca.explained_variance_[0]


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


def get_state_tracts(usps, year = 2015, epsg = 2163):

    tr_gdf = gpd.GeoDataFrame.from_postgis("""SELECT rn, pop, 
                                                     ST_Area(ST_Transform(tr.geom, epsg)) area,
                                                     ST_Centroid(ST_Transform(tr.geom, epsg)) geometry,
                                                     ST_AsText(ST_Transform(tr.geom, epsg)) shape,
                                                     ST_X(ST_Centroid(ST_Transform(tr.geom, epsg))) x, 
                                                     ST_Y(ST_Centroid(ST_Transform(tr.geom, epsg))) y
                                              FROM census_tracts_{} tr, states,
                                                   (SELECT
                                                      state, county, tract,
                                                      row_number() over (PARTITION BY state ORDER BY county, tract NULLS LAST) - 1 AS rn
                                                    FROM census_tracts_2015) rn
                                              WHERE 
                                                rn.state = tr.state AND rn.county = tr.county AND rn.tract = tr.tract AND
                                                tr.state = fips AND states.usps = UPPER('{}');""".format(year, usps),
                                           con = psycopg2.connect(database = "census", user = user, password = passwd,
                                                                  host = "saxon.harris.uchicago.edu", port = 5432),
                                           geom_col = "geometry", index_col = "rn", crs = from_epsg(epsg))

    tr_gdf["shape"] = gpd.GeoSeries(tr_gdf["shape"].apply(wkt.loads))

    return tr_gdf


def evaluate_state(cd_gdf, tr_gdf):
    
    tr_gdf = gpd.tools.sjoin(tr_gdf, cd_gdf.reset_index(), op = "within")
    tr_gdf = tr_gdf.reset_index()[["cd", "pop", "area", "x", "y", "geometry", "shape"]]

    tr_gdf = tr_gdf.merge(pd.DataFrame({"ctr"  : tr_gdf.groupby("cd").apply(ctr, w = "area"),
                                        "pctr" : tr_gdf.groupby("cd").apply(ctr, w = "pop")}).reset_index(), on = "cd")

    tr_gdf["dctr"] = tr_gdf.distance(tr_gdf.set_geometry("ctr"))
    tr_gdf["dpctr"] = tr_gdf.distance(tr_gdf.set_geometry("pctr"))

    boundaries = shapely.ops.unary_union(cd_gdf.boundary)
    tr_gdf["dperim"] = tr_gdf.distance(boundaries)

    tr_cd_group = tr_gdf.groupby("cd")

    cd_gdf["ctr"]        = gpd.GeoSeries(tr_cd_group.apply(ctr, w = "area"))
    cd_gdf["pctr"]       = gpd.GeoSeries(tr_cd_group.apply(ctr, w = "pop"))

    cd_gdf["R_isoa"]     = np.sqrt(cd_gdf.area / math.pi)
    cd_gdf["R_isop"]     = cd_gdf.length / (2 * math.pi)
    cd_gdf["R_LIC"]      = cd_gdf.geometry.apply(R_inscribed)
    cd_gdf["R_SCC"]      = tr_cd_group.apply(R_circumscribing)
    cd_gdf["R_mean"]     = tr_cd_group.apply(R_mean)
    cd_gdf["R_harm"]     = tr_cd_group.apply(R_harm)
    cd_gdf["R_dyn"]      = tr_cd_group.apply(R_dyn)

    dist_pop = gpd.tools.sjoin(cd_gdf[["geometry"]].reset_index(),
                               tr_gdf[["geometry", "pop"]], op = "contains")[["cd", "pop"]]

    hull_pop = gpd.tools.sjoin(cd_gdf.set_geometry(cd_gdf.convex_hull).reset_index(), 
                               tr_gdf[["geometry", "pop"]], op = "contains")[["cd", "pop"]]
    

    cd_gdf["pop"] = dist_pop.groupby("cd").sum()
    cd_gdf["pop_hull"]   = hull_pop.groupby("cd").sum()

    state_mp = shapely.ops.unary_union(list(cd_gdf.geometry))
    cd_gdf["A_hull"]     = cd_gdf.convex_hull.intersection(state_mp).area

    cd_gdf["IP_d"]       = tr_cd_group.apply(avg_interperson_distance)

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
    cd_gdf["obj_hull_pop"]    = cd_gdf["pop"] / cd_gdf["pop_hull"]
    cd_gdf["obj_hull_area"]   = cd_gdf.area / cd_gdf["A_hull"]
    cd_gdf["obj_rohrbach"]    = tr_cd_group.apply(rohrbach)

    cd_gdf["pop_ratio"]       = cd_gdf["pop"] * cd_gdf.shape[0] / cd_gdf["pop"].sum()

    cd_gdf["a_sq_mi"]         = cd_gdf.area * 3.8610216e-7



def efficiency_gap_rep(cd_gdf, tr_gdf):
    
    tr = gpd.tools.sjoin(tr_gdf, cd_gdf[["geometry"]].reset_index(), op = "within")
    cd = tr.filter(regex = "cd|^[DR][901][02468]$").groupby("cd").sum()
    
    years = [x[1:] for x in list(cd.filter(regex = "^D").columns)]
    
    EG = {}
    for y in years:
        cd["D" + y + "_waste"] = np.where(cd["D" + y] > cd["R" + y], cd["D" + y] - 0.5 * (cd["R" + y] + cd["D" + y]), cd["D" + y])
        cd["R" + y + "_waste"] = np.where(cd["R" + y] > cd["D" + y], cd["R" + y] - 0.5 * (cd["R" + y] + cd["D" + y]), cd["R" + y])
        egy = cd.filter(regex = y + "_waste").sum()/(cd["D" + y].sum() + cd["R" + y ].sum())
        EG[int(y)] = egy["R" + y + "_waste"] - egy["D" + y + "_waste"]
    
    return sum(EG.values()) / len(years)


min_rep = pd.read_csv("min_rep.csv")
min_rep = min_rep[min_rep.Session == 115].copy(deep = True)
min_rep.reset_index(drop = True, inplace = True)

probit_b = sm_dm.Probit(min_rep["BRep"], min_rep[["BFrac", "const"]]).fit()
probit_h = sm_dm.Probit(min_rep["HRep"], min_rep[["HFrac", "const"]]).fit()

def minority_seats(cd_gdf, tr_gdf):
    
    tr = gpd.tools.sjoin(tr_gdf, cd_gdf[["geometry"]].reset_index(), op = "within")
    cd = tr.filter(regex = "cd|vap").groupby("cd").sum()
    
    cd["BFrac"]  = cd["black_vap"]    / cd["total_vap"]
    cd["BSeats"] = norm.sf(- (cd["BFrac"] * probit_b.params["BFrac"] + probit_b.params.const))

    cd["HFrac"]  = cd["hispanic_vap"] / cd["total_vap"]
    cd["HSeats"] = norm.sf(- (cd["HFrac"] * probit_h.params["HFrac"] + probit_h.params.const))

    # print(cd["BSeats"].sum(), cd["HSeats"].sum())

    return cd["BSeats"].sum() + cd["HSeats"].sum()
    

dec = pd.read_csv("decennial_compactness.csv")

pca = PCA(n_components = 2)
pca.fit(dec.filter(regex = "^obj"))
cols = list(dec.filter(regex = "^obj").columns)


