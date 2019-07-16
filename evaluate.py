#!/usr/bin/env python

from dist_tools import *
from shutil import copyfile as cp


pol_states = ["fl", "il", "la", "md", "mn", "nc", "pa", "tn", "tx", "va", "wi"]
states = ["ca", "ny", "al", "oh"]
states = ["nc", "tx", "fl"]

# files = sorted(glob.glob("cd/*.geojson") + \
#                glob.glob("../chalk/s3/res/*/power/s26*/c*/*geojson") + \
#                glob.glob("../chalk/s3/res/*/split/s001/*geojson") + \
#                sum([glob.glob("../chalk/s3/res/res/{}/*/s26*/c*/*geojson".format(s)) for s in states], []))
# 
# files = glob.glob("../chalk/s3/res/*/path_frac/s26*/c*/*geojson") + \
#         glob.glob("../chalk/s3/res/*/hull_*/s26*/c*/*geojson") + \
#         glob.glob("../chalk/s3/res/*/split/s001/*geojson") 
# 
# files = glob.glob("maps/*split*geojson")

# files = glob.glob("cd/*115*geojson")

# files = glob.glob("../cluscious/res/wi_as/power/s27*/c*/*geojson")
# files = glob.glob("cd/wi_as_2016L.geojson")
files = glob.glob("extras/*geojson")

files.sort()

print("Running", len(files), "maps.")

last_state, map_json = None, {}
cols = ["state", "method", "seed", "cycle", "file",
        "spatial", "dseats", "rseats", "eg", "bseats", "hseats"]

# with open("map_listing.csv", "w") as out: pass
# print(files)

processed = pd.read_csv("map_listing.csv", names = cols)\
              .drop_duplicates("file", keep = "last")["file"].values

for fi, f in enumerate(files): 
    
    fo = f.replace("../chalk/s3/res/", "")
    fo = fo.replace("../cluscious/res/", "")
    fo = fo.replace("res/", "")
    fo = fo.replace("maps/", "")
    fo = fo.replace("extras/", "")
    fo = re.sub(r"cd/([a-z][a-z])_(1[01][1457])", r"\1/cd_\2", fo)
    fo = re.sub(r"cd/wi_as_2016L", r"wi_as/2016L", fo)
    fo = fo.replace("/final", "")

    if "maps/" in f:
        fl = re.sub(r"([a-z][a-z])_([a-z_]+)_(s[0-9]{3})_(c[0-9]{3}).geojson", r"\1+\2+\3+\4", fo)
        fl = re.sub(r"([a-z][a-z])_(cd_1[01][1457]).geojson", r"\1+\2", fl)
        fl = re.sub(r"([a-z][a-z])_split_s001.geojson", r"\1+split+s001+c000", fl)
        file_list = fl.split("+")
    elif "extras" in f:
        fl = re.sub(r"([a-z][a-z])_(cd_1[01][14567]).geojson", r"\1+\2", fo)
        file_list = fl.split("+")
    else:
        file_list = fo.replace(".geojson", "").split("/")

    fo = fo.replace("/", "_")
    fo = "maps/" + fo

    if fo in processed: continue
    if not "maps/" in f: cp(f, fo)

    with open("evaluate.sub", "a") as out: out.write(fo + "\n")
    
    if len(file_list) < 3: file_list.append("s001")
    if len(file_list) < 4: file_list.append("c000")
    print(file_list)
    
    state_tag, method, seed, cycle = file_list
    state = state_tag.split("_")[0]

    if state_tag != last_state:
        epsg, fips, seats = get_state_info(state)
        trdf = get_state_cells(state, 2015, epsg, "tract")
        last_state = state_tag

    
    plan_dict = {"state" : state_tag, "method" : method, "file" : fo,
                 "seed" : int(seed[1:]), "cycle" : int(cycle[1:])}

    # gdf = gpd.read_file(fo).rename(columns = {"id" : "cd"}).set_index("cd").to_crs(epsg = epsg)
    gdf = gpd.read_file(fo).to_crs(epsg = epsg)

    if state_tag != state:
        seats = gdf.shape[0]
    

    if "cd" in method or "201" in method:

        plan_dict["spatial"] = gdf["PCA1"].mean()

    else:

        try:    evaluate_spatial(gdf, trdf)
        except: continue
        plan_dict["spatial"] = gdf["obj_pca1"].mean()


    if state in pol_states:
        plan_dict["eg"] = efficiency_gap_rep(gdf)
        plan_dict["dseats"], plan_dict["rseats"] = party_seats(gdf, seats)
    else:
        plan_dict["eg"] = -1
        plan_dict["dseats"] = -1
        plan_dict["rseats"] = -1

    plan_dict["bseats"], plan_dict["hseats"] = minority_seats(gdf, seats)
    
    with open("map_listing.csv", "a") as out:
        pd.DataFrame([plan_dict])[cols].to_csv(out, float_format='%.4g', header = False, index = False)


df = pd.read_csv("map_listing.csv", names = cols).drop_duplicates("file", keep = "last")\
       .sort_values(by = ["state", "method", "seed", "cycle"])

for s in df.state.unique():
    sdf = df.loc[df.state == s, cols].to_json("{}_dir.json".format(s), orient = "records")



