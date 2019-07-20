#!/usr/bin/env python 

from dist_tools import *

# Calculate all spatial, demographic, and vote scores of
# districts enacted for the 107th, 11th, and 114th Congresses.

columns = {'cd': "ID", "pop_ratio" : "Pop./Target", "a_sq_mi" : 'Area [sq mi]', 
           "obj_dyn_radius": 'DynamicRadius', "obj_polsby" : 'IPQ', 
           "obj_axis" : 'AxisRatio', "obj_exchange" : 'Exchange', "obj_lic" : 'InscrCircle', 
           "obj_harm_radius" : 'HarmonicRadius', "obj_inertia_a" : 'InertiaArea',
           "obj_mean_radius" : 'MeanRadius', "obj_scc" : 'CircCircle', "obj_hull_pop" : 'HullPop',
           "obj_inertia_p" : 'InertiaPop', "obj_ip_dist" : "InterpersonalDistance",
           "obj_rohrbach" : 'DistPerimeter', "obj_hull_area" : 'HullArea',
           "obj_pca1" : "PCA1", "obj_pca2" : "PCA2",
           "pop" : "Population"}


metrics_cols = ['usps', 'congress', 'year', 'cd', 'pop', 'pop_ratio',
		'R_isoa', 'R_isop', 'R_LIC', 'R_SCC', 'R_mean', 'R_harm', 'R_dyn', 'pop_hull', 'A_hull', 'IP_d', 
		'obj_ip_dist', 'obj_polsby', 'obj_lic', 'obj_scc', 'obj_inertia_a', 'obj_inertia_p', 'obj_mean_radius',
		'obj_harm_radius', 'obj_dyn_radius', 'obj_axis', 'obj_exchange', 'obj_hull_pop', 'obj_hull_area', 'obj_rohrbach']

os.makedirs("cd", exist_ok = True)

us_states = ["al", "ak", "az", "ar", "ca", "co", "ct", "de", "dc", "fl",
             "ga", "hi", "id", "il", "in", "ia", "ks", "ky", "la", "me",
             "md", "ma", "mi", "mn", "ms", "mo", "mt", "ne", "nv", "nh",
             "nj", "nm", "ny", "nc", "nd", "oh", "ok", "or", "pa", "ri",
             "sc", "sd", "tn", "tx", "ut", "vt", "va", "wa", "wv", "wi", "wy"]

# first = ["pa", "fl", "wi", "nc", "md", "mn", "il", "wi", "tx", "ca", "va"]

for usps in us_states:

    for year, sessn in [[1990, 107], [2000, 111], [2010, 114]]:
        
        epsg, fips, seats = get_state_info(usps)
        if seats == 1: continue

        # Logging, since this is super slow.
        print(usps, year, sessn, seats, epsg, fips)
        with open("run_historic.sub", "a") as out:
            out.write("{} {} {} {} {} {}\n".format(usps, year, sessn, seats, epsg, fips))

        # CD boundaries; block geometries for 2000, 2010 and block groups for 1990.
        cd_gdf = get_state_cd_map(usps, sessn, year, epsg)
        pt_gdf = get_state_cells(usps, year, epsg, "bg" if (year == 1990) else "blocks")

        try:
            
            # Evaluate the spatial indicies.
            evaluate_spatial(cd_gdf, pt_gdf)

        except MemoryError:
            with open("run_historic.sub", "a") as out:
            	out.write(" >> Recovering from MemoryError -- switching to block groups\n")

            # If that fails due to memory, downshift to block groups.
            cd_gdf = get_state_cd_map(usps, sessn, year, epsg)
            bl_gdf = get_state_cells(usps, year, epsg, "bg")
            evaluate_spatial(cd_gdf, bl_gdf)
        
        # Same thing for demographics and votes.
        evaluate_demographics(cd_gdf, pt_gdf, usps, epsg)
        vote_list = sorted(list(cd_gdf.filter(regex = r"Vote|Share|Frac|Pop").columns))
                
        # Save it.
        output_geojson(cd_gdf.reset_index().rename(columns = columns)\
                             [list(columns.values()) + vote_list + ["geometry"]],
                       "cd/{}_{}.geojson".format(usps, sessn))
        
        with open("metrics.csv", "a") as out:
            
            cd_gdf.reset_index()[metrics_cols].to_csv(out, header = True, index = False)



# Other legislative districts.
for usps, tag, house, leg_year, data_year in []: # [["wi", "wi_as", "L", 2016, 2010]]

    epsg, fips, _ = get_state_info(usps)

    with open("run_historic.sub", "a") as out:
        out.write("{} {} {} {} {}\n".format(usps, leg_year, house, epsg, fips))
    print(usps, house, leg_year)

    cd_gdf = get_state_ld_map(usps, house, leg_year, epsg)
    pt_gdf = get_state_cells(usps, data_year, epsg, "blocks")

    cd_gdf.index.name = "cd"
    cd_gdf.rename(columns = {"house" : "congress"}, inplace = True)

    try:
        evaluate_spatial(cd_gdf, pt_gdf)

    except MemoryError:
        with open("run_historic.sub", "a") as out:
            out.write(" >> Recovering from MemoryError -- switching to block groups\n")

        cd_gdf = get_state_ld_map(usps, house, leg_year, epsg)
        bl_gdf = get_state_cells(usps, data_year, epsg, "bg")
        evaluate_spatial(cd_gdf, bl_gdf)
    
    evaluate_demographics(cd_gdf, pt_gdf, usps, epsg)
    vote_list = sorted(list(cd_gdf.filter(regex = r"Vote|Share|Frac|Pop").columns))

    cd_gdf["usps"] = tag
            
    output_geojson(cd_gdf.reset_index().rename(columns = columns)\
                         [list(columns.values()) + vote_list + ["geometry"]],
                   "cd/{}_{}{}.geojson".format(tag, leg_year, house))
    
    with open("metrics.csv", "a") as out:
        
        cd_gdf.reset_index()[metrics_cols].to_csv(out, header = True, index = False)




