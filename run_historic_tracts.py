#!/usr/bin/env python3

from dist_tools import *

metrics = pd.DataFrame()

for year, sessn in [[1990, 107], [2000, 111], [2010, 114]]:
    for usps in us_states:

        epsg, fips, seats = get_state_info(usps)
        if seats == 1: continue

        print(usps, year, sessn, seats, epsg, fips)

        cd_gdf = get_state_cd_map(usps, sessn, year, epsg)
        tr_gdf = get_state_cells(usps, year, epsg, cell = "tract")
        evaluate_spatial(cd_gdf, tr_gdf)
        metrics = metrics.append(cd_gdf.to_crs(epsg = 2163))
        
metrics = metrics.reset_index()

cols = ['usps', 'congress', 'year', 'cd', 
        'R_isoa', 'R_isop', 'R_LIC', 'R_SCC', 'R_mean', 'R_harm', 'R_dyn',
        'pop', 'pop_hull', 'A_hull', 'IP_d', 
        'obj_ip_dist', 'obj_polsby', 'obj_lic', 'obj_scc', 'obj_inertia_a', 'obj_inertia_p', 
        'obj_mean_radius', 'obj_harm_radius', 'obj_dyn_radius', 'obj_axis', 'obj_exchange', 
        'obj_hull_pop', 'obj_hull_area', 'obj_rohrbach']

# metrics[cols + ['geometry']].to_file("data/decennial_compactness.shp")
# metrics[cols].to_csv("data/decennial_compactness.csv", index = False)

