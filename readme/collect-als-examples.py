"""Retrieve measured 3DEP nodes for the new rGEDI-selected forest footprints."""
from pathlib import Path
import importlib.util,csv,time,os
import numpy as np
root=Path(os.environ.get('RGEDI_FIGURE_DATA','D:/rGEDI_manuscript/Manuscript_v07/forest_examples_v15'))
source=Path(__file__).with_name('als_3dep_helpers.py')
spec=importlib.util.spec_from_file_location('als',source)
als=importlib.util.module_from_spec(spec);spec.loader.exec_module(als)
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(p,keys,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(keys);w.writerows(rows)
for folder,project,landscape in [('sierra_nevada','CA_SierraNevada_6_2022','Sierra Nevada'),('palm_coast','FL_2017PalmCoast_1_C22','Palm Coast')]:
    out=root/folder
    while not (out/'candidates.csv').exists():
        print('Waiting for rGEDI selection:',landscape,flush=True);time.sleep(15)
    candidate=read(out/'candidates.csv')[0]
    shot=candidate['shot_number'];lon=float(candidate['lon_lowestmode']);lat=float(candidate['lat_lowestmode'])
    url='https://s3-us-west-2.amazonaws.com/usgs-lidar-public/'+project
    cloud,original,terrain=als.profile(als.around_shot(url,lon,lat,cache=out/'ept_nodes'))
    print(landscape,'returns',len(original),'ALS height 95th percentile',np.quantile(original[:,3],.95),flush=True)
    if np.quantile(original[:,3],.95)<8:raise ValueError('Selected footprint lacks a forest ALS profile')
    save(out/'points.csv',['shot','landscape','relative_easting_m','relative_northing_m','relative_elevation_m','height_m','class'],[(shot,landscape,*r) for r in original])
    save(out/'terrain.csv',['shot','landscape','relative_easting_m','relative_northing_m','relative_elevation_m'],[(shot,landscape,*r) for r in terrain])
    print('Saved measured ALS data for',landscape,flush=True)
