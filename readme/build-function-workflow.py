"""Original vector diagram of the rGEDI public workflow and supporting functions."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
fig, ax = plt.subplots(figsize=(11, 9))
ax.set(xlim=(0, 12), ylim=(0, 10.3)); ax.axis('off')
fig.subplots_adjust(left=.015, right=.985, top=.99, bottom=.015)
colors = ['#285F89', '#197F78', '#6C5092', '#AD6630']
def card(x, y, w, h, title, text, color, size=9.5):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.025,rounding_size=0.08',
                              facecolor='white',edgecolor=color,linewidth=1.1))
    ax.add_patch(FancyBboxPatch((x,y+h-.43),w,.43,boxstyle='round,pad=0.025,rounding_size=0.06',
                              facecolor=color,edgecolor=color))
    ax.text(x+.12,y+h-.20,title,color='white',fontsize=11,weight='bold',va='center')
    ax.text(x+.13,y+(h-.43)/2,text,fontsize=size,va='center',linespacing=1.5,color='#24353E')
def arrow(start,end,style='-'):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=12,
                               linewidth=1.15,color='#78868C',linestyle=style))
for i,(title,text) in enumerate([
    ('Local files','HDF5 and GeoTIFF\nDisk and shared storage'),
    ('NASA Earthdata','Original GEDI granules\nAuthenticated cloud access'),
    ('Google Earth Engine','GEDI and ancillary imagery\nCloud processing'),
    ('GFZ gediDB','Analysis-ready footprints\nPublic TileDB service')]):
    card(.15+3*i,8.85,2.7,1.15,title,text,colors[0],9.2)
ax.text(6,8.53,'gedifinder(source = ...)  Ã¢â‚¬Â¢  product readers',ha='center',fontsize=9,color=colors[0],weight='bold')
for x in (1.5,4.5,7.5,10.5):
    ax.plot([x,x],[8.83,8.25],color='#78868C',lw=1)
ax.plot([1.5,10.5],[8.25,8.25],color='#78868C',lw=1)
arrow((2.05,8.25),(2.05,7.99))
card(.15,6.30,3.7,1.65,'1  Read and extract',
     'readLevel1B/2A/2B/3/4A/4B()\ngetLevel1BWF()  getLevel2AM()\ngetLevel2BVPM()  getLevel4A()\ngetLevel2BPAI/PAVDProfile()',colors[0],9.5)
card(4.15,6.30,3.7,1.65,'2  Filter and clip',
     'Quality, beam and solar-elevation filters\nclip()\nclipLevel*Geometry()\nPreserve shot identifiers and spatial support',colors[1],9.5)
card(8.15,6.30,3.7,1.65,'3  Grid and summarize',
     'rasterize()  polyStats()\ngridStatsLevel2AM/2BVPM/4A()\naccumulateGEDIGranules()\nassembleGEDIArchiveYears()',colors[1],9.5)
arrow((3.88,7.13),(4.12,7.13)); arrow((7.88,7.13),(8.12,7.13))
card(.15,3.95,3.7,1.65,'4  Simulate waveforms',
     'ALS point clouds and plot locations\ngediWFSimulator()\ngediWFMetrics()\nSimulated waveforms and structural metrics',colors[3],9.5)
card(4.15,3.95,3.7,1.65,'5  Sample and model',
     'sample()  extractGEDIAncillary()\nextractEE()  varSel()\nfit_model()  fit_metrics()\nTraining and validation',colors[2],9.5)
card(8.15,3.95,3.7,1.65,'6  Predict and export',
     'predict()  predictGEDIH5()\nmap_create()  map_view()  map_download()\nee_image_to_drive()\nee_check_task_status()',colors[2],9.5)
arrow((6,6.27),(6,5.63)); arrow((3.88,4.78),(4.12,4.78)); arrow((7.88,4.78),(8.12,4.78))
ax.text(6,3.56,'Supporting functions connect each stage',ha='center',weight='bold',fontsize=10,color='#24353E')
card(.15,1.58,3.7,1.65,'Access and configuration',
     'rGEDI_configure()  earthdata_login()\nee_initialize()  gee_available()\ngedidb_available()\ngedidb_disconnect()',colors[0],9.5)
card(4.15,1.58,3.7,1.65,'Spatial and sampling helpers',
     'mosaic/project/resample/aggregateLevel*()\nrandom/spaced/grid/stratifiedSampling()\ngeomSampling()  rasterSampling()\nget_catalog_id()  search_datasets()',colors[1],9.3)
card(8.15,1.58,3.7,1.65,'Inspect and visualize',
     'getTrack()  plot_gedi_orbit_animation()\nplotWFMetrics()  plotPAI/PAVDProfile()\nplotLevel3/4A/4B()\nclose()',colors[3],9.5)
for i,(label,col) in enumerate(zip(['Access and extraction','Spatial processing','Modelling and mapping','Simulation and visualization'],colors)):
    ax.scatter(.3+i*3,.97,s=60,color=col,marker='s')
    ax.text(.47+i*3,.97,label,fontsize=9,va='center',color='#24353E')
ax.text(.15,.38,'Grouped labels abbreviate related product-specific functions. The complete callable inventory is available in the package reference manual.',fontsize=8,color='#59676E')
for ext in ('png','pdf'):
    fig.savefig(ROOT / f'readme/fig-package-workflow.{ext}',dpi=400,facecolor='white')
