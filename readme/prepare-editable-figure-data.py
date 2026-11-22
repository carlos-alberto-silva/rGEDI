from pathlib import Path
import csv,json,sys,os
import numpy as np
import matplotlib.pyplot as plt
repo=Path.cwd()
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
sites=read(repo/'readme/data/us_forests/sites.csv')
bins=read(repo/'readme/data/us_forests/bins.csv')
for s in sites:
    s['waveform']=sorted([{'height':float(r['height_m']),'amplitude':float(r['amplitude_pct'])} for r in bins if r['shot']==s['shot']],key=lambda r:r['height'])
model=Path(sys.argv[1] if len(sys.argv)>1 else os.environ.get('RGEDI_FOREST_OUTPUT','D:/rGEDI_manuscript/Manuscript_v07/forest_model'))
predictions=read(model/'predictions.csv')
stats=read(model/'validation.csv')
models=[]
for st in stats:
    rows=[r for r in predictions if r['stage']==st['stage']]
    x=np.array([float(r['observed']) for r in rows]); y=np.array([float(r['predicted']) for r in rows])
    coef=np.polyfit(x,y,1)
    models.append({'stage':st['stage'],'observed':x.tolist(),'predicted':y.tolist(),'stats':st,'slope':float(coef[0]),'intercept':float(coef[1])})
colors=['#%02x%02x%02x'%tuple(int(255*v) for v in plt.cm.viridis(i/63)[:3]) for i in range(64)]
profile=read(repo/'readme/data/appalachian_inset_waveform.csv')
profile_wave={'shot':profile[0]['shot_number'],'date':'2019-05-31','rh98':float(profile[0]['rh98']),'title':'Observed L1B waveform','waveform':sorted([{'height':float(r['height']),'amplitude':float(r['amplitude'])} for r in profile],key=lambda r:r['height'])}
(repo/'readme/data/editable_figure_values.json').write_text(json.dumps({'sites':sites,'models':models,'viridis':colors,'profileWave':profile_wave}),encoding='utf-8')
print('Saved exact waveform and model values for editable PowerPoint graphics.')
