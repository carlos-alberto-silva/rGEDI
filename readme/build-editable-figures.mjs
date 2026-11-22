// Original editable scientific figures. Run from the rGEDI repository root.
import fs from 'node:fs/promises';
import {readFileSync} from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import {pathToFileURL} from 'node:url';
const modules=process.env.ARTIFACT_NODE_MODULES || path.resolve(path.dirname(process.execPath),'..','node_modules');
process.env.RUNTIME_NODE_MODULES=modules;
const {Presentation,PresentationFile}=await import(pathToFileURL(path.join(modules,'@oai/artifact-tool/dist/artifact_tool.mjs')).href);
const skill=process.env.PRESENTATION_SKILL || path.join(os.homedir(),'.codex/plugins/cache/openai-primary-runtime/presentations/26.1112.11217/skills/presentations');
const {finalizePresentation,applyPresentationChartFont}=await import(pathToFileURL(path.join(skill,'container_tools/artifact_tool_utils.mjs')).href);
const repo=process.cwd();
const workspace=process.env.FIGURE_WORKSPACE || repo;
const build=path.join(workspace,'07_others/fire_revision/editable_build');
const out=path.join(workspace,'03_figures/fire_revision/editable');
const python=process.env.ARTIFACT_PYTHON || path.resolve(modules,'../../python/python.exe');
await fs.mkdir(build,{recursive:true}); await fs.mkdir(out,{recursive:true});
const data=JSON.parse(await fs.readFile(path.join(repo,'readme/data/editable_figure_values.json'),'utf8'));
const font='Arial', ink='#24353E';
function text(s,value,x,y,w,h,size=20,color=ink,bold=false){
  const a=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  a.text=value;a.text.style={fontFamily:font,fontSize:size,color,bold};a.text.insets=0;return a;
}
function rect(s,x,y,w,h,fill,line=fill){return s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill,line:{fill:line,width:1}})}
function line(s,x1,y1,x2,y2,color='#718086',width=1.3,dash='solid'){
  return s.shapes.add({geometry:'custom',position:{left:Math.min(x1,x2),top:Math.min(y1,y2),width:Math.max(1,Math.abs(x2-x1)),height:Math.max(1,Math.abs(y2-y1))},fill:'none',line:{fill:color,width,style:dash},customPaths:[{width:Math.max(1,Math.abs(x2-x1)),height:Math.max(1,Math.abs(y2-y1)),commands:[{moveTo:{x:x1<=x2?0:Math.abs(x2-x1),y:y1<=y2?0:Math.abs(y2-y1)}},{lineTo:{x:x1<=x2?Math.abs(x2-x1):0,y:y1<=y2?Math.abs(y2-y1):0}}]}]});
}
function arrow(s,x,y,w,h,down=false){s.shapes.add({geometry:down?'downArrow':'rightArrow',position:{left:x,top:y,width:w,height:h},fill:'#7D8B90',line:{fill:'none',width:0}})}
function image(s,name,x,y,w,h){return s.images.add({blob:new Uint8Array(readFileSync(path.join(repo,'readme',name))),contentType:'image/png',alt:name,position:{left:x,top:y,width:w,height:h},fit:'contain'})}
function card(s,x,y,w,h,title,body,color){
  rect(s,x,y,w,h,'#FFFFFF',color);rect(s,x,y,w,42,color);
  text(s,title,x+13,y+9,w-26,27,21,'#FFFFFF',true);
  text(s,body,x+13,y+55,w-26,h-63,17.5);
}
async function finish(p,name,w,h,requirements={}){
 name=name.replace('_v12','_v22');
 const candidate=path.join(build,name+'.draft.pptx');
  await (await PresentationFile.exportPptx(p)).save(candidate);
  const result=await finalizePresentation({workspaceDir:workspace,candidatePath:candidate,finalPath:path.join(out,name+'.pptx'),
    pythonExecutable:python,integrityValidatorPath:path.join(skill,'container_tools/inspect_presentation_package_integrity.py'),
    layoutValidatorPath:path.join(skill,'container_tools/inspect_presentation_layout_geometry.py'),
    layoutArgs:['--expected-slide-size-emu',`${w*9525},${h*9525}`,'--validate-heading-fit'],
    fontPolicy:{basis:'design',families:[font]},verifyArtifactToolImport:true,
    receiptPath:path.join(build,name+'.validation.json'),...requirements});
  for(let i=0;i<p.slides.items.length;i++){
    const png=await p.export({slide:p.slides.items[i],format:'png',scale:1});
    await fs.writeFile(path.join(build,`${name}-${i+1}.png`),new Uint8Array(await png.arrayBuffer()));
  }
  console.log('Saved',name,result.finalPath||'');
}

function wave(s,site,x,y,pw,ph,cumulative=true,scale=1){
 const yy=z=>y+(60-z)/70*ph;
 line(s,x,y,x,y+ph,ink);line(s,x,y+ph,x+pw,y+ph,ink);
 [0,20,40,60].forEach(z=>{line(s,x,yy(z),x+pw,yy(z),'#E2E8E8',.7);text(s,`${z} m`,x-45,yy(z)-9,43,22,14*scale)});
 [0,50,100].forEach(v=>text(s,String(v),x+pw*v/100-12,y+ph+6,40,22,15*scale));
 const pts=site.waveform.filter(r=>r.height>=-10&&r.height<=60);
 const commands=pts.map((r,i)=>({[i?'lineTo':'moveTo']:{x:r.amplitude/100*pw,y:(60-r.height)/70*ph}}));
 s.shapes.add({geometry:'custom',position:{left:x,top:y,width:pw,height:ph},fill:'#DDE9ED',line:{fill:'none',width:0},customPaths:[{width:pw,height:ph,commands:[{moveTo:{x:0,y:ph}},...commands.map(o=>({lineTo:Object.values(o)[0]})),{lineTo:{x:0,y:0}},{close:{}}]}]});
 s.shapes.add({geometry:'custom',name:`GEDI waveform ${site.shot}`,position:{left:x,top:y,width:pw,height:ph},fill:'none',line:{fill:'#225C78',width:1.5},customPaths:[{width:pw,height:ph,commands}]});
 if(cumulative){
  let total=0;const energy=[0];for(let j=1;j<pts.length;j++){total+=(pts[j].height-pts[j-1].height)*(pts[j].amplitude+pts[j-1].amplitude)/2;energy.push(total)}
  s.shapes.add({geometry:'custom',name:'Normalized cumulative return energy',position:{left:x,top:y,width:pw,height:ph},fill:'none',line:{fill:'#303C3B',width:1.4},customPaths:[{width:pw,height:ph,commands:pts.map((r,i)=>({[i?'lineTo':'moveTo']:{x:energy[i]/total*pw,y:(60-r.height)/70*ph}}))}]});
  line(s,x,y,x+pw,y,ink);[0,.5,1].forEach(v=>text(s,String(v),x+pw*v-10,y-24,35,22,14*scale));
  text(s,'Cumulative energy',x-15,y-49,pw+40,24,16*scale);
 }
 line(s,x,yy(0),x+pw,yy(0),'#636B65',1.1,'dash');
 line(s,x,yy(+site.rh98),x+pw,yy(+site.rh98),'#A12B74',1.2,'dash');
 text(s,`RH98 ${(+site.rh98).toFixed(1)} m`,x+10,yy(+site.rh98)-23,pw-5,24,16*scale,'#A12B74');
 text(s,'Ground 0 m',x+10,yy(0)+5,pw-5,24,15*scale,'#636B65');
 text(s,site.title||`GEDI ${site.date}`,x-5,y-(cumulative?76:29),pw+40,25,17*scale);
 text(s,'Peak return (%)',x-8,y+ph+31,pw+40,27,17*scale);
}

// The complete function diagram uses native editable shapes and text.
{
 const w=1280,h=1160,p=Presentation.create({slideSize:{width:w,height:h}}),s=p.slides.add();s.background.fill='#FFFFFF';
 const c=['#285F89','#197F78','#6C5092','#AD6630'];
 text(s,'rGEDI functions and workflows',30,15,1220,48,34,ink,true);
 [['Local files','HDF5 and GeoTIFF\nDisk and shared storage'],['NASA Earthdata','Original granules\nAuthenticated cloud access'],['Google Earth Engine','GEDI and ancillary imagery\nCloud processing'],['GFZ gediDB','Analysis-ready footprints\nPublic TileDB service']].forEach((r,i)=>card(s,30+315*i,90,295,125,...r,c[0]));
 text(s,'gedifinder(source = ...) and product readers',280,288,950,30,20,c[0],true);
 line(s,170,280,1110,280);[170,485,800,1110].forEach(x=>line(s,x,216,x,280));arrow(s,211,283,18,35,true);
 const cards=[
 ['1  Read and extract','readLevel1B/2A/2B/3/4A/4B()\ngetLevel1BWF()  getLevel2AM()\ngetLevel2BVPM()  getLevel4A()\ngetLevel2BPAI/PAVDProfile()',c[0]],
 ['2  Filter and clip','Quality, beam and solar-elevation filters\nclip()\nclipLevel*Geometry()\nPreserve shot IDs and spatial support',c[1]],
 ['3  Grid and summarize','rasterize()  polyStats()\ngridStatsLevel2AM/2BVPM/4A()\naccumulateGEDIGranules()\nassembleGEDIArchiveYears()',c[1]],
 ['4  Simulate waveforms','ALS point clouds and plot locations\ngediWFSimulator()\ngediWFMetrics()\nSimulated returns and structural metrics',c[3]],
 ['5  Sample and model','sample()  extractGEDIAncillary()\nextractEE()  varSel()\nfit_model()  fit_metrics()\nTraining and validation',c[2]],
 ['6  Predict and export','predict()  predictGEDIH5()\nmap_create/view/download()\nee_image_to_drive()\nee_check_task_status()',c[2]]];
 cards.forEach((r,i)=>card(s,30+(i%3)*420,325+Math.floor(i/3)*235,390,195,...r));
 [325,560].forEach(y=>{arrow(s,423,y+85,24,17);arrow(s,843,y+85,24,17)});arrow(s,629,524,18,32,true);
 text(s,'Supporting functions',30,790,1220,35,26,ink,true);
 [['Access and configuration','rGEDI_configure()  earthdata_login()\nee_initialize()  gee_available()\ngedidb_available/disconnect()',c[0]],
 ['Spatial and sampling helpers','mosaic/project/resample/aggregateLevel*()\nrandom/spaced/grid/stratifiedSampling()\ngeomSampling()  rasterSampling()',c[1]],
 ['Inspect and visualize','getTrack()  plot_gedi_orbit_animation()\nplotWFMetrics()  plotPAI/PAVDProfile()\nplotLevel3/4A/4B()  close()',c[3]]].forEach((r,i)=>card(s,30+420*i,840,390,160,...r));
 ['Access and extraction','Spatial processing','Modelling and mapping','Simulation and visualization'].forEach((r,i)=>{rect(s,30+i*315,1040,15,15,c[i]);text(s,r,53+i*315,1037,285,32,17)});
 text(s,'Grouped labels abbreviate related functions. See the package reference manual for full names and arguments.',30,1100,1220,35,17);
 s.speakerNotes.text='Original rGEDI workflow diagram. All boxes, labels, arrows and colours are editable PowerPoint shapes.';
 await finish(p,'rGEDI_function_workflow_v12',w,h);
}

// The full six-location plate retains separate clouds and editable waveform paths.
{
 const w=1440,h=2150,p=Presentation.create({slideSize:{width:w,height:h}}),s=p.slides.add();s.background.fill='#FFFFFF';
 image(s,'fig-us-forest-sites.png',0,0,w,610);

 for(let i=0;i<data.sites.length;i++){
   const site=data.sites[i],col=i%2,row=Math.floor(i/2),x=col*720,y=720+row*450;
   image(s,`data/us_forests/cloud_${site.id}.png`,x+5,y-5,465,365);
   wave(s,site,x+490,y+5,202,305);
 }
 data.viridis.forEach((c,i)=>rect(s,430+i*8.5,2020,8.6,16,c));
 [0,10,20,30,40,50,60].forEach(v=>text(s,String(v),425+v/60*544,2044,45,25,18));
 text(s,'ALS return height above local ground (m)',430,2080,650,27,20);
 text(s,'All ALS clips have a 25 m diameter. Waveform height is relative to GEDI ground.',260,2117,1080,26,17);
 s.speakerNotes.text='Measured USGS 3DEP point clouds and archived NASA GEDI L1B returns. Shot IDs, granules and dates are in readme/data/us_forests/sites.csv. Each cloud is a separate high-resolution image. Waveform paths, reference lines, ticks and labels are editable shapes. ALS years differ from GEDI at five sites. Cumulative energy is integrated over the displayed -10 to 60 m interval and does not reproduce the mission L2A processing. Shared axes do not imply simultaneous acquisition.';
 await finish(p,'rGEDI_US_ecosystem_waveforms_v12',w,h);
}

// The waveform inset remains native and independently editable.
{
 const w=1400,h=1310,p=Presentation.create({slideSize:{width:w,height:h}}),s=p.slides.add();s.background.fill='#FFFFFF';
 image(s,'data/profile_map.png',0,0,790,440);image(s,'data/profile_zoom.png',805,0,595,440);
 image(s,'data/profile_profile_base.png',0,445,1400,425);
 rect(s,1110,455,280,380,'#FFFFFF','#91A69C');
 wave(s,data.profileWave,1160,525,185,230,false,.73);
 image(s,'data/profile_vert.png',0,880,690,420);image(s,'data/profile_box.png',705,880,695,420);
 s.speakerNotes.text='NPS Great Smoky Mountains National Park AOI. A uses quality-screened GEE L2A V2 footprints from May through September 2019. B uses geolocation of all eight tracks from one Earthdata L2A V3 granule, orbit 02629 part 03. C uses screened L2B V3 PAVD and matched L2A V3 RH98, with the actual L1B waveform from shot 26290500300273355. The PAVD display cells meet between sample centres. D and E compare mapped NLCD 2019 forest classes from extractEE(). Source data and native-bin summaries are archived. Map and dense profile panels are independent images. The inset waveform, RH98 line, ground line, axes and labels are editable shapes. Source scripts control internal plotting details. https://www.nps.gov/grsm/planyourvisit/access_data.htm';
 await finish(p,'rGEDI_forest_characterization_v12',w,h);
}

// Native XY charts preserve the observations and fitted relationships for editing.
{
 const w=1400,h=1000,p=Presentation.create({slideSize:{width:w,height:h}}),s=p.slides.add();s.background.fill='#FFFFFF';
 text(s,'Forest biomass model training and validation',35,22,1330,55,34,ink,true);
 for(let i=0;i<data.models.length;i++){
  const m=data.models[i],left=20+i*690;
  const chart=s.charts.add('scatter',{position:{left,top:115,width:670,height:735},title:m.stage,
   titleTextStyle:{fontSize:24,typeface:font,bold:true},hasLegend:true,legend:{position:'bottom',textStyle:{fontSize:17,typeface:font}},
   scatterOptions:{style:'lineWithMarkers'},
   series:[{name:'Observations',xValues:m.observed.map(x=>+x.toPrecision(14)),values:m.predicted.map(x=>+x.toPrecision(14)),marker:{symbol:'circle',size:3,fill:'#268A83'},line:{fill:'none',width:0}},
    {name:'1:1',xValues:[0,750],values:[0,750],marker:{symbol:'none'},line:{fill:'#7D8587',width:1.2,style:'dash'}},
    {name:'Linear fit',xValues:[0,750],values:[m.intercept,m.intercept+750*m.slope].map(x=>+x.toPrecision(14)),marker:{symbol:'none'},line:{fill:'#A12B74',width:1.8}}],
   xAxis:{min:0,max:750,majorUnit:150,title:{text:'Observed GEDI AGBD (Mg/ha)',textStyle:{fontSize:20,typeface:font}},textStyle:{fontSize:18,typeface:font},majorGridlines:{fill:'#E8EDED',width:.7}},
   yAxis:{min:0,max:750,majorUnit:150,title:{text:'Predicted AGBD (Mg/ha)',textStyle:{fontSize:20,typeface:font}},textStyle:{fontSize:18,typeface:font},majorGridlines:{fill:'#E8EDED',width:.7}}});
  applyPresentationChartFont(chart,{fontFamily:font});
  text(s,`n = ${m.stats.n}\nR\u00b2 = ${(+m.stats.r2).toFixed(2)}\nRMSE = ${(+m.stats.rmse).toFixed(1)} Mg/ha\nBias = ${(+m.stats.bias).toFixed(1)} Mg/ha`,left+90,205,320,145,21);
 }
 text(s,'Statistics include all observations. One training observation exceeds the displayed range of 0 to 750 Mg/ha.',35,910,1330,50,20);
 s.speakerNotes.text='Native editable scatter charts contain every saved prediction. The final model pools 392 samples from 2019, 499 from 2021 and 486 from 2022. Each validation fold trains on 2019/2021 outside withheld 2.5 km blocks and predicts 2022 inside withheld blocks. GEDI L4A biomass is the reference response. The PPTX uses a common dot colour to keep native chart editing simple. The repository PNG uses density colours.';
 const maps=p.slides.add();maps.background.fill='#FFFFFF';
 image(maps,'fig-forest-landscape.png',0,30,700,900);image(maps,'fig-gedi-wall-to-wall.png',700,30,700,900);
 maps.speakerNotes.text='Independent high-resolution map objects show Landsat 8 natural colour and the 2019 biomass estimate. Both include the 2020 Mangum Fire boundary, north arrow, scale and legend. Raster source and rendering code are saved in readme/plot-forest-maps.R. Landsat and predictor access use rGEDI GEE interfaces.';
 await finish(p,'rGEDI_biomass_model_and_maps_v12',w,h,{requiredNativeChartOwnerSlides:[1],materializeLiteralChartWorkbooks:true});
}
