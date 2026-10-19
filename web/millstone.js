import {instrumentSpec} from './instrument-plot.js';
import init,{acf_diagnostic} from './pkg/guisdap_web.js';
import {fieldPanel,curve,iqr} from './scientific-plots.js';
import {EISCAT_LIMITS} from './plot-style.js';
import {mathAsset} from './math-labels.js';
const el=id=>document.getElementById(id),query=new URLSearchParams(location.search),antenna=query.get('antenna')||'ZENITH',date=query.get('date');
const clock=t=>new Date(Math.round(t)*1000).toISOString().slice(11,19),json=async u=>{const r=await fetch(u,{cache:u.startsWith('api/')?'no-store':'default'});if(!r.ok)throw Error(await r.text());return r.json();},floats=async u=>{const r=await fetch(u);if(!r.ok)throw Error(`${r.status}: ${u}`);return new Float32Array(await r.arrayBuffer());};
const nearest=(a,x)=>a.reduce((best,v,i)=>Math.abs(v-x)<Math.abs(a[best]-x)?i:best,0);
const median=a=>{const b=[...a].sort((a,b)=>a-b);return b.length?(b[(b.length-1)>>1]+b[b.length>>1])/2:NaN;};
let rawSourceData,rawReceiverGains=null,displayGain=null,fitAzimuth=[],fitElevation=[];
let m,root,parameterData,rawData,plasmaData,dayStart,selected=0,selectedTime,selectedRange=0,plasmaIndex=0,points=[],revision=0,key,saveChain=Promise.resolve(),detailRequest=0,detailTimer,plasmaField,individualField,rawField,plasmaRawField,maps=[],profiles=[],lastPlasmaId;
let estimateTimer,estimateAbort,estimateRequest=0;let magicEstimate=null;let reviewGeneration=0;let saveError=null;let activeView='plasma';const allFields=[];
const records=()=>m.plasma_records[plasmaIndex];
const isNative=()=>m?.analysis_kind==='native';
function selectAnalysis(id){const params={analysis:id,time:selectedTime,range_m:m.range_m[selectedRange]};if(parent!==window)parent.postMessage({type:'ruisdap-analysis',...params},location.origin);else{const u=new URL(location.href);for(const [k,v] of Object.entries(params))u.searchParams.set(k,v);u.searchParams.delete('start');u.searchParams.delete('end');u.searchParams.set('view',activeView);location.assign(u);}}
function profileCard(parent,id,title){const card=document.createElement('article');card.className='plot';const heading=document.createElement('h3');heading.textContent=title;heading.className='scientific-heading';const asset=mathAsset(title);if(asset){const img=document.createElement('img');img.src=asset;img.alt=title;img.className='profile-math-heading';heading.replaceChildren(img);}card.append(heading);const plot=document.createElement('div');plot.id=id;card.append(plot);el(parent).append(card);return id;}
function drawProfile(id,values,ranges,{log=false,title='',caption='',errors=null,limits=null,compact=false}={}){
 const valid=values.filter(v=>Number.isFinite(v)&&(!log||v>0)),lo=valid.length?Math.min(...valid):0,hi=valid.length?Math.max(...valid):1;
 const xy=values.flatMap((v,i)=>[Number.isFinite(v)&&(!log||v>0)?(log?Math.log10(v):v):null,ranges[i]/1000]);
 if(errors)for(let i=0;i<values.length;i++){const e=errors[i],v=values[i];if(Number.isFinite(v)&&Number.isFinite(e)&&e>=0){valid.push(v+e);if(!log||v-e>0)valid.push(v-e);}}
 let x0=log?Math.log10(Math.max(lo,1e-30)):lo,x1=log?Math.log10(Math.max(hi,1e-29)):hi;const pad=Math.max((x1-x0)*.06,Math.abs(x1)*.001,.01);x0=Math.min(x0,...valid.map(v=>log?Math.log10(Math.max(v,1e-30)):v));x1=Math.max(x1,...valid.map(v=>log?Math.log10(Math.max(v,1e-30)):v));x0-=pad;x1+=pad;
 if(limits){x0=log?Math.log10(limits[0]):limits[0];x1=log?Math.log10(limits[1]):limits[1];}
 const r=m.range_m[selectedRange]/1000,errorSegments=[];if(errors)for(let i=0;i<values.length;i++){const e=errors[i],v=values[i];if(e!==null&&Number.isFinite(e)&&e>=0&&Number.isFinite(v)){const a=log?(v-e>0?Math.log10(v-e):x0):v-e,b=log?(v+e>0?Math.log10(v+e):x0):v+e;const y=ranges[i]/1000,cap=(ranges.at(-1)-ranges[0])/1000*.003;errorSegments.push(a,y,b,y,a,y-cap,a,y+cap,b,y-cap,b,y+cap);}}
 curve(id,{title,bounds:[x0,x1,...(compact?densityRange():[ranges[0]/1000,ranges.at(-1)/1000])],margins:compact?[12,22,12,58]:undefined,xlabel:title,xformat:log?'log10':'',ylabel:compact?'':'Slant range / km',yformat:compact?'none':'',pick:(x,y)=>chooseRange(y*1000),series:[{xy:errorSegments,segments:true,color:'#4d6170'},{xy,line:true,radius:1.8},{xy:[x0,r,x1,r],color:'#986327'}],caption});
}
function densityRange(){
 if(!individualField)return [m.range_m[0]/1000,m.range_m.at(-1)/1000];
 const b=individualField.viewport.bounds(),a=individualField.item.spec.bounds;
 return [a[2]+b[1]*(a[3]-a[2]),a[2]+b[3]*(a[3]-a[2])];
}
function drawDensityProfile(){
 const f=m.fits[selected];drawProfile('mh-density-profile',fitDensities(f),m.range_m,{log:true,compact:true,title:displayGain?'Ne / m⁻³':'Ne proxy',caption:'',errors:fitDensityErrors(f)});
}
function calibratedDensity(proxy,te,ti){
 if(!(proxy>0&&te>0&&ti>0))return NaN;
 const k=4*Math.PI*m.radar_frequency_hz/299792458,a=k*k*8.8541878188e-12*1.380649e-23*te/(1.602176634e-19**2),ratio=1+te/ti;
 const forward=ne=>ne*ratio/((1+a/ne)*(ratio+a/ne));let lo=0,hi=Math.max(proxy,1);
 while(forward(hi)<proxy)hi*=2;
 for(let i=0;i<65;i++){const mid=(lo+hi)/2;if(forward(mid)<proxy)lo=mid;else hi=mid;}return (lo+hi)/2;
}
function fitDensities(f){if(isNative())return f.parameters.map(row=>row[0]);return f.parameters.map((row,i)=>displayGain?calibratedDensity((f.uncalibrated_density_proxy?.[i]??row[0])*displayGain,row[2],row[1]):row[0]);}
function fitDensityErrors(f){if(isNative())return f.standard_error?.map(row=>row[0]);const densityValues=fitDensities(f);return f.parameters.map((row,i)=>{
 const sigma=(f.uncalibrated_standard_error||f.standard_error)?.[i]?.[0];if(!(Number.isFinite(sigma)&&sigma>=0))return null;
 if(!displayGain)return sigma;
 const ne=densityValues[i],te=row[2],ti=row[1],k=4*Math.PI*m.radar_frequency_hz/299792458,a=k*k*8.8541878188e-12*1.380649e-23*te/(1.602176634e-19**2),ratio=1+te/ti;
 const proxy=ne*ratio/((1+a/ne)*(ratio+a/ne)),derivative=proxy*(3/ne-1/(ne+a)-ratio/(ratio*ne+a));
 return displayGain*sigma/derivative;
});}
function fitCalibrationWindow(){
 if(!m||isNative()||!points.length||!maps.length)return;
 const a=Math.max(0,Math.min(15*3600,...points.map(p=>p.time_unix_s-dayStart-600))),b=Math.min(86400,Math.max(21*3600,...points.map(p=>p.time_unix_s-dayStart+600)));
 setTimeWindow([a,b]);
}
function sharedWindow(){const f=maps[0],a=f.item.spec.bounds,b=f.viewport.bounds();return [a[0]+(a[1]-a[0])*b[0],a[0]+(a[1]-a[0])*b[2]];}
function setTimeWindow(bounds){for(const f of allFields){const a=f.item.spec.bounds,v=f.viewport.bounds(),span=a[1]-a[0];f.viewport.reset();f.viewport.zoom_box((bounds[0]-a[0])/span,v[1],(bounds[1]-a[0])/span,v[3]);f.redraw();}drawInstruments();}
function drawInstruments(){
 if(!m||!maps.length)return;const time=m.fit_times.map(t=>t-dayStart),bounds=sharedWindow();
 const records=m.fits.map((f,i)=>({start:f.start_unix_s-dayStart,end:f.end_unix_s-dayStart,azimuth_deg:fitAzimuth[i]??f.azimuth_deg,elevation_deg:fitElevation[i]??f.elevation_deg,system_temperature_k:f.system_temperature_k,transmitter_w:f.transmitter_w}));
 curve('mh-instruments',instrumentSpec(records,bounds,{resetViewport:true,selected:selectedTime-dayStart,pick:x=>chooseTime(x+dayStart),onView:setTimeWindow}));
}
function applyDisplayCalibration(gain){
 if(!m||!maps.length)return;displayGain=isNative()?m.calibration.gain:Number.isFinite(gain)&&gain>0?gain:(m.calibration?.gain||null);
 const title=displayGain?'Ne / m⁻³':'Ne calibration unavailable',nt=m.fit_times.length,nr=m.range_m.length,values=new Float32Array(nt*nr);
 for(let t=0;t<nt;t++){const ne=fitDensities(m.fits[t]);for(let r=0;r<nr;r++)values[r*nt+t]=displayGain?ne[r]:NaN;}
 maps[0]=maps[0].update(values,{title,limits:EISCAT_LIMITS[0],mode:1});
 const physicalRaw=!!displayGain&&(m.raw_units==='m^-3'||rawReceiverGains);
 if(physicalRaw){const n=m.raw_times.length;rawData=new Float32Array(rawSourceData.length);for(let i=0;i<rawData.length;i++)rawData[i]=m.raw_units==='m^-3'?rawSourceData[i]*displayGain/m.calibration.gain:rawReceiverGains[i%n]>0?2*displayGain*rawSourceData[i]/rawReceiverGains[i%n]:NaN;}
 else rawData=rawSourceData;
 const rawTitle=physicalRaw?'Raw Ne / m⁻³':'Raw Ne calibration unavailable',settings={title:rawTitle,mode:1,colormap:'plasma',limits:EISCAT_LIMITS[4]};
 const shown=physicalRaw?rawData:new Float32Array(rawData.length).fill(NaN);rawField=rawField.update(shown,settings);plasmaRawField=plasmaRawField.update(shown,settings);
 el('mh-density-note').textContent=displayGain?`Ne uses C = ${displayGain.toExponential(5)}. Raw Ne assumes Te/Ti = 1 and a short 10 µs lag; recorded receiver gain and transmit power are removed. `:'Absolute Ne requires a saved calibration or included plasma-line calibration points. Select plasma-line features to estimate C.';
 if(isNative())el('mh-density-note').textContent='Native fits include pulse-energy/range-response normalization and recorded transmit power. Velocity toward/away sign remains unverified.';
 if(Number.isFinite(selectedTime))renderSelected();
}
async function loadRawCalibration(){if(isNative())return;try{const q=await json(`api/catalog/millstone/raw-calibration?date=${m.date}&antenna=${m.antenna}`);if(q.time_unix_s.length!==m.raw_times.length||q.time_unix_s.some((t,i)=>Math.abs(t-m.raw_times[i])>.01))throw Error('Receiver-gain times disagree with raw profiles.');rawReceiverGains=q.receiver_gain;fitAzimuth=q.azimuth_deg;fitElevation=q.elevation_deg;applyDisplayCalibration(magicEstimate?.magic_constant||m.calibration?.gain);drawInstruments();}catch(e){el('mh-density-note').textContent=e.message;}}
function chooseTime(t,{plasma=false}={}){
 plasmaIndex=nearest(m.plasma_times,t);selectedTime=plasma?m.plasma_times[plasmaIndex]:t;selected=nearest(m.fit_times,selectedTime);if(isNative()&&!plasma)selectedTime=m.fit_times[selected];el('millstone-time').value=selected;el('millstone-plasma-time').value=plasmaIndex;renderSelected();
}
function chooseRange(range){selectedRange=nearest(m.range_m,range);el('millstone-range').value=selectedRange;el('mh-range').value=(m.range_m[selectedRange]/1000).toFixed(2);if(Number.isFinite(selectedTime))renderSelected();}
function renderSelected(){
 if(!m)return;const f=m.fits[selected],r=m.range_m[selectedRange];el('selected-time').textContent=`${m.date} · ${clock(selectedTime)} UTC · ${m.antenna}`;
 el('selected-fit-context').textContent=`Ion fit ${clock(f.start_unix_s)}–${clock(f.end_unix_s)} UTC · Ptx ${(f.transmitter_w/1000).toFixed(0)} kW · Tsys ${f.system_temperature_k.toFixed(0)} K`;
 if(isNative()){el('selected-fit-context').textContent+=` · ${f.source_integrations} sampled integrations, ${f.exposure_s.toFixed(1)} s exposure · ${f.fit_status[selectedRange]}`;el('mh-fit-download').href=new URL(f.hdf5_url,new URL(root,location.href));}
 el('millstone-range-value').textContent=`${(r/1000).toFixed(1)} km`;el('mh-range').value=(r/1000).toFixed(2);
 el('plasma-integration-context').textContent=`Plasma line ${clock(records().start)}–${clock(records().end)} UTC · ${plasmaIndex+1} / ${m.plasma_times.length}`;
 for(let p=0;p<4;p++)drawProfile(profiles[p],p===0?fitDensities(f):f.parameters.map(v=>v[p]??NaN),m.range_m,{log:p===0,title:p===0?(displayGain?'Ne / m⁻³':'Ne calibration unavailable'):m.parameter_names[p],errors:p===0?fitDensityErrors(f):f.standard_error?.map(row=>row[p]),limits:p>0?EISCAT_LIMITS[p]:undefined,caption:isNative()?'Rust ACF fit with 1σ parameter errors; calibration, composition and shared-lag covariance excluded.':p===0?(m.calibration?'Calibrated density at fixed source Ti/Te.':'Source fitted density proxy; calibration is reviewed below.'):p===2?'Te uncertainty unavailable: source stores Te/Ti uncertainty without its covariance with Ti.':'Stored fit with 1σ source error bars.'});
 const rawIndex=nearest(m.raw_times,selectedTime),values=m.raw_range_m.map((_,r)=>rawData[r*m.raw_times.length+rawIndex]);
 drawProfile(profiles[4],values,m.raw_range_m,{title:displayGain&&(m.raw_units==='m^-3'||rawReceiverGains)?'Raw Ne / m⁻³':'Raw Ne calibration unavailable',caption:`Native 10 µs lag · ${clock(m.raw_times[rawIndex])} UTC`});
 drawDensityProfile();
 for(const f of allFields)f.setCursor(selectedTime-dayStart,r/1000);
 calculation();markers();individualField?.setCursor(Number(el('mh-frequency').value),r/1000);detailRequest++;if(individualField&&lastPlasmaId!==records().id){individualField.item.canvas.style.pointerEvents='none';individualField.article.style.opacity='.45';el('mh-save').disabled=true;}clearTimeout(detailTimer);detailTimer=setTimeout(loadDiagnostics,120);
}
function candidate(){
 if(isNative())throw Error('Calibration editing uses the historical analysis. Select Historical to edit the reviewed points.');
 const f=m.fits[selected],source=records(),t=(source.start+source.end)/2,frequency=Number(el('mh-frequency').value)*1e6,gate=selectedRange,range=m.range_m[gate],[stored,Ti,Te]=f.parameters[gate],proxy=f.uncalibrated_density_proxy?.[gate]??stored,sigma=Number(el('mh-frequency-error').value)*1000;
 if(!Number.isFinite(frequency)||Math.abs(frequency)<1e6||Math.abs(frequency)>12.5e6||!(proxy>0&&Ti>0&&Te>0))throw Error('Select a plasma-line frequency and a valid ion-line range.');
 if(!(Number.isFinite(f.transmitter_w)&&f.transmitter_w>0))throw Error('Missing or invalid recorded transmit power in ion fit.');
 if(!(sigma>0)||!Number.isFinite(sigma))throw Error('Frequency uncertainty must be positive.');if(Math.abs(m.fit_times[selected]-t)>300)throw Error('No ion-line fit within 300 seconds of this plasma integration.');
 const e=1.602176634e-19,me=9.1093837139e-31,eps=8.8541878188e-12,kb=1.380649e-23,k=4*Math.PI*m.radar_frequency_hz/299792458,omega2=(2*Math.PI*frequency)**2-3*k*k*kb*Te/me;
 if(!(omega2>0))throw Error('Thermal correction exceeds the selected frequency squared.');const ne=eps*me*omega2/e**2,q=k*k*eps*kb*Te/(ne*e*e),ratio=Te/Ti,gain=ne*(1+ratio)/((1+q)*(1+q+ratio))/proxy;
 const flags=[];if(q>=.1)flags.push('large Debye correction');if(f.transmitter_w<400e3)flags.push('low transmit power');if(f.system_temperature_k>2000)flags.push('high system temperature');if(f.space_object_count?.[gate]>0)flags.push('space object in fit gate');
 return {id:[...crypto.getRandomValues(new Uint8Array(16))].map(v=>v.toString(16).padStart(2,'0')).join(''),date:m.date,time_unix_s:t,observation_id:source.id,fit_id:f.id,station:m.station,antenna:m.antenna,mode:'lpi-30',frequency_hz:frequency,frequency_standard_error_hz:sigma,range_m:range,gate,plasma_density_m3:ne,k2_debye2:q,multiplier:gain,included:true,flags};
}
function calculation(){try{const p=candidate();el('mh-calculation').textContent=`Plasma-line Ne ${p.plasma_density_m3.toExponential(3)} m⁻³ · ${(p.range_m/1000).toFixed(1)} km · magic constant C ${p.multiplier.toExponential(4)} W m⁻³ K⁻¹ km⁻²${p.flags.length?' · '+p.flags.join('; '):''}`;}catch(e){el('mh-calculation').textContent=e.message;}}
function addPoint(){try{const p=candidate();if(points.some(a=>a.observation_id===p.observation_id&&a.gate===p.gate&&Math.abs(a.frequency_hz-p.frequency_hz)<1000))throw Error('This feature is already saved.');points.push(p);saved(true);}catch(e){el('mh-calculation').textContent=e.message;}}
function removePoint(index){if(isNative()){el('mh-detail-status').textContent='Select Historical to edit these calibration points.';return;}points.splice(index,1);saved(true);}
function markers(){
 if(!m)return;const marker=(p,i,x,y)=>({x,y,number:i+1,included:p.included,label:`Remove calibration point ${i+1}: ${clock(p.time_unix_s)} UTC, ${(p.frequency_hz/1e6).toFixed(3)} MHz, ${(p.range_m/1000).toFixed(1)} km`,remove:()=>removePoint(i)});
 plasmaField?.markers(points.map((p,i)=>marker(p,i,p.time_unix_s-dayStart,p.frequency_hz/1e6)));
 plasmaRawField?.markers(points.map((p,i)=>marker(p,i,p.time_unix_s-dayStart,p.range_m/1000)));
 individualField?.markers(points.flatMap((p,i)=>p.observation_id===records().id?[marker(p,i,p.frequency_hz/1e6,p.range_m/1000)]:[]));
}
function showCalibrationPoint(point){
 el('mh-frequency').value=(point.frequency_hz/1e6).toFixed(5);
 chooseRange(point.range_m);chooseTime(point.time_unix_s,{plasma:true});
 // Keep a selected point visible even outside the current time window.
 for(const field of allFields){
  const b=field.viewport.bounds(),span=field.xe.at(-1)-field.xe[0],x=(point.time_unix_s-dayStart-field.xe[0])/span;
  if(x>=b[0]&&x<=b[2])continue;
  const width=b[2]-b[0],lo=clip(x-width/2,0,1-width);
  field.viewport.reset();field.viewport.zoom_box(lo,b[1],lo+width,b[3]);field.redraw();
 }
 el('plasma-individual').scrollIntoView({block:'center'});
}
function saved(persist=false,updateEstimate=true){
 if(persist){reviewGeneration++;magicEstimate=null;}
 const table=el('cal-points');table.replaceChildren();const head=document.createElement('tr');for(const label of ['Use','UTC','Frequency / MHz','Range / km','Magic C','']){const th=document.createElement('th');th.textContent=label;head.append(th);}table.append(head);
 for(const [i,p]of points.entries()){const row=document.createElement('tr'),use=document.createElement('input');use.type='checkbox';use.checked=p.included;use.setAttribute('aria-label',`Include calibration point ${i+1}`);use.onchange=()=>{p.included=use.checked;saved(true);};const del=document.createElement('button');del.textContent=`Remove ${i+1}`;del.onclick=()=>removePoint(i);const view=document.createElement('button');view.className='cal-point-time';view.textContent=clock(p.time_unix_s);view.setAttribute('aria-label',`Show calibration point ${i+1}`);view.onclick=()=>showCalibrationPoint(p);for(const value of [use,view, (p.frequency_hz/1e6).toFixed(4),(p.range_m/1000).toFixed(1),p.multiplier.toExponential(4),del]){const td=document.createElement('td');if(value instanceof Node)td.append(value);else td.textContent=value;row.append(td);}table.append(row);}
 const groups=new Map();for(const p of points.filter(p=>p.included)){const a=groups.get(p.fit_id)||[];a.push(p.multiplier);groups.set(p.fit_id,a);}const gains=[...groups.values()].map(median),gain=median(gains);
 el('mh-calibration-summary').textContent=gains.length?`${points.filter(p=>p.included).length} included point${points.filter(p=>p.included).length===1?'':'s'} · ${gains.length} independent ion fit${gains.length===1?'':'s'}.`:'No included calibration points. Tick Use for the points to include.';
 const a=el('mh-export');if(a.dataset.url)URL.revokeObjectURL(a.dataset.url);a.dataset.url=URL.createObjectURL(new Blob([JSON.stringify({schema:'ruisdap-calibration-review-v1',revision,points,magic_constant_units:'W m^-3 K^-1 km^-2',magic_constant_estimate:magicEstimate},null,2)],{type:'application/json'}));a.href=a.dataset.url;
 localStorage.setItem(`ruisdap-millstone-points-${m.date}-${m.antenna}`,JSON.stringify(points));markers();
 el('mh-calibration-gains').hidden=!points.length;const good=points.filter(p=>p.multiplier>0),times=good.map(p=>p.time_unix_s-dayStart),v=good.map(p=>Math.log10(p.multiplier));if(good.length)curve('mh-calibration-gains',{title:'Magic constant C',bounds:[Math.min(...times)-60,Math.max(...times)+60,Math.min(...v)-.15,Math.max(...v)+.15],xlabel:'Time / UTC',xformat:'utc',ylabel:'C / (W m⁻³ K⁻¹ km⁻²)',yformat:'log10',series:[... [true,false].map(included=>({xy:good.filter(p=>p.included===included).flatMap(p=>[p.time_unix_s-dayStart,Math.log10(p.multiplier)]),style:included?0:1,radius:4,color:included?'#236b8e':'#758b98'}))],caption:'Per-point C previews. The estimate uses source-verified values and groups points sharing an ion fit.'});
 if(persist){const snapshot=structuredClone(points);saveChain=saveChain.then(async()=>{el('mh-save-status').textContent='Saving calibration…';const r=await fetch(`api/calibrations/${key}`,{method:'PUT',headers:{'Content-Type':'application/json'},body:JSON.stringify({revision,points:snapshot})});if(!r.ok)throw Error(`Calibration save ${r.status}; reload before editing if another tab changed this review.`);revision=(await r.json()).revision;saveError=null;el('mh-save-status').textContent=`Saved on avaruus · revision ${revision}`;}).catch(e=>{saveError=e.message;el('mh-save-status').textContent=e.message;});}
 if(persist)fitCalibrationWindow();
 if(updateEstimate)scheduleEstimate();
}
function scheduleEstimate(){
 clearTimeout(estimateTimer);estimateAbort?.abort();const request=++estimateRequest;
 magicEstimate=null;
 if(!points.some(p=>p.included)){el('mh-magic-constant').textContent=points.length?'No included points. Tick Use or choose Include all.':'Select plasma-line features; C is estimated automatically.';return;}
 el('mh-magic-constant').textContent='Updating magic constant…';
 estimateTimer=setTimeout(()=>estimateMagicConstant(request),120);
}
async function estimateMagicConstant(request){
 const generation=reviewGeneration,output=el('mh-magic-constant');
 try{
  await saveChain;
  if(request!==estimateRequest||generation!==reviewGeneration)return;
  if(saveError)throw Error(saveError);
  const requestedRevision=revision,controller=new AbortController();estimateAbort=controller;
  const response=await fetch('api/catalog/millstone/estimate',{method:'POST',headers:{'Content-Type':'application/json'},signal:controller.signal,body:JSON.stringify({notebook_key:key,revision:requestedRevision,date:m.date,antenna:m.antenna})});
  if(!response.ok)throw Error(await response.text());
  const result=await response.json();
  if(request!==estimateRequest||generation!==reviewGeneration||revision!==requestedRevision)return;
  magicEstimate=result;applyDisplayCalibration(result.magic_constant);
  for(const checked of result.points){const index=points.findIndex(p=>p.id===checked.id);if(index>=0)points[index]=checked;}
  saved(false,false);
  output.textContent=`C = ${result.magic_constant.toExponential(6)} W m⁻³ K⁻¹ km⁻² · ${result.independent_fit_count} independent ion fit${result.independent_fit_count===1?'':'s'} · ${result.point_count} point${result.point_count===1?'':'s'}${result.between_fit_scatter===null?' · between-fit scatter unavailable with one fit':` · between-fit scatter ${result.between_fit_scatter.toExponential(3)}`}${m.calibration?` · current C ${m.calibration.gain.toExponential(6)}, correction ×${(result.magic_constant/m.calibration.gain).toFixed(4)}`:''}`;
 }catch(error){if(request===estimateRequest&&error.name!=='AbortError')output.textContent=error.message;}
 finally{if(request===estimateRequest)estimateAbort=null;}
}
async function loadNotebook(){
 const name=`ruisdap-millstone-notebook-${m.date}-${m.antenna}`;key=localStorage.getItem(name);if(!key){key=[...crypto.getRandomValues(new Uint8Array(32))].map(v=>v.toString(16).padStart(2,'0')).join('');localStorage.setItem(name,key);}
 const n=await json(`api/calibrations/${key}`);revision=n.revision;points=n.points.filter(p=>p.antenna===m.antenna&&(!p.date||p.date===m.date));const inherited=!points.length&&m.calibration?.points?.length;if(inherited)points=structuredClone(m.calibration.points);saved(!!inherited);el('mh-save-status').textContent=`Loaded avaruus calibration · revision ${revision}`;
}
async function loadDiagnostics(){
 const request=++detailRequest,record=records();el('mh-detail-status').textContent='Loading selected native measurements…';
 try{
  const work=[];
  if(activeView==='plasma'&&lastPlasmaId!==record.id)work.push((async()=>{const p=await json(`api/catalog/millstone/plasma/${record.id}`),values=await floats(p.values_url);if(request!==detailRequest)return;individualField=fieldPanel('plasma-individual','mh-plasma-individual','Individual plasma-line power / source units',p.frequency_hz.map(f=>f/1e6),p.range_m.map(r=>r/1000),values,{time:false,palette:1,colormap:'plasma',mode:2,onView:drawDensityProfile,ylabel:'Slant range / km',xlabel:'Signed frequency offset / MHz',pick:(f,r)=>{chooseRange(r*1000);el('mh-frequency').value=f.toFixed(5);selectedTime=(record.start+record.end)/2;selected=nearest(m.fit_times,selectedTime);if(isNative())renderSelected();else{calculation();addPoint();}}});drawDensityProfile();individualField.setCursor(Number(el('mh-frequency').value),m.range_m[selectedRange]/1000);lastPlasmaId=record.id;individualField.item.canvas.style.pointerEvents='auto';individualField.article.style.opacity='1';el('mh-save').disabled=false;markers();})());
  if(activeView==='ion'&&isNative())work.push(loadFitDiagnostics(request));
  if(activeView==='ion'&&!isNative())work.push((async()=>{const acf=await json(`api/catalog/millstone/acf?date=${m.date}&antenna=${m.antenna}&time=${selectedTime}&range_m=${m.range_m[selectedRange]}`);if(request!==detailRequest)return;const plots=JSON.parse(acf_diagnostic(JSON.stringify(acf)));for(const [id,spec]of [['mh-acf-real',plots.real],['mh-acf-imag',plots.imaginary],['mh-acf-spectrum',plots.spectrum]])curve(id,spec);el('mh-acf-context').textContent=`Native ACF ${clock(acf.time)} UTC · ${(acf.range_m/1000).toFixed(1)} km · receiver-gain corrected. These are individual source integrations, not the full ${clock(m.fits[selected].start_unix_s)}–${clock(m.fits[selected].end_unix_s)} ion fit average.`;})());
  await Promise.all(work);if(request===detailRequest)el('mh-detail-status').textContent=`Individual plasma line ${clock(record.start)}–${clock(record.end)} UTC. Click a feature to add a calibration point; click a marker to remove it.`;
 }catch(e){if(request===detailRequest){el('mh-detail-status').textContent=e.message;if(activeView==='ion'){el('mh-acf-context').textContent=e.message;for(const id of ['mh-acf-real','mh-acf-imag','mh-acf-spectrum'])curve(id,{bounds:[0,1,-1,1],series:[],caption:e.message});}}}
}

const nativeDiagnostics=new Map();
async function loadFitDiagnostics(request){
 const f=m.fits[selected],url=new URL(f.diagnostics_url,new URL(root,location.href)).href;
 if(!nativeDiagnostics.has(url)){nativeDiagnostics.set(url,await json(url));if(nativeDiagnostics.size>4)nativeDiagnostics.delete(nativeDiagnostics.keys().next().value);}
 const data=nativeDiagnostics.get(url),g=data.gates[selectedRange];if(request!==detailRequest)return;
 const input={lag_s:data.lag_s,observed:g.observed,sigma:g.sigma};if(!g.observed.some((z,i)=>z.every(Number.isFinite)&&Number.isFinite(g.sigma[i])))throw Error('No retained ACF measurements at this selected range and integration.');const plots=JSON.parse(acf_diagnostic(JSON.stringify(input)));
 const spectralInput={...input,observed:input.observed.map(z=>[z[0],z[1]===null?null:-z[1]])};
 plots.spectrum=JSON.parse(acf_diagnostic(JSON.stringify(spectralInput))).spectrum;
 const modelInput={lag_s:[],observed:[],sigma:[]};
 for(let i=0;i<data.lag_s.length;i++)if(g.model[i].every(Number.isFinite)&&g.observed[i].every(Number.isFinite)&&Number.isFinite(g.sigma[i])){modelInput.lag_s.push(data.lag_s[i]);modelInput.observed.push(g.model[i]);modelInput.sigma.push(0);}
 if(modelInput.lag_s.length){
  const model=JSON.parse(acf_diagnostic(JSON.stringify(modelInput)));
  model.spectrum=JSON.parse(acf_diagnostic(JSON.stringify({...modelInput,observed:modelInput.observed.map(z=>[z[0],-z[1]])}))).spectrum;
  for(const key of ['real','imaginary','spectrum']){const series=structuredClone(model[key].series.at(-1));series.color='#c24334';series.radius=0;series.line=true;plots[key].series.push(series);}
 }
 for(const [part,key]of ['real','imaginary'].entries()){
  const values=[];for(let i=0;i<data.lag_s.length;i++){if(Number.isFinite(g.observed[i][part])&&Number.isFinite(g.sigma[i]))values.push(g.observed[i][part]-g.sigma[i],g.observed[i][part]+g.sigma[i]);if(Number.isFinite(g.model[i][part]))values.push(g.model[i][part]);}
  const lo=Math.min(0,...values),hi=Math.max(0,...values),pad=Math.max((hi-lo)*.08,1e-12);plots[key].bounds[2]=lo-pad;plots[key].bounds[3]=hi+pad;plots[key].ylabel='ACF / ('+data.units+')';plots[key].caption='Measured ACF ±1σ (blue); fitted model (red). Diagonal lag variance; shared-lag covariance unavailable.';
 }
 const spectrum=plots.spectrum;
 for(const series of spectrum.series)for(let i=1;i<series.xy.length;i+=2)series.xy[i]*=1000;
 for(const band of spectrum.bands||[])for(let i=0;i<band.values.length;i++)if(i%3)band.values[i]*=1000;
 const powers=spectrum.series.flatMap(s=>s.xy.filter((_,i)=>i%2)).filter(Number.isFinite),bands=(spectrum.bands||[]).flatMap(b=>b.values.filter((_,i)=>i%3)).filter(Number.isFinite),lo=Math.min(0,...powers,...bands),hi=Math.max(0,...powers,...bands),pad=Math.max((hi-lo)*.08,1e-15);
 spectrum.bounds[2]=lo-pad;spectrum.bounds[3]=hi+pad;spectrum.ylabel='ACF DFT / ('+data.units+' kHz⁻¹)';spectrum.caption='DFT of the same measured lags: blue data and red fitted model. Missing zero lag leaves the baseline undetermined. Native negative-lag convention retained.';
 for(const [id,key]of [['mh-acf-real','real'],['mh-acf-imag','imaginary'],['mh-acf-spectrum','spectrum']])curve(id,plots[key]);
 el('mh-acf-context').textContent=`Rust fit ${clock(f.start_unix_s)}–${clock(f.end_unix_s)} UTC · ${(g.range_m/1000).toFixed(1)} km · ${f.source_integrations} sampled integrations (${f.exposure_s.toFixed(1)} s exposure) · ${g.status}${Number.isFinite(g.reduced_chi_square)?' · reduced χ² '+g.reduced_chi_square.toFixed(2):''}. Observed and model ACFs use the same fit average, receiver response and per-source transmit-power normalization.`;
}

function setView(kind){activeView=kind;el('millstone-plasma').hidden=kind==='ion';el('millstone-ion').hidden=kind!=='ion';el('mh-acf-section').hidden=kind!=='ion';if(m){for(const f of allFields)f.redraw();requestAnimationFrame(renderSelected);}}
window.addEventListener('message',e=>{if(e.origin===location.origin&&e.source===parent&&e.data?.type==='ruisdap-view')setView(e.data.view);});
activeView=query.get('view')==='ion'?'ion':'plasma';
try{
 await init();await document.fonts.load('18px RustplotlibPagella');const config=await json('api/catalog/millstone/datasets'),d=config.datasets.find(d=>d.antenna===antenna&&(!date||d.date===date));if(!d)throw Error('Dataset unavailable');const choices=[{id:'historical',label:'Historical · full day',data_root:d.data_root},...(d.analyses||[])],analysis=query.get('analysis')||'historical',choice=choices.find(a=>a.id===analysis);if(!choice)throw Error('Requested analysis is unavailable');for(const a of choices){const o=document.createElement('option');o.value=a.id;o.textContent=a.label;el('mh-analysis').append(o);}el('mh-analysis').value=analysis;el('mh-analysis').onchange=()=>selectAnalysis(el('mh-analysis').value);root=choice.data_root;m=await json(root+'manifest.json');el('mh-analysis-note').textContent=isNative()?'Experimental Rust · sparse 180 s fits · O+ · 200–540 km. Raw Ne and plasma-line panels use historical context.':'Historical full-day analysis with reviewed plasma-line density calibration.';el('mh-fit-download').hidden=!isNative();el('mh-calibration-editor').hidden=isNative();[parameterData,rawData,plasmaData]=await Promise.all([m.parameters_url,m.raw_url,m.plasma_url].map(u=>floats(new URL(u,new URL(root,location.href)))));rawSourceData=rawData.slice();dayStart=Date.parse(m.date+'T00:00:00Z')/1000;
 const time=m.fit_times.map(t=>t-dayStart),ranges=m.range_m.map(r=>r/1000),rawTime=m.raw_times.map(t=>t-dayStart),rawRange=m.raw_range_m.map(r=>r/1000),fitIntervals=m.fits.map(f=>[f.start_unix_s-dayStart,f.end_unix_s-dayStart]),nt=time.length,nr=ranges.length;
 plasmaField=fieldPanel('plasma-timeline','mh-plasma-field','Collapsed plasma-line power / source units',m.plasma_times.map(t=>t-dayStart),m.plasma_frequency_hz.map(f=>f/1e6),plasmaData,{onView:drawInstruments,intervals:m.plasma_records.map(r=>[r.start-dayStart,r.end-dayStart]),colormap:'plasma',mode:2,ylabel:'Plasma frequency / MHz',pick:(t,f)=>{el('mh-frequency').value=f.toFixed(5);chooseTime(t+dayStart,{plasma:true});}});allFields.push(plasmaField);
 for(let p=0;p<4;p++){const f=fieldPanel('parameter-maps','mh-map-'+p,m.parameter_names[p],time,ranges,parameterData.subarray(p*nt*nr,(p+1)*nt*nr),{onView:drawInstruments,intervals:fitIntervals,mode:p===0?1:0,palette:p===3?0:1,colormap:p===3?'seismic':'eiscat',sharedTime:true,limits:p===0&&!m.calibration?undefined:EISCAT_LIMITS[p],pick:(t,r)=>{chooseRange(r*1000);chooseTime(t+dayStart);}});maps.push(f);allFields.push(f);profiles.push(profileCard('parameter-profiles','mh-profile-'+p,m.parameter_names[p]));}
 const rawOptions={onView:drawInstruments,mode:2,palette:1,colormap:'plasma',pick:(t,r)=>{chooseRange(r*1000);chooseTime(t+dayStart,{plasma:true});}};
 rawField=fieldPanel('parameter-maps','mh-raw-field','Raw Ne / m⁻³',rawTime,rawRange,rawData,{...rawOptions,sharedTime:true,sharedTimeLast:true});plasmaRawField=fieldPanel('plasma-raw-ne','mh-plasma-raw-field','Raw Ne / m⁻³',rawTime,rawRange,rawData,rawOptions);allFields.push(rawField,plasmaRawField);profiles.push(profileCard('parameter-profiles','mh-profile-raw','Raw Ne'));
 profileCard('plasma-density-profile','mh-density-profile','Ion-line Ne');
 for(const [id,max,fn]of [['millstone-time',m.fit_times.length-1,i=>chooseTime(m.fit_times[i])],['millstone-plasma-time',m.plasma_times.length-1,i=>chooseTime(m.plasma_times[i],{plasma:true})],['millstone-range',m.range_m.length-1,i=>chooseRange(m.range_m[i])]]){el(id).max=max;el(id).oninput=()=>fn(Number(el(id).value));}
 for(const [id,lo,hi]of [['mh-full-day',0,86400],['mh-eclipse-window',15*3600,21*3600]])el(id).onclick=()=>{for(const f of allFields){f.viewport.reset();f.viewport.zoom_box(clip((lo-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),0,clip((hi-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),1);f.redraw();}drawInstruments();};
 for(const f of allFields){f.viewport.zoom_box(clip((15*3600-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),0,clip((21*3600-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),1);f.redraw();}
 el('mh-frequency').oninput=calculation;el('mh-frequency-error').oninput=calculation;el('mh-range').onchange=()=>chooseRange(Number(el('mh-range').value)*1000);el('mh-save').onclick=addPoint;
 el('mh-include-all').onclick=()=>{if(points.some(p=>!p.included)){for(const p of points)p.included=true;saved(true);}};
 el('mh-recalculate').onclick=async()=>{try{await saveChain;if(saveError)throw Error(saveError);const token=el('mh-processing-key').value.trim();if(!token)throw Error('Enter the processing key under Processing access.');el('mh-recalculate').disabled=true;el('mh-recalculate-status').textContent='Recalculating density from the saved review…';const r=await fetch('api/catalog/millstone/calibrate',{method:'POST',headers:{'Content-Type':'application/json',Authorization:'Bearer '+token},body:JSON.stringify({notebook_key:key,revision,date:m.date,antenna:m.antenna})});if(!r.ok)throw Error(await r.text());const result=await r.json();el('mh-recalculate-status').textContent=`${result.message} Magic constant C ${result.gain.toExponential(5)} W m⁻³ K⁻¹ km⁻². Reload the dataset to open the new density product.`;const link=document.createElement('a');link.href=result.hdf5_url;link.textContent='Download calibrated HDF5';el('mh-recalculate-status').append(' ',link);const reload=document.createElement('button');reload.textContent='Open recalculated data';reload.onclick=()=>location.reload();el('mh-recalculate-status').append(' ',reload);}catch(e){el('mh-recalculate-status').textContent=e.message;}finally{el('mh-recalculate').disabled=false;}};
 if(isNative()){points=structuredClone(m.calibration?.points||[]);markers();}else{await loadNotebook();fitCalibrationWindow();}applyDisplayCalibration(m.calibration?.gain);loadRawCalibration();drawInstruments();chooseRange(300e3);chooseTime(Number(query.get('time'))|| (isNative()?m.fit_times[0]:m.plasma_times.find(t=>t>=dayStart+18*3600)||m.fit_times[0]),{plasma:!isNative()});if(query.has('range_m'))chooseRange(Number(query.get('range_m')));const windowStart=Number(query.get('start')),windowEnd=Number(query.get('end'));if(isNative()){const [lo,hi]=m.analysis_default_window,a=Math.max(lo,windowStart),b=Math.min(hi,windowEnd);setTimeWindow((windowEnd>windowStart&&b>a?[a,b]:[lo,hi]).map(t=>t-dayStart));}else setTimeWindow([0,86400]);setView(activeView);
 if(m.calibration&&!isNative())el('mh-density-note').textContent=`Plasma-line calibrated Ne · C ${m.calibration.gain.toExponential(5)} W m⁻³ K⁻¹ km⁻² · revision ${m.calibration.revision}. Source Ti, Te, velocity and their fit remain unchanged. Raw Ne uses Te/Ti = 1.`;
 el('millstone-status').textContent=`${choice.label} · ${m.antenna} · ${m.fits.length} fitted integrations · ${m.raw_times.length} ${isNative()?'historical raw profiles':'native ACF profiles'} · ${m.plasma_times.length} plasma-line integrations`;
 parent.postMessage({type:'ruisdap-ready'},location.origin);
}catch(e){el('millstone-status').textContent='Unable to load Millstone measurements: '+e.message;console.error(e);}
function clip(v,a,b){return Math.max(a,Math.min(b,v));}
