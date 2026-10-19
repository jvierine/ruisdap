import init,{acf_diagnostic} from './pkg/guisdap_web.js';
import {fieldPanel,curve,iqr} from './scientific-plots.js';
const el=id=>document.getElementById(id),query=new URLSearchParams(location.search),antenna=query.get('antenna')||'ZENITH',date=query.get('date');
const clock=t=>new Date(Math.round(t)*1000).toISOString().slice(11,19),json=async u=>{const r=await fetch(u);if(!r.ok)throw Error(await r.text());return r.json();},floats=async u=>{const r=await fetch(u);if(!r.ok)throw Error(`${r.status}: ${u}`);return new Float32Array(await r.arrayBuffer());};
const nearest=(a,x)=>a.reduce((best,v,i)=>Math.abs(v-x)<Math.abs(a[best]-x)?i:best,0);
const median=a=>{const b=[...a].sort((a,b)=>a-b);return b.length?(b[(b.length-1)>>1]+b[b.length>>1])/2:NaN;};
let m,root,parameterData,rawData,plasmaData,dayStart,selected=0,selectedTime,selectedRange=0,plasmaIndex=0,points=[],revision=0,key,saveChain=Promise.resolve(),detailRequest=0,detailTimer,plasmaField,individualField,rawField,plasmaRawField,maps=[],profiles=[],lastPlasmaId;
let saveError=null;let activeView='plasma';const allFields=[];
const records=()=>m.plasma_records[plasmaIndex];
function profileCard(parent,id,title){const card=document.createElement('article');card.className='plot';const heading=document.createElement('h3');heading.textContent=title;heading.className='scientific-heading';card.append(heading);const plot=document.createElement('div');plot.id=id;card.append(plot);el(parent).append(card);return id;}
function drawProfile(id,values,ranges,{log=false,title='',caption='',errors=null}={}){
 const valid=values.filter(v=>Number.isFinite(v)&&(!log||v>0)),lo=valid.length?Math.min(...valid):0,hi=valid.length?Math.max(...valid):1;
 const xy=values.flatMap((v,i)=>[Number.isFinite(v)&&(!log||v>0)?(log?Math.log10(v):v):null,ranges[i]/1000]);
 let x0=log?Math.log10(Math.max(lo,1e-30)):lo,x1=log?Math.log10(Math.max(hi,1e-29)):hi;const pad=Math.max((x1-x0)*.06,Math.abs(x1)*.001,.01);x0-=pad;x1+=pad;
 const r=m.range_m[selectedRange]/1000,errorSegments=[];if(errors)for(let i=0;i<values.length;i++){const e=errors[i],v=values[i];if(e!==null&&Number.isFinite(e)&&e>=0&&Number.isFinite(v)){const a=log?(v-e>0?Math.log10(v-e):x0):v-e,b=log?(v+e>0?Math.log10(v+e):x0):v+e;errorSegments.push(a,ranges[i]/1000,b,ranges[i]/1000);}}
 curve(id,{title,bounds:[x0,x1,ranges[0]/1000,ranges.at(-1)/1000],xlabel:title,xformat:log?'log10':'',ylabel:'Slant range / km',pick:(x,y)=>chooseRange(y*1000),series:[{xy:errorSegments,segments:true,color:'#7b929e'},{xy,line:true,radius:1.2},{xy:[x0,r,x1,r],color:'#986327'}],caption});
}
function chooseTime(t,{plasma=false}={}){
 plasmaIndex=nearest(m.plasma_times,t);selectedTime=plasma?m.plasma_times[plasmaIndex]:t;selected=nearest(m.fit_times,selectedTime);el('millstone-time').value=selected;el('millstone-plasma-time').value=plasmaIndex;renderSelected();
}
function chooseRange(range){selectedRange=nearest(m.range_m,range);el('millstone-range').value=selectedRange;el('mh-range').value=(m.range_m[selectedRange]/1000).toFixed(2);if(Number.isFinite(selectedTime))renderSelected();}
function renderSelected(){
 if(!m)return;const f=m.fits[selected],r=m.range_m[selectedRange];el('selected-time').textContent=`${m.date} · ${clock(selectedTime)} UTC · ${m.antenna}`;
 el('selected-fit-context').textContent=`Ion fit ${clock(f.start_unix_s)}–${clock(f.end_unix_s)} UTC · Ptx ${(f.transmitter_w/1000).toFixed(0)} kW · Tsys ${f.system_temperature_k.toFixed(0)} K`;
 el('millstone-range-value').textContent=`${(r/1000).toFixed(1)} km`;el('mh-range').value=(r/1000).toFixed(2);
 el('plasma-integration-context').textContent=`Plasma line ${clock(records().start)}–${clock(records().end)} UTC · ${plasmaIndex+1} / ${m.plasma_times.length}`;
 for(let p=0;p<4;p++)drawProfile(profiles[p],f.parameters.map(v=>v[p]??NaN),m.range_m,{log:p===0,title:m.parameter_names[p],errors:f.standard_error?.map(row=>row[p]),caption:p===0?(m.calibration?'Calibrated density at fixed source Ti/Te.':'Source fitted density proxy; calibration is reviewed below.'):'Stored fit; source uncertainty retained.'});
 const rawIndex=nearest(m.raw_times,selectedTime),values=m.raw_range_m.map((_,r)=>rawData[r*m.raw_times.length+rawIndex]);
 drawProfile(profiles[4],values,m.raw_range_m,{title:'Signed raw Ne proxy',caption:`Native 10 µs lag · ${clock(m.raw_times[rawIndex])} UTC`});
 drawProfile('mh-density-profile',f.parameters.map(v=>v[0]??NaN),m.range_m,{log:true,title:m.parameter_names[0],caption:'Range selection is shared with the ACF and individual plasma-line plot.'});
 for(const f of allFields)f.setCursor(selectedTime-dayStart,r/1000);
 calculation();markers();detailRequest++;if(individualField&&lastPlasmaId!==records().id){individualField.item.canvas.style.pointerEvents='none';individualField.article.style.opacity='.45';el('mh-save').disabled=true;}clearTimeout(detailTimer);detailTimer=setTimeout(loadDiagnostics,120);
}
function candidate(){
 const f=m.fits[selected],source=records(),t=(source.start+source.end)/2,frequency=Number(el('mh-frequency').value)*1e6,gate=selectedRange,range=m.range_m[gate],[stored,Ti,Te]=f.parameters[gate],proxy=f.uncalibrated_density_proxy?.[gate]??stored,sigma=Number(el('mh-frequency-error').value)*1000;
 if(!Number.isFinite(frequency)||Math.abs(frequency)<1e6||Math.abs(frequency)>12.5e6||!(proxy>0&&Ti>0&&Te>0))throw Error('Select a plasma-line frequency and a valid ion-line range.');
 if(!(sigma>0)||!Number.isFinite(sigma))throw Error('Frequency uncertainty must be positive.');if(Math.abs(m.fit_times[selected]-t)>300)throw Error('No ion-line fit within 300 seconds of this plasma integration.');
 const e=1.602176634e-19,me=9.1093837139e-31,eps=8.8541878188e-12,kb=1.380649e-23,k=4*Math.PI*m.radar_frequency_hz/299792458,omega2=(2*Math.PI*frequency)**2-3*k*k*kb*Te/me;
 if(!(omega2>0))throw Error('Thermal correction exceeds the selected frequency squared.');const ne=eps*me*omega2/e**2,q=k*k*eps*kb*Te/(ne*e*e),ratio=Te/Ti,gain=ne*(1+ratio)/((1+q)*(1+q+ratio))/proxy;
 const flags=[];if(q>=.1)flags.push('large Debye correction');if(f.transmitter_w<400e3)flags.push('low transmit power');if(f.system_temperature_k>2000)flags.push('high system temperature');if(f.space_object_count?.[gate]>0)flags.push('space object in fit gate');
 return {id:[...crypto.getRandomValues(new Uint8Array(16))].map(v=>v.toString(16).padStart(2,'0')).join(''),date:m.date,time_unix_s:t,observation_id:source.id,fit_id:f.id,station:m.station,antenna:m.antenna,mode:'lpi-30',frequency_hz:frequency,frequency_standard_error_hz:sigma,range_m:range,gate,plasma_density_m3:ne,k2_debye2:q,multiplier:gain,included:false,flags};
}
function calculation(){try{const p=candidate();el('mh-calculation').textContent=`Plasma-line Ne ${p.plasma_density_m3.toExponential(3)} m⁻³ · ${(p.range_m/1000).toFixed(1)} km · gain ${p.multiplier.toExponential(4)}${p.flags.length?' · '+p.flags.join('; '):''}`;}catch(e){el('mh-calculation').textContent=e.message;}}
function addPoint(){try{const p=candidate();if(points.some(a=>a.observation_id===p.observation_id&&a.gate===p.gate&&Math.abs(a.frequency_hz-p.frequency_hz)<1000))throw Error('This feature is already saved.');points.push(p);saved(true);}catch(e){el('mh-calculation').textContent=e.message;}}
function removePoint(index){points.splice(index,1);saved(true);}
function markers(){
 if(!m)return;const marker=(p,i,x,y)=>({x,y,included:p.included,label:`Remove calibration point ${i+1}: ${clock(p.time_unix_s)} UTC, ${(p.frequency_hz/1e6).toFixed(3)} MHz, ${(p.range_m/1000).toFixed(1)} km`,remove:()=>removePoint(i)});
 plasmaField?.markers(points.map((p,i)=>marker(p,i,p.time_unix_s-dayStart,p.frequency_hz/1e6)));
 plasmaRawField?.markers(points.map((p,i)=>marker(p,i,p.time_unix_s-dayStart,p.range_m/1000)));
 individualField?.markers(points.flatMap((p,i)=>p.observation_id===records().id?[marker(p,i,p.frequency_hz/1e6,p.range_m/1000)]:[]));
}
function saved(persist=false){
 const table=el('cal-points');table.replaceChildren();const head=document.createElement('tr');for(const label of ['Use','UTC','Frequency / MHz','Range / km','Gain','']){const th=document.createElement('th');th.textContent=label;head.append(th);}table.append(head);
 for(const [i,p]of points.entries()){const row=document.createElement('tr'),use=document.createElement('input');use.type='checkbox';use.checked=p.included;use.setAttribute('aria-label',`Include calibration point ${i+1}`);use.onchange=()=>{p.included=use.checked;saved(true);};const del=document.createElement('button');del.textContent=`Remove ${i+1}`;del.onclick=()=>removePoint(i);for(const value of [use,clock(p.time_unix_s), (p.frequency_hz/1e6).toFixed(4),(p.range_m/1000).toFixed(1),p.multiplier.toExponential(4),del]){const td=document.createElement('td');if(value instanceof Node)td.append(value);else td.textContent=value;row.append(td);}table.append(row);}
 const groups=new Map();for(const p of points.filter(p=>p.included)){const a=groups.get(p.fit_id)||[];a.push(p.multiplier);groups.set(p.fit_id,a);}const gains=[...groups.values()].map(median),gain=median(gains);
 el('mh-calibration-summary').textContent=gains.length?`Reviewed gain ${gain.toExponential(5)} · ${gains.length} independent fit integrations · ${points.filter(p=>p.included).length} included points. Median within each fit, then across fits.`:'No included calibration points. New points are excluded until reviewed.';
 const a=el('mh-export');if(a.dataset.url)URL.revokeObjectURL(a.dataset.url);a.dataset.url=URL.createObjectURL(new Blob([JSON.stringify({schema:'ruisdap-calibration-review-v1',revision,points},null,2)],{type:'application/json'}));a.href=a.dataset.url;
 localStorage.setItem(`ruisdap-millstone-points-${m.date}-${m.antenna}`,JSON.stringify(points));markers();
 el('mh-calibration-gains').hidden=!points.length;const good=points.filter(p=>p.multiplier>0),times=good.map(p=>p.time_unix_s-dayStart),v=good.map(p=>Math.log10(p.multiplier));if(good.length)curve('mh-calibration-gains',{title:'Calibration gain',bounds:[Math.min(...times)-60,Math.max(...times)+60,Math.min(...v)-.15,Math.max(...v)+.15],xlabel:'Time / UTC',xformat:'utc',ylabel:'Gain',yformat:'log10',series:[{xy:good.flatMap(p=>[p.time_unix_s-dayStart,Math.log10(p.multiplier)]),style:1,radius:4}],caption:'Repeated points from one ion integration are grouped before combining gains.'});
 if(persist){const snapshot=structuredClone(points);saveChain=saveChain.then(async()=>{el('mh-save-status').textContent='Saving calibration…';const r=await fetch(`api/calibrations/${key}`,{method:'PUT',headers:{'Content-Type':'application/json'},body:JSON.stringify({revision,points:snapshot})});if(!r.ok)throw Error(`Calibration save ${r.status}; reload before editing if another tab changed this review.`);revision=(await r.json()).revision;saveError=null;el('mh-save-status').textContent=`Saved on avaruus · revision ${revision}`;}).catch(e=>{saveError=e.message;el('mh-save-status').textContent=e.message;});}
}
async function loadNotebook(){
 const name=`ruisdap-millstone-notebook-${m.date}-${m.antenna}`;key=localStorage.getItem(name);if(!key){key=[...crypto.getRandomValues(new Uint8Array(32))].map(v=>v.toString(16).padStart(2,'0')).join('');localStorage.setItem(name,key);}
 const n=await json(`api/calibrations/${key}`);revision=n.revision;points=n.points.filter(p=>p.antenna===m.antenna&&(!p.date||p.date===m.date));saved();el('mh-save-status').textContent=`Loaded avaruus calibration · revision ${revision}`;
}
async function loadDiagnostics(){
 const request=++detailRequest,record=records();el('mh-detail-status').textContent='Loading selected native measurements…';
 try{
  const work=[];
  if(lastPlasmaId!==record.id)work.push((async()=>{const p=await json(`api/catalog/millstone/plasma/${record.id}`),values=await floats(p.values_url);if(request!==detailRequest)return;individualField=fieldPanel('plasma-individual','mh-plasma-individual','Individual plasma-line power / source units',p.frequency_hz.map(f=>f/1e6),p.range_m.map(r=>r/1000),values,{time:false,palette:1,mode:2,ylabel:'Slant range / km',xlabel:'Signed frequency offset / MHz',pick:(f,r)=>{chooseRange(r*1000);el('mh-frequency').value=f.toFixed(5);selectedTime=(record.start+record.end)/2;selected=nearest(m.fit_times,selectedTime);calculation();addPoint();}});lastPlasmaId=record.id;individualField.item.canvas.style.pointerEvents='auto';individualField.article.style.opacity='1';el('mh-save').disabled=false;markers();})());
  work.push((async()=>{const acf=await json(`api/catalog/millstone/acf?date=${m.date}&antenna=${m.antenna}&time=${selectedTime}&range_m=${m.range_m[selectedRange]}`);if(request!==detailRequest)return;const plots=JSON.parse(acf_diagnostic(JSON.stringify(acf)));for(const [id,spec]of [['mh-acf-real',plots.real],['mh-acf-imag',plots.imaginary],['mh-acf-spectrum',plots.spectrum]])curve(id,spec);el('mh-acf-context').textContent=`Native ACF ${clock(acf.time)} UTC · ${(acf.range_m/1000).toFixed(1)} km · receiver-gain corrected. These are individual source integrations, not the full ${clock(m.fits[selected].start_unix_s)}–${clock(m.fits[selected].end_unix_s)} ion fit average.`;})());
  await Promise.all(work);if(request===detailRequest)el('mh-detail-status').textContent=`Individual plasma line ${clock(record.start)}–${clock(record.end)} UTC. Click a feature to add a calibration point; click a circle to remove it.`;
 }catch(e){if(request===detailRequest)el('mh-detail-status').textContent=e.message;}
}
function setView(kind){activeView=kind;el('millstone-plasma').hidden=kind==='ion';el('millstone-ion').hidden=kind!=='ion';if(m){for(const f of allFields)f.redraw();requestAnimationFrame(renderSelected);}}
window.addEventListener('message',e=>{if(e.origin===location.origin&&e.source===parent&&e.data?.type==='ruisdap-view')setView(e.data.view);});
try{
 await init();await document.fonts.load('18px RustplotlibPagella');const config=await json('catalog/millstone-config.json').catch(()=>json('millstone-config.json')),d=config.datasets.find(d=>d.antenna===antenna&&(!date||d.date===date));if(!d)throw Error('Dataset unavailable');root=d.data_root;m=await json(root+'manifest.json');[parameterData,rawData,plasmaData]=await Promise.all([m.parameters_url,m.raw_url,m.plasma_url].map(u=>floats(new URL(u,new URL(root,location.href)))));dayStart=Date.parse(m.date+'T00:00:00Z')/1000;
 const time=m.fit_times.map(t=>t-dayStart),ranges=m.range_m.map(r=>r/1000),rawTime=m.raw_times.map(t=>t-dayStart),rawRange=m.raw_range_m.map(r=>r/1000),fitIntervals=m.fits.map(f=>[f.start_unix_s-dayStart,f.end_unix_s-dayStart]),nt=time.length,nr=ranges.length;
 plasmaField=fieldPanel('plasma-timeline','mh-plasma-field','Collapsed plasma-line power / source units',m.plasma_times.map(t=>t-dayStart),m.plasma_frequency_hz.map(f=>f/1e6),plasmaData,{intervals:m.plasma_records.map(r=>[r.start-dayStart,r.end-dayStart]),mode:2,ylabel:'Plasma frequency / MHz',pick:(t,f)=>{el('mh-frequency').value=f.toFixed(5);chooseTime(t+dayStart,{plasma:true});}});allFields.push(plasmaField);
 for(let p=0;p<4;p++){const f=fieldPanel('parameter-maps','mh-map-'+p,m.parameter_names[p],time,ranges,parameterData.subarray(p*nt*nr,(p+1)*nt*nr),{intervals:fitIntervals,mode:p===0?1:0,palette:p===3?0:1,pick:(t,r)=>{chooseRange(r*1000);chooseTime(t+dayStart);}});maps.push(f);allFields.push(f);profiles.push(profileCard('parameter-profiles','mh-profile-'+p,m.parameter_names[p]));}
 const rawOptions={mode:2,palette:0,pick:(t,r)=>{chooseRange(r*1000);chooseTime(t+dayStart,{plasma:true});}};
 rawField=fieldPanel('parameter-maps','mh-raw-field','Signed raw Ne proxy · 10 µs lag',rawTime,rawRange,rawData,rawOptions);plasmaRawField=fieldPanel('plasma-raw-ne','mh-plasma-raw-field','Signed raw Ne proxy · 10 µs lag',rawTime,rawRange,rawData,rawOptions);allFields.push(rawField,plasmaRawField);profiles.push(profileCard('parameter-profiles','mh-profile-raw','Signed raw Ne proxy'));
 profileCard('plasma-density-profile','mh-density-profile','Ion-line density counterpart');
 for(const [id,max,fn]of [['millstone-time',m.fit_times.length-1,i=>chooseTime(m.fit_times[i])],['millstone-plasma-time',m.plasma_times.length-1,i=>chooseTime(m.plasma_times[i],{plasma:true})],['millstone-range',m.range_m.length-1,i=>chooseRange(m.range_m[i])]]){el(id).max=max;el(id).oninput=()=>fn(Number(el(id).value));}
 for(const [id,lo,hi]of [['mh-full-day',0,86400],['mh-eclipse-window',15*3600,21*3600]])el(id).onclick=()=>{for(const f of allFields){f.viewport.reset();f.viewport.zoom_box(clip((lo-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),0,clip((hi-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),1);f.redraw();}};
 for(const f of allFields){f.viewport.zoom_box(clip((15*3600-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),0,clip((21*3600-f.xe[0])/(f.xe.at(-1)-f.xe[0]),0,1),1);f.redraw();}
 el('mh-frequency').oninput=calculation;el('mh-frequency-error').oninput=calculation;el('mh-range').onchange=()=>chooseRange(Number(el('mh-range').value)*1000);el('mh-save').onclick=addPoint;
 el('mh-recalculate').onclick=async()=>{try{await saveChain;if(saveError)throw Error(saveError);const token=el('mh-processing-key').value.trim();if(!token)throw Error('Enter the processing key under Processing access.');el('mh-recalculate').disabled=true;el('mh-recalculate-status').textContent='Recalculating density from the saved review…';const r=await fetch('api/catalog/millstone/calibrate',{method:'POST',headers:{'Content-Type':'application/json',Authorization:'Bearer '+token},body:JSON.stringify({notebook_key:key,revision,date:m.date,antenna:m.antenna})});if(!r.ok)throw Error(await r.text());const result=await r.json();el('mh-recalculate-status').textContent=`${result.message} Gain ${result.gain.toExponential(5)}. Reload the dataset to open the new density product.`;const link=document.createElement('a');link.href=result.hdf5_url;link.textContent='Download calibrated HDF5';el('mh-recalculate-status').append(' ',link);const reload=document.createElement('button');reload.textContent='Open recalculated data';reload.onclick=()=>location.reload();el('mh-recalculate-status').append(' ',reload);}catch(e){el('mh-recalculate-status').textContent=e.message;}finally{el('mh-recalculate').disabled=false;}};
 await loadNotebook();chooseRange(300e3);chooseTime(m.plasma_times.find(t=>t>=dayStart+18*3600)||m.fit_times[0],{plasma:true});setView(activeView);
 if(m.calibration)el('mh-density-note').textContent=`Plasma-line calibrated Ne · gain ${m.calibration.gain.toExponential(5)} · revision ${m.calibration.revision}. Source Ti, Te, velocity and their fit remain unchanged. Raw 10 µs lag remains a signed source proxy.`;
 el('millstone-status').textContent=`${m.antenna} · ${m.fits.length} fitted integrations · ${m.raw_times.length} native ACF profiles · ${m.plasma_times.length} plasma-line integrations`;
 parent.postMessage({type:'ruisdap-ready'},location.origin);
}catch(e){el('millstone-status').textContent='Unable to load Millstone measurements: '+e.message;console.error(e);}
function clip(v,a,b){return Math.max(a,Math.min(b,v));}
