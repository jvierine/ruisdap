// One dataset owns its two scientific views. Loaded views stay alive when hidden.
const el=id=>document.getElementById(id), query=new URLSearchParams(location.search);
const buttons=[...document.querySelectorAll('[data-view]')], frames=new Map();
let routes={}, active='plasma', canonical=false;
async function json(url){const r=await fetch(url);if(!r.ok)throw Error(`Could not load ${url}: ${r.status}`);return r.json();}
function frame(kind,url){
 const node=document.createElement('iframe');node.className='dataset-frame';node.title=kind==='ion'?'Ion-line measurements':'Plasma-line measurements and calibration';
 node.src=url;node.addEventListener('load',()=>{node.contentWindow.postMessage({type:'ruisdap-view',view:active},location.origin);if(!new URL(url,location.href).pathname.endsWith('/plasma.html')){node.dataset.ready='true';if(!node.hidden)el('dataset-status').textContent='';}});
 el('dataset-views').append(node);frames.set(kind,node);return node;
}
function show(kind){
 active=kind;const url=new URL(location.href);url.searchParams.set('view',kind);history.replaceState(null,'',url);
 for(const b of buttons)b.setAttribute('aria-pressed',String(b.dataset.view===kind));
 for(const f of frames.values())f.hidden=true;
 if(canonical){const f=frames.get('paired')||frame('paired',routes.plasma);f.hidden=false;f.contentWindow.postMessage({type:'ruisdap-view',view:kind},location.origin);el('dataset-status').textContent='';return;}
 if(!routes[kind]){el('dataset-status').textContent=`No ${kind}-line product has been published for this dataset.`;return;}
 const f=frames.get(kind)||frame(kind,routes[kind]);f.hidden=false;el('dataset-status').textContent=f.dataset.ready?'':`Loading ${kind}-line measurements…`;
}
for(const b of buttons)b.addEventListener('click',()=>show(b.dataset.view));
window.addEventListener('message',event=>{
 if(event.origin!==location.origin)return;
 const target=[...frames.values()].find(f=>f.contentWindow===event.source);if(!target)return;
 if(event.data?.type==='ruisdap-ready'){target.dataset.ready='true';document.body.dataset.readyMs=performance.now().toFixed(1);if(!target.hidden)el('dataset-status').textContent='';}
 if(event.data?.type==='ruisdap-error'&&!target.hidden)el('dataset-status').textContent='Unable to load measurements: '+event.data.message;
 if(event.data?.type==='ruisdap-height')target.style.height=Math.max(700,Math.min(100000,event.data.height))+'px';
 if(event.data?.type==='ruisdap-analysis'&&canonical){const u=new URL(location.href);for(const k of ['analysis','time','range_m'])u.searchParams.set(k,event.data[k]);u.searchParams.delete('start');u.searchParams.delete('end');u.searchParams.set('view',active);location.assign(u);return;}
 if(event.data?.type==='ruisdap-time')for(const f of frames.values())if(f!==target)f.contentWindow.postMessage(event.data,location.origin);
});
try{
 if(query.has('millstone')){
  const antenna=query.get('millstone');if(!['ZENITH','MISA'].includes(antenna))throw Error('Unknown Millstone antenna');
  canonical=true;routes.plasma='millstone.html?antenna='+antenna+'&date='+encodeURIComponent(query.get('date')||'2024-04-08');
  el('dataset-title').textContent='Millstone Hill · '+antenna;el('dataset-context').textContent=(query.get('date')||'2024-04-08')+' · UTC';
 }else if(query.has('observation')){
  const config=await json('observations-config.json'), records=[];
  for(const index of config.indexes){const root=index.slice(0,index.lastIndexOf('/')+1);for(const r of(await json(index)).records)records.push({...r,root});}
  const seed=records.find(r=>r.id===query.get('observation'));if(!seed)throw Error('Observation not found');
  const date=new Date(seed.start_unix_s*1000).toISOString().slice(0,10);
  const paired=records.filter(r=>r.station===seed.station&&r.antenna===seed.antenna&&new Date(r.start_unix_s*1000).toISOString().startsWith(date));
  const plasma=paired.find(r=>r.has_spectra)||seed;
  canonical=true;routes.plasma=seed.station==="Millstone Hill" ? "millstone.html?antenna="+encodeURIComponent(seed.antenna)+"&date="+date : "observations.html?record="+encodeURIComponent(plasma.id);
  el('dataset-title').textContent=seed.station+' · '+seed.antenna;el('dataset-context').textContent=date+' · UTC';
 }else{
  const inventory=await json('catalog/catalog.json'), e=inventory.experiments.find(e=>e.id===query.get('experiment'));
  if(!e)throw Error('Experiment not found in the current inventory');
  el('dataset-title').textContent=e.experiment;el('dataset-context').textContent=e.date+' · '+e.radar+' · UTC';
  for(const kind of ['plasma','ion']){
   const p=e.products.find(p=>p.complete&&p.label.toLowerCase().startsWith(kind+'-line'));
   if(p){const link=p.link.startsWith('./?')?p.link.replace('./?','plasma.html?'):p.link;const u=new URL(link,location.href);u.searchParams.delete('original');u.searchParams.set('experiment',e.id);routes[kind]=u.pathname+u.search;}
  }
  if(!Object.keys(routes).length){el('dataset-processing').hidden=false;el('dataset-process').disabled=!e.supported;el('dataset-status').textContent=e.supported?'A validated processing recipe is available.':'No validated processing recipe is configured for this experiment and date.';}
  el('dataset-process').addEventListener('click',async()=>{
   try{const r=await fetch(`api/catalog/jobs/${e.id}`,{method:'POST',headers:{Authorization:'Bearer '+el('dataset-key').value.trim()}});if(!r.ok)throw Error(await r.text());el('dataset-process-status').textContent='Queued. See Processing status.';}catch(e){el('dataset-process-status').textContent=e.message;}
  });
 }
 if(canonical&&routes.plasma?.startsWith('millstone.html')){const u=new URL(routes.plasma,location.href);for(const key of ['analysis','time','range_m','start','end'])if(query.has(key))u.searchParams.set(key,query.get(key));routes.plasma=u.pathname+u.search;}
 const requestedView=(query.get('view')||'').replace(/\.+$/,'');show(['ion','plasma'].includes(requestedView)?requestedView:(routes.plasma?'plasma':'ion'));
}catch(e){el('dataset-status').textContent=e.message;el('dataset-title').textContent='Dataset unavailable';}
