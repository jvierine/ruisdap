// Shared instrument axes and units for every ISR dataset.
export function instrumentSpec(records,bounds,options={}){
 const axes=[['azimuth_deg','Azimuth / °','#236b8e'],['elevation_deg','Elevation / °','#814aa1'],['system_temperature_k','Tsys / K','#b46c19'],['transmitter_w','Ptx / MW','#be3838']];
 const yaxes=[],series=[];
 for(const [key,label,color]of axes){
  const scale=key==='transmitter_w'?1e-6:1,visible=records.filter(r=>r.end>bounds[0]&&r.start<bounds[1]);
  const vals=visible.map(r=>r[key]).filter(Number.isFinite).map(v=>v*scale),lo=vals.length?Math.min(...vals):0,hi=vals.length?Math.max(...vals):1,pad=Math.max((hi-lo)*.08,Math.abs(hi)*.01,.01),limits=[lo-pad,hi+pad];
  yaxes.push({label,color,limits});let previous=null;const xy=[];
  for(const r of visible){if(previous!==null&&r.start-previous>1)xy.push(null,null);const v=r[key];if(Number.isFinite(v)){const y=(v*scale-limits[0])/(limits[1]-limits[0]);xy.push(Math.max(bounds[0],r.start),y,Math.min(bounds[1],r.end),y);}else xy.push(null,null);previous=r.end;}
  series.push({color,xy});
 }
 if(Number.isFinite(options.selected))series.push({color:'#986327',xy:[options.selected,0,options.selected,1]});
 return {title:'Instrument diagnostics',height:340,margins:[160,32,180,58],bounds:[...bounds,0,1],xlabel:'Time / UTC',xformat:'utc',yformat:'none',yaxes,series,...options};
}
