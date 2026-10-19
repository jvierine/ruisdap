// Fetch only the temporal resolution needed by the current screen window.
// Tiles contain arithmetic means in physical density units, never log averages.
export class BeataPower {
 constructor(url){this.url=url;this.cache=new Map();this.generation=0;this.manifest=fetch(url,{cache:'no-store'}).then(r=>{if(!r.ok)throw Error('BEATA power manifest '+r.status);return r.json();});}
 async tile(url){if(this.cache.has(url)){const v=this.cache.get(url);this.cache.delete(url);this.cache.set(url,v);return v;}const r=await fetch(url);if(!r.ok)throw Error('BEATA power tile '+r.status);const v=new Float32Array(await r.arrayBuffer());this.cache.set(url,v);while(this.cache.size>12)this.cache.delete(this.cache.keys().next().value);return v;}
 async field(spec,pixels,apply){
  const generation=++this.generation,m=await this.manifest,[start,end]=spec.bounds;
  const target=Math.max(256,Math.min(1600,pixels)),wanted=(end-start)/target;
  const level=m.levels.filter(l=>l.cadence_s<=wanted).at(-1)||m.levels[0];
  const tiles=level.tiles.filter(t=>t.end>start&&t.start<end);
  if(!tiles.length)return;
  const loaded=await Promise.all(tiles.map(async t=>[t,await this.tile(new URL(t.url,new URL(this.url,location.href)).href)]));
  if(generation!==this.generation)return;
  const cols=[],xe=[];
  for(const [t,data]of loaded)for(let j=0;j<t.columns;j++){
   const lo=Math.max(start,t.x_edges[j]),hi=Math.min(end,t.x_edges[j+1]);if(hi-lo<=.02)continue;
   if(!xe.length){xe.push(start);if(lo-start>.02){cols.push(null);xe.push(lo);}}else if(lo>xe.at(-1)+.02){cols.push(null);xe.push(lo);}
   cols.push([data,j,t.columns]);xe.push(hi);
  }
  if(!cols.length)return;if(end-xe.at(-1)>.02){cols.push(null);xe.push(end);}else xe[xe.length-1]=end;
  const y=m.altitude_km,ye=[y[0]-(y[1]-y[0])/2,...y.slice(1).map((v,i)=>(v+y[i])/2),y.at(-1)+(y.at(-1)-y.at(-2))/2];
  const values=new Float32Array(cols.length*y.length);values.fill(NaN);const gain=(spec.magic_constant??m.definition.base_magic_constant)/m.definition.base_magic_constant;
  for(let r=0;r<y.length;r++)for(let j=0;j<cols.length;j++){const col=cols[j];if(col)values[r*cols.length+j]=col[0][r*col[2]+col[1]]*gain;}
  apply(values,cols.length,y.length,{...spec,x_edges:new Float32Array(xe.map(t=>t-Math.floor(start/86400)*86400)),dataOriginX:Math.floor(start/86400)*86400,y_edges:new Float32Array(ye),powerLod:level.factor,powerCadence:level.cadence_s});
 }
}
