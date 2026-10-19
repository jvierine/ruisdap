// Shared cell boundaries and robust display limits; input samples are unchanged.
export function timeCells(times,values,ny,intervals,gap=3600){
 const steps=times.slice(1).map((v,i)=>v-times[i]).filter(v=>v>0).sort((a,b)=>a-b),dt=steps[Math.floor(steps.length/2)]||1;
 const bounds=[],cols=[];
 for(let i=0;i<times.length;i++){
  let lo=intervals?.[i]?.[0]??times[i]-dt/2,hi=intervals?.[i]?.[1]??times[i]+dt/2;
  if(i&&lo-bounds.at(-1)<gap){lo=(times[i-1]+times[i])/2;bounds[bounds.length-1]=lo;}
  else if(i&&lo>bounds.at(-1)){cols.push(-1);bounds.push(lo);}
  if(!bounds.length)bounds.push(lo);cols.push(i);bounds.push(hi);
 }
 const data=new Float32Array(cols.length*ny);data.fill(NaN);
 for(let r=0;r<ny;r++)for(let i=0;i<cols.length;i++)if(cols[i]>=0)data[r*cols.length+i]=values[r*times.length+cols[i]];
 return [new Float32Array(bounds),data];
}
export function edges(x){return new Float32Array(x.length===1?[x[0]-.5,x[0]+.5]:[x[0]-(x[1]-x[0])/2,...x.slice(1).map((v,i)=>(v+x[i])/2),x.at(-1)+(x.at(-1)-x.at(-2))/2]);}
export function iqr(values,log=false){
 const v=[],stride=Math.max(1,Math.floor(values.length/50000));for(let i=0;i<values.length;i+=stride)if(Number.isFinite(values[i])&&(!log||values[i]>0))v.push(log?Math.log10(values[i]):values[i]);v.sort((a,b)=>a-b);
 if(!v.length)return log?[1,10]:[0,1];const q=p=>v[Math.floor(p*(v.length-1))],d=Math.max(q(.75)-q(.25),Math.abs(q(.75))*1e-6,1e-12),a=Math.max(q(.01),q(.25)-1.5*d),b=Math.max(a+d*1e-6,Math.min(q(.99),q(.75)+1.5*d));return log?[10**a,10**b]:[a,b];
}
