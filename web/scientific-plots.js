import {FieldView,SeriesView,PlotViewport} from './pkg/guisdap_web.js';
import {bindPlotNavigation} from './vendor/rustplotlib/navigation.js';
const charts=new Map(),curves=new Map();
const clip=(v,a,b)=>Math.max(a,Math.min(b,v));
import {timeCells,edges,iqr} from './plot-data.js';
export {timeCells,edges,iqr};
export function fieldPanel(parent,id,title,x,y,values,options={}){
 const old=document.getElementById(id);const article=old?.closest('article')||document.createElement('article');article.className='plot parameter-row scientific-map';
 if(!old)article.innerHTML=`<h3>${title}</h3><div class="ion-raster"><canvas id="${id}" aria-label="${title}"></canvas><div class="ion-axes" id="${id}-axes"></div><div class="plot-markers"></div></div><div class="shared-colour-controls"><label>Colour minimum <input type="range" aria-label="${title} colour minimum"><output></output></label><label>Colour maximum <input type="range" aria-label="${title} colour maximum"><output></output></label><button>Auto IQR</button></div>`;
 if(!old)document.getElementById(parent).append(article);
 let xe,filled;if(options.time!==false)[xe,filled]=timeCells(x,values,y.length,options.intervals);else{xe=edges(x);filled=values;}
 const ye=edges(y),initial=options.limits||iqr(values,options.mode===1),spec={...options,title,bounds:[xe[0],xe.at(-1),ye[0],ye.at(-1)],limits:initial,x_edges:xe,y_edges:ye};
 const item=renderField(id,filled,xe.length-1,ye.length-1,spec);
 return {article,item,view:item.view,viewport:item.viewport,redraw:()=>item.redraw(),xe,markers:points=>item.markers(points),setCursor:(x,y)=>item.setCursor(x,y)};
}
export function renderField(id,values,nx,ny,spec){
 let p=charts.get(id);
 if(!p){
  const canvas=document.getElementById(id),box=canvas.parentElement,axesHolder=box.querySelector('.ion-axes')||box.appendChild(document.createElement('div'));
  box.classList.add('ion-raster');axesHolder.classList.add('ion-axes');axesHolder.replaceChildren();
  const axes=document.createElement('canvas');axes.className='scientific-axes';axesHolder.append(axes);
  const bar=document.createElement('canvas');bar.className='scientific-bar';box.append(bar);
  const barAxes=document.createElement('canvas');barAxes.className='scientific-bar-axes';box.append(barAxes);
  const markerLayer=box.querySelector('.plot-markers')||box.appendChild(document.createElement('div'));markerLayer.className='plot-markers';
  const cursor=document.createElement('div');cursor.className='scientific-cursor';box.append(cursor);
  p={canvas,box,view:new FieldView(canvas),axes:new SeriesView(axes),axesCanvas:axes,bar:new FieldView(bar),barCanvas:bar,barAxes:new SeriesView(barAxes),barAxesCanvas:barAxes,viewport:new PlotViewport(),spec,markerLayer,cursor,points:[],selected:null,queued:false};charts.set(id,p);
  p.setCursor=(x,y)=>{p.selected=[x,y];p.redraw();};
  p.markers=points=>{p.points=points;p.redraw();};
  p.redraw=()=>{if(p.queued)return;p.queued=true;requestAnimationFrame(()=>{p.queued=false;draw(p);});};
  if(!spec.external_navigation)bindPlotNavigation(canvas,p.viewport,p.redraw);
  let down=null,moved=false;canvas.addEventListener('pointerdown',e=>{down=[e.clientX,e.clientY];moved=false;});canvas.addEventListener('pointermove',e=>{if(down&&Math.hypot(e.clientX-down[0],e.clientY-down[1])>4)moved=true;});
  canvas.addEventListener('click',e=>{if(moved)return;const r=canvas.getBoundingClientRect(),b=visibleBounds(p);p.spec.pick?.(b[0]+(e.clientX-r.left)/r.width*(b[1]-b[0]),b[3]-(e.clientY-r.top)/r.height*(b[3]-b[2]));});
  new ResizeObserver(p.redraw).observe(box);
 }
 for(const layer of p.box.querySelectorAll('[id^=cal-broad-],#cal-power'))Object.assign(layer.style,{inset:'auto',left:'76px',top:'22px',width:'calc(100% - 208px)',height:'calc(100% - 80px)'});
 p.values=values.slice();p.nx=nx;p.ny=ny;p.spec=spec;
 upload(p);installControls(p);p.redraw();return p;
}
function upload(p){
 const s=p.spec,normalized=s.normalized,mode=s.mode||0,palette=s.palette??1;
 if(s.x_edges&&s.y_edges)p.view.set_rectilinear_data(s.x_edges,s.y_edges,p.values);else p.view.set_data(p.nx,p.ny,p.values);
 p.view.set_scale(normalized?0:s.limits[0],normalized?1:s.limits[1],normalized?0:mode,s.normWidth||Math.max(Math.abs(s.limits[1]-s.limits[0])/10,1e-20),palette);
 p.bar.set_data(1,256,new Float32Array(Array.from({length:256},(_,i)=>i/255)));p.bar.set_scale(0,1,0,1,palette);
}
function installControls(p){
 const controls=p.box.parentElement.querySelector('.shared-colour-controls');if(!controls)return;
 const log=p.spec.mode===1,range=iqr(p.values,log),inputs=controls.querySelectorAll('input'),outputs=controls.querySelectorAll('output');
 for(let i=0;i<2;i++){let lo=log?Math.log10(Math.max(range[0]*.01,1e-30)):Math.min(range[0],p.spec.limits[0]),hi=log?Math.log10(range[1]*100):Math.max(range[1],p.spec.limits[1]);if(!log){const d=hi-lo||1;lo-=2*d;hi+=4*d;}inputs[i].min=lo;inputs[i].max=hi;inputs[i].step=(hi-lo)/10000;inputs[i].value=log?Math.log10(p.spec.limits[i]):p.spec.limits[i];outputs[i].textContent=p.spec.limits[i].toPrecision(4);
  inputs[i].oninput=()=>{const v=log?10**Number(inputs[i].value):Number(inputs[i].value),limits=[...p.spec.limits];limits[i]=v;if(limits[1]<=limits[0])return;if(p.spec.palette===0&&p.spec.symmetric)limits[1-i]=-v;p.spec.limits=limits;p.spec.onLimits?.(...limits);upload(p);installControls(p);p.redraw();};}
 controls.querySelector('button').onclick=()=>{p.spec.limits=range;p.spec.onLimits?.(...range);upload(p);installControls(p);p.redraw();};
}
function visibleBounds(p){const a=p.spec.bounds,b=p.viewport.bounds();return [a[0]+(a[1]-a[0])*b[0],a[0]+(a[1]-a[0])*b[2],a[2]+(a[3]-a[2])*b[1],a[2]+(a[3]-a[2])*b[3]];}
function draw(p){
 if(!p.box.clientWidth||p.box.closest('[hidden]'))return;
 const started=performance.now(),s=p.spec,w=p.box.clientWidth,h=p.box.clientHeight,dpr=devicePixelRatio||1,b=visibleBounds(p),margins=[76,22,132,58];
 p.axes.axes({width:w,height:h,dpr,bounds:b,margins,x:{label:s.time===false?(s.xlabel||'Frequency offset / MHz'):'Time / UTC',format:s.time===false?'':'utc'},y:{label:s.ylabel||'Slant range / km',format:''}});
 const [l,t,pw,ph]=p.axes.axes_rect();Object.assign(p.canvas.style,{left:l+'px',top:t+'px',width:pw+'px',height:ph+'px'});if(p.canvas.width!==Math.round(pw*dpr))p.canvas.width=Math.round(pw*dpr);if(p.canvas.height!==Math.round(ph*dpr))p.canvas.height=Math.round(ph*dpr);p.view.set_view(...p.viewport.bounds());
 const barLeft=w-108;for(const c of [p.barCanvas,p.barAxesCanvas]){c.style.left=barLeft+'px';c.style.top='0';}
 Object.assign(p.barCanvas.style,{top:t+'px',width:'16px',height:ph+'px'});if(p.barCanvas.width!==Math.round(16*dpr))p.barCanvas.width=Math.round(16*dpr);if(p.barCanvas.height!==Math.round(ph*dpr))p.barCanvas.height=Math.round(ph*dpr);p.bar.draw();
 p.barAxesCanvas.style.width='108px';p.barAxesCanvas.style.height=h+'px';
 const log=s.mode===1;p.barAxes.axes({width:108,height:h,dpr,bounds:[0,1,log?Math.log10(s.limits[0]):s.limits[0],log?Math.log10(s.limits[1]):s.limits[1]],margins:[0,t,92,58],right_y:true,no_frame:true,x:{format:'none'},y:{format:log?'log10':''}});
 p.barCanvas.style.zIndex=2;p.barAxesCanvas.style.zIndex=1;
 p.markerLayer.replaceChildren();for(const point of p.points){const x=l+(point.x-b[0])/(b[1]-b[0])*pw,y=t+(b[3]-point.y)/(b[3]-b[2])*ph;if(x<l||x>l+pw||y<t||y>t+ph)continue;const marker=document.createElement('button');marker.className='cal-marker';marker.style.left=x+'px';marker.style.top=y+'px';marker.type='button';marker.dataset.included=String(!!point.included);marker.title=(point.included?'Included':'Excluded')+' · '+(point.label||'Remove calibration point');marker.setAttribute('aria-label',point.label||'Remove calibration point');marker.onclick=e=>{e.stopPropagation();point.remove?.();};p.markerLayer.append(marker);}
 const c=p.selected;p.cursor.hidden=!c||c[0]<b[0]||c[0]>b[1];if(c){Object.assign(p.cursor.style,{left:(l+(c[0]-b[0])/(b[1]-b[0])*pw)+'px',top:t+'px',height:ph+'px'});}
 p.canvas.dataset.drawMs=(performance.now()-started).toFixed(3);p.canvas.dataset.renderer='rustplotlib';
}
export function curve(id,spec){
 const host=document.getElementById(id);if(!host)return;let p=curves.get(id);
 if(!p){host.replaceChildren();const canvas=document.createElement('canvas');canvas.setAttribute('aria-label',spec.title||id);canvas.className='scientific-curve';host.append(canvas);const note=document.createElement('p');note.className='caption';host.append(note);p={host,canvas,note,view:new SeriesView(canvas),viewport:new PlotViewport(),pending:false};curves.set(id,p);new ResizeObserver(()=>curve(id,p.spec)).observe(host);bindPlotNavigation(canvas,p.viewport,()=>curve(id,p.spec),()=>{const r=p.view.axes_rect();return {x:r[0],y:r[1],width:r[2],height:r[3]};});let down=null,moved=false;canvas.addEventListener('pointerdown',e=>{down=[e.clientX,e.clientY];moved=false;});canvas.addEventListener('pointermove',e=>{if(down&&Math.hypot(e.clientX-down[0],e.clientY-down[1])>4)moved=true;});canvas.addEventListener('click',e=>{if(moved||!p.spec.pick)return;const r=canvas.getBoundingClientRect(),a=p.view.axes_rect(),b=p.spec.bounds,v=p.viewport.bounds(),x=(e.clientX-r.left-a[0])/a[2],y=1-(e.clientY-r.top-a[1])/a[3];p.spec.pick(b[0]+(b[1]-b[0])*(v[0]+x*(v[2]-v[0])),b[2]+(b[3]-b[2])*(v[1]+y*(v[3]-v[1])));});}
 p.spec=spec;p.canvas.style.height=(spec.height||300)+'px';if(!host.clientWidth||host.closest('[hidden]')||p.pending)return;p.pending=true;requestAnimationFrame(()=>{p.pending=false;drawCurve(p);});
}
function drawCurve(p){const spec=p.spec,host=p.host,b=spec.bounds,v=p.viewport.bounds(),bounds=[b[0]+(b[1]-b[0])*v[0],b[0]+(b[1]-b[0])*v[2],b[2]+(b[3]-b[2])*v[1],b[2]+(b[3]-b[2])*v[3]];
 p.view.axes({width:host.clientWidth,height:spec.height||300,dpr:devicePixelRatio||1,bounds,x:{label:spec.xlabel||'',format:spec.xformat||''},y:{label:spec.ylabel||'',format:spec.yformat||''}});
 for(const b of spec.bands||[])p.view.band(new Float64Array(b.values.map(v=>v??NaN)),b.color);
 for(const series of spec.series||[]){const xy=new Float64Array(series.xy.map(v=>v??NaN));if(series.segments)p.view.segments(xy,series.color||'#236b8e');else if(series.style!==undefined)p.view.markers(xy,series.color||'#236b8e',series.radius||3,series.style);else p.view.series(xy,series.color||'#236b8e',series.line!==false,series.radius||0);}
 p.note.textContent=spec.caption||'';p.canvas.dataset.renderer='rustplotlib';
}
window.ruisdapField=(id,values,nx,ny,spec)=>{
 const match=id.match(/^fit-raster-(\d)$/);
 if(match)spec.onLimits=(lo,hi)=>{for(const [edge,value]of [[0,lo],[1,hi]])document.getElementById(`fit-scale-${match[1]}-${edge}`).value=value;document.getElementById(`fit-scale-${match[1]}-0`).dispatchEvent(new Event('input',{bubbles:true}));};
 return renderField(id,values,nx,ny,spec);
};
window.ruisdapCurve=curve;
