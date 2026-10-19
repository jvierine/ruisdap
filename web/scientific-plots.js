import {instrumentSpec} from './instrument-plot.js';
import {BeataPower} from './beata-power.js';
import {FieldView,SeriesView,PlotViewport} from './pkg/guisdap_web.js';
import {bindPlotNavigation} from './vendor/rustplotlib/navigation.js';
import {EISCAT_COLORMAP} from './plot-style.js';
import {PLASMA_COLORMAP} from './plasma-colormap.js';
import {mathAsset,mathLabel} from './math-labels.js';
const charts=new Map(),curves=new Map();
// Keep the common time-axis edges valid in a narrow embedded viewer as well
// as on a full desktop. Shared maps, debris and diagnostics use the same edges.
function fitMargins(width,margins){
 if(width>=900)return margins;
 const [left,top,right,bottom]=margins;
 return [Math.min(left,Math.max(52,width*160/900)),top,Math.min(right,Math.max(120,width*180/900)),bottom];
}
let openColours=null;
document.addEventListener('keydown',e=>{if(e.key==='Escape'&&openColours){openColours.hidden=true;openColours.parentElement.querySelector('.scientific-colour-hit')?.setAttribute('aria-expanded','false');openColours=null;}});
const clip=(v,a,b)=>Math.max(a,Math.min(b,v));
import {timeCells,edges,iqr} from './plot-data.js';
export {timeCells,edges,iqr};
export function fieldPanel(parent,id,title,x,y,values,options={}){
 const old=document.getElementById(id);const article=old?.closest('article')||document.createElement('article');article.className='plot parameter-row scientific-map';
 if(!old)article.innerHTML=`<h3>${title}</h3><div class="ion-raster"><canvas id="${id}" aria-label="${title}"></canvas><div class="ion-axes" id="${id}-axes"></div><div class="plot-markers"></div></div><div class="shared-colour-controls" hidden><label>Colour minimum <input type="range" aria-label="${title} colour minimum"><output></output></label><label>Colour maximum <input type="range" aria-label="${title} colour maximum"><output></output></label><button class="auto-iqr">Auto IQR</button><button class="standard-limits">Standard limits</button><button class="close-colours" aria-label="Close colour limits">Close</button></div>`;
 if(!old)document.getElementById(parent).append(article);
 article.querySelector('h3').textContent=title;
 article.classList.toggle('shared-time-row',!!options.sharedTime);article.classList.toggle('shared-time-last',!!options.sharedTimeLast);
 let xe,filled;if(options.time!==false)[xe,filled]=timeCells(x,values,y.length,options.intervals);else{xe=edges(x);filled=values;}
 const ye=edges(y),initial=options.limits||iqr(values,options.mode===1),spec={timeGroup:options.time===false?null:'dataset-time',...options,title,bounds:[...(options.timeBounds||[xe[0],xe.at(-1)]),ye[0],ye.at(-1)],limits:initial,defaultLimits:options.limits,x_edges:xe,y_edges:ye};
 const item=renderField(id,filled,xe.length-1,ye.length-1,spec);
 return {article,item,view:item.view,viewport:item.viewport,redraw:()=>item.redraw(),xe,update:(data,settings={})=>fieldPanel(parent,id,settings.title||item.spec.title,x,y,data,{...options,...settings}),markers:points=>item.markers(points),setCursor:(x,y)=>item.setCursor(x,y)};
}
export function renderField(id,values,nx,ny,spec){
 let p=charts.get(id);
 if(!p){
  const canvas=document.getElementById(id),box=canvas.parentElement,axesHolder=box.querySelector('.ion-axes')||box.appendChild(document.createElement('div'));
  canvas.style.zIndex='1';box.classList.add('ion-raster');axesHolder.classList.add('ion-axes');axesHolder.replaceChildren();
  const axes=document.createElement('canvas');axes.className='scientific-axes';axesHolder.append(axes);
  const bar=document.createElement('canvas');bar.className='scientific-bar';box.append(bar);
  const barAxes=document.createElement('canvas');barAxes.className='scientific-bar-axes';box.append(barAxes);
  const markerLayer=box.querySelector('.plot-markers')||box.appendChild(document.createElement('div'));markerLayer.className='plot-markers';
  const cursor=document.createElement('div');cursor.className='scientific-cursor';box.append(cursor);
  const rangeCursor=document.createElement('div');rangeCursor.className='scientific-range-cursor';box.append(rangeCursor);
  p={canvas,box,view:new FieldView(canvas),axes:new SeriesView(axes),axesCanvas:axes,bar:new FieldView(bar),barCanvas:bar,barAxes:new SeriesView(barAxes),barAxesCanvas:barAxes,viewport:new PlotViewport(),spec,markerLayer,cursor,rangeCursor,points:[],selected:null,queued:false};charts.set(id,p);
  p.setCursor=(x,y)=>{p.selected=[x,y];p.redraw();};
  p.markers=points=>{p.points=points;p.redraw();};
  p.redraw=()=>{if(p.queued)return;p.queued=true;requestAnimationFrame(()=>{p.queued=false;draw(p);});};
  if(!spec.external_navigation)bindPlotNavigation(canvas,p.viewport,()=>{syncTime(p);p.spec.onView?.();p.redraw();});
  let down=null,moved=false;canvas.addEventListener('pointerdown',e=>{down=[e.clientX,e.clientY];moved=false;});canvas.addEventListener('pointermove',e=>{if(down&&Math.hypot(e.clientX-down[0],e.clientY-down[1])>4)moved=true;});
  canvas.addEventListener('click',e=>{if(moved)return;const r=canvas.getBoundingClientRect(),b=visibleBounds(p);p.spec.pick?.(b[0]+(e.clientX-r.left)/r.width*(b[1]-b[0]),b[3]-(e.clientY-r.top)/r.height*(b[3]-b[2]));});
  new ResizeObserver(p.redraw).observe(box);
 }
 for(const layer of p.box.querySelectorAll('[id^=cal-broad-],#cal-power'))Object.assign(layer.style,{inset:'auto',left:'76px',top:'22px',width:'calc(100% - 208px)',height:'calc(100% - 80px)'});
 p.values=values.slice();p.nx=nx;p.ny=ny;p.spec=spec;
 p.box.parentElement.classList.toggle('shared-time-row',!!spec.sharedTime);p.box.parentElement.classList.toggle('shared-time-last',!!spec.sharedTimeLast);
 upload(p);installControls(p);p.redraw();return p;
}
function upload(p){
 const s=p.spec,normalized=s.normalized,mode=s.mode||0,palette=s.palette??1;
 if(s.x_edges&&s.y_edges)p.view.set_rectilinear_data(s.x_edges,s.y_edges,p.values);else p.view.set_data(p.nx,p.ny,p.values);
 p.view.set_scale(normalized?0:s.limits[0],normalized?1:s.limits[1],normalized?0:mode,s.normWidth||Math.max(Math.abs(s.limits[1]-s.limits[0])/10,1e-20),palette);
 p.bar.set_data(1,256,new Float32Array(Array.from({length:256},(_,i)=>i/255)));p.bar.set_scale(0,1,0,1,palette);
 if(s.colormap!=='seismic'){const rgb=s.colormap==='plasma'?PLASMA_COLORMAP:EISCAT_COLORMAP;p.view.set_colormap(rgb);p.bar.set_colormap(rgb);}
}
function installControls(p){
 const controls=p.box.parentElement.querySelector('.shared-colour-controls');if(!controls)return;
 if(!p.colourButton){
  const button=document.createElement('button');button.type='button';button.className='scientific-colour-hit';button.setAttribute('aria-label',`${p.spec.title} colour limits`);button.title='Adjust colour limits';controls.id=p.canvas.id+'-colour-controls';button.setAttribute('aria-controls',controls.id);button.setAttribute('aria-expanded','false');
  button.onclick=()=>{const show=controls.hidden;if(openColours){openColours.hidden=true;const previous=openColours.parentElement.querySelector('.scientific-colour-hit');previous?.setAttribute('aria-expanded','false');}controls.hidden=!show;openColours=show?controls:null;button.setAttribute('aria-expanded',String(show));};
  controls.querySelector('.close-colours').onclick=()=>{controls.hidden=true;openColours=null;button.setAttribute('aria-expanded','false');button.focus();};p.box.append(button);p.colourButton=button;
 }
 const log=p.spec.mode===1,range=iqr(p.values,log),inputs=controls.querySelectorAll('input'),outputs=controls.querySelectorAll('output');
 for(let i=0;i<2;i++){inputs[i].setAttribute('aria-label',`${p.spec.title} colour ${i?'maximum':'minimum'}`);let lo=log?Math.log10(Math.max(range[0]*.01,1e-30)):Math.min(range[0],p.spec.limits[0]),hi=log?Math.log10(range[1]*100):Math.max(range[1],p.spec.limits[1]);if(!log){const d=hi-lo||1;lo-=2*d;hi+=4*d;}inputs[i].min=lo;inputs[i].max=hi;inputs[i].step=(hi-lo)/10000;inputs[i].value=log?Math.log10(p.spec.limits[i]):p.spec.limits[i];outputs[i].textContent=p.spec.limits[i].toPrecision(4);
  inputs[i].oninput=()=>{const v=log?10**Number(inputs[i].value):Number(inputs[i].value),limits=[...p.spec.limits];limits[i]=v;if(limits[1]<=limits[0])return;if(p.spec.palette===0&&p.spec.symmetric)limits[1-i]=-v;p.spec.limits=limits;p.spec.onLimits?.(...limits);upload(p);installControls(p);p.redraw();};}
 controls.querySelector('.auto-iqr').onclick=()=>{p.spec.limits=range;p.spec.onLimits?.(...range);upload(p);installControls(p);p.redraw();};
 const standard=controls.querySelector('.standard-limits');standard.hidden=!p.spec.defaultLimits;standard.onclick=()=>{p.spec.limits=[...p.spec.defaultLimits];upload(p);installControls(p);p.redraw();};
}
function syncTime(source){
 const group=source.spec.timeGroup||source.box.closest('.parameter-stack');if(!group)return;
 const b=visibleBounds(source);
 for(const target of charts.values()){
  if(target===source||(target.spec.timeGroup||target.box.closest('.parameter-stack'))!==group)continue;
  const v=target.viewport.bounds(),a=target.spec.bounds,span=a[1]-a[0];
  target.viewport.reset();target.viewport.zoom_box((b[0]-a[0])/span,v[1],(b[1]-a[0])/span,v[3]);target.redraw();
 }
}
function visibleBounds(p){const a=p.spec.bounds,b=p.viewport.bounds();return [a[0]+(a[1]-a[0])*b[0],a[0]+(a[1]-a[0])*b[2],a[2]+(a[3]-a[2])*b[1],a[2]+(a[3]-a[2])*b[3]];}
function draw(p){
 if(!p.box.clientWidth||p.box.closest('[hidden]'))return;
 const started=performance.now(),s=p.spec,w=p.box.clientWidth,h=p.box.clientHeight,dpr=devicePixelRatio||1,b=visibleBounds(p),bottom=s.sharedTime&&!s.sharedTimeLast?10:58,top=s.sharedTime?10:22,margins=fitMargins(w,[s.sharedTime?160:76,top,s.sharedTime?180:132,bottom]);
 p.axes.axes({width:w,height:h,dpr,bounds:b,margins,x:{label:s.time===false?(s.xlabel||'Frequency offset / MHz'):'Time / UTC',format:s.sharedTime&&!s.sharedTimeLast?'none':s.time===false?'':'utc'},y:{label:s.ylabel||'Slant range / km',format:''}});
 const [l,t,pw,ph]=p.axes.axes_rect();for(const overlay of p.box.querySelectorAll('canvas[id$="-quality"]'))Object.assign(overlay.style,{left:l+'px',top:t+'px',width:pw+'px',height:ph+'px'});Object.assign(p.canvas.style,{left:l+'px',top:t+'px',width:pw+'px',height:ph+'px'});if(p.canvas.width!==Math.round(pw*dpr))p.canvas.width=Math.round(pw*dpr);if(p.canvas.height!==Math.round(ph*dpr))p.canvas.height=Math.round(ph*dpr);if(s.x_edges&&s.y_edges){const x=s.x_edges,y=s.y_edges;p.view.set_view((b[0]-(s.dataOriginX||0)-x[0])/(x.at(-1)-x[0]),(b[2]-y[0])/(y.at(-1)-y[0]),(b[1]-(s.dataOriginX||0)-x[0])/(x.at(-1)-x[0]),(b[3]-y[0])/(y.at(-1)-y[0]));}else p.view.set_view(...p.viewport.bounds());
 const barLeft=w-margins[2]+24;for(const c of [p.barCanvas,p.barAxesCanvas]){c.style.left=barLeft+'px';c.style.top='0';}
 Object.assign(p.barCanvas.style,{top:t+'px',width:'16px',height:ph+'px'});if(p.barCanvas.width!==Math.round(16*dpr))p.barCanvas.width=Math.round(16*dpr);if(p.barCanvas.height!==Math.round(ph*dpr))p.barCanvas.height=Math.round(ph*dpr);p.bar.draw();
 p.barAxesCanvas.style.width='108px';p.barAxesCanvas.style.height=h+'px';
 const log=s.mode===1;p.barAxes.axes({width:108,height:h,dpr,bounds:[0,1,log?Math.log10(s.limits[0]):s.limits[0],log?Math.log10(s.limits[1]):s.limits[1]],margins:[0,t,92,bottom],right_y:true,no_frame:true,x:{format:'none'},y:{label:mathAsset(s.colorbarTitle||s.title)?'':s.colorbarTitle||s.title||'',format:log?'log10':''}});
 mathLabel(p,'bar',mathAsset(s.colorbarTitle||s.title),s.colorbarTitle||s.title,barLeft+93,t+ph/2,true);
 p.canvas.setAttribute('aria-label',s.title||p.canvas.id);
 if(p.colourButton){p.colourButton.setAttribute('aria-label',`${s.title} colour limits`);}
 if(p.colourButton)Object.assign(p.colourButton.style,{left:barLeft+'px',top:t+'px',width:'108px',height:ph+'px'});
 p.canvas.dataset.colormap=s.colormap||'eiscat';p.canvas.dataset.limits=JSON.stringify(s.limits);p.canvas.dataset.timeAxis=String(!(s.sharedTime&&!s.sharedTimeLast));
 p.barCanvas.style.zIndex=2;p.barAxesCanvas.style.zIndex=1;
 p.markerLayer.replaceChildren();for(const point of p.points){const x=l+(point.x-b[0])/(b[1]-b[0])*pw,y=t+(b[3]-point.y)/(b[3]-b[2])*ph;if(x<l||x>l+pw||y<t||y>t+ph)continue;const marker=document.createElement('button');marker.className='cal-marker';marker.style.left=x+'px';marker.style.top=y+'px';marker.innerHTML='<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M2 12h6m8 0h6M12 2v6m0 8v6" class="marker-outline"/><path d="M2 12h6m8 0h6M12 2v6m0 8v6" class="marker-ink"/></svg>';if(point.number){const label=document.createElement('span');label.className='cal-marker-label';label.textContent=point.number;marker.append(label);}marker.type='button';marker.dataset.included=String(!!point.included);marker.title=(point.included?'Included':'Excluded')+' · '+(point.label||'Remove calibration point');marker.setAttribute('aria-label',point.label||'Remove calibration point');marker.onclick=e=>{e.stopPropagation();point.remove?.();};p.markerLayer.append(marker);}
 const c=p.selected;p.rangeCursor.hidden=!c||c[1]<b[2]||c[1]>b[3];if(c)Object.assign(p.rangeCursor.style,{left:l+'px',top:(t+(b[3]-c[1])/(b[3]-b[2])*ph)+'px',width:pw+'px'});p.cursor.hidden=!c||c[0]<b[0]||c[0]>b[1];if(c){Object.assign(p.cursor.style,{left:(l+(c[0]-b[0])/(b[1]-b[0])*pw)+'px',top:t+'px',height:ph+'px'});}
 p.canvas.dataset.visibleBounds=JSON.stringify(b);p.canvas.dataset.axesRect=JSON.stringify([l,t,pw,ph]);
 p.canvas.dataset.drawMs=(performance.now()-started).toFixed(3);p.canvas.dataset.renderer='rustplotlib';
}
export function curve(id,spec){
 if(spec.instrument_records)spec=instrumentSpec(spec.instrument_records,spec.bounds.slice(0,2),{selected:spec.selected});
 const host=document.getElementById(id);if(!host)return;host.style.position='relative';let p=curves.get(id);
 if(!p){host.replaceChildren();const canvas=document.createElement('canvas');canvas.setAttribute('aria-label',spec.title||id);canvas.className='scientific-curve';host.append(canvas);const note=document.createElement('p');note.className='caption';host.append(note);p={host,canvas,note,view:new SeriesView(canvas),viewport:new PlotViewport(),pending:false};curves.set(id,p);new ResizeObserver(()=>curve(id,p.spec)).observe(host);bindPlotNavigation(canvas,p.viewport,()=>{const a=p.spec.bounds,v=p.viewport.bounds();p.spec.onView?.([a[0]+(a[1]-a[0])*v[0],a[0]+(a[1]-a[0])*v[2]]);curve(id,p.spec);},()=>{const r=p.view.axes_rect();return {x:r[0],y:r[1],width:r[2],height:r[3]};});let down=null,moved=false;canvas.addEventListener('pointerdown',e=>{down=[e.clientX,e.clientY];moved=false;});canvas.addEventListener('pointermove',e=>{if(down&&Math.hypot(e.clientX-down[0],e.clientY-down[1])>4)moved=true;});canvas.addEventListener('click',e=>{if(moved||!p.spec.pick)return;const r=canvas.getBoundingClientRect(),a=p.view.axes_rect(),b=p.spec.bounds,v=p.viewport.bounds(),x=(e.clientX-r.left-a[0])/a[2],y=1-(e.clientY-r.top-a[1])/a[3];p.spec.pick(b[0]+(b[1]-b[0])*(v[0]+x*(v[2]-v[0])),b[2]+(b[3]-b[2])*(v[1]+y*(v[3]-v[1])));});}
 if(/^fit-instrument-/.test(id))spec.margins=[160,32,180,58];
 p.spec=spec;if(spec.resetViewport)p.viewport.reset();p.canvas.style.height=(spec.height||(host.closest('.plasma-profile-pair')?520:300))+'px';if(!host.clientWidth||host.closest('[hidden]')||p.pending)return;p.pending=true;requestAnimationFrame(()=>{p.pending=false;drawCurve(p);});
}
function drawCurve(p){const spec=p.spec,host=p.host,b=spec.bounds,v=p.viewport.bounds(),bounds=[b[0]+(b[1]-b[0])*v[0],b[0]+(b[1]-b[0])*v[2],b[2]+(b[3]-b[2])*v[1],b[2]+(b[3]-b[2])*v[3]];
 p.view.axes({width:host.clientWidth,height:spec.height||(host.closest('.plasma-profile-pair')?520:300),dpr:devicePixelRatio||1,bounds,...(spec.margins?{margins:fitMargins(host.clientWidth,spec.margins)}:{}),x:{label:mathAsset(spec.xlabel)?'':spec.xlabel||'',format:spec.xformat||''},y:{label:mathAsset(spec.ylabel)?'':spec.ylabel||'',format:spec.yformat||''}});
 const rect=p.view.axes_rect(),height=spec.height||(host.closest('.plasma-profile-pair')?520:300);
 mathLabel(p,'x',mathAsset(spec.xlabel),spec.xlabel,rect[0]+rect[2]/2,height-17);
 mathLabel(p,'y',mathAsset(spec.ylabel),spec.ylabel,18,rect[1]+rect[3]/2,true);
 if(spec.yaxes){
  p.extraAxes??=[];host.style.position='relative';
  const legend=p.legend||host.appendChild(document.createElement('div'));p.legend=legend;legend.className='multi-axis-legend';legend.replaceChildren();Object.assign(legend.style,{left:rect[0]+'px',top:'0px',width:rect[2]+'px',zIndex:'2',rowGap:'0px',lineHeight:'16px'});
  spec.yaxes.forEach((axis,i)=>{
   let a=p.extraAxes[i];if(!a){const canvas=document.createElement('canvas');canvas.className='scientific-extra-axis';host.append(canvas);a={canvas,view:new SeriesView(canvas)};p.extraAxes[i]=a;}
   const right=i>=2,offset=right?(i-2)*80:-i*80,left=right?rect[0]+rect[2]+offset-16:rect[0]+offset-80;
   Object.assign(a.canvas.style,{left:left+'px',top:'0',width:'100px',height:height+'px'});
   a.view.axes({width:100,height,dpr:devicePixelRatio||1,bounds:[0,1,...axis.limits],margins:right?[0,rect[1],84,height-rect[1]-rect[3]]:[80,rect[1],10,height-rect[1]-rect[3]],right_y:right,no_frame:true,x:{format:'none'},y:{label:mathAsset(axis.label)?'':axis.label}});
   a.view.series(new Float64Array([right?1:0,axis.limits[0],right?1:0,axis.limits[1]]),axis.color,true,0);
   mathLabel(p,'axis-'+i,mathAsset(axis.label),axis.label,right?left+85:left+18,rect[1]+rect[3]/2,true);
   const label=document.createElement('span');label.style.color=axis.color;label.textContent=axis.label;legend.append(label);
  });
 }
 for(const b of spec.bands||[])p.view.band(new Float64Array(b.values.map(v=>v??NaN)),b.color);
 for(const series of spec.series||[]){const xy=new Float64Array(series.xy.map(v=>v??NaN));if(series.segments)p.view.segments(xy,series.color||'#236b8e');else if(series.style!==undefined)p.view.markers(xy,series.color||'#236b8e',series.radius||3,series.style);else p.view.series(xy,series.color||'#236b8e',series.line!==false,series.radius||0);}
 p.canvas.dataset.visibleBounds=JSON.stringify(bounds);p.canvas.dataset.axesRect=JSON.stringify(rect);p.note.textContent=spec.caption||'';p.canvas.dataset.renderer='rustplotlib';
}
let beataPower=null;
window.ruisdapField=(id,values,nx,ny,spec)=>{
 const match=id.match(/^fit-raster-(\d)$/);
 if(match){spec.colormap=match[1]==='3'?'seismic':match[1]==='4'?'plasma':'eiscat';spec.sharedTime=true;spec.sharedTimeLast=match[1]==='4';}
 if(match)spec.onLimits=(lo,hi)=>{for(const [edge,value]of [[0,lo],[1,hi]])document.getElementById(`fit-scale-${match[1]}-${edge}`).value=value;document.getElementById(`fit-scale-${match[1]}-0`).dispatchEvent(new Event('input',{bubbles:true}));};
 const p=renderField(id,values,nx,ny,spec);
 if(id==='fit-raster-4'&&spec.power_tiles_url){
  if(!beataPower||beataPower.url!==spec.power_tiles_url)beataPower=new BeataPower(spec.power_tiles_url);
  beataPower.field(spec,p.canvas.clientWidth,(data,nx,ny,s)=>{renderField(id,data,nx,ny,s);p.canvas.dataset.powerLod=s.powerLod;p.canvas.dataset.powerCadence=s.powerCadence;const note=document.getElementById('beata-power-detail');if(note)note.textContent=s.powerLod===1?'Native power profiles · approximately 0.357 s per profile':'Power averaged to '+s.powerCadence.toFixed(2)+' s · zoom in for native profiles';}).catch(e=>{const note=document.getElementById('beata-power-detail');if(note)note.textContent='High-cadence power unavailable: '+(e.message||String(e));console.error(e);});
 }
 return p;
};
window.ruisdapCurve=curve;
