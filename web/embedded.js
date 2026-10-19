// Embedded views retain their WASM state, GPU textures and selected integration.
if(window.parent!==window){
 document.documentElement.classList.add('dataset-embedded');
 let pending=false;
 const height=()=>{if(pending)return;pending=true;requestAnimationFrame(()=>{pending=false;parent.postMessage({type:'ruisdap-height',height:document.body.getBoundingClientRect().height+2},location.origin);});};
 new ResizeObserver(height).observe(document.body);height();
 window.addEventListener('message',event=>{
  if(event.origin!==location.origin||event.source!==parent)return;
  if(event.data?.type==='ruisdap-view'&&document.getElementById('observation-plasma-panel')){
   const ion=event.data.view==='ion';
   document.getElementById('observation-plasma-panel').parentElement.classList.toggle('ion-only',ion);
   for(const id of ['observation-calibration','observation-band-control','observation-display-control'])document.getElementById(id).hidden=ion;
   document.getElementById('observation-heading').textContent=ion?'Ion-line profiles':'Plasma-line measurements and calibration';
   height();
  }
 });
}
