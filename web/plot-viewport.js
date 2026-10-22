export function setViewportBounds(viewport,bounds){
 // Rust's zoom_box rejects spans below 0.001. Round outward before reset,
 // otherwise floating-point roundoff at the zoom floor resets the whole view.
 const b=[...bounds],minimum=.001000001;
 for(const i of [0,1]){const size=Math.max(minimum,Math.min(1,b[i+2]-b[i])),center=(b[i]+b[i+2])/2;b[i]=Math.max(0,Math.min(1-size,center-size/2));b[i+2]=b[i]+size;}
 const current=viewport.bounds();if(b.every((v,i)=>Math.abs(v-current[i])<1e-10))return;
 viewport.reset();viewport.zoom_box(...b);
}
