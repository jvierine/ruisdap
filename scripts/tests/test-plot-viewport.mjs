import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import init,{PlotViewport} from '../../web/pkg/guisdap_web.js';
import {setViewportBounds} from '../../web/plot-viewport.js';
await init({module_or_path:await readFile(new URL('../../web/pkg/guisdap_web_bg.wasm',import.meta.url))});
const source=new PlotViewport(),raw=new PlotViewport(),map=new PlotViewport();
for(let i=0;i<40;i++){
 source.zoom(.5,.5,.5);
 const b=[...source.bounds()];
 setViewportBounds(raw,b);
 setViewportBounds(map,[b[0],0,b[2],1]);
 setViewportBounds(source,[raw.bounds()[0],b[1],raw.bounds()[2],b[3]]);
 const actual=source.bounds();
 assert.ok(actual[2]-actual[0]<.51,'Repeated zoom reset the time window');
 assert.ok(actual[3]-actual[1]<.51,'Repeated zoom reset the range window');
 assert.ok(Math.abs(raw.bounds()[0]-actual[0])<1e-8);
 assert.ok(Math.abs(map.bounds()[2]-actual[2])<1e-8);
}
const tiny=[.5,.4,.50000000000001,.40000000000001];
setViewportBounds(source,tiny);
for(let i=0;i<100;i++)setViewportBounds(source,[...source.bounds()]);
assert.ok(source.bounds()[2]-source.bounds()[0]<.002,'Synchronization reset the minimum time window');
assert.ok(source.bounds()[3]-source.bounds()[1]<.002,'Synchronization reset the minimum range window');
for(const bounds of [[0,0,1,1],[-.1,.2,.2,.8],[.99999,.99999,1,1]]){
 setViewportBounds(source,bounds);
 assert.ok([...source.bounds()].every(v=>Number.isFinite(v)&&v>=0&&v<=1));
}
source.free();raw.free();map.free();
console.log('Rust viewport: repeated zoom and linked time/range synchronization retain the zoom floor.');
