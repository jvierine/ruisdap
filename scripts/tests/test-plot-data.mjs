import assert from 'node:assert/strict';
import {timeCells,edges,iqr} from '../../web/plot-data.js';
const times=[0,120,240],values=new Float32Array([1,2,3,10,20,30]);
const [boundary,data]=timeCells(times,values,2,[[-30,30],[90,150],[210,270]]);
assert.deepEqual([...boundary],[-30,60,180,270]);
assert.deepEqual([...data],[1,2,3,10,20,30]);
assert.deepEqual(times,[0,120,240]);assert.deepEqual([...values],[1,2,3,10,20,30]);
const [gap,missing]=timeCells([0,7200],new Float32Array([5,10]),1,[[-30,30],[7170,7230]]);
assert.deepEqual([...gap],[-30,30,7170,7230]);assert.equal(missing[0],5);assert.ok(Number.isNaN(missing[1]));assert.equal(missing[2],10);
const [hour,hourData]=timeCells([30,3690],new Float32Array([1,2]),1,[[0,60],[3660,3720]]);
assert.equal(hour.length,4);assert.ok(Number.isNaN(hourData[1]));
assert.deepEqual([...edges([1,3,7])],[0,2,5,9]);
for(const log of [false,true]){
 for(const samples of [[],[NaN,-Infinity],[4,4,4]]){
  const limits=iqr(samples,log);assert.ok(limits.every(Number.isFinite));
  assert.ok(limits[1]>limits[0]);if(log)assert.ok(limits[0]>0);
 }
 const sample=Array.from({length:1000},(_,i)=>100+i%100);sample.push(1e20);
 const limits=iqr(sample,log);assert.ok(limits[0]>=100&&limits[1]<=199);
}
console.log('Plot data: short gaps extended, hour gaps blank, samples and ordering preserved.');
