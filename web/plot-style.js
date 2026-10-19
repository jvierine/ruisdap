// GUISDAP vizu.m defaults and its 78-colour myb palette (white endpoint cut).
// https://github.com/ingemarh/guisdap ; units are physical, never proxy units.
export const EISCAT_LIMITS = [[1e10,1e12],[0,3000],[0,4000],[-200,200],[1e9,1e12]];
export const EISCAT_COLORMAP = (()=>{
 const stops=[[0,0,0],[0,0,.5],[0,0,1],[0,.5,.5],[0,1,0],[.5,1,0],[1,1,0],[1,.5,0],[1,0,0],[1,0,.5],[1,0,1],[1,.5,1]];
 const knots=stops.map((_,i)=>Math.round(i*77/11)),rgb=new Float32Array(78*3);
 for(let i=0,j=0;i<78;i++){
  while(j<knots.length-2&&i>knots[j+1])j++;
  const t=(i-knots[j])/(knots[j+1]-knots[j]);
  for(let c=0;c<3;c++)rgb[i*3+c]=Math.tanh(stops[j][c]*(1-t)+stops[j+1][c]*t)/Math.tanh(1);
 }
 return rgb;
})();
