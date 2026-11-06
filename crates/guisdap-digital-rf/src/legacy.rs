//! Digital RF 1.0 compound-i16 recordings (Arecibo 2016 archive).
//! File-local reads preserve absolute sample indices and indexed gaps.
use anyhow::{ensure, Result};
use hdf5::{File, Dataset, H5Type};
use num_complex::Complex32;
use std::{path::Path,ops::Range};
#[repr(C)]
#[derive(Clone,Copy)]
struct IQ { r:i16, i:i16 }
unsafe impl H5Type for IQ {
 fn type_descriptor()->hdf5::types::TypeDescriptor {
  use hdf5::types::{CompoundField,CompoundType,TypeDescriptor};
  TypeDescriptor::Compound(CompoundType{fields:vec![CompoundField::typed::<i16>("r",0,0),CompoundField::typed::<i16>("i",2,1)],size:4})
 }
}
pub struct LegacyFile { _file:File, data:Dataset, pub sample_rate_hz:u64, index:Vec<(u64,usize)> }
pub struct LegacyBlock { pub start_sample:u64, pub samples:Vec<Complex32> }
impl LegacyFile {
 pub fn open(path:impl AsRef<Path>)->Result<Self>{
  let file=File::open(path)?;let data=file.dataset("rf_data")?;
  ensure!(data.shape().len()==2 && data.shape()[1]==1,"expected one complex channel");
  ensure!(data.dtype()?.to_descriptor()?==IQ::type_descriptor(),"unsupported legacy sample type (expected compound i16 r/i)");
  let sr=data.attr("sample_rate")?.read_raw::<f64>()?;
  ensure!(sr.len()==1 && sr[0].is_finite() && sr[0]>0. && sr[0].fract()==0. && sr[0]<1e12,"invalid legacy sample rate");
  let idx=file.dataset("rf_data_index")?.read_2d::<u64>()?;
  ensure!(idx.ncols()==2 && idx.nrows()>0,"invalid legacy RF index");
  let mut index=vec![];
  for row in idx.rows(){index.push((row[0],usize::try_from(row[1])?));}
  ensure!(index[0].1==0 && index.iter().all(|x|x.1<data.shape()[0]),"invalid legacy offsets");
  for w in index.windows(2){ensure!(w[1].1>w[0].1 && w[0].0.checked_add((w[1].1-w[0].1)as u64).is_some_and(|end|end<=w[1].0),"overlapping legacy index");}
  Ok(Self{_file:file,data,sample_rate_hz:sr[0]as u64,index})
 }
 pub fn rows(&self)->usize{self.data.shape()[0]}
 pub fn read_rows(&self,rows:Range<usize>)->Result<Vec<LegacyBlock>>{
  ensure!(rows.start<=rows.end && rows.end<=self.rows(),"legacy read outside file");
  let mut out=vec![];
  for (k,&(absolute,offset)) in self.index.iter().enumerate(){
   let end=self.index.get(k+1).map_or(self.rows(),|x|x.1);
   let a=offset.max(rows.start);let b=end.min(rows.end);if a>=b{continue;}
   let raw=self.data.read_slice_1d::<IQ,_>((a..b,0))?;
   out.push(LegacyBlock{start_sample:absolute.checked_add((a-offset)as u64).ok_or_else(||anyhow::anyhow!("absolute sample overflow"))?,samples:raw.iter().map(|x|Complex32::new(x.r as f32,x.i as f32)).collect()});
  }
  Ok(out)
 }
}
#[cfg(test)]
mod tests {
 use super::*;
 #[test]fn gaps_signed_iq_and_boundaries()->Result<()>{
  let dir=tempfile::tempdir()?;let path=dir.path().join("v1.h5");
  {
   let f=File::create(&path)?;
   let d=f.new_dataset::<IQ>().shape([4,1]).create("rf_data")?;
   d.write_raw(&[IQ{r:-32768,i:32767},IQ{r:3,i:-4},IQ{r:5,i:6},IQ{r:7,i:8}])?;
   d.new_attr::<f64>().shape([1]).create("sample_rate")?.write_raw(&[25e6])?;
   f.new_dataset::<u64>().shape([2,2]).create("rf_data_index")?.write_raw(&[100,0,200,2])?;
  }
  let f=LegacyFile::open(&path)?;let b=f.read_rows(0..4)?;
  assert_eq!(b.len(),2);assert_eq!(b[0].start_sample,100);assert_eq!(b[1].start_sample,200);
  assert_eq!(b[0].samples[0],Complex32::new(-32768.,32767.));
  let b=f.read_rows(1..3)?;assert_eq!(b[0].start_sample,101);assert_eq!(b[1].start_sample,200);
  assert!(f.read_rows(0..5).is_err());Ok(())
 }
}
