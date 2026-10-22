//! Independent peak-power forward equation versus the unit-energy LPI design.
//! Receiver spectral filtering remains the fitter's frequency response, not a
//! scalar baud correction. This test isolates pulse energy and range measure.
use guisdap_millstone::{lpi,signal,C};
#[test]
fn thirty_and_hundred_microsecond_codes_and_mixed_sequence_recover_same_density_scale(){
 let truth=[C::new(2.,-0.3),C::new(4.,0.5),C::new(3.,-0.2)];
 for decimation in [10usize,30,100] {
 for sequence in [vec![30;8],vec![100;8],vec![30,100,30,100,30,100,30,100]] {
  let mut rows=Vec::new();let mut measurements=Vec::new();
  for (code,baud)in sequence.into_iter().enumerate(){
   let dr=decimation as f64*299792.458/(2.*1e6); // one fixed range basis, including mixed codes
   let n=1800;let mut peak=vec![C::default();n];
   for bit in 0..16{let sign=if ((bit*13+code*7)^(bit>>1)^code).count_ones()%2==0{1.}else{-1.};for j in 0..baud{peak[bit*baud+j]=C::new(sign,0.);}}
   let energy=peak.iter().map(|z|z.norm_sqr()).sum::<f64>();
   let unit:Vec<_>=peak.iter().map(|z|*z/energy.sqrt()).collect();
   let ambiguity=signal::lag_products(&unit,10,decimation).unwrap();
   for m in 2..ambiguity.len(){
    // Forward voltage covariance in physical peak-power units. Integrate
    // directly over individual measurement samples and range basis cells.
    let mut y=C::default();let mut row=Vec::new();
    for r in 0..truth.len(){
     let delay=(r+2)*decimation;
     for t in m*decimation..(m+1)*decimation{if t>=delay&&t+10-delay<n{y+=truth[r]*dr*peak[t-delay]*peak[t+10-delay].conj();}}
     if m>=r+2{let a=ambiguity[m-r-2]*energy*dr;if a.norm_sqr()>1e-30{row.push((r,a));}}
    }
    if !row.is_empty(){rows.push(row);measurements.push(y);}
   }
  }
  let (recovered,_,_)=lpi::solve(&rows,&measurements,&vec![1.;rows.len()],3).unwrap();
  for (got,want)in recovered.iter().zip(truth){assert!((*got-want).norm()<1e-10,"{got} != {want}");}
 }
 }
}
