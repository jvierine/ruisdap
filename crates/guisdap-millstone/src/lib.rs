//! Millstone processing kernels. Voltage cancellation precedes filtering and
//! lag products; failed fits become explicit masks, never repaired voltages.
pub mod debris;
pub mod fitting;
pub mod lpi;
pub mod product;
pub mod signal;
pub use num_complex::Complex64 as C;

pub mod archive;
