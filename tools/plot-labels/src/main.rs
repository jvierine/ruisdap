// Pre-render once with Rustplotlib's native mathpazo LaTeX typesetter.
// Browser plots reuse these glyphs without a TeX runtime or per-frame typesetting.
use anyhow::Result;
use std::{fs,path::PathBuf};
fn main()->Result<()>{
 let out=PathBuf::from(std::env::args().nth(1).expect("output directory"));fs::create_dir_all(&out)?;
 let labels=[("ne",r"$n_e\,/\,\mathrm{m}^{-3}$"),("raw-ne",r"$\mathrm{Raw}\ n_e\,/\,\mathrm{m}^{-3}$"),("ti",r"$T_i\,/\,\mathrm{K}$"),("te",r"$T_e\,/\,\mathrm{K}$"),("vi",r"$v_i\,/\,(\mathrm{m}\,\mathrm{s}^{-1})$"),("tsys",r"$T_\mathrm{sys}\,/\,\mathrm{K}$"),("ptx",r"$P_\mathrm{tx}\,/\,\mathrm{MW}$"),("ne-proxy",r"$n_e\ \mathrm{proxy}$")];
 let latex=rustplotlib::Latex{cache:out.join(".latex-cache"),font_pt:18.,..Default::default()};
 let glyphs=latex.typeset(labels.iter().map(|(_,label)|label.to_string()))?;
 for (name,label) in labels{let g=&glyphs[label];let v=g.view_box;fs::write(out.join(format!("{name}.svg")),format!("<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"{}pt\" height=\"{}pt\" viewBox=\"{} {} {} {}\" fill=\"#314858\">{}</svg>",g.width_pt,g.height_pt,v[0],v[1],v[2],v[3],g.body))?;}
 Ok(())
}
