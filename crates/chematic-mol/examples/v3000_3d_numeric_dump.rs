//! Regenerated 3D V3000 output from both native conformer writers.
use chematic_core::{Coords3D, Point3};
use chematic_mol::{
    parse_mol_v3000_with_coords_3d, write_mol_v3000_with_conformer,
    write_mol_v3000_with_conformer_checked,
};
use serde_json::{Value, json};
fn main() {
    let path = std::env::args().nth(1).expect("fixture path");
    let packet: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let rows:Vec<_>=packet["cases"].as_array().unwrap().iter().map(|c| {
  let (mol,meta,coords)=parse_mol_v3000_with_coords_3d(c["input"].as_str().unwrap()).unwrap();
  let conformer=Coords3D{points:coords.into_iter().map(|(x,y,z)|Point3::new(x,y,z)).collect()};
  json!({"id":c["id"],"plain":write_mol_v3000_with_conformer(&mol,&meta,&conformer),"checked":write_mol_v3000_with_conformer_checked(&mol,&meta,&conformer).unwrap()})
 }).collect();
    println!("{}", json!({"rows":rows}));
}
