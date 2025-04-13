import path from "path";

import { readFile } from "./utils";

// TODO: use @mat3ra/standata
export { default as Silicon } from "../fixtures/si-standata.json";

export const TOLERANCE = 1e-3;

export { default as FeO } from "../fixtures/FeO.json";
export { default as Na4Cl4 } from "../fixtures/Na4Cl4.json";
export { default as Na4Cl4Cartesian } from "../fixtures/Na4Cl4-cartesian.json";
export { default as SiSupercell } from "../fixtures/Si-supercell.json";
export { default as C2H4 } from "../fixtures/C2H4.json";
export { default as C2H4Translated } from "../fixtures/C2H4-translated.json";
export { default as Na } from "../fixtures/Na.json";
export { default as OSiBasis } from "../fixtures/OSi-basis.json";
export { default as Si2Basis } from "../fixtures/Si2-basis.json";
export { default as Ge2Basis } from "../fixtures/Ge2-basis.json";
export { default as AsGeBasis } from "../fixtures/AsGe-basis.json";
export { default as FeLiSiBasis } from "../fixtures/FeLiSi-basis.json";
export { default as LiFeSiBasis } from "../fixtures/LiFeSi-basis.json";
export { default as atomicConstraints } from "../fixtures/atomic-constraints.json";
export { default as Si2BasisRepeated } from "../fixtures/Si2-basis-repeated.json";
export { default as H2HInitial } from "../fixtures/H2+H-initial.json";
export { default as H2HFinal } from "../fixtures/H2+H-final.json";
export { default as H2HImage } from "../fixtures/H2+H-image.json";
export { default as SiSlab } from "../fixtures/Si-slab.json";
export { default as Zr1H23Zr1H1 } from "../fixtures/Zr1H23Zr1H1.json";
export { default as Graphene } from "../fixtures/Graphene.json";
export { default as NiHex } from "../fixtures/Ni-hex.json";

// NON-JSON
export const FIXTURES_DIR = path.resolve(__dirname, "../fixtures");

// NON-JSON
export const Na4Cl4Poscar = readFile(path.join(FIXTURES_DIR, "Na4Cl4.poscar"));
export const SiPWSCFInput = readFile(path.join(FIXTURES_DIR, "Si-pwscf.in"));
export const Zr1H23Zr1H1Poscar = readFile(path.join(FIXTURES_DIR, "Zr1H23Zr1H1.poscar"));
export const H2OPoscar = readFile(path.join(FIXTURES_DIR, "H2O.poscar"));
export const GraphenePoscar = readFile(path.join(FIXTURES_DIR, "Graphene.poscar"));
export const NiHexPoscar = readFile(path.join(FIXTURES_DIR, "Ni-hex.poscar"));
