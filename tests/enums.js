import path from "path";

import { readFile, readJSONFile } from "./utils";

export const TOLERANCE = 1e-3;

export const FIXTURES_DIR = path.resolve(__dirname, "./fixtures");
export const Si = readJSONFile(path.join(FIXTURES_DIR, "Si.json"));
export const Na4Cl4 = readJSONFile(path.join(FIXTURES_DIR, "Na4Cl4.json"));
export const Na4Cl4Cartesian = readJSONFile(path.join(FIXTURES_DIR, "Na4Cl4-cartesian.json"));
export const Na4Cl4Poscar = readFile(path.join(FIXTURES_DIR, "Na4Cl4.poscar"));
export const SiSupercell = readJSONFile(path.join(FIXTURES_DIR, "Si-supercell.json"));
export const C2H4 = readJSONFile(path.join(FIXTURES_DIR, "C2H4.json"));
export const C2H4Translated = readJSONFile(path.join(FIXTURES_DIR, "C2H4-translated.json"));
export const Na = readJSONFile(path.join(FIXTURES_DIR, "Na.json"));

export const OSiBasis = readJSONFile(path.join(FIXTURES_DIR, "OSi-basis.json"));
export const Si2Basis = readJSONFile(path.join(FIXTURES_DIR, "Si2-basis.json"));
export const Ge2Basis = readJSONFile(path.join(FIXTURES_DIR, "Ge2-basis.json"));
export const AsGeBasis = readJSONFile(path.join(FIXTURES_DIR, "AsGe-basis.json"));
export const FeLiSiBasis = readJSONFile(path.join(FIXTURES_DIR, "FeLiSi-basis.json"));
export const LiFeSiBasis = readJSONFile(path.join(FIXTURES_DIR, "LiFeSi-basis.json"));
export const atomicConstraints = readJSONFile(path.join(FIXTURES_DIR, "atomic-constraints.json"));
export const Si2BasisRepeated = readJSONFile(path.join(FIXTURES_DIR, "Si2-basis-repeated.json"));
export const H2HInitial = readJSONFile(path.join(FIXTURES_DIR, "H2+H-initial.json"));
export const H2HFinal = readJSONFile(path.join(FIXTURES_DIR, "H2+H-final.json"));
export const H2HImage = readJSONFile(path.join(FIXTURES_DIR, "H2+H-image.json"));
export const SiSlab = readJSONFile(path.join(FIXTURES_DIR, "Si-slab.json"));
export const SiPWSCFInput = readFile(path.join(FIXTURES_DIR, "Si-pwscf.in"));

export const Zr1H23Zr1H1 = readJSONFile(path.join(FIXTURES_DIR, "Zr1H23Zr1H1.json"));
export const Zr1H23Zr1H1Poscar = readFile(path.join(FIXTURES_DIR, "Zr1H23Zr1H1.poscar"));
export const H2O = readFile(path.join(FIXTURES_DIR, "H2O.poscar"));

export const Graphene = readJSONFile(path.join(FIXTURES_DIR, "Graphene.json"));
export const GraphenePoscar = readFile(path.join(FIXTURES_DIR, "Graphene.poscar"));
export const GraphenePWSCFInput = readFile(path.join(FIXTURES_DIR, "Graphene-pwscf.in"));
export const NiHex = readJSONFile(path.join(FIXTURES_DIR, "Ni-hex.json"));
export const NiHexPoscar = readFile(path.join(FIXTURES_DIR, "Ni-hex.poscar"));
export const NiCub = readJSONFile(path.join(FIXTURES_DIR, "Ni-cub.json"));
export const NiCubPWSCFInput = readFile(path.join(FIXTURES_DIR, "Ni-cub-pwscf.in"));
export const Sb2S3Orc = readJSONFile(path.join(FIXTURES_DIR, "Sb2S3-ORC.json"));
export const Sb2S3OrcPWSCFInput = readFile(path.join(FIXTURES_DIR, "Sb2S3-ORC-pwscf.in"));
export const BNHex = readJSONFile(path.join(FIXTURES_DIR, "BN-hex.json"));
export const BNHexPWSCFInput = readFile(path.join(FIXTURES_DIR, "BN-hex-pwscf.in"));

export const FortranFile1 = readFile(path.join(FIXTURES_DIR, "fortran-file-1.txt"));
export const FortranFile1JSON = readJSONFile(path.join(FIXTURES_DIR, "fortran-file-1.json"));
export const FortranFileInvalid = readFile(path.join(FIXTURES_DIR, "fortran-file-invalid.txt"));
export const FortranFileNoCards = readFile(path.join(FIXTURES_DIR, "fortran-file-no-cards.txt"));
