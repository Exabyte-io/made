import path from "path";
import {readFile, readJSONFile} from "./utils";

export const FIXTURES_DIR = path.resolve(__dirname, "./fixtures");
export const Si = readJSONFile(path.join(FIXTURES_DIR, "Si.json"));
export const Na4Cl4 = readJSONFile(path.join(FIXTURES_DIR, "Na4Cl4.json"));
export const Na4Cl4Poscar = readFile(path.join(FIXTURES_DIR, "Na4Cl4.poscar"));
export const SiSupercell = readJSONFile(path.join(FIXTURES_DIR, "Si-supercell.json"));

export const OSiBasis = readJSONFile(path.join(FIXTURES_DIR, "OSi-basis.json"));
export const Si2Basis = readJSONFile(path.join(FIXTURES_DIR, "Si2-basis.json"));
export const Ge2Basis = readJSONFile(path.join(FIXTURES_DIR, "Ge2-basis.json"));
export const AsGeBasis = readJSONFile(path.join(FIXTURES_DIR, "AsGe-basis.json"));
export const FeLiSiBasis = readJSONFile(path.join(FIXTURES_DIR, "FeLiSi-basis.json"));
export const LiFeSiBasis = readJSONFile(path.join(FIXTURES_DIR, "LiFeSi-basis.json"));
