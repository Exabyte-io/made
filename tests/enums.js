import path from "path";
import {readFile, readJSONFile} from "./utils";

export const FIXTURES_DIR = path.resolve(__dirname, "./fixtures");
export const Si = readJSONFile(path.join(FIXTURES_DIR, "Si.json"));
export const Na4Cl4 = readJSONFile(path.join(FIXTURES_DIR, "Na4Cl4.json"));
export const Na4Cl4Poscar = readFile(path.join(FIXTURES_DIR, "Na4Cl4.poscar"));
export const SiSupercell = readJSONFile(path.join(FIXTURES_DIR, "Si-supercell.json"));
export const FeLiSiBasis = readJSONFile(path.join(FIXTURES_DIR, "FeLiSi-basis.json"));
export const LiFeSiBasis = readJSONFile(path.join(FIXTURES_DIR, "LiFeSi-basis.json"));
