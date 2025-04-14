import fs from "fs";
import path from "path";

// TODO: use @mat3ra/standata

import Silicon from "../fixtures/si-standata.json";
import FeO from "../fixtures/FeO.json";
import Na4Cl4 from "../fixtures/Na4Cl4.json";
import Na4Cl4Cartesian from "../fixtures/Na4Cl4-cartesian.json";
import SiSupercell from "../fixtures/Si-supercell.json";
import C2H4 from "../fixtures/C2H4.json";
import C2H4Translated from "../fixtures/C2H4-translated.json";
import Na from "../fixtures/Na.json";
import OSiBasis from "../fixtures/OSi-basis.json";
import Si2Basis from "../fixtures/Si2-basis.json";
import Ge2Basis from "../fixtures/Ge2-basis.json";
import AsGeBasis from "../fixtures/AsGe-basis.json";
import FeLiSiBasis from "../fixtures/FeLiSi-basis.json";
import LiFeSiBasis from "../fixtures/LiFeSi-basis.json";
import atomicConstraints from "../fixtures/atomic-constraints.json";
import Si2BasisRepeated from "../fixtures/Si2-basis-repeated.json";
import H2HInitial from "../fixtures/H2+H-initial.json";
import H2HFinal from "../fixtures/H2+H-final.json";
import H2HImage from "../fixtures/H2+H-image.json";
import SiSlab from "../fixtures/Si-slab.json";
import Zr1H23Zr1H1 from "../fixtures/Zr1H23Zr1H1.json";
import Graphene from "../fixtures/Graphene.json";
import NiHex from "../fixtures/Ni-hex.json";
import SiSlab100 from "../fixtures/si-slab-100.json";
import SiSlab111 from "../fixtures/si-slab-111-0.5-vacuum-ratio.json";
import SiSlab111NoVacuum from "../fixtures/si-slab-111-0-vacuum.json";


export {
    // Silicon
    Silicon,
    SiSupercell,
    SiSlab,
    SiSlab100,
    SiSlab111,
    SiSlab111NoVacuum,
    Si2Basis,
    Si2BasisRepeated,
    // Other
    FeO,
    Na4Cl4,
    Na4Cl4Cartesian,
    C2H4,
    C2H4Translated,
    Na,
    OSiBasis,
    Ge2Basis,
    AsGeBasis,
    FeLiSiBasis,
    LiFeSiBasis,
    atomicConstraints,
    H2HInitial,
    H2HFinal,
    H2HImage,
    Zr1H23Zr1H1,
    Graphene,
    NiHex,
};

// NON-JSON
export const FIXTURES_DIR = path.resolve(__dirname, "../fixtures");

// NON-JSON
export const Na4Cl4Poscar = fs.readFileSync(path.join(FIXTURES_DIR, "Na4Cl4.poscar"), "utf8");
export const SiPWSCFInput = fs.readFileSync(path.join(FIXTURES_DIR, "Si-pwscf.in"), "utf8");
export const Zr1H23Zr1H1Poscar = fs.readFileSync(
    path.join(FIXTURES_DIR, "Zr1H23Zr1H1.poscar"),
    "utf8",
);
export const H2OPoscar = fs.readFileSync(path.join(FIXTURES_DIR, "H2O.poscar"), "utf8");
export const GraphenePoscar = fs.readFileSync(path.join(FIXTURES_DIR, "Graphene.poscar"), "utf8");
export const NiHexPoscar = fs.readFileSync(path.join(FIXTURES_DIR, "Ni-hex.poscar"), "utf8");
