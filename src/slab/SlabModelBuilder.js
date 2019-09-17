
//package burai.app.project.viewer.modeler.slabmodel;

//import burai.atoms.model.Cell;

import {SlabModelStem} from './SlabModelStem';
import {Atom} from './Atom';
import {Cell} from './Cell';
import {Lattice} from './Lattice';
import {Constants} from './Constants';
import {Material} from "../material";
import {Made} from "@exabyte-io/made.js";

export class SlabModelBuilder {
    constructor(material) {

//        if (cell == null) {
//            throw new IllegalArgumentException("cell is null.");
//        }


        const lattice = [[material.lattice.a, 0.0, 0.0], [0.0, material.lattice.b, 0.0], [0.0, 0.0, material.lattice.c]];

//        Cell cell = null;
//        try {
            this.cell = new Cell(lattice);
//        } catch (ZeroVolumCellException e) {
//            throw new IOException(e);
//        }

        //this.cell.stopResolving();
        for (const i in material.basis.coordinates) {
            const x = material.basis.coordinates[i].value[0];
            const y = material.basis.coordinates[i].value[1];
            const z = material.basis.coordinates[i].value[2];
            this.cell.addAtom1(material.basis.elements[i].value, x, y, z);
        }

        //this.cell.restartResolving();



//        this.cell = new Cell(lattice);
    }

    build(h, k, l) {
//        SlabModel slabModel = null;
//        try {
            const slabModel = new SlabModelStem(this.cell, h, k, l);
//        } catch (MillerIndexException e) {
            //e.printStackTrace();
//            return null;
//        }

        const slabModels = slabModel.getSlabModels();
        if (slabModels == null) return null;
        if (slabModels.length < 1) return null;

        const retSlabModels = [];

        for (const slabModel2 of slabModels) {
            const cell = Cell.getEmptyCell();

            if (cell === null) {
                return [];
            }

            if (!slabModel2.putOnCell(cell)) {
                return [];
            }

            const ibrav = Lattice.getBravais(cell.lattice);
            const lattConst = Lattice.getLatticeConstants(14, cell.lattice, false);

            for (const atom of cell.atoms) {
                const pos = cell.convertToLatticePosition(atom.x, atom.y, atom.z);
                atom.x = pos[0];
                atom.y = pos[1];
                atom.z = pos[2];
            }

            retSlabModels.push({slabModel:slabModel2, cell:cell, latticeType:SlabModelBuilder.convertIbravToLatticeType(ibrav), lattConst:lattConst});
        }

        return retSlabModels;
    }

    static convertIbravToLatticeType(ibrav) {
        switch (ibrav) {
        case 1: // "Cubic P (sc)"
            return Made.LATTICE_TYPE.CUB;
        case 2: // "Cubic F (fcc)"
            return Made.LATTICE_TYPE.FCC;
        case 3: // "Cubic I (bcc)"
            return Made.LATTICE_TYPE.BCC;
        case -3: // "Cubic I (bcc), more symmetric"
            return Made.LATTICE_TYPE.BCC;
        case 4: // "Hexagonal and Trigonal P"
            return Made.LATTICE_TYPE.HEX;

        case 5: // "Trigonal R, 3fold axis c"
            return Made.LATTICE_TYPE.TRI;
        case -5: // "Trigonal R, 3fold axis 111"
            return Made.LATTICE_TYPE.TRI;

        case 6: // "Tetragonal P (st)"
            return Made.LATTICE_TYPE.TET;
        case 7: // "Tetragonal I (bct)"
            return Made.LATTICE_TYPE.BCT;
        case 8: // "Orthorhombic P"
            return Made.LATTICE_TYPE.ORC;
        case 9: // "Orthorhombic base-c. #1"
            return Made.LATTICE_TYPE.ORCC;
        case -9: // "Orthorhombic base-c. #2"
            return Made.LATTICE_TYPE.ORCC;
    //        case 91:
    //            break;
        case 10: // "Orthorhombic face-c."
            return Made.LATTICE_TYPE.ORCF;
        case 11: // "Orthorhombic body-c."
            return Made.LATTICE_TYPE.ORCI;
        case 12: // "Monoclinic P, unique axis c"
            return Made.LATTICE_TYPE.MCL;
        case -12: // "Monoclinic P, unique axis b"
            return Made.LATTICE_TYPE.MCL;
        case 13: // "Monoclinic base-centered"
            return Made.LATTICE_TYPE.MCLC;
    //        case -13:
    //            break;
        case 14: // "Triclinic"
            return Made.LATTICE_TYPE.TRI;
        default:
            return '';
        }
    }
}

