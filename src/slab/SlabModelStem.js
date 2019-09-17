
//package burai.app.project.viewer.modeler.slabmodel;

//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.HashMap;
//import java.util.LinkedHashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//import java.util.UUID;

//import burai.app.project.viewer.modeler.ModelerBase;
//import burai.atoms.model.Atom;
//import burai.atoms.model.Cell;
//import burai.atoms.model.exception.ZeroVolumCellException;
//import burai.atoms.model.property.AtomProperty;
//import burai.atoms.model.property.CellProperty;
//import burai.com.math.Lattice;
//import burai.com.math.Matrix3D;

import {SlabModel} from './SlabModel';
import {SlabGenom} from './SlabGenom';
import {SlabModelLeaf} from './SlabModelLeaf';
import {Atom} from './Atom';
import {Matrix3D} from './Matrix3D';
import {Lattice} from './Lattice';
import {Constants} from './Constants';
import {ConstantAtoms} from './ConstantAtoms';
import {CellProperty} from './CellProperty';

export class SlabModelStem extends SlabModel {
    static get DET_THR() { return 1.0e-8; }
    static get PACK_THR() { return 1.0e-6; } // internal coordinate
    static get OFFSET_THR() { return 1.0e-12; } // angstrom
    static get THICK_THR() { return 1.0e-12; } // angstrom
    static get POSIT_THR() { return 1.0e-4; } // angstrom
    static get VALUE_THR() { return 1.0e-12; }
    static get STOIX_THR() { return 1.0e-6; }

    static get STOIX_NATOM() { return "NAtom"; }

    static get STEP_FOR_GENOMS() { return 0.50; } // angstrom
    static get STEP_FOR_CTHICK() { return 0.05; } // internal coordinate
    static get MAX_FOR_CTHICK() { return 20; }

    static get SLAB_FIX_THR() { return 0.1; } // angstrom
    static get SLAB_FIX_RATE() { return 0.5; } // internal coordinate

    #miller1 = null;
    #miller2 = null;
    #miller3 = null;

    #intercept1 = null;
    #intercept2 = null;
    #intercept3 = null;
    #numIntercept = null;
    #hasIntercept1 = null;
    #hasIntercept2 = null;
    #hasIntercept3 = null;

    #vector1 = null;
    #vector2 = null;
    #vector3 = null;

    #boundBox = null;

    #lattConst = null;
    #lattUnit = null;
    #lattAuxi = null;
    #lattSlab = null;

    #entryUnit = null;
    #entryAuxi = null;
    #entrySlab = null;

    #stoixUnit = null;

    #codeAuxi = null;
    #codeSlab = null;

    #currOffset = null;
    #currThickness = null;
    #currVacuum = null;
    #currScaleA = null;
    #currScaleB = null;

    /**
     * cell is not changed.
     */
    constructor(cell, h, k, l) {
        super();

        if (cell == null) {
            throw new IllegalArgumentException("cell is null.");
        }

        if (!this.setupMillers(cell, h, k, l)) {
            throw new MillerIndexException();
        }

        if (!this.setupUnitAtomsInCell(cell)) {
            throw new MillerIndexException();
        }

        if (!this.setupUnitAtomsInSlab()) {
            throw new MillerIndexException();
        }

        this.lattAuxi = null;
        this.lattSlab = null;
        this.entryAuxi = null;
        this.entrySlab = null;
        this.codeAuxi = null;
        this.codeSlab = null;

        this.currOffset = null;
        this.currThickness = null;
        this.currVacuum = null;
        this.currScaleA = null;
        this.currScaleB = null;
    }

    getSlabModels() {
        if (this.lattUnit === null || this.lattUnit.length < 3) {
            return [ new SlabModelLeaf(this, this.offset) ];
        }
        if (this.lattUnit[2] === null || this.lattUnit[2].length < 3) {
            return [ new SlabModelLeaf(this, this.offset) ];
        }

        const nstep = parseInt(this.lattUnit[2][2] / SlabModelStem.STEP_FOR_GENOMS);

        if (nstep < 2) {
            return new [ new SlabModelLeaf(this, this.offset) ];
        }

        const slabGenoms = new Map();

        for (let i = 0; i < nstep; i++) {
            const offset = i / nstep;

            const slabGenom = this.getSlabGenom(offset);
            if (slabGenom !== null) {
                let include = false;
                for (const [mapSlabGenom, mapOffset] of slabGenoms) {
                  if (mapSlabGenom.equals(slabGenom) === true) {
                    include = true;
                    break;
                  }
                }

                if (include === false) {
                  slabGenoms.set(slabGenom, offset);
                }
            }
        }

        if (slabGenoms.size === 0) {
            return [ new SlabModelLeaf(this, this.offset) ];
        }

        const slabModels = [];
        for (const [key, offset] of slabGenoms) {
            slabModels.push(new SlabModelLeaf(this, offset));
        }

        return slabModels;
    }

    getSlabGenom(offset) {
        if (this.lattUnit === null || this.lattUnit.length < 3) {
            return null;
        }
        if (this.lattUnit[2] === null || this.lattUnit[2].length < 3) {
            return null;
        }

        if (this.entryUnit === null || this.entryUnit.length === 0) {
            return null;
        }

        const natom = this.entryUnit.length;
        let iatom = natom;

        const names = [];
        const coords = [];

        for (let i = 0; i < natom; i++) {
            let entry = this.entryUnit[i];
            if (entry === null) {
                return null;
            }

            const c1 = entry.c + offset;
            let c2 = c1;
            c2 -= Math.floor(c2);

            const dc = Math.abs(c2 - 1.0);
            let dz = dc * this.lattUnit[2][2];
            if (dz < SlabModelStem.POSIT_THR) {
                c2 -= 1.0;
            }

            dz = Math.abs(c1 - c2) * this.lattUnit[2][2];
            if (iatom >= natom && dz < SlabModelStem.POSIT_THR) {
                iatom = i;
            }

            names.push(entry.name);
            coords.push(c2 * this.lattUnit[2][2]);
        }

        const names2 = [];
        const coords2 = [];
        for (let i = 0; i < natom; i++) {
            names2.push('');
            coords2.push(0.0);
        }

        let index = 0;
        for (let i = iatom; i < natom; i++) {
            names2[index] = names[i];
            coords2[index] = coords[i];
            index++;
        }

        for (let i = 0; i < iatom; i++) {
            names2[index] = names[i];
            coords2[index] = coords[i];
            index++;
        }

        return new SlabGenom(names2, coords2);
    }

//    updateCell(cell) {
//        return this.updateCell(cell, this);
//    }

    updateCell(cell, slabModel) {
        if (cell === null) {
            return false;
        }

        if (slabModel === null) {
            return false;
        }

        if (!this.setupAuxiAtoms(slabModel)) {
            return false;
        }

        if (!this.setupSlabAtoms(slabModel)) {
            return false;
        }

        if (this.lattSlab === null || this.lattSlab.length < 3) {
            return false;
        }

        if (this.entrySlab === null || this.entrySlab.length === 0) {
            return false;
        }

        //cell.stopResolving();

        //try {
            cell.moveLattice(this.lattSlab);
        //} catch (e) {
        //    if (e instanceof ZeroVolumCellException) {
        //      cell.restartResolving();
        //      return false;
        //    }
        //}

        //cell.setProperty(CellProperty.AXIS, "z");

        const natom = this.entrySlab.length;
        const natom2 = cell.numAtoms(true);
        let refAtoms = null;
        if (natom === natom2) {
            refAtoms = cell.listAtoms(true);
        }

        let zMax = -Number.MAX_VALUE;
        let zMin = +Number.MAX_VALUE;
        for (const entry of this.entrySlab) {
            if (entry !== null) {
                zMax = Math.max(zMax, entry.z);
                zMin = Math.min(zMin, entry.z);
            }
        }

        if (refAtoms !== null && refAtoms.length == natom) {
            for (let i = 0; i < natom; i++) {
                const entry = this.entrySlab.get(i);
                if (entry === null) {
                    continue;
                }

                const name = entry.name;
                if (name === null || name.length === 0) {
                    continue;
                }

                let atom = refAtoms[i];
                if (atom === null) {
                    atom = new Atom(name, entry.x, entry.y, entry.z);
                    cell.addAtom2(atom);

                } else {
                    const name2 = atom.getName();
                    if (name !== name2) {
                        atom.setName(name);
                    }

                    const x2 = atom.getX();
                    const y2 = atom.getY();
                    const z2 = atom.getZ();
                    if (Math.abs(entry.x - x2) > 0.0 || Math.abs(entry.y - y2) > 0.0 || Math.abs(entry.z - z2) > 0.0) {
                        atom.moveTo(entry.x, entry.y, entry.z);
                    }
                }

                let toFix = false;
                if ((zMax - zMin) < SlabModelStem.SLAB_FIX_THR) {
                    toFix = true;
                } else if ((entry.z - zMin - SlabModelStem.SLAB_FIX_THR) / (zMax - zMin) < SlabModelStem.SLAB_FIX_RATE) {
                    toFix = true;
                }

                //const xFix = atom.booleanProperty(AtomProperty.FIXED_X);
                //const yFix = atom.booleanProperty(AtomProperty.FIXED_Y);
                //const zFix = atom.booleanProperty(AtomProperty.FIXED_Z);

                //if (xFix !== toFix) {
                //    atom.setProperty(AtomProperty.FIXED_X, toFix);
                //}
                //if (yFix !== toFix) {
                //    atom.setProperty(AtomProperty.FIXED_Y, toFix);
                //}
                //if (zFix !== toFix) {
                //    atom.setProperty(AtomProperty.FIXED_Z, toFix);
                //}
            }

        } else {
            cell.removeAllAtoms();

            for (const entry of this.entrySlab) {
                if (entry === null) {
                    continue;
                }

                const name = entry.name;
                if (name === null || name.length === 0) {
                    continue;
                }

                let toFix = false;
                if ((zMax - zMin) < SlabModelStem.SLAB_FIX_THR) {
                    toFix = true;
                } else if ((entry.z - zMin - SlabModelStem.SLAB_FIX_THR) / (zMax - zMin) < SlabModelStem.SLAB_FIX_RATE) {
                    toFix = true;
                }

                const atom = new Atom(name, entry.x, entry.y, entry.z);
                //atom.setProperty(AtomProperty.FIXED_X, toFix);
                //atom.setProperty(AtomProperty.FIXED_Y, toFix);
                //atom.setProperty(AtomProperty.FIXED_Z, toFix);

                cell.addAtom2(atom);
            }
        }

//        cell.restartResolving();
        return true;
    }

    setupSlabAtoms(slabModel) {
        if (slabModel === null) {
            return false;
        }

        // check status
        let sameCondition = true;
        let vacuumOnlyChanged = true;

        if (this.codeSlab === null || this.codeSlab !== this.codeAuxi) {
            sameCondition = false;
            vacuumOnlyChanged = false;
        }

        if (this.currVacuum === null || Math.abs(slabModel.vacuum - this.currVacuum) > SlabModelStem.VALUE_THR) {
            sameCondition = false;
        }

        if (this.currScaleA === null || slabModel.scaleA !== parseInt(this.currScaleA)) {
            sameCondition = false;
            vacuumOnlyChanged = false;
        }

        if (this.currScaleB === null || slabModel.scaleB !== parseInt(this.currScaleB)) {
            sameCondition = false;
            vacuumOnlyChanged = false;
        }

        if (sameCondition) {
            return true;
        }

        this.codeSlab = null;
        this.currVacuum = null;
        this.currScaleA = null;
        this.currScaleB = null;

        // lattice
        if (this.lattUnit === null) {
            return false;
        }
        if (this.lattAuxi === null || this.lattAuxi.length < 3) {
            return false;
        }
        if (this.lattAuxi[2] === null || this.lattAuxi[2].length < 3) {
            return false;
        }

        let lattSlab2 = null;
        if (vacuumOnlyChanged) {
            lattSlab2 = this.lattSlab === null ? null : Matrix3D.copy2D(this.lattSlab);
        }

        this.lattSlab = Matrix3D.copy2D(this.lattUnit);
        if (this.lattSlab === null || this.lattSlab.length < 3) {
            return false;
        }

        const aScale = Math.max(1, slabModel.scaleA);
        this.lattSlab[0] = Matrix3D.multAlpha1D(aScale, this.lattSlab[0]);
        if (this.lattSlab[0] === null || this.lattSlab[0].length < 3) {
            return false;
        }

        const bScale = Math.max(1, slabModel.scaleB);
        this.lattSlab[1] = Matrix3D.multAlpha1D(bScale, this.lattSlab[1]);
        if (this.lattSlab[1] === null || this.lattSlab[1].length < 3) {
            return false;
        }

        const zSlab = this.lattSlab[2][2];
        const zTotal = this.lattAuxi[2][2] + 2.0 * Math.max(0.0, slabModel.vacuum);
        const zScale = zSlab === 0.0 ? 1.0 : (zTotal / zSlab);
        this.lattSlab[2] = Matrix3D.multAlpha1D(zScale, this.lattSlab[2]);
        if (this.lattSlab[2] === null || this.lattSlab[2].length < 3) {
            return false;
        }

        // atoms
        if (vacuumOnlyChanged) {
            // shift atoms
            if (this.entrySlab == null || this.entrySlab.length === 0) {
                return false;
            }

            if (lattSlab2 == null || lattSlab2.length < 3) {
                return false;
            }
            if (lattSlab2[2] == null || lattSlab2[2].length < 3) {
                return false;
            }

            const tx = 0.5 * (this.lattSlab[2][0] - lattSlab2[2][0]);
            const ty = 0.5 * (this.lattSlab[2][1] - lattSlab2[2][1]);
            const tz = 0.5 * (this.lattSlab[2][2] - lattSlab2[2][2]);

            for (const entry of this.entrySlab) {
                if (entry === null) {
                    continue;
                }

                entry.x += tx;
                entry.y += ty;
                entry.z += tz;
            }

        } else {
            // create atoms
            const natom = this.entryAuxi.length;
            if (this.entryAuxi === null || natom < 1) {
                return false;
            }

            if ((aScale * bScale * natom) >= ConstantAtoms.MAX_NUM_ATOMS) {
                return false;
            }

            this.entrySlab = [];

            const tx = 0.5 * (this.lattSlab[2][0] - this.lattAuxi[2][0]);
            const ty = 0.5 * (this.lattSlab[2][1] - this.lattAuxi[2][1]);
            const tz = 0.5 * (this.lattSlab[2][2] - this.lattAuxi[2][2]);

            for (let ia = 0; ia < aScale; ia++) {
                const ra = ia / aScale;
                for (let ib = 0; ib < bScale; ib++) {
                    const rb = ib / bScale;

                    const vx = ra * this.lattSlab[0][0] + rb * this.lattSlab[1][0];
                    const vy = ra * this.lattSlab[0][1] + rb * this.lattSlab[1][1];

                    for (const entry of this.entryAuxi) {
                        if (entry === null) {
                            continue;
                        }

                        const entry2 = new AtomEntry(this.lattSlab);
                        entry2.name = entry.name;
                        entry2.x = entry.x + tx + vx;
                        entry2.y = entry.y + ty + vy;
                        entry2.z = entry.z + tz;

                        this.entrySlab.push(entry2);
                    }
                }
            }
        }

        // store status
        this.codeSlab = this.codeAuxi;
        this.currVacuum = slabModel.vacuum;
        this.currScaleA = slabModel.scaleA;
        this.currScaleB = slabModel.scaleB;

        return true;
    }

    setupAuxiAtoms(slabModel) {
        if (slabModel === null) {
            return false;
        }

        // check status
        let sameCondition = true;

        if (this.currOffset == null || Math.abs(slabModel.offset - this.currOffset) > SlabModelStem.VALUE_THR) {
            sameCondition = false;
        }

        if (this.currThickness == null || Math.abs(slabModel.thickness - this.currThickness) > SlabModelStem.VALUE_THR) {
            sameCondition = false;
        }

        if (sameCondition) {
            return true;
        }

        this.codeAuxi = null;
        this.currOffset = null;
        this.currThickness = null;

        // init lattice
        if (this.lattUnit === null) {
            return false;
        }

        this.lattAuxi = Matrix3D.copy2D(this.lattUnit);
        if (this.lattAuxi === null || this.lattAuxi.length < 3) {
            return false;
        }
        if (this.lattAuxi[0] === null || this.lattAuxi[0].length < 3) {
            return false;
        }
        if (this.lattAuxi[1] === null || this.lattAuxi[1].length < 3) {
            return false;
        }
        if (this.lattAuxi[2] === null || this.lattAuxi[2].length < 3) {
            return false;
        }

        // init atoms
        if (this.entryUnit === null || this.entryUnit.length === 0) {
            return false;
        }

        const cOffset = slabModel.offset - Math.floor(slabModel.offset);
        let cThick = Math.max(0.0, slabModel.thickness);

        for (let istep = 0; istep < SlabModelStem.MAX_FOR_CTHICK; istep++) {

            const hasOffset = this.lattAuxi[2][2] * Math.abs(Math.min(cOffset, 1.0 - cOffset)) > SlabModelStem.OFFSET_THR;
            const hasThick = this.lattAuxi[2][2] * Math.abs(cThick - 1.0) > SlabModelStem.THICK_THR;
            if (hasOffset || hasThick) {
                this.entryAuxi = [];
            } else {
                this.entryAuxi = this.entryUnit;
            }

            const nThick = parseInt(Math.ceil(cThick) + 0.1);
            for (let iThick = 0; iThick < nThick; iThick++) {
                let entryBuffer = null;
                let stoixBuffer = null;
                if (iThick === (nThick - 1)) {
                    entryBuffer = [];
                    stoixBuffer = new Map();
                }

                for (const phase of [ true, false ]) {

                    for (const entry of this.entryUnit) {
                        if (entry === null) {
                            continue;
                        }

                        const name = entry.name;
                        if (name === null || name.length === 0) {
                            continue;
                        }

                        const a = entry.a;
                        const b = entry.b;
                        let c = entry.c + cOffset;
                        const c_ = c;
                        c -= Math.floor(c);

                        let dc = Math.abs(c - 1.0);
                        let dz = dc * this.lattAuxi[2][2];
                        if (dz < SlabModelStem.POSIT_THR) {
                            c -= 1.0;
                        }

                        const zshift = Math.abs(c - c_) * this.lattAuxi[2][2];
                        if (zshift < SlabModelStem.POSIT_THR) {
                            if (!phase) {
                                continue;
                            }
                        } else {
                            if (phase) {
                                continue;
                            }
                        }

                        c -= iThick;

                        dc = c - (1.0 - cThick);
                        dz = dc * this.lattAuxi[2][2];
                        if (dz < (-2.0 * SlabModelStem.POSIT_THR)) {
                            continue;
                        }

                        let entry2 = null;
                        if (this.entryAuxi !== this.entryUnit) {
                            entry2 = new AtomEntry(this.lattAuxi);
                        } else {
                            entry2 = entry;
                        }

                        entry2.name = name;
                        entry2.x = a * this.lattAuxi[0][0] + b * this.lattAuxi[1][0] + c * this.lattAuxi[2][0];
                        entry2.y = a * this.lattAuxi[0][1] + b * this.lattAuxi[1][1] + c * this.lattAuxi[2][1];
                        entry2.z = a * this.lattAuxi[0][2] + b * this.lattAuxi[1][2] + c * this.lattAuxi[2][2];

                        if (this.entryAuxi !== this.entryUnit) {
                            if (iThick < (nThick - 1)) {
                                this.entryAuxi.push(entry2);

                            } else {
                                entryBuffer.push(entry2);

                                stoixBuffer.set(SlabModelStem.STOIX_NATOM, entryBuffer.length);
                                if (stoixBuffer.has(name)) {
                                    const value = stoixBuffer.get(name);
                                    stoixBuffer.set(name, value + 1.0);
                                } else {
                                    stoixBuffer.set(name, 1.0);
                                }

                                if (this.stoixUnit !== null && this.equalsStoichiometry(this.stoixUnit, stoixBuffer)) {
                                    this.entryAuxi = this.entryAuxi.concat(entryBuffer);
                                    entryBuffer = [];
                                    stoixBuffer.clear();
                                }
                            }
                        }
                    }
                }
            }

            if (this.entryAuxi.length > 0) {
                break;
            }

            cThick += SlabModelStem.STEP_FOR_CTHICK;
        }

        if (this.entryAuxi.length >= ConstantAtoms.MAX_NUM_ATOMS) {
            return false;
        }

        // modify lattice
        let zMax = 0.0;
        let zMin = 0.0;
        let zFirst = true;
        for (const entry of this.entryAuxi) {
            if (entry === null) {
                continue;
            }

            if (zFirst) {
                zFirst = false;
                zMax = entry.z;
                zMin = entry.z;
            } else {
                zMax = Math.max(zMax, entry.z);
                zMin = Math.min(zMin, entry.z);
            }
        }

        const zSlab = this.lattAuxi[2][2];
        const zDelta = Math.max(zMax - zMin, 0.0);
        const zScale = zSlab === 0.0 ? 1.0 : (zDelta / zSlab);
        this.lattAuxi[2] = Matrix3D.multAlpha1D(zScale, this.lattAuxi[2]);
        if (this.lattAuxi[2] === null || this.lattAuxi[2].length < 3) {
            return false;
        }

        // modify atoms
        let xMin = 0.0;
        let yMin = 0.0;
        if (this.lattAuxi[2][2] > SlabModelStem.THICK_THR) {
            xMin = zMin * this.lattAuxi[2][0] / this.lattAuxi[2][2];
            yMin = zMin * this.lattAuxi[2][1] / this.lattAuxi[2][2];
        }

        for (const entry of this.entryAuxi) {
            if (entry === null) {
                continue;
            }

            entry.x -= xMin;
            entry.y -= yMin;
            entry.z -= zMin;
        }

        // store status
//        const uuid = UUID.randomUUID();
//        if (uuid !== null) {
//            this.codeAuxi = uuid.toString();
//        }

        this.codeAuxi = new Date().getTime().toString(16) + Math.floor(1000 * Math.random()).toString(16)

        this.currOffset = slabModel.offset;
        this.currThickness = slabModel.thickness;

        return true;
    }

    setupUnitAtomsInCell(cell) {
        if (cell === null) {
            return false;
        }

        const atoms = cell.listAtoms(true);
        if (atoms === null || atoms.length < 1) {
            return false;
        }

        this.entryUnit = [];

        const lattice = cell.copyLattice();

        for (const atom of atoms) {
            if (atom === null) {
                return false;
            }

            const name = atom.getName();
            if (name === null || name.length === 0) {
                return false;
            }

            const x = atom.getX();
            const y = atom.getY();
            const z = atom.getZ();

            const coord = cell.convertToLatticePosition(x, y, z);
            if (coord === null || coord.length < 3) {
                return false;
            }

            if (lattice !== null) {
                const entry = new AtomEntry(lattice);
                entry.name = name;
                entry.a = coord[0];
                entry.b = coord[1];
                entry.c = coord[2];
                this.entryUnit.push(entry);
            }
        }

        return true;
    }

    setupUnitAtomsInSlab() {
        const entryUnit_ = this.entryUnit;
        const natom = entryUnit_.length;
        if (natom < 1) {
            return false;
        }

        const lattInt = [];
        lattInt.push([ this.vector1[0], this.vector1[1], this.vector1[2] ]);
        lattInt.push([ this.vector2[0], this.vector2[1], this.vector2[2] ]);
        lattInt.push([ this.vector3[0], this.vector3[1], this.vector3[2] ]);

        const detLatt = Matrix3D.determinant(lattInt);
        if (Math.abs(detLatt) < SlabModelStem.DET_THR) {
            console.error("volume is zero.");
            return false;
        }
        if (detLatt < 0.0) {
            console.error("parity is incorrect.");
            return false;
        }

        const detLatt2 = Math.round(detLatt);
        if (Math.abs(detLatt - detLatt2) >= SlabModelStem.DET_THR) {
            console.error("not integer volume of lattice.");
            return false;
        }

        const nsize = parseInt(detLatt2 + 0.1);
        if ((nsize * natom) >= Constants.MAX_NUM_ATOMS) {
            return false;
        }

        const invLatt = Matrix3D.inverse(lattInt);
        if (invLatt === null || invLatt.length < 3) {
            return false;
        }
        if (invLatt[0] === null || invLatt[0].length < 3) {
            return false;
        }
        if (invLatt[1] === null || invLatt[1].length < 3) {
            return false;
        }
        if (invLatt[2] === null || invLatt[2].length < 3) {
            return false;
        }

        this.entryUnit = [];

        for (let ia = this.boundBox[0][0]; ia <= this.boundBox[0][1]; ia++) {
            const ta = ia;
            for (let ib = this.boundBox[1][0]; ib <= this.boundBox[1][1]; ib++) {
                const tb = ib;
                for (let ic = this.boundBox[2][0]; ic <= this.boundBox[2][1]; ic++) {
                    const tc = ic;

                    for (const entry of entryUnit_) {
                        const a = entry.a + ta;
                        const b = entry.b + tb;
                        const c = entry.c + tc;
                        const a2 = a * invLatt[0][0] + b * invLatt[1][0] + c * invLatt[2][0];
                        const b2 = a * invLatt[0][1] + b * invLatt[1][1] + c * invLatt[2][1];
                        const c2 = a * invLatt[0][2] + b * invLatt[1][2] + c * invLatt[2][2];

                        if (-SlabModelStem.PACK_THR <= a2 && a2 < (1.0 + SlabModelStem.PACK_THR) &&
                                -SlabModelStem.PACK_THR <= b2 && b2 < (1.0 + SlabModelStem.PACK_THR) &&
                                -SlabModelStem.PACK_THR <= c2 && c2 < (1.0 + SlabModelStem.PACK_THR)) {

                            const entry2 = new AtomEntry(this.lattUnit);
                            entry2.name = entry.name;
                            entry2.a = a2;
                            entry2.b = b2;
                            entry2.c = c2;

                            let include = false;
                            for (const entry3 of this.entryUnit) {
                              if (entry3.equals(entry2)) {
                                include = true;
                                break;
                              }
                            }

                            if (!include) {
                                this.entryUnit.push(entry2);
                            }
                        }
                    }
                }
            }
        }

        if (this.entryUnit.length === 0) {
            return false;
        }

        for (const entry of this.entryUnit) {
            entry.a -= Math.floor(entry.a);
            entry.b -= Math.floor(entry.b);
            entry.c -= Math.floor(entry.c);

            const dc = Math.abs(entry.c - 1.0);
            const dz = dc * this.lattUnit[2][2];
            if (dz < SlabModelStem.POSIT_THR) {
                entry.c -= 1.0;
            }
        }

        this.entryUnit.sort(function(entry1, entry2) {
            if (entry2 === null) {
                return -1;
            }

            const dc = entry1.c - entry2.c;
            let dx = dc * entry1.lattice[2][0];
            let dy = dc * entry1.lattice[2][1];
            let dz = dc * entry1.lattice[2][2];
            let rr = dx * dx + dy * dy + dz * dz;
            if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
                if (dc > 0.0) {
                    return -1;
                } else {
                    return 1;
                }
            }

            const db = entry1.b - entry2.b;
            dx = db * entry1.lattice[1][0];
            dy = db * entry1.lattice[1][1];
            dz = db * entry1.lattice[1][2];
            rr = dx * dx + dy * dy + dz * dz;
            if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
                if (db > 0.0) {
                    return 1;
                } else {
                    return -1;
                }
            }

            const da = entry1.a - entry2.a;
            dx = da * entry1.lattice[0][0];
            dy = da * entry1.lattice[0][1];
            dz = da * entry1.lattice[0][2];
            rr = dx * dx + dy * dy + dz * dz;
            if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
                if (da > 0.0) {
                    return 1;
                } else {
                    return -1;
                }
            }

            if (entry1.name === null) {
                if (entry2.name === null) {
                    return 0;
                } else {
                    return 1;
                }
            }

            return entry1.name.localeCompare(entry2.name);
        } );

        // keep stoichiometry
        this.stoixUnit = new Map();
        if (!this.setupStoichiometry(this.entryUnit, this.stoixUnit)) {
            return false;
        }

        return true;
    }

    setupStoichiometry(entryList, stoixMap) {
        if (entryList === null || entryList.length === 0) {
            return false;
        }

        if (stoixMap === null) {
            return false;
        }
        if (stoixMap.size > 0) {
            stoixMap.clear();
        }

        let natom = 0;

        for (const entry of entryList) {
            const name = entry == null ? null : entry.name;
            if (name === null || name.length === 0) {
                continue;
            }

            natom++;

            if (stoixMap.has(name)) {
                const value = stoixMap.get(name);
                stoixMap.set(name, value + 1.0);
            } else {
                stoixMap.set(name, 1.0);
            }
        }

        if (natom < 1) {
            return false;
        }

        stoixMap.set(SlabModelStem.STOIX_NATOM, natom);

        return true;
    }

    equalsStoichiometry(stoixMap1, stoixMap2) {
        if (stoixMap1 === null || stoixMap1.size === 0) {
            return false;
        }

        if (stoixMap2 === null || stoixMap2.size === 0) {
            return false;
        }

        const names1 = stoixMap1.keys();
        if (names1 === null) {
            return false;
        }

        const names2 = stoixMap2.keys();
        if (names2 === null) {
            return false;
        }

        if (names1.length !== names2.length) {
            return false;
        }

        let natom1 = 0.0;
        if (stoixMap1.has(SlabModelStem.STOIX_NATOM)) {
            natom1 = stoixMap1.get(SlabModelStem.STOIX_NATOM);
        }

        if (natom1 <= 0.0) {
            return false;
        }

        let natom2 = 0.0;
        if (stoixMap2.has(SlabModelStem.STOIX_NATOM)) {
            natom2 = stoixMap2.get(SlabModelStem.STOIX_NATOM);
        }

        if (natom2 <= 0.0) {
            return false;
        }

        for (const name of names1) {
            if (!stoixMap2.has(name)) {
                return false;
            }

            const value1 = stoixMap1.get(name) / natom1;
            const value2 = stoixMap2.get(name) / natom2;
            if (Math.abs(value1 - value2) > SlabModelStem.STOIX_THR) {
                return false;
            }
        }

        return true;
    }

    setupMillers(cell, h, k, l) {
        if (cell === null) {
            return false;
        }

        if (h === 0 && k === 0 && l === 0) {
            return false;
        }

        this.miller1 = h;
        this.miller2 = k;
        this.miller3 = l;

        if (!this.setupIntercepts()) {
            return false;
        }

        if (!this.setupVectors()) {
            return false;
        }

        if (!this.setupBoundaryBox()) {
            return false;
        }

        if (!this.setupLattice(cell)) {
            return false;
        }

        return true;
    }

    setupIntercepts() {
        let scaleMin = 1;
        let scaleMax = 1;
        this.numIntercept = 0;

        if (this.miller1 != 0) {
            scaleMin = Math.max(scaleMin, Math.abs(this.miller1));
            scaleMax *= Math.abs(this.miller1);
            this.numIntercept++;
            this.hasIntercept1 = true;
        } else {
            this.hasIntercept1 = false;
        }

        if (this.miller2 != 0) {
            scaleMin = Math.max(scaleMin, Math.abs(this.miller2));
            scaleMax *= Math.abs(this.miller2);
            this.numIntercept++;
            this.hasIntercept2 = true;
        } else {
            this.hasIntercept2 = false;
        }

        if (this.miller3 != 0) {
            scaleMin = Math.max(scaleMin, Math.abs(this.miller3));
            scaleMax *= Math.abs(this.miller3);
            this.numIntercept++;
            this.hasIntercept3 = true;
        } else {
            this.hasIntercept3 = false;
        }

        if (scaleMin < 1) {
            console.error("scaleMin is not positive.");
            return false;
        }

        if (scaleMax < scaleMin) {
            console.error("scaleMax < scaleMin.");
            return false;
        }

        if (this.numIntercept < 1) {
            console.error("there are no intercepts.");
            return false;
        }

        let scale = 0;
        for (let i = scaleMin; i <= scaleMax; i++) {
            if (this.hasIntercept1 && (i % this.miller1) != 0) {
                continue;
            }
            if (this.hasIntercept2 && (i % this.miller2) != 0) {
                continue;
            }
            if (this.hasIntercept3 && (i % this.miller3) != 0) {
                continue;
            }

            scale = i;
            break;
        }

        if (scale < 1) {
            console.error("cannot detect scale.");
            return false;
        }

        this.intercept1 = this.hasIntercept1 ? (scale / this.miller1) : 0;
        this.intercept2 = this.hasIntercept2 ? (scale / this.miller2) : 0;
        this.intercept3 = this.hasIntercept3 ? (scale / this.miller3) : 0;

        return true;
    }

    setupVectors() {
        this.vector1 = [ 0, 0, 0 ];
        this.vector2 = [ 0, 0, 0 ];
        this.vector3 = [ 0, 0, 0 ];

        if (this.numIntercept <= 1) {
            this.setupVectors1();
        } else if (this.numIntercept <= 2) {
            this.setupVectors2();
        } else {
            this.setupVectors3();
        }

        return true;
    }

    setupVectors1() {
        if (this.hasIntercept1) {
            if (this.intercept1 > 0) {
                this.vector1[1] = 1;
                this.vector2[2] = 1;
                this.vector3[0] = 1;
            } else {
                this.vector1[2] = 1;
                this.vector2[1] = 1;
                this.vector3[0] = -1;
            }

        } else if (this.hasIntercept2) {
            if (this.intercept2 > 0) {
                this.vector1[2] = 1;
                this.vector2[0] = 1;
                this.vector3[1] = 1;
            } else {
                this.vector1[0] = 1;
                this.vector2[2] = 1;
                this.vector3[1] = -1;
            }

        } else if (this.hasIntercept3) {
            if (this.intercept3 > 0) {
                this.vector1[0] = 1;
                this.vector2[1] = 1;
                this.vector3[2] = 1;
            } else {
                this.vector1[1] = 1;
                this.vector2[0] = 1;
                this.vector3[2] = -1;
            }
        }
    }

    setupVectors2() {
        if (!this.hasIntercept3) { // cat in A-B plane
            const sign1 = Math.sign(this.intercept1);
            const sign2 = Math.sign(this.intercept2);
            this.vector1[2] = sign1 * sign2;
            this.vector2[0] = +this.intercept1;
            this.vector2[1] = -this.intercept2;
            this.vector3[0] = sign1;
            this.vector3[1] = sign2;

        } else if (!this.hasIntercept2) { // cat in A-C plane
            const sign1 = Math.sign(this.intercept1);
            const sign3 = Math.sign(this.intercept3);
            this.vector1[1] = sign1 * sign3;
            this.vector2[0] = -this.intercept1;
            this.vector2[2] = +this.intercept3;
            this.vector3[0] = sign1;
            this.vector3[2] = sign3;

        } else if (!this.hasIntercept1) { // cat in B-C plane
            const sign2 = Math.sign(this.intercept2);
            const sign3 = Math.sign(this.intercept3);
            this.vector1[0] = sign2 * sign3;
            this.vector2[1] = +this.intercept2;
            this.vector2[2] = -this.intercept3;
            this.vector3[1] = sign2;
            this.vector3[2] = sign3;
        }
    }

    setupVectors3() {
        const sign1 = Math.sign(this.intercept1);
        const sign2 = Math.sign(this.intercept2);
        const sign3 = Math.sign(this.intercept3);
        if (sign3 > 0) {
            this.vector1[1] = +sign1 * this.intercept2;
            this.vector1[2] = -sign1 * this.intercept3;
            this.vector2[0] = -sign2 * this.intercept1;
            this.vector2[2] = +sign2 * this.intercept3;
        } else {
            this.vector1[0] = -sign1 * this.intercept1;
            this.vector1[2] = +sign1 * this.intercept3;
            this.vector2[1] = +sign2 * this.intercept2;
            this.vector2[2] = -sign2 * this.intercept3;
        }
        this.vector3[0] = sign1;
        this.vector3[1] = sign2;
        this.vector3[2] = sign3;
    }

    setupBoundaryBox() {

        this.boundBox = [[0, 0], [0, 0], [0, 0]];

        for (let i = 0; i < 3; i++) {
            this.boundBox[i][0] = 0;
            this.boundBox[i][1] = 0;

            if (this.vector1[i] < 0) {
                this.boundBox[i][0] += this.vector1[i];
            } else {
                this.boundBox[i][1] += this.vector1[i];
            }

            if (this.vector2[i] < 0) {
                this.boundBox[i][0] += this.vector2[i];
            } else {
                this.boundBox[i][1] += this.vector2[i];
            }

            if (this.vector3[i] < 0) {
                this.boundBox[i][0] += this.vector3[i];
            } else {
                this.boundBox[i][1] += this.vector3[i];
            }
        }

        return true;
    }

    setupLattice(cell) {
        if (cell === null) {
            return false;
        }

        const lattInt = [];
        lattInt.push([ this.vector1[0], this.vector1[1], this.vector1[2] ]);
        lattInt.push([ this.vector2[0], this.vector2[1], this.vector2[2] ]);
        lattInt.push([ this.vector3[0], this.vector3[1], this.vector3[2] ]);

        const lattUnit0 = Matrix3D.mult2D2D(lattInt, cell.copyLattice());
        this.lattConst = lattUnit0 === null ? null : Lattice.getCellDm(14, lattUnit0);
        if (this.lattConst === null || this.lattConst.length < 6) {
            return false;
        }

        this.lattUnit = Lattice.getCell(14, this.lattConst);
        if (this.lattUnit === null || this.lattUnit.length < 3) {
            return false;
        }
        if (this.lattUnit[0] === null || this.lattUnit[0].length < 3) {
            return false;
        }
        if (this.lattUnit[1] === null || this.lattUnit[1].length < 3) {
            return false;
        }
        if (this.lattUnit[2] === null || this.lattUnit[2].length < 3) {
            return false;
        }

        return true;
    }

}

class AtomEntry {
    #lattice = null;

    name = null;
    a = null;
    b = null;
    c = null;
    x = null;
    y = null;
    z = null;

    constructor(lattice) {
        if (lattice === null) {
            throw new IllegalArgumentException("lattice is null.");
        }

        this.lattice = lattice;
        this.name = null;
        this.a = 0.0;
        this.b = 0.0;
        this.c = 0.0;
        this.x = 0.0;
        this.y = 0.0;
        this.z = 0.0;
    }

    compareTo(other) {
        return compareToStatic(this, other);
    }

    static compareToStatic(entry1, entry2) {
        if (entry2 === null) {
            return -1;
        }

        let dx;
        let dy;
        let dz;
        let rr;

        const dc = entry1.c - entry2.c;
        dx = dc * entry1.lattice[2][0];
        dy = dc * entry1.lattice[2][1];
        dz = dc * entry1.lattice[2][2];
        rr = dx * dx + dy * dy + dz * dz;
        if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
            if (dc > 0.0) {
                return -1;
            } else {
                return 1;
            }
        }

        const db = entry1.b - entry2.b;
        dx = db * entry1.lattice[1][0];
        dy = db * entry1.lattice[1][1];
        dz = db * entry1.lattice[1][2];
        rr = dx * dx + dy * dy + dz * dz;
        if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
            if (db > 0.0) {
                return 1;
            } else {
                return -1;
            }
        }

        const da = entry1.a - entry2.a;
        dx = da * entry1.lattice[0][0];
        dy = da * entry1.lattice[0][1];
        dz = da * entry1.lattice[0][2];
        rr = dx * dx + dy * dy + dz * dz;
        if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
            if (da > 0.0) {
                return 1;
            } else {
                return -1;
            }
        }

        if (entry1.name == null) {
            if (entry2.name == null) {
                return 0;
            } else {
                return 1;
            }
        }

        return entry1.name.compareTo(entry2.name);
    }

    hashCode() {
        return this.name == null ? 0 : this.name.hashCode();
    }

    equals(obj) {
        return AtomEntry.equalsStatic(this, obj);
    }

    static equalsStatic(entry, obj) {
        if (entry === obj) {
            return true;
        }
        if (obj === null) {
            return false;
        }
        //if (entry.getClass() !== obj.getClass()) {
        //    return false;
        //}

        const other = obj;
        if (entry.name === null) {
            if (entry.name !== other.name) {
                return false;
            }
        } else {
            if (entry.name !== other.name) {
                return false;
            }
        }

        let da = Math.abs(entry.a - other.a);
        da = Math.min(da, Math.abs(entry.a - other.a + Math.sign(0.5 - entry.a)));

        let db = Math.abs(entry.b - other.b);
        db = Math.min(db, Math.abs(entry.b - other.b + Math.sign(0.5 - entry.b)));

        let dc = Math.abs(entry.c - other.c);
        dc = Math.min(dc, Math.abs(entry.c - other.c + Math.sign(0.5 - entry.c)));

        const dx = da * entry.lattice[0][0] + db * entry.lattice[1][0] + dc * entry.lattice[2][0];
        const dy = da * entry.lattice[0][1] + db * entry.lattice[1][1] + dc * entry.lattice[2][1];
        const dz = da * entry.lattice[0][2] + db * entry.lattice[1][2] + dc * entry.lattice[2][2];
        const rr = dx * dx + dy * dy + dz * dz;
        if (rr > SlabModelStem.POSIT_THR * SlabModelStem.POSIT_THR) {
            return false;
        }

        return true;
    }
}

