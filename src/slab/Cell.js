
//package burai.atoms.model;

//import java.util.ArrayList;
//import java.util.List;

//import burai.atoms.model.event.CellEvent;
//import burai.atoms.model.event.CellEventListener;
//import burai.atoms.model.exception.ZeroVolumCellException;
//import burai.com.consts.ConstantAtoms;
//import burai.com.math.Lattice;
//import burai.com.math.Matrix3D;

import {Matrix3D} from './Matrix3D';
import {ConstantAtoms} from './ConstantAtoms';
import {Atom} from './Atom';

//export class Cell extends Model<CellEvent, CellEventListener> {
export class Cell {
    static get ATOMS_POSITION_WITH_LATTICE() { return 0; }
    static get ATOMS_POSITION_SCALED() { return 1; }
    static get ATOMS_POSITION_LEFT() { return 2; }

    static get MAX_ATOMS_TO_RESOLVE() { return ConstantAtoms.MAX_NUM_ATOMS; }

    static get MIN_VOLUME() { return 1.0e-6; }

    static get MIN_BOUNDARY() { return 3.0; }

    static get THR_LATTICE() { return 1.0e-4; }

    #lattice = null;

    #volume = null;

    #normLattice = null;

    #recLattice = null;

    #atoms = null;

    #bonds = null;

    #maxAtomsToResolve = null;

    #resolverStopping = null;

    #atomsResolver = null;

    #bondsResolver = null;

    static getEmptyCell() {
        try {
            return new Cell(Matrix3D.unitOne());
        } catch (e) {
            if(e /*ZeroVolumCellException*/) {
                e.printStackTrace();
            }
        }

        return null;
    }

    constructor(lattice, maxAtomsToResolve = Cell.MAX_ATOMS_TO_RESOLVE) {
        this.checkLattice(lattice);
        this.setupLattice(lattice);

        this.atoms = null;
        this.bonds = null;

        //this.maxAtomsToResolve = Math.max(0, maxAtomsToResolve);

        //this.resolverStopping = false;

//        this.atomsResolver = new AtomsResolver(this);

//        this.bondsResolver = null;
//        if (this.maxAtomsToResolve > 0) {
//            this.bondsResolver = new BondsResolver(this);
//        }
    }

    checkLattice(lattice) {
        if (lattice == null || lattice.length < 3) {
            throw new IllegalArgumentException("lattice is null or too short.");
        }

        for (let i = 0; i < 3; i++) {
            if (lattice[i] == null || lattice.length < 3) {
                throw new IllegalArgumentException("lattice[" + i + "] is null or too short.");
            }
        }

        const volume = this.calcVolumeByLattice(lattice);
        if (volume < Cell.MIN_VOLUME) {
            throw new ZeroVolumCellException();
        }
    }

    setupLattice(lattice) {
        this.lattice = Matrix3D.copy2D(lattice);
        this.calcVolume();
        this.calcNormLattice();
        this.calcRecLattice();
    }

    calcVolumeByLattice(lattice) {
        return Math.abs(Matrix3D.determinant(lattice));
    }

    calcVolume() {
        this.volume = this.calcVolumeByLattice(this.lattice);
    }

    calcNormLattice() {
        this.normLattice = [0.0, 0.0, 0.0];
        this.normLattice[0] = Matrix3D.norm(this.lattice[0]);
        this.normLattice[1] = Matrix3D.norm(this.lattice[1]);
        this.normLattice[2] = Matrix3D.norm(this.lattice[2]);
    }

    calcRecLattice() {
        this.recLattice = Matrix3D.inverse(this.lattice);
    }

/*
    @Override
    protected CellEvent createEvent() {
        return new CellEvent(this);
    }
*/

/*
    protected double[][] getLattice() {
        return this.lattice;
    }
*/

/*
    protected double[][] getRecLattice() {
        return this.recLattice;
    }
*/

/*
    protected double[] getNormLattice() {
        return this.normLattice;
    }
*/

/*
    public double getVolume() {
        return this.volume;
    }
*/

/*
    protected double getBoundaryVolume() {
        Atom atom0 = null;
        if (this.atoms != null && !(this.atoms.isEmpty())) {
            atom0 = this.atoms.get(0);
        }

        double x0 = atom0 == null ? 0.0 : atom0.getX();
        double y0 = atom0 == null ? 0.0 : atom0.getY();
        double z0 = atom0 == null ? 0.0 : atom0.getZ();
        double a0 = x0 * this.recLattice[0][0] + y0 * this.recLattice[1][0] + z0 * this.recLattice[2][0];
        double b0 = x0 * this.recLattice[0][1] + y0 * this.recLattice[1][1] + z0 * this.recLattice[2][1];
        double c0 = x0 * this.recLattice[0][2] + y0 * this.recLattice[1][2] + z0 * this.recLattice[2][2];
        double aMin = a0;
        double aMax = a0;
        double bMin = b0;
        double bMax = b0;
        double cMin = c0;
        double cMax = c0;

        int natom = this.atoms == null ? 0 : this.atoms.size();
        for (int i = 1; i < natom; i++) {
            Atom atom = this.atoms.get(i);
            double x = atom.getX();
            double y = atom.getY();
            double z = atom.getZ();
            double a = x * this.recLattice[0][0] + y * this.recLattice[1][0] + z * this.recLattice[2][0];
            double b = x * this.recLattice[0][1] + y * this.recLattice[1][1] + z * this.recLattice[2][1];
            double c = x * this.recLattice[0][2] + y * this.recLattice[1][2] + z * this.recLattice[2][2];
            aMin = Math.min(aMin, a);
            aMax = Math.max(aMax, a);
            bMin = Math.min(bMin, b);
            bMax = Math.max(bMax, b);
            cMin = Math.min(cMin, c);
            cMax = Math.max(cMax, c);
        }

        double volume = this.volume;
        volume *= Math.max((aMax - aMin), Math.min(MIN_BOUNDARY / this.normLattice[0], 1.0));
        volume *= Math.max((bMax - bMin), Math.min(MIN_BOUNDARY / this.normLattice[1], 1.0));
        volume *= Math.max((cMax - cMin), Math.min(MIN_BOUNDARY / this.normLattice[2], 1.0));

        return volume;
    }
*/

    copyLattice() {
        return Matrix3D.copy2D(this.lattice);
    }

/*
    private double getLatticeConstant(int ibrav, int i, boolean asCos) {
        double[] lattConst = Lattice.getLatticeConstants(ibrav, this.lattice, asCos);
        if (lattConst == null || lattConst.length < 6) {
            return 0.0;
        }

        return lattConst[i];
    }
*/

/*
    public double getA() {
        return Lattice.getA(this.lattice);
    }
*/

/*
    public double getA(int ibrav) {
        return this.getLatticeConstant(ibrav, 0, true);
    }
*/

/*
    public double getB() {
        return Lattice.getB(this.lattice);
    }
*/

/*
    public double getB(int ibrav) {
        return this.getLatticeConstant(ibrav, 1, true);
    }
*/

/*
    public double getC() {
        return Lattice.getC(this.lattice);
    }
*/

/*
    public double getC(int ibrav) {
        return this.getLatticeConstant(ibrav, 2, true);
    }
*/

/*
    public double getAlpha() {
        return Lattice.getAlpha(this.lattice);
    }
*/

/*
    public double getAlpha(int ibrav) {
        return this.getLatticeConstant(ibrav, 3, false);
    }
*/

/*
    public double getBeta() {
        return Lattice.getBeta(this.lattice);
    }
*/

/*
    public double getBeta(int ibrav) {
        return this.getLatticeConstant(ibrav, 4, false);
    }
*/

/*
    public double getGamma() {
        return Lattice.getGamma(this.lattice);
    }
*/

/*
    public double getGamma(int ibrav) {
        return this.getLatticeConstant(ibrav, 5, false);
    }
*/

/*
    public double getCosAlpha() {
        return Lattice.getCosAlpha(this.lattice);
    }
*/

/*
    public double getCosAlpha(int ibrav) {
        return this.getLatticeConstant(ibrav, 3, true);
    }
*/

/*
    public double getCosBeta() {
        return Lattice.getCosBeta(this.lattice);
    }
*/

/*
    public double getCosBeta(int ibrav) {
        return this.getLatticeConstant(ibrav, 4, true);
    }
*/

/*
    public double getCosGamma() {
        return Lattice.getCosGamma(this.lattice);
    }
*/

/*
    public double getCosGamma(int ibrav) {
        return this.getLatticeConstant(ibrav, 5, true);
    }
*/

    convertToCartesianPosition(a, b, c, lattice) {
        return Matrix3D.mult1D2D([ a, b, c ], lattice);
    }

/*
    public double[] convertToCartesianPosition(double a, double b, double c) {
        return this.convertToCartesianPosition(a, b, c, this.lattice);
    }
*/

    convertToLatticePositionRecLattice(x, y, z, recLattice) {
        return Matrix3D.mult1D2D([ x, y, z ], recLattice);
    }

    convertToLatticePosition(x, y, z) {
        return this.convertToLatticePositionRecLattice(x, y, z, this.recLattice);
    }

/*
    public boolean isInCell(double x, double y, double z) {
        double[] position = this.convertToLatticePosition(x, y, z);
        double a = position[0];
        double b = position[1];
        double c = position[2];

        boolean inCell = true;
        inCell = inCell && (0.0 <= a) && (a < 1.0);
        inCell = inCell && (0.0 <= b) && (b < 1.0);
        inCell = inCell && (0.0 <= c) && (c < 1.0);

        return inCell;
    }
*/

/*
    protected List<Atom> getAtoms() {
        return this.atoms;
    }
*/

    numAtoms(masterOnly) {
        if (this.atoms === null || this.atoms.length === 0) {
            return 0;
        }

        let natom = this.atoms.length;

        if (masterOnly) {
            for (const atom of this.atoms) {
                if (atom === null || atom.isSlaveAtom()) {
                    natom--;
                }
            }
        }

        return natom;
    }

/*
    public int numAtoms() {
        return this.numAtoms(false);
    }
*/

    listAtoms(masterOnly) {
        if (this.atoms === null) {
            return null;
        }

        let atoms2 = this.atoms;

        if (masterOnly) {
            atoms2 = [];
            for (const atom of this.atoms) {
                if (!atom.isSlaveAtom()) {
                    atoms2.push(atom);
                }
            }
        }

        return atoms2;
    }

/*
    public Atom[] listAtoms() {
        return this.listAtoms(false);
    }
*/

/*
    protected List<Bond> getBonds() {
        return this.bonds;
    }
*/

/*
    public int numBonds() {
        if (this.bonds == null || this.bonds.isEmpty()) {
            return 0;
        }

        return this.bonds.size();
    }
*/

/*
    public Bond[] listBonds() {
        if (this.bonds == null) {
            return null;
        }

        return this.bonds.toArray(new Bond[this.bonds.size()]);
    }
*/

/*
    public boolean isResolving() {
        return (this.numAtoms(true) <= this.maxAtomsToResolve);
    }
*/

/*
    public void stopResolving() {
        this.resolverStopping = true;

        if (this.atomsResolver != null) {
            this.atomsResolver.setAuto(false);
        }

        if (this.bondsResolver != null) {
            this.bondsResolver.setAuto(false);
        }
    }
*/

/*
    public void restartResolving() {
        this.resolverStopping = false;

        if (this.atomsResolver != null) {
            this.atomsResolver.setAuto(true);
            this.atomsResolver.resolve();
        }

        if (this.bondsResolver != null && this.isResolving()) {
            this.bondsResolver.setAuto(true);
            this.bondsResolver.resolve();
        }
    }
*/

    equalsLattice(lattice) {
        this.checkLattice(lattice);

        return Matrix3D.equals2D2DThreshold(this.lattice, lattice, Cell.THR_LATTICE);
    }

    moveLattice(lattice) {
        return this.moveLatticeAtomsPositionRefAtoms(lattice, Cell.ATOMS_POSITION_WITH_LATTICE, null);
    }

/*
    public boolean moveLattice(double lattice[][], int atomsPosition) throws ZeroVolumCellException {
        return this.moveLattice(lattice, atomsPosition, null);
    }
*/

    moveLatticeAtomsPositionRefAtoms(lattice, atomsPosition, refAtoms) {
        if (this.equalsLattice(lattice)) {
            return false;
        }

//        boolean orgAutoAtoms = false;
//        if (this.atomsResolver != null) {
//            orgAutoAtoms = this.atomsResolver.isAuto();
//            this.atomsResolver.setAuto(false);
//        }

//        boolean orgAutoBonds = false;
//        if (this.bondsResolver != null) {
//            orgAutoBonds = this.bondsResolver.isAuto();
//            this.bondsResolver.setAuto(false);
//        }

        if (this.atoms !== null) {
            if (atomsPosition === Cell.ATOMS_POSITION_WITH_LATTICE) {
                const atomList = this.listAtoms(true);
                const withRef = (refAtoms !== null && atomList.length <= refAtoms.size);

                for (let i = 0; i < atomList.length; i++) {
                    let atom = null;
                    if (withRef) {
                        atom = refAtoms.get(i);
                        atom = (atom === null) ? atomList[i] : atom;
                    } else {
                        atom = atomList[i];
                    }

                    let position = null;
                    let x = atom.getX();
                    let y = atom.getY();
                    let z = atom.getZ();
                    position = this.convertToLatticePosition(x, y, z);
                    const a = position[0];
                    const b = position[1];
                    const c = position[2];
                    position = this.convertToCartesianPosition(a, b, c, lattice);
                    x = position[0];
                    y = position[1];
                    z = position[2];
                    atom.moveTo(x, y, z);
                }

            } else if (atomsPosition === Cell.ATOMS_POSITION_SCALED) {
                const oldScale = this.normLattice[0];
                const newScale = Matrix3D.norm(lattice[0]);
                const atomList = this.listAtoms(true);
                const withRef = (refAtoms !== null && atomList.length <= refAtoms.size);

                for (let i = 0; i < atomList.length; i++) {
                    let atom = null;
                    if (withRef) {
                        atom = refAtoms.get(i);
                        atom = (atom === null) ? atomList[i] : atom;
                    } else {
                        atom = atomList[i];
                    }

                    const x = (newScale / oldScale) * atom.getX();
                    const y = (newScale / oldScale) * atom.getY();
                    const z = (newScale / oldScale) * atom.getZ();
                    atomList[i].moveTo(x, y, z);
                }

            } else if (atomsPosition === Cell.ATOMS_POSITION_LEFT) {
                if (refAtoms !== null && refAtoms.size !== 0) {
                    const atomList = this.listAtoms(true);

                    let natom = 0;
                    if (atomList.length <= refAtoms.size) {
                        natom = atomList.length;
                    }

                    for (let i = 0; i < natom; i++) {
                        let atom = refAtoms.get(i);
                        atom = (atom === null) ? atomList[i] : atom;
                        const x = atom.getX();
                        const y = atom.getY();
                        const z = atom.getZ();
                        atomList[i].moveTo(x, y, z);
                    }
                }
            }
        }

        this.setupLattice(lattice);

//        if (this.atomsResolver != null) {
//            this.atomsResolver.setAuto(orgAutoAtoms);
//        }

//        if (this.bondsResolver != null) {
//            this.bondsResolver.setAuto(orgAutoBonds);
//        }

//        if (this.listeners != null) {
//            CellEvent event = new CellEvent(this);
//            event.setLattice(lattice);
//            for (CellEventListener listener : this.listeners) {
//                listener.onLatticeMoved(event);
//            }
//        }

        return true;
    }

/*
    public boolean hasAtomAt(double a, double b, double c) {
        double[] position = this.convertToCartesianPosition(a, b, c);
        double x = position[0];
        double y = position[1];
        double z = position[2];
        return this.hasAtomAt(new Atom(null, x, y, z));
    }
*/

/*
    public boolean hasAtomAt(Atom atom) {
        if (atom == null) {
            return false;
        }

        if (this.atoms == null || this.atoms.isEmpty()) {
            return false;
        }

        if (this.atomsResolver != null) {
            this.atomsResolver.packAtomIntoCell(atom);
        }

        for (Atom atom2 : this.atoms) {
            if (atom.equalsPosition(atom2)) {
                return true;
            }
        }

        return false;
    }
*/

/*
    public int indexOfAtom(Atom atom) {
        if (atom == null) {
            return -1;
        }

        if (this.atoms != null) {
            return this.atoms.indexOf(atom);
        }

        return -1;
    }
*/

    addAtom1(name, a, b, c) {
        const position = this.convertToCartesianPosition(a, b, c, this.lattice);
        const x = position[0];
        const y = position[1];
        const z = position[2];
        return this.addAtom2(new Atom(name, x, y, z));
    }

    addAtom2(atom) {
        if (atom === null) {
            return false;
        }

        if (this.atoms === null) {
            this.atoms = [];
        }

        //if (!atom.isSlaveAtom()) {
        //    if (this.atomsResolver != null) {
        //        this.atomsResolver.packAtomIntoCell(atom);
        //    }
        //}

        if (this.atoms.includes(atom)) {
            return false;
        }

        this.atoms.push(atom);
        //if (!status) {
        //    return false;
        //}

        //if (this.bondsResolver != null && (!this.resolverStopping)) {
        //    boolean auto1 = this.bondsResolver.isAuto();
        //    boolean auto2 = this.isResolving();
        //    if (auto1 && (!auto2)) {
        //        this.removeAllBonds();
        //        this.bondsResolver.setAuto(false);
        //    }
        //}

        //if (this.listeners != null) {
        //    CellEvent event = new CellEvent(this);
        //    event.setAtom(atom);
        //    for (CellEventListener listener : this.listeners) {
        //        listener.onAtomAdded(event);
        //    }
        //}

        return true;
    }

    removeAtom(atom) {
        if (atom === null) {
            return false;
        }

        if (this.atoms === null) {
            return false;
        }

        const index = this.atoms.indexOf(atom);
        if (index < 0) {
            return false;
        }

        const atom2 = this.atoms[index];
        this.atoms.splice(index, 1);
//        atom2.notDisplay();
        if (atom2.isSlaveAtom()) {
            atom2.setMasterAtom(null);
        }

//        if (this.listeners != null) {
//            CellEvent event = new CellEvent(this);
//            event.setAtom(atom2);
//            for (CellEventListener listener : this.listeners) {
//                listener.onAtomRemoved(event);
//            }
//        }

//        if (this.bondsResolver != null && (!this.resolverStopping)) {
//            boolean auto1 = this.bondsResolver.isAuto();
//            boolean auto2 = this.isResolving();
//            if ((!auto1) && auto2) {
//                this.bondsResolver.setAuto(true);
//                this.bondsResolver.resolve();
//            }
//        }

        return true;
    }

    removeAllAtoms() {
        const atomList = this.listAtoms();
        if (atomList === null || atomList.length < 1) {
            return;
        }

//        boolean orgAuto = false;
//        if (this.bondsResolver != null) {
//            orgAuto = this.bondsResolver.isAuto();
//            this.bondsResolver.setAuto(false);
//        }

        for (const atom of atomList) {
            this.removeAtom(atom);
        }

//        if (this.bondsResolver != null) {
//            this.bondsResolver.setAuto(orgAuto);
//            this.bondsResolver.resolve();
//        }
    }

/*
    protected Bond pickBond(Atom atom1, Atom atom2) {
        return this.pickBond(atom1, atom2, this.bonds);
    }
*/

/*
    protected Bond pickBond(Atom atom1, Atom atom2, List<Bond> bonds) {
        if (bonds == null || bonds.isEmpty()) {
            return null;
        }

        for (Bond bond : bonds) {
            Atom refAtom1 = bond.getAtom1();
            Atom refAtom2 = bond.getAtom2();
            if (refAtom1 == atom1 && refAtom2 == atom2) {
                return bond;
            }
            if (refAtom1 == atom2 && refAtom2 == atom1) {
                return bond;
            }
        }

        return null;
    }
*/

/*
    protected List<Bond> pickBonds(Atom atom1) {
        if (this.bonds == null || this.bonds.isEmpty()) {
            return null;
        }

        List<Bond> bonds = new ArrayList<Bond>();

        for (Bond bond : this.bonds) {
            Atom refAtom1 = bond.getAtom1();
            Atom refAtom2 = bond.getAtom2();
            if (refAtom1 == atom1 || refAtom2 == atom1) {
                bonds.add(bond);
            }
        }

        return bonds;
    }
*/

/*
    protected boolean addBond(Bond bond) {
        if (bond == null) {
            return false;
        }

        if (this.bonds == null) {
            this.bonds = new ArrayList<Bond>();
        }

        if (this.bonds.contains(bond)) {
            return false;
        }

        boolean status = this.bonds.add(bond);
        if (!status) {
            return false;
        }

        if (this.listeners != null) {
            CellEvent event = new CellEvent(this);
            event.setBond(bond);
            for (CellEventListener listener : this.listeners) {
                listener.onBondAdded(event);
            }
        }

        return true;
    }
*/

/*
    protected boolean removeBond(Bond bond) {
        if (bond == null) {
            return false;
        }

        if (this.bonds == null) {
            return false;
        }

        int index = this.bonds.indexOf(bond);
        if (index < 0) {
            return false;
        }

        Bond bond2 = this.bonds.remove(index);
        bond2.notDisplay();
        bond2.detachFromAtoms();

        if (this.listeners != null) {
            CellEvent event = new CellEvent(this);
            event.setBond(bond2);
            for (CellEventListener listener : this.listeners) {
                listener.onBondRemoved(event);
            }
        }

        return true;
    }
*/

/*
    protected void removeAllBonds() {
        Bond[] bondList = this.listBonds();

        for (Bond bond : bondList) {
            this.removeBond(bond);
        }
    }
*/

/*
    @Override
    public void flushListeners() {
        super.flushListeners();

        if (this.atoms != null) {
            for (Atom atom : this.atoms) {
                atom.flushListeners();
            }
        }

        if (this.bonds != null) {
            for (Bond bond : this.bonds) {
                bond.flushListeners();
            }
        }
    }
*/

/*
    @Override
    public void display() {
        super.display();

        if (this.atoms != null) {
            for (Atom atom : this.atoms) {
                atom.display();
            }
        }

        if (this.bonds != null) {
            for (Bond bond : this.bonds) {
                bond.display();
            }
        }
    }
*/

/*
    @Override
    public void notDisplay() {
        super.notDisplay();

        if (this.atoms != null) {
            for (Atom atom : this.atoms) {
                atom.notDisplay();
            }
        }

        if (this.bonds != null) {
            for (Bond bond : this.bonds) {
                bond.notDisplay();
            }
        }
    }
*/
}
