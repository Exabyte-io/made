import { expect } from "chai";

import { Basis } from "../../../src/js/basis/basis";
import {
    AsGeBasis,
    C2H4,
    C2H4Translated,
    FeLiSiBasis,
    LiFeSiBasis,
    Na,
    Na4Cl4,
    Na4Cl4Cartesian,
} from "../enums";
import { assertDeepAlmostEqual } from "../utils";

describe("Basis", () => {
    it("should return true if basises are equal", () => {
        const basis1 = new Basis(FeLiSiBasis);
        const basis2 = new Basis(LiFeSiBasis);
        expect(basis1.isEqualTo(basis2)).to.be.equal(true);
    });

    it("should return true when basis is compared to its clone", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.isEqualTo(basis.clone())).to.be.equal(true);
    });

    it("should return jsonified basis", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.toJSON()).to.be.deep.almost.equal(Na4Cl4.basis);
    });

    it("should return true if cells are equal", () => {
        const basis1 = new Basis(FeLiSiBasis);
        const basis2 = new Basis(LiFeSiBasis);
        expect(basis1.hasEquivalentCellTo(basis2)).to.be.equal(true);
    });

    /**
     * Elements
     */

    it("should return Na4Cl4 as unitCellFormula", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.unitCellFormula).to.be.equal("Na4Cl4");
    });

    it("should return NaCl as formula", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.formula).to.be.equal("NaCl");
    });

    it("should return elements", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.elements).to.be.deep.almost.equal(Na4Cl4.basis.elements);
    });

    it("should return unique elements", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.uniqueElements).to.be.deep.equal(["Na", "Cl"]);
    });

    it("should return elements count", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.elementCounts).to.be.deep.equal([
            {
                value: "Na",
                count: 4,
            },
            {
                value: "Cl",
                count: 4,
            },
        ]);
    });

    it("should set elements", () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.elements = AsGeBasis.elements;
        expect(basis.elements).to.be.deep.equal(AsGeBasis.elements);
    });

    /**
     * Coordinates
     */

    it("should return coordinates", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.coordinates).to.be.deep.almost.equal(Na4Cl4.basis.coordinates);
    });

    it("should set coordinates", () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.coordinates = AsGeBasis.coordinates;
        expect(basis.coordinates).to.be.deep.almost.equal(AsGeBasis.coordinates);
    });

    /**
     * Units
     */

    it("should return true if basis is in crystal units", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.isInCrystalUnits).to.be.equal(true);
    });

    it("should convert crystal to cartesian", () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.toCartesian();
        expect(basis.isInCartesianUnits).to.be.equal(true);
        expect(basis.coordinates).to.be.deep.almost.equal(Na4Cl4Cartesian.basis.coordinates);
    });

    /**
     * Atoms
     */

    it("should add a new atom", () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.addAtom({
            element: "Ge",
            coordinate: [0, 0.25, 0.25],
        });
        expect(basis.elements).to.be.deep.equal([
            ...Na4Cl4.basis.elements,
            {
                id: 8,
                value: "Ge",
            },
        ]);
        expect(basis.coordinates).to.be.deep.almost.equal([
            ...Na4Cl4.basis.coordinates,
            {
                id: 8,
                value: [0, 0.25, 0.25],
            },
        ]);
    });

    it("should remove an atom", () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.removeAtom({
            id: 3,
            element: "Na",
            coordinate: [0.5, 0.5, 0],
        });
        expect(basis.elements.length).to.be.equal(7);
        expect(basis.coordinates.length).to.be.equal(7);
    });

    it("should return number of atoms", () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.nAtoms).to.be.equal(8);
    });

    /**
     * Hash
     */

    it("should return the string", () => {
        const basis = new Basis(Na4Cl4.basis);
        const string =
            "Cl 0,0,0.5;Cl 0,0.5,0;Cl 0.5,0,0;Cl 0.5,0.5,0.5;Na 0,0,0;Na 0,0.5,0.5;Na 0.5,0,0.5;Na 0.5,0.5,0;";
        expect(basis.getAsSortedString()).to.be.equal(string);
    });

    /**
     * Representation
     */

    it("should return standard representation", () => {
        const basis = new Basis(Na4Cl4Cartesian.basis);
        expect(basis.standardRepresentation).to.be.deep.almost.equal(Na4Cl4.basis);
    });

    /**
     * Minimum lattice size generated for a molecule with more than one atom.
     * The minimumLatticeSize is the maximum pairwise distance between two atoms
     * of the structure, C2H4 in this case.
     */
    it("should return minimum lattice size for a molecule", () => {
        const basis = new Basis(C2H4.basis);
        const minimumLatticeSize = 3.162;
        const latticeSize = basis.getMinimumLatticeSize();
        expect(latticeSize).to.be.equal(minimumLatticeSize);
    });

    /**
     * Minimum lattice size generated for a structure made up of a single atom.
     * The minimumLatticeSize is the atomic radii of the atom, Na in this case, in units of angstroms.
     */
    it("should return minimum lattice size for an atom", () => {
        const basis = new Basis(Na.basis);
        const minimumLatticeSize = 1.8;
        const latticeSize = basis.getMinimumLatticeSize();
        expect(latticeSize).to.be.equal(minimumLatticeSize);
    });

    //* * Pairwise Distance */
    it("should return max distance", () => {
        const basis = new Basis(C2H4.basis);
        const maxDistance = 1.581;
        expect(basis.maxPairwiseDistance).to.be.equal(maxDistance);
    });

    //* * Center of Coordinates */
    it("should return center of coordinates", () => {
        const basis = new Basis(C2H4.basis);
        const centerOfCoordinatesArray = [0.3333, -0.0417, 0.0917];
        const basisCenterOfCoordinatesPoint = basis.centerOfCoordinatesPoint;
        assertDeepAlmostEqual(basisCenterOfCoordinatesPoint, centerOfCoordinatesArray);
    });

    //* * Translation by Vector */
    it("should return the updated basis coordinates", () => {
        const basis = new Basis(C2H4.basis);
        const translationVector = [1.6917, 1.9667, 2];
        basis.translateByVector(translationVector);
        assertDeepAlmostEqual(basis.coordinates, C2H4Translated.basis.coordinates);
    });

    it("should strip label and return atomic symbol", () => {
        expect(Basis.stripLabelToGetElementSymbol("Fe1")).to.be.equal("Fe");
        expect(Basis.stripLabelToGetElementSymbol("Fe11")).to.be.equal("Fe");
        expect(Basis.stripLabelToGetElementSymbol("fe_a")).to.be.equal("Fe");
        expect(Basis.stripLabelToGetElementSymbol("Fe-a1")).to.be.equal("Fe");
        expect(Basis.stripLabelToGetElementSymbol("FE1")).to.be.equal("Fe");
        expect(Basis.stripLabelToGetElementSymbol("c_a")).to.be.equal("C");
    });
});
