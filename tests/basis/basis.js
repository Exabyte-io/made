import { expect } from 'chai';

import { Basis } from '../../src/basis/basis';
import { AsGeBasis, FeLiSiBasis, LiFeSiBasis, Na4Cl4, Na4Cl4Cartesian } from '../enums';

describe('Basis', () => {
    it('should return true if basises are equal', () => {
        const basis1 = new Basis(FeLiSiBasis);
        const basis2 = new Basis(LiFeSiBasis);
        expect(basis1.isEqualTo(basis2)).to.be.equal(true);
    });

    it('should return true when basis is compared to its clone', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.isEqualTo(basis.clone())).to.be.equal(true);
    });

    it('should return jsonified basis', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.toJSON()).to.be.deep.almost.equal(Na4Cl4.basis);
    });

    it('should return true if cells are equal', () => {
        const basis1 = new Basis(FeLiSiBasis);
        const basis2 = new Basis(LiFeSiBasis);
        expect(basis1.hasEquivalentCellTo(basis2)).to.be.equal(true);
    });

    /**
     * Elements
     */

    it('should return Na4Cl4 as unitCellFormula', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.unitCellFormula).to.be.equal('Na4Cl4');
    });

    it('should return NaCl as formula', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.formula).to.be.equal('NaCl');
    });

    it('should return elements', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.elements).to.be.deep.almost.equal(Na4Cl4.basis.elements);
    });

    it('should return unique elements', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.uniqueElements).to.be.deep.equal(['Na', 'Cl']);
    });

    it('should return elements count', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.elementCounts).to.be.deep.equal([
            {
                value: 'Na',
                count: 4,
            },
            {
                value: 'Cl',
                count: 4,
            },
        ]);
    });

    it('should set elements', () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.elements = AsGeBasis.elements;
        expect(basis.elements).to.be.deep.equal(AsGeBasis.elements);
    });

    /**
     * Coordinates
     */

    it('should return coordinates', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.coordinates).to.be.deep.almost.equal(Na4Cl4.basis.coordinates);
    });

    it('should set coordinates', () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.coordinates = AsGeBasis.coordinates;
        expect(basis.coordinates).to.be.deep.almost.equal(AsGeBasis.coordinates);
    });

    /**
     * Units
     */

    it('should return true if basis is in crystal units', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.isInCrystalUnits).to.be.equal(true);
    });

    it('should convert crystal to cartesian', () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.toCartesian();
        expect(basis.isInCartesianUnits).to.be.equal(true);
        expect(basis.coordinates).to.be.deep.almost.equal(Na4Cl4Cartesian.basis.coordinates);
    });

    /**
     * Atoms
     */

    it('should add a new atom', () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.addAtom({
            element: 'Ge',
            coordinate: [0, 0.25, 0.25],
        });
        expect(basis.elements).to.be.deep.equal([
            ...Na4Cl4.basis.elements,
            {
                id: 8,
                value: 'Ge',
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

    it('should remove an atom', () => {
        const basis = new Basis(Na4Cl4.basis);
        basis.removeAtom({
            id: 3,
            element: 'Na',
            coordinate: [0.5, 0.5, 0],
        });
        expect(basis.elements.length).to.be.equal(7);
        expect(basis.coordinates.length).to.be.equal(7);
    });

    it('should return number of atoms', () => {
        const basis = new Basis(Na4Cl4.basis);
        expect(basis.nAtoms).to.be.equal(8);
    });

    /**
     * Hash
     */

    it('should return the string', () => {
        const basis = new Basis(Na4Cl4.basis);
        const string =
            'Cl 0,0,0.5;Cl 0,0.5,0;Cl 0.5,0,0;Cl 0.5,0.5,0.5;Na 0,0,0;Na 0,0.5,0.5;Na 0.5,0,0.5;Na 0.5,0.5,0;';
        expect(basis.getAsSortedString()).to.be.equal(string);
    });

    /**
     * Representation
     */

    it('should return standard representation', () => {
        const basis = new Basis(Na4Cl4Cartesian.basis);
        expect(basis.standardRepresentation).to.be.deep.almost.equal(Na4Cl4.basis);
    });
});
