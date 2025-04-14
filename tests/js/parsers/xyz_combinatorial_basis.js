import { expect } from "chai";

import { Basis } from "../../../src/js/basis/basis";
import { CombinatorialBasis } from "../../../src/js/parsers/xyz_combinatorial_basis";
import { AsGeBasis, FeLiSiBasis, Ge2Basis, OSiBasis, Si2Basis } from "../fixtures";

describe("Parsers:CombinatorialBasis", () => {
    it("toBasisConfig", () => {
        // eslint-disable-next-line new-cap
        const basisConfig2 = CombinatorialBasis.toBasisConfig([
            {
                element: "Si",
                coordinates: [0, 0, 0],
            },
            {
                element: "Li",
                coordinates: [0.3, 0.3, 0.3],
            },
            {
                element: "Fe",
                coordinates: [0.5, 0.5, 0.5],
            },
        ]);
        const [basis1, basis2] = [new Basis(FeLiSiBasis), new Basis(basisConfig2)];

        expect(basis1.isEqualTo(basis2)).to.equal(true);
    });

    it("Regular XYZ", () => {
        const xyz = `
            Si 0.0 0.0 0.0
            O 0.5 0.5 0.5
        `;
        const basis = new CombinatorialBasis(xyz);
        expect(basis.uniqueElements).deep.equal(["O", "Si"]);
        const basis1 = new Basis(basis.allBasisConfigs[0]);
        const basis2 = new Basis(OSiBasis);
        expect(basis1.isEqualTo(basis2)).to.equal(true);
    });

    it("Permutation XYZ", () => {
        const permutation = `
            Si/Ge/As 0.0 0.0 0.0
            Si/Ge 0.5 0.5 0.5
        `;
        const basis = new CombinatorialBasis(permutation);
        expect(basis.uniqueElements).deep.equal(["As", "Ge", "Si"]);
        const basisList1 = basis.allBasisConfigs.map((c) => new Basis(c));
        const basisList2 = [Si2Basis, Ge2Basis, AsGeBasis].map((c) => new Basis(c));

        basisList1.forEach((basis1) => {
            const condition = basisList2
                .map((basis2) => basis2.isEqualTo(basis1))
                .reduce((a, b) => a || b);
            expect(condition).to.equal(true);
        });
    });

    it("Combination XYZ", () => {
        const combination = `
            Si,Ge,As 0.0 0.0 0.0
            Si,Ge 0.5 0.5 0.5
        `;
        const basis = new CombinatorialBasis(combination);
        expect(basis.uniqueElements).deep.equal(["As", "Ge", "Si"]);

        const basisList1 = basis.allBasisConfigs.map((c) => new Basis(c));
        expect(basisList1.length).to.equal(6);

        const basisList2 = [
            CombinatorialBasis.toBasisConfig([
                {
                    element: "Si",
                    coordinates: [0, 0, 0],
                },
                {
                    element: "Si",
                    coordinates: [0.5, 0.5, 0.5],
                },
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: "Si",
                    coordinates: [0, 0, 0],
                },
                {
                    element: "Ge",
                    coordinates: [0.5, 0.5, 0.5],
                },
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: "Ge",
                    coordinates: [0, 0, 0],
                },
                {
                    element: "Si",
                    coordinates: [0.5, 0.5, 0.5],
                },
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: "Ge",
                    coordinates: [0, 0, 0],
                },
                {
                    element: "Ge",
                    coordinates: [0.5, 0.5, 0.5],
                },
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: "As",
                    coordinates: [0, 0, 0],
                },
                {
                    element: "Si",
                    coordinates: [0.5, 0.5, 0.5],
                },
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: "As",
                    coordinates: [0, 0, 0],
                },
                {
                    element: "Ge",
                    coordinates: [0.5, 0.5, 0.5],
                },
            ]),
        ].map((c) => new Basis(c));

        basisList1.forEach((basis1) => {
            const condition = basisList2
                .map((basis2) => basis2.isEqualTo(basis1))
                .reduce((a, b) => a || b);
            expect(condition).to.equal(true);
        });
    });

    it("Mixing combination and permutation should fail", () => {
        const mixedInMultipleLines = `
            Si/Ge/As 0.0 0.0 0.0\n
            Si,Ge 0.5 0.5 0.5
        `;
        const mixedInOneLine = `
            Si/Ge,As 0.0 0.0 0.0
        `;
        // Wrapping into function as recommended in: http://chaijs.com/api/bdd/
        expect(() => new CombinatorialBasis(mixedInMultipleLines)).to.throw();
        expect(() => new CombinatorialBasis(mixedInOneLine)).to.throw();
    });
});
