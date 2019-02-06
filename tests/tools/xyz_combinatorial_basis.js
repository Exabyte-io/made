import {expect} from "chai";

import {Basis} from "../../src/basis/basis";
import {CombinatorialBasis} from "../../src/parsers/xyz_combinatorial_basis";

describe('CombinatorialBasis', function () {
    it('toBasisConfig', function () {
        const basisConfig1 = {
            elements: [
                {
                    id: 1,
                    value: 'Fe'
                },
                {
                    id: 2,
                    value: 'Li'
                },
                {
                    id: 3,
                    value: 'Si'
                }
            ],
            coordinates: [
                {
                    id: 1,
                    value: [0.5, 0.5, 0.5]
                },
                {
                    id: 2,
                    value: [0.3, 0.3, 0.3]
                },
                {
                    id: 3,
                    value: [0, 0, 0]
                }
            ]
        };

        const basisConfig2 = new CombinatorialBasis.toBasisConfig([
            {
                element: 'Si',
                coordinates: [0, 0, 0]
            },
            {
                element: 'Li',
                coordinates: [0.3, 0.3, 0.3]
            },
            {
                element: 'Fe',
                coordinates: [0.5, 0.5, 0.5]
            }
        ]);

        const [basis1, basis2] = [new Basis(basisConfig1), new Basis(basisConfig2)];

        expect(basis1.isEqualTo(basis2)).to.equal(true);
    });

    it('Regular XYZ', function () {
        const xyz = `
            Si 0.0 0.0 0.0\n
            O 0.5 0.5 0.5
        `;

        const basis = new CombinatorialBasis(xyz);
        expect(basis.uniqueElements).deep.equal(['O', 'Si']);
        const basis1 = new Basis(basis.allBasisConfigs[0]);
        const basis2 = new Basis({
            elements: [
                {
                    id: 1,
                    value: 'O'
                },
                {
                    id: 2,
                    value: 'Si'
                }
            ],
            coordinates: [
                {
                    id: 1,
                    value: [0.5, 0.5, 0.5]
                },
                {
                    id: 2,
                    value: [0, 0, 0]
                }
            ]
        });
        expect(basis1.isEqualTo(basis2)).to.equal(true);
    });

    it('Permutation XYZ', function () {
        const permutation = `
            Si/Ge/As 0.0 0.0 0.0\n
            Si/Ge 0.5 0.5 0.5
        `;
        const basis = new CombinatorialBasis(permutation);
        expect(basis.uniqueElements).deep.equal(['As', 'Ge', 'Si']);

        const basisList1 = basis.allBasisConfigs.map(c => new Basis(c));

        const basisList2 = [
            {
                elements: [
                    {
                        id: 1,
                        value: 'Si'
                    },
                    {
                        id: 2,
                        value: 'Si'
                    }
                ],
                coordinates: [
                    {
                        id: 1,
                        value: [0, 0, 0]
                    },
                    {
                        id: 2,
                        value: [0.5, 0.5, 0.5]
                    }
                ]
            },
            {
                elements: [
                    {
                        id: 1,
                        value: 'Ge'
                    },
                    {
                        id: 2,
                        value: 'Ge'
                    }
                ],
                coordinates: [
                    {
                        id: 1,
                        value: [0, 0, 0]
                    },
                    {
                        id: 2,
                        value: [0.5, 0.5, 0.5]
                    }
                ]
            },
            {
                elements: [
                    {
                        id: 1,
                        value: 'As'
                    },
                    {
                        id: 2,
                        value: 'Ge'
                    }
                ],
                coordinates: [
                    {
                        id: 1,
                        value: [0, 0, 0]
                    },
                    {
                        id: 2,
                        value: [0.5, 0.5, 0.5]
                    }
                ]
            }
        ].map(c => new Basis(c));

        basisList1.forEach(basis1 => {
            const condition = basisList2.map(basis2 => basis2.isEqualTo(basis1)).reduce((a, b) => a || b);
            expect(condition).to.equal(true);
        })

    });

    it('Combination XYZ', function () {
        const combination = `
            Si,Ge,As 0.0 0.0 0.0\n
            Si,Ge 0.5 0.5 0.5
        `;
        const basis = new CombinatorialBasis(combination);
        expect(basis.uniqueElements).deep.equal(['As', 'Ge', 'Si']);

        const basisList1 = basis.allBasisConfigs.map(c => new Basis(c));
        expect(basisList1.length).to.equal(6);

        const basisList2 = [
            CombinatorialBasis.toBasisConfig([
                {
                    element: 'Si',
                    coordinates: [0, 0, 0]
                },
                {
                    element: 'Si',
                    coordinates: [0.5, 0.5, 0.5]
                }
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: 'Si',
                    coordinates: [0, 0, 0]
                },
                {
                    element: 'Ge',
                    coordinates: [0.5, 0.5, 0.5]
                }
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: 'Ge',
                    coordinates: [0, 0, 0]
                },
                {
                    element: 'Si',
                    coordinates: [0.5, 0.5, 0.5]
                }
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: 'Ge',
                    coordinates: [0, 0, 0]
                },
                {
                    element: 'Ge',
                    coordinates: [0.5, 0.5, 0.5]
                }
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: 'As',
                    coordinates: [0, 0, 0]
                },
                {
                    element: 'Si',
                    coordinates: [0.5, 0.5, 0.5]
                }
            ]),
            CombinatorialBasis.toBasisConfig([
                {
                    element: 'As',
                    coordinates: [0, 0, 0]
                },
                {
                    element: 'Ge',
                    coordinates: [0.5, 0.5, 0.5]
                }
            ])
        ].map(c => new Basis(c));

        basisList1.forEach(basis1 => {
            const condition = basisList2.map(basis2 => basis2.isEqualTo(basis1)).reduce((a, b) => a || b);
            expect(condition).to.equal(true);
        })
    });

    it('Mixing combination and permutation should fail', function () {
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
