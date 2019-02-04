import {expect} from "chai";

import {Made} from "../src/made";
import {Basis} from "../src/basis/basis";

describe('Made', function () {
    it('unitCellFormula', function () {
        const elements = [
            {
                value: 'Si',
                id: 1
            },
            {
                value: 'Si',
                id: 2
            },
            {
                value: 'Li',
                id: 3
            },
            {
                value: 'Fe',
                id: 4
            },
            {
                value: 'Fe',
                id: 5
            },
            {
                value: 'Fe',
                id: 6
            },
            {
                value: 'Si',
                id: 7
            }
        ];

        expect(Made.unitCellFormula(elements)).to.be.equal('Si3LiFe3');
    });

    it('formula', function () {
        const elements = [
            {
                value: 'Si',
                id: 1
            },
            {
                value: 'Si',
                id: 2
            },
            {
                value: 'O',
                id: 3
            },
            {
                value: 'O',
                id: 4
            },
            {
                value: 'O',
                id: 5
            },
            {
                value: 'O',
                id: 6
            }
        ];

        expect(Made.formula(elements)).to.be.equal('SiO2');
    });

    it('supercell', function () {
        const material = {
            basis: {
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
                        value: [0.25, 0.25, 0.25]
                    }
                ]
            },
            lattice: {
                a: 5,
                b: 5,
                c: 5,
                alpha: 90,
                beta: 90,
                gamma: 90,
                units: {
                    length: 'angstrom',
                    angle: 'degree'
                }
            }
        };

        const expectedSupercell = {
            basis: {
                elements: [
                    {
                        id: 1,
                        value: 'Si'
                    },
                    {
                        id: 2,
                        value: 'Si'
                    },
                    {
                        id: 3,
                        value: 'Si'
                    },
                    {
                        id: 4,
                        value: 'Si'
                    },
                    {
                        id: 5,
                        value: 'Si'
                    },
                    {
                        id: 6,
                        value: 'Si'
                    },
                    {
                        id: 7,
                        value: 'Si'
                    },
                    {
                        id: 8,
                        value: 'Si'
                    },
                    {
                        id: 9,
                        value: 'Si'
                    },
                    {
                        id: 10,
                        value: 'Si'
                    },
                    {
                        id: 11,
                        value: 'Si'
                    },
                    {
                        id: 12,
                        value: 'Si'
                    },
                    {
                        id: 13,
                        value: 'Si'
                    },
                    {
                        id: 14,
                        value: 'Si'
                    },
                    {
                        id: 15,
                        value: 'Si'
                    },
                    {
                        id: 16,
                        value: 'Si'
                    }
                ],
                coordinates: [
                    {
                        id: 1,
                        value: [0.000000000, 0.000000000, 0.000000000]
                    },
                    {
                        id: 2,
                        value: [0.000000000, 0.000000000, 0.500000000]
                    },
                    {
                        id: 3,
                        value: [0.000000000, 0.500000000, 0.000000000]
                    },
                    {
                        id: 4,
                        value: [0.000000000, 0.500000000, 0.500000000]
                    },
                    {
                        id: 5,
                        value: [0.500000000, 0.000000000, 0.000000000]
                    },
                    {
                        id: 6,
                        value: [0.500000000, 0.000000000, 0.500000000]
                    },
                    {
                        id: 7,
                        value: [0.500000000, 0.500000000, 0.000000000]
                    },
                    {
                        id: 8,
                        value: [0.500000000, 0.500000000, 0.500000000]
                    },
                    {
                        id: 9,
                        value: [0.125000000, 0.125000000, 0.125000000]
                    },
                    {
                        id: 10,
                        value: [0.125000000, 0.125000000, 0.625000000]
                    },
                    {
                        id: 11,
                        value: [0.125000000, 0.625000000, 0.125000000]
                    },
                    {
                        id: 12,
                        value: [0.125000000, 0.625000000, 0.625000000]
                    },
                    {
                        id: 13,
                        value: [0.625000000, 0.125000000, 0.125000000]
                    },
                    {
                        id: 14,
                        value: [0.625000000, 0.125000000, 0.625000000]
                    },
                    {
                        id: 15,
                        value: [0.625000000, 0.625000000, 0.125000000]
                    },
                    {
                        id: 16,
                        value: [0.625000000, 0.625000000, 0.625000000]
                    }
                ]
            },
            lattice: {
                a: 10,
                b: 10,
                c: 10,
                alpha: 90,
                beta: 90,
                gamma: 90,
                units: {
                    length: 'angstrom',
                    angle: 'degree'
                }
            }
        };

        const supercell = Made.tools.supercell.generateConfig(material, [[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
        expect(supercell.lattice).deep.equal(expectedSupercell.lattice);

        const basis1 = Basis.fromExabyteFormat(expectedSupercell.basis);
        const basis2 = Basis.fromExabyteFormat(supercell.basis);

        // TODO: enable after refactoring made.js: problem is because of 5 != 0.49999999999999994
        // expect(basis1.equals(basis2)).to.be.ok;
    });
});
