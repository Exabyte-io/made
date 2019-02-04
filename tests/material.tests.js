import {expect} from "chai";

import {Material} from "../src/material";

describe('Material#getAsPOSCAR', function () {
    it('should return poscar with sorted basis items', function () {
        const material = new Material({
            name: 'Material #1',
            basis: {
                elements: [
                    {
                        id: 1,
                        value: "Li"
                    },
                    {
                        id: 2,
                        value: "Si"
                    },
                    {
                        id: 3,
                        value: "Fe"
                    },
                    {
                        id: 4,
                        value: "Li"
                    },
                    {
                        id: 5,
                        value: "Si"
                    },
                    {
                        id: 6,
                        value: "Si"
                    }
                ],
                coordinates: [
                    {
                        id: 1,
                        value: [0, 0, 0]
                    },
                    {
                        id: 2,
                        value: [0.25, 0, 0]
                    },
                    {
                        id: 3,
                        value: [0, 0.25, 0]
                    },
                    {
                        id: 4,
                        value: [0, 0, 0.25]
                    },
                    {
                        id: 5,
                        value: [0.25, 0.5, 0]
                    },
                    {
                        id: 6,
                        value: [0.5, 0.25, 0.25]
                    }
                ],
                name: "basis",
                units: "crystal"
            },
            lattice: {
                a: 5,
                b: 5,
                c: 5,
                alpha: 60,
                beta: 60,
                gamma: 60,
                units: {
                    length: "angstrom",
                    angle: "degree"
                },
                type: "FCC",
                vectors: {
                    a: [5, 0, 0],
                    b: [0, 5, 0],
                    c: [0, 0, 5],
                    name: "lattice vectors",
                    alat: 1,
                    units: "angstrom"
                }
            }
        });

        const poscar = `Material #1
1.0
   4.330127000	   0.000000000	   2.500000000
   1.443376000	   4.082483000	   2.500000000
   0.000000000	   0.000000000	   5.000000000
Fe Li Si
1 2 3
direct
   0.000000000    0.250000000    0.000000000 Fe
   0.000000000    0.000000000    0.000000000 Li
   0.000000000    0.000000000    0.250000000 Li
   0.250000000    0.000000000    0.000000000 Si
   0.250000000    0.500000000    0.000000000 Si
   0.500000000    0.250000000    0.250000000 Si`;

        expect(material.getAsPOSCAR()).to.be.equal(poscar);
    });
});
