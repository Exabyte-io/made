import {coordinatesGetCenterOfSpaceAsVector, coordinatesGetMaxPairwiseDistance} from "../../src/tools/coordinates";
import {Basis} from "../../src/basis/basis";
import {assertDeepAlmostEqual} from "../utils";
import {ArrayWithIds} from "../../src/abstract/array_with_ids"
import {C2H4} from "../enums";
import {expect} from "chai";


describe('Tools:Coordinates', function () {

    it('should return center of coordinates array', function() {
        const basis = new Basis(C2H4.basis);
        const basisCoordinatesArray = new ArrayWithIds(basis.coordinates).array;
        const nAtoms = new ArrayWithIds(basis.elements).array.length;
        const centerOfCoordinates = [ 0.3333, -0.0417, 0.0917 ];
        assertDeepAlmostEqual(coordinatesGetCenterOfSpaceAsVector(basisCoordinatesArray, nAtoms),centerOfCoordinates);
    })

    it('should return max pairwise distance between atoms of a structure', function() {
        const basis = new Basis(C2H4.basis);
        const basisCoordinates = new ArrayWithIds((basis.coordinates));
        const nAtoms = new ArrayWithIds(basis.elements).array.length;
        const maxDistance = 1.5811388300841898;
        expect(coordinatesGetMaxPairwiseDistance(basisCoordinates, nAtoms)).to.be.equal(maxDistance);
        }
    )
});
