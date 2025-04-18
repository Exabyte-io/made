import { ArrayWithIds } from "@mat3ra/code";
import { PointSchema as Coordinate } from "@mat3ra/esse/dist/js/types";

/**
 * Specialized class for handling atomic coordinates in a basis
 * Extends ArrayWithIds to provide coordinate-specific functionality
 */
export class Coordinates extends ArrayWithIds<Coordinate> {}
