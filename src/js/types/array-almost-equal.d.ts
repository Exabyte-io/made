/**
 * Declaration file for the array-almost-equal package
 * Used for comparing arrays of numbers with a specified tolerance
 */

declare module "array-almost-equal" {
    /**
     * Compares two arrays of numbers and returns true if all elements are
     * equal within the specified tolerance
     *
     * @param a - First array to compare
     * @param b - Second array to compare
     * @param tolerance - The maximum allowed difference between array elements (default: 1e-6)
     * @returns boolean - True if arrays are almost equal, false otherwise
     */
    function almostEqual(a: number[], b: number[], tolerance?: number): boolean;

    export default almostEqual;
}
