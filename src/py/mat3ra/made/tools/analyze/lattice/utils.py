"""
Vector-based lattice type detection functions.

Based on primitive cell definitions from Setyawan & Curtarolo (2010)
and the JavaScript implementation in src/js/cell/primitive_cell.ts
"""
import math
from typing import List

import numpy as np

from ....lattice import LatticeTypeEnum


def detect_lattice_type_from_vectors(
    vectors: List[List[float]], 
    tolerance: float = 0.1, 
    angle_tolerance: float = 5.0
) -> LatticeTypeEnum:
    """
    Detect lattice type from lattice vectors using pattern matching.
    
    Based on the primitive cell definitions from Setyawan & Curtarolo (2010).
    
    Args:
        vectors: 3x3 array of lattice vectors
        tolerance: Tolerance for length comparisons (Angstroms)
        angle_tolerance: Tolerance for angle comparisons (degrees)
        
    Returns:
        Detected lattice type
    """
    if len(vectors) != 3 or any(len(v) != 3 for v in vectors):
        raise ValueError("Expected 3x3 array of lattice vectors")
    
    # Convert to numpy for easier calculations
    v = np.array(vectors)
    
    # Calculate lengths
    a = np.linalg.norm(v[0])
    b = np.linalg.norm(v[1])
    c = np.linalg.norm(v[2])
    
    # Calculate angles between vectors (in degrees)
    def angle_between_vectors(v1, v2):
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
        return np.degrees(np.arccos(cos_angle))
    
    alpha = angle_between_vectors(v[1], v[2])  # angle between b and c
    beta = angle_between_vectors(v[0], v[2])   # angle between a and c
    gamma = angle_between_vectors(v[0], v[1])  # angle between a and b
    
    # Helper functions for comparisons
    def is_equal(val1, val2, tol=tolerance):
        return abs(val1 - val2) <= tol
    
    def is_angle_equal(angle1, angle2, tol=angle_tolerance):
        return abs(angle1 - angle2) <= tol
    
    def is_right_angle(angle, tol=angle_tolerance):
        return is_angle_equal(angle, 90.0, tol)
    
    def is_120_angle(angle, tol=angle_tolerance):
        return is_angle_equal(angle, 120.0, tol)
    
    def is_60_angle(angle, tol=angle_tolerance):
        return is_angle_equal(angle, 60.0, tol)
    
    # Check for specific lattice patterns
    
    # 1. Cubic patterns (CUB, FCC, BCC)
    if is_equal(a, b) and is_equal(b, c):
        # All lengths equal
        if is_right_angle(alpha) and is_right_angle(beta) and is_right_angle(gamma):
            return LatticeTypeEnum.CUB
        elif is_60_angle(alpha) and is_60_angle(beta) and is_60_angle(gamma):
            # FCC primitive cell: equal lengths, 60° angles
            return LatticeTypeEnum.FCC
        elif is_angle_equal(alpha, beta) and is_angle_equal(beta, gamma):
            # Could be BCC or RHL - check vector patterns
            # BCC primitive vectors have specific pattern
            if _is_bcc_pattern(v, tolerance):
                return LatticeTypeEnum.BCC
            else:
                return LatticeTypeEnum.RHL
    
    # 2. Tetragonal patterns (TET, BCT)
    if is_equal(a, b) and not is_equal(a, c):
        if is_right_angle(alpha) and is_right_angle(beta) and is_right_angle(gamma):
            return LatticeTypeEnum.TET
        elif _is_bct_pattern(v, tolerance):
            return LatticeTypeEnum.BCT
    
    # 3. Hexagonal pattern
    if is_equal(a, b) and not is_equal(a, c):
        if is_right_angle(alpha) and is_right_angle(beta) and is_120_angle(gamma):
            return LatticeTypeEnum.HEX
        elif _is_hex_pattern(v, tolerance):
            return LatticeTypeEnum.HEX
    
    # 4. Orthorhombic patterns (ORC, ORCF, ORCI, ORCC)
    if not is_equal(a, b) and not is_equal(b, c) and not is_equal(a, c):
        if is_right_angle(alpha) and is_right_angle(beta) and is_right_angle(gamma):
            return LatticeTypeEnum.ORC
        elif _is_orcf_pattern(v, tolerance):
            return LatticeTypeEnum.ORCF
        elif _is_orci_pattern(v, tolerance):
            return LatticeTypeEnum.ORCI
        elif _is_orcc_pattern(v, tolerance):
            return LatticeTypeEnum.ORCC
    
    # 5. Monoclinic patterns (MCL, MCLC)
    if is_right_angle(alpha) and is_right_angle(gamma) and not is_right_angle(beta):
        if _is_mclc_pattern(v, tolerance):
            return LatticeTypeEnum.MCLC
        else:
            return LatticeTypeEnum.MCL
    
    # 6. Rhombohedral (already checked above in cubic section)
    if is_equal(a, b) and is_equal(b, c) and is_angle_equal(alpha, beta) and is_angle_equal(beta, gamma):
        return LatticeTypeEnum.RHL
    
    # Default to triclinic if no pattern matches
    return LatticeTypeEnum.TRI


def _is_bcc_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match BCC primitive cell pattern."""
    # BCC primitive: [-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]
    v = vectors
    a_est = np.linalg.norm(v[0])
    
    # Expected BCC pattern (normalized)
    expected = np.array([
        [-0.5, 0.5, 0.5],
        [0.5, -0.5, 0.5],
        [0.5, 0.5, -0.5]
    ]) * a_est
    
    # Check if current vectors match expected pattern (allowing for rotations)
    return _vectors_match_pattern(v, expected, tolerance)


def _is_bct_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match BCT primitive cell pattern."""
    # BCT primitive: [-a/2, a/2, c/2], [a/2, -a/2, c/2], [a/2, a/2, -c/2]
    v = vectors
    lengths = [np.linalg.norm(vec) for vec in v]
    a_est = min(lengths)
    c_est = max(lengths)
    
    expected = np.array([
        [-a_est/2, a_est/2, c_est/2],
        [a_est/2, -a_est/2, c_est/2],
        [a_est/2, a_est/2, -c_est/2]
    ])
    
    return _vectors_match_pattern(v, expected, tolerance)


def _is_hex_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match HEX primitive cell pattern."""
    # HEX primitive: [a/2, -a*sqrt(3)/2, 0], [a/2, a*sqrt(3)/2, 0], [0, 0, c]
    v = vectors
    
    # Find the vector along z-axis (should be [0, 0, c])
    z_vector_idx = None
    for i, vec in enumerate(v):
        if abs(vec[0]) < tolerance and abs(vec[1]) < tolerance and abs(vec[2]) > tolerance:
            z_vector_idx = i
            break
    
    if z_vector_idx is None:
        return False
    
    c_est = abs(v[z_vector_idx][2])
    
    # Check other two vectors
    other_indices = [i for i in range(3) if i != z_vector_idx]
    v1, v2 = v[other_indices[0]], v[other_indices[1]]
    
    # They should have equal lengths and specific pattern
    if not abs(np.linalg.norm(v1) - np.linalg.norm(v2)) < tolerance:
        return False
    
    a_est = np.linalg.norm(v1)
    
    # Check if they match hexagonal pattern
    expected_1 = np.array([a_est/2, -a_est*math.sqrt(3)/2, 0])
    expected_2 = np.array([a_est/2, a_est*math.sqrt(3)/2, 0])
    
    return (_vector_matches(v1, expected_1, tolerance) and _vector_matches(v2, expected_2, tolerance)) or \
           (_vector_matches(v1, expected_2, tolerance) and _vector_matches(v2, expected_1, tolerance))


def _is_orcf_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match ORCF primitive cell pattern."""
    # ORCF primitive: [0, b/2, c/2], [a/2, 0, c/2], [a/2, b/2, 0]
    v = vectors
    
    # Each vector should have one zero component
    zero_components = []
    for vec in v:
        zero_count = sum(1 for x in vec if abs(x) < tolerance)
        if zero_count != 1:
            return False
        zero_components.append([i for i, x in enumerate(vec) if abs(x) < tolerance][0])
    
    # Should have one zero in each dimension
    return set(zero_components) == {0, 1, 2}


def _is_orci_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match ORCI primitive cell pattern."""
    # ORCI primitive: [-a/2, b/2, c/2], [a/2, -b/2, c/2], [a/2, b/2, -c/2]
    # Similar to BCC but with different lengths
    v = vectors
    
    # Check if it has the body-centered pattern with different lengths
    lengths = [np.linalg.norm(vec) for vec in v]
    if len(set(np.round(lengths, 3))) == 1:  # All equal lengths -> not ORCI
        return False
    
    # Check sign pattern
    sign_patterns = []
    for vec in v:
        signs = [1 if x > tolerance else (-1 if x < -tolerance else 0) for x in vec]
        sign_patterns.append(signs)
    
    # Should have body-centered sign pattern
    expected_patterns = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
    return _sign_patterns_match(sign_patterns, expected_patterns)


def _is_orcc_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match ORCC primitive cell pattern."""
    # ORCC primitive: [a/2, b/2, 0], [-a/2, b/2, 0], [0, 0, c]
    v = vectors
    
    # One vector should be along z-axis
    z_vector_count = sum(1 for vec in v if abs(vec[0]) < tolerance and abs(vec[1]) < tolerance)
    if z_vector_count != 1:
        return False
    
    # Other two should be in xy-plane with specific pattern
    xy_vectors = [vec for vec in v if not (abs(vec[0]) < tolerance and abs(vec[1]) < tolerance)]
    if len(xy_vectors) != 2:
        return False
    
    # Check if they have opposite x-components and same y-components
    v1, v2 = xy_vectors
    return (abs(v1[0] + v2[0]) < tolerance and abs(v1[1] - v2[1]) < tolerance and 
            abs(v1[2]) < tolerance and abs(v2[2]) < tolerance)


def _is_mclc_pattern(vectors: np.ndarray, tolerance: float) -> bool:
    """Check if vectors match MCLC primitive cell pattern."""
    # MCLC primitive: [a/2, b/2, 0], [-a/2, b/2, 0], [0, c*cos(alpha), c*sin(alpha)]
    v = vectors
    
    # Similar to ORCC but with one vector at an angle
    xy_vectors = []
    angled_vectors = []
    
    for vec in v:
        if abs(vec[2]) < tolerance:
            xy_vectors.append(vec)
        else:
            angled_vectors.append(vec)
    
    return len(xy_vectors) == 2 and len(angled_vectors) == 1


def _vectors_match_pattern(v1: np.ndarray, v2: np.ndarray, tolerance: float) -> bool:
    """Check if two sets of vectors match within tolerance (allowing rotations)."""
    # Simple check: compare sorted lengths and angles
    lengths1 = sorted([np.linalg.norm(vec) for vec in v1])
    lengths2 = sorted([np.linalg.norm(vec) for vec in v2])
    
    if not all(abs(l1 - l2) < tolerance for l1, l2 in zip(lengths1, lengths2)):
        return False
    
    # Check angles between vectors
    def get_angles(vectors):
        angles = []
        for i in range(3):
            for j in range(i+1, 3):
                cos_angle = np.dot(vectors[i], vectors[j]) / (np.linalg.norm(vectors[i]) * np.linalg.norm(vectors[j]))
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angles.append(np.degrees(np.arccos(cos_angle)))
        return sorted(angles)
    
    angles1 = get_angles(v1)
    angles2 = get_angles(v2)
    
    return all(abs(a1 - a2) < tolerance * 10 for a1, a2 in zip(angles1, angles2))  # More lenient for angles


def _vector_matches(v1: np.ndarray, v2: np.ndarray, tolerance: float) -> bool:
    """Check if two vectors match within tolerance."""
    return np.allclose(v1, v2, atol=tolerance)


def _sign_patterns_match(patterns1: List[List[int]], patterns2: List[List[int]]) -> bool:
    """Check if sign patterns match (allowing permutations)."""
    from itertools import permutations
    
    for perm in permutations(patterns2):
        if patterns1 == list(perm):
            return True
    return False
