import numpy as np


def are_vectors_colinear(v1, v2, tol=1e-3):
    v1 = np.array(v1)
    v2 = np.array(v2)
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 < tol or norm2 < tol:
        return False
    cos_angle = np.dot(v1, v2) / (norm1 * norm2)
    # Colinear if |cos_angle| close to 1
    return np.abs(np.abs(cos_angle) - 1) < tol
