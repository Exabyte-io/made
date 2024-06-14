import numpy as np
from constants import ATOMIC_COORD_UNITS, UNITS
from lattice import Lattice  # Assuming you have a similar Lattice class in Python
from paths import paths  # Assuming this imports a dictionary named paths
from symmetry_points import symmetry_points  # Assuming this imports a function
from scipy.spatial.distance import cdist


class ReciprocalLattice(Lattice):
    @property
    def reciprocal_vectors(self):
        vectors_ = self.vectors.vector_arrays
        a = np.linalg.norm(vectors_[0])
        divider = np.dot(vectors_[0], np.cross(vectors_[1], vectors_[2])) / a
        return [
            np.cross(vectors_[1], vectors_[2]) / divider,
            np.cross(vectors_[2], vectors_[0]) / divider,
            np.cross(vectors_[0], vectors_[1]) / divider,
        ]

    @property
    def reciprocal_vector_norms(self):
        return [np.linalg.norm(vec) for vec in self.reciprocal_vectors]

    @property
    def reciprocal_vector_ratios(self):
        norms = self.reciprocal_vector_norms
        max_norm = max(norms)
        return [n / max_norm for n in norms]

    def get_cartesian_coordinates(self, point):
        return np.dot(point, self.reciprocal_vectors)

    @property
    def symmetry_points(self):
        return symmetry_points(self)

    @property
    def default_kpoint_path(self):
        return paths.get(self.type_extended, paths.get(self.type))

    def extract_kpoint_path(self, data_points=[]):
        kpoint_path = []
        symm_points = self.symmetry_points

        for index, point in enumerate(data_points):
            # Find the symmetry point that almost equals the data point
            for symm_point in symm_points:
                if np.allclose(symm_point['coordinates'], point, atol=1e-4):
                    kpoint_path.append({
                        'point': symm_point['point'],
                        'steps': index,
                        'coordinates': symm_point['coordinates'],
                    })
                    break
        return kpoint_path

    def calculate_dimension(self, n_points, index):
        norms = self.reciprocal_vector_norms
        other_indices = [0, 1, 2]
        other_indices.remove(index)
        j, k = other_indices
        N = np.cbrt((n_points * norms[index] ** 2) / (norms[j] * norms[k]))
        return max(1, np.ceil(N))

    def get_dimensions_from_points_count(self, n_kpoints):
        return [self.calculate_dimension(n_kpoints, i) for i in range(3)]

    @property
    def conversion_table(self):
        a = self.a
        return {
            ATOMIC_COORD_UNITS['cartesian']: {
                UNITS['angstrom']: (2 * np.pi) / a,
            },
            UNITS['angstrom']: {
                ATOMIC_COORD_UNITS['cartesian']: a / (2 * np.pi),
            },
        }

    def get_dimensions_from_spacing(self, spacing, units=ATOMIC_COORD_UNITS['cartesian']):
        factor = self.conversion_table[units][ATOMIC_COORD_UNITS['cartesian']] or 1
        return [max(1, np.ceil(np.round(norm / (spacing * factor), 4))) for norm in self.reciprocal_vector_norms]

    def get_spacing_from_dimensions(self, dimensions, units=ATOMIC_COORD_UNITS['cartesian']):
        factor = self.conversion_table[ATOMIC_COORD_UNITS['cartesian']][units] or 1
        norms = self.reciprocal_vector_norms
        return factor * np.mean([norms[i] / max(1, dim) for i, dim in enumerate(dimensions)])
