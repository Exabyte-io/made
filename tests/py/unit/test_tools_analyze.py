import numpy as np
from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.other import get_average_interlayer_distance, get_surface_area
from mat3ra.made.tools.analyze.rdf import RadialDistributionFunction
from unit.fixtures.generated.fixtures import INTERFACE_ATOMS

from .fixtures.nanoribbon import GRAPHENE_ZIGZAG_NANORIBBON


def test_calculate_average_interlayer_distance():
    distance = get_average_interlayer_distance(INTERFACE_ATOMS, 1, 2)
    assert np.isclose(distance, 4.0725)


def test_calculate_surface_area():
    atoms = bulk("Si", cubic=False)
    area = get_surface_area(atoms)
    assert np.isclose(area, 12.7673)


def test_radial_distribution_function():
    material = Material.create(GRAPHENE_ZIGZAG_NANORIBBON)

    rdf = RadialDistributionFunction.from_material(material, cutoff=10.0, bin_size=0.1)

    # Test that RDF is non-negative
    assert np.all(rdf.rdf >= 0), "RDF contains negative values."

    # Test the first peak properties
    assert rdf.first_peak_index > 0, "First peak index should be greater than 0."
    assert rdf.first_peak_value > 0, "First peak value should be greater than 0."
    assert rdf.first_peak_width > 0, "First peak width should be greater than 0."
    assert rdf.first_peak_distance > 0, "First peak distance should be greater than 0."

    # Test that `is_within_first_peak` works as expected
    assert rdf.is_within_first_peak(rdf.first_peak_distance), "First peak distance should be within the first peak."
    assert not rdf.is_within_first_peak(0), "Zero distance should not be within the first peak."
    assert not rdf.is_within_first_peak(10.0), "Cutoff distance should not be within the first peak."

    # Test that RDF drops to zero near the cutoff
    assert np.isclose(rdf.rdf[-1], 0, atol=1e-2), "RDF should approach zero near the cutoff."

    # Specific material related tests
    assert np.isclose(rdf.first_peak_distance, 1.42, atol=0.1), "First peak distance should be close to 1.42."
