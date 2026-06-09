import json
import os

_CONSTANTS_FILE = os.path.join(os.path.dirname(__file__), "../../../../constants.json")
with open(_CONSTANTS_FILE) as f:
    _CONSTANTS = json.load(f)

DEFAULT_NON_PERIODIC_MIN_LATTICE_SIZE = _CONSTANTS["molecule"]["nonPeriodicMinimumLatticeSize"]
DIATOMIC_LATTICE_PADDING_FACTOR = _CONSTANTS["molecule"]["diatomicLatticePaddingFactor"]
MOLECULAR_LATTICE_PADDING_FACTOR = _CONSTANTS["molecule"]["molecularLatticePaddingFactor"]
PRECISION_MAP = _CONSTANTS["precision"]
