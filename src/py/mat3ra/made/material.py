import hashlib
from typing import Any, List, Dict, Union

from @mat3ra.code.dist.js.entity import HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity
from @mat3ra.code.dist.js.entity.in_memory import AnyObject
from @mat3ra.code.dist.js.types import ConsistencyCheck, DerivedPropertiesSchema, FileSourceSchema, InChIRepresentationSchema, MaterialSchema
from crypto-js import CryptoJS
from .basis.constrained_basis import ConstrainedBasis
from .cell.conventional_cell import isConventionalCellSameAsPrimitiveForLatticeType, PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES, PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS
from .constants import ATOMIC_COORD_UNITS, units
from .constraints.constraints import Constraint
from .lattice.lattice import Lattice
from .lattice.lattice_vectors import BravaisConfigProps
from .parsers.parsers import parsers
from .tools.supercell import supercellTools
from .types import MaterialJSON

defaultMaterialConfig = {
    "name": "Silicon FCC",
    "basis": {
        "elements": [
            {
                "id": 1,
                "value": "Si",
            },
            {
                "id": 2,
                "value": "Si",
            },
        ],
        "coordinates": [
            {
                "id": 1,
                "value": [0.0, 0.0, 0.0],
            },
            {
                "id": 2,
                "value": [0.25, 0.25, 0.25],
            },
        ],
        "units": ATOMIC_COORD_UNITS.crystal,
    },
    "lattice": {
        "type": "FCC",
        "a": 3.867,
        "b": 3.867,
        "c": 3.867,
        "alpha": 60,
        "beta": 60,
        "gamma": 60,
        "units": {
            "length": units.angstrom,
            "angle": units.degree,
        },
    },
}

MaterialSchemaJSON = Dict[str, Union[MaterialSchema, AnyObject]]
MaterialBaseEntity = Type[HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity]
MaterialBaseEntityConstructor = Type[MaterialBaseEntity]

def MaterialMixin(superclass: MaterialBaseEntityConstructor) -> MaterialBaseEntityConstructor:
    class MadeMaterial(superclass):
        _json: MaterialSchemaJSON

        def __init__(self, *config: Any) -> None:
            super().__init__(*config)
            self.name = super().name or self.formula

        def toJSON(self) -> MaterialJSON:
            return {
                **super().toJSON(),
                "lattice": self.Lattice.toJSON(),
                "basis": self.Basis.toJSON(),
                "name": self.name,
                "isNonPeriodic": self.isNonPeriodic,
            }

        @property
        def defaultConfig(self) -> Dict[str, Any]:
            return defaultMaterialConfig

        @property
        def src(self) -> FileSourceSchema:
            return self.prop("src") as FileSourceSchema

        @src.setter
        def src(self, src: FileSourceSchema) -> None:
            self.setProp("src", src)

        def updateFormula(self) -> None:
            self.setProp("formula", self.Basis.formula)
            self.setProp("unitCellFormula", self.Basis.unitCellFormula)

        @property
        def isNonPeriodic(self) -> bool:
            return self.prop("isNonPeriodic", False)

        @isNonPeriodic.setter
        def isNonPeriodic(self, bool: bool) -> None:
            self.setProp("isNonPeriodic", bool)

        def getDerivedPropertyByName(self, name: str) -> DerivedPropertiesSchema:
            return next((x for x in self.getDerivedProperties() if x.name == name), None)

        def getDerivedProperties(self) -> DerivedPropertiesSchema:
            return self.prop("derivedProperties", [])

        @property
        def formula(self) -> str:
            return self.prop("formula") or self.Basis.formula

        @property
        def unitCellFormula(self) -> str:
            return self.prop("unitCellFormula") or self.Basis.unitCellFormula

        def setBasis(self, textOrObject: Union[str, BasisConfig], format: str = None, unitz: str = None) -> None:
            basis = None
            if format == "xyz":
                basis = parsers.xyz.toBasisConfig(textOrObject, unitz)
            else:
                basis = textOrObject
            self.setProp("basis", basis)
            self.updateFormula()

        def setBasisConstraints(self, constraints: List[Constraint]) -> None:
            self.setBasis({ **self.basis, "constraints": constraints })

        @property
        def basis(self) -> BasisConfig:
            return self.prop("basis") as BasisConfig

        @property
        def Basis(self) -> ConstrainedBasis:
            return ConstrainedBasis({ **self.basis, "cell": self.Lattice.vectorArrays })

        @property
        def uniqueElements(self) -> List[str]:
            return self.Basis.uniqueElements

        @property
        def lattice(self) -> Union[BravaisConfigProps, None]:
            return self.prop("lattice")

        @lattice.setter
        def lattice(self, config: Union[BravaisConfigProps, None]) -> None:
            self.setProp("lattice", config)

        @property
        def Lattice(self) -> Lattice:
            return Lattice(self.lattice)

        def getInchiStringForHash(self) -> str:
            inchi = self.getDerivedPropertyByName("inchi")
            if inchi:
                return inchi.value
            raise Exception("Hash cannot be created. Missing InChI string in derivedProperties")

        def calculateHash(self, salt: str = "", isScaled: bool = False, bypassNonPeriodicCheck: bool = False) -> str:
            message = ""
            if not self.isNonPeriodic or bypassNonPeriodicCheck:
                message = f"{self.Basis.hashString}#{self.Lattice.getHashString(isScaled)}#{salt}"
            else:
                message = self.getInchiStringForHash()
            return hashlib.md5(message.encode()).hexdigest()

        @property
        def hash(self) -> str:
            return self.prop("hash") as str

        @hash.setter
        def hash(self, hash: str) -> None:
            self.setProp("hash", hash)

        @property
        def scaledHash(self) -> str:
            return self.calculateHash("", True)

        def toCrystal(self) -> None:
            basis = self.Basis
            basis.toCrystal()
            self.setProp("basis", basis.toJSON())

        def toCartesian(self) -> None:
            basis = self.Basis
            basis.toCartesian()
            self.setProp("basis", basis.toJSON())

        def getBasisAsXyz(self, fractional: bool = False) -> str:
            return parsers.xyz.fromMaterial(self.toJSON(), fractional)

        def getAsQEFormat(self) -> str:
            return parsers.espresso.toEspressoFormat(self.toJSON())

        def getAsPOSCAR(self, ignoreOriginal: bool = False, omitConstraints: bool = False) -> str:
            src = self.src
            if src and src.extension == "poscar" and not ignoreOriginal:
                return self.src.text
            return parsers.poscar.toPoscar(self.toJSON(), omitConstraints)

        def getACopyWithConventionalCell(self) -> MadeMaterial:
            material = self.clone()
            if isConventionalCellSameAsPrimitiveForLatticeType(self.Lattice.type):
                return material
            conventionalSupercellMatrix = PRIMITIVE_TO_CONVENTIONAL_CELL_MULTIPLIERS[self.Lattice.type]
            conventionalLatticeType = PRIMITIVE_TO_CONVENTIONAL_CELL_LATTICE_TYPES[self.Lattice.type]
            config = supercellTools.generateConfig(material, conventionalSupercellMatrix)
            config.lattice.type = conventionalLatticeType
            config.name = f"{material.name} - conventional cell"
            return self.constructor(config)

        def getConsistencyChecks(self) -> List[ConsistencyCheck]:
            basisChecks = self.getBasisConsistencyChecks()
            return basisChecks

        def getBasisConsistencyChecks(self) -> List[ConsistencyCheck]:
            checks = []
            limit = 1000
            basis = self.Basis
            if len(self.Basis.elements) < limit:
                overlappingAtomsGroups = basis.getOverlappingAtoms()
                for group in overlappingAtomsGroups:
                    id1, id2, element1, element2 = group
                    checks.extend([
                        {
                            "key": f"basis.coordinates.{id1}",
                            "name": "atomsOverlap",
                            "severity": "warning",
                            "message": f"Atom {element1} is too close to {element2} at position {id2 + 1}",
                        },
                        {
                            "key": f"basis.coordinates.{id2}",
                            "name": "atomsOverlap",
                            "severity": "warning",
                            "message": f"Atom {element2} is too close to {element1} at position {id1 + 1}",
                        },
                    ])
            return checks

    return MadeMaterial

Material = MaterialMixin(HasConsistencyChecksHasMetadataNamedDefaultableInMemoryEntity)


