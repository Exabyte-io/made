from typing import Tuple, List
from pydantic import BaseModel

from ..slab.termination import Termination


class TerminationPair(BaseModel):
    film_termination: Termination
    substrate_termination: Termination

    def __str__(self):
        # return str(self.to_tuple_of_str())
        return f"({self.film_termination}, {self.substrate_termination})"

    def __repr__(self):
        return self.__str__()

    def __init__(self, film_termination: Termination, substrate_termination: Termination):
        super().__init__(film_termination=film_termination, substrate_termination=substrate_termination)

    @classmethod
    def from_tuple_of_str(cls, termination_pair: Tuple[str, str]):
        film_termination = Termination.from_string(termination_pair[0])
        substrate_termination = Termination.from_string(termination_pair[1])
        return cls(film_termination=film_termination, substrate_termination=substrate_termination)

    def to_tuple_of_str(self) -> Tuple[str, str]:
        return str(self.film_termination), str(self.substrate_termination)

    @classmethod
    def from_pymatgen(cls, termination_pair: Tuple[str, str]):
        return cls.from_tuple_of_str(termination_pair)

    def to_pymatgen(self) -> Tuple[str, str]:
        return self.to_tuple_of_str()


def safely_select_termination_pair(
    provided_termination_pair: TerminationPair, generated_termination_pairs: List[TerminationPair]
) -> TerminationPair:
    """
    Attempt finding provided in generated terminations to find a complete match,
    if match isn't found, get terminations with equivalent chemical elements,
    if that fails, return the first generated termination pair.
    """
    provided_film_termination = provided_termination_pair.film_termination
    provided_substrate_termination = provided_termination_pair.substrate_termination
    hotfix_termination_pair = provided_termination_pair
    if provided_termination_pair not in generated_termination_pairs:
        for termination_pair in generated_termination_pairs:
            generated_film_termination = termination_pair.film_termination
            generated_substrate_termination = termination_pair.substrate_termination
            if (
                generated_film_termination.chemical_elements == provided_film_termination.chemical_elements
                and generated_substrate_termination.chemical_elements
                == provided_substrate_termination.chemical_elements
            ):
                hotfix_termination_pair = termination_pair
            else:
                hotfix_termination_pair = generated_termination_pairs[0]
            print("Interface will be built with terminations: ", hotfix_termination_pair)

    return hotfix_termination_pair
