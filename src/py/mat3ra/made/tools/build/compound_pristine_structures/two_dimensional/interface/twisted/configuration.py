from ..base.configuration import InterfaceConfiguration


class TwistedNanoribbonsInterfaceConfiguration(InterfaceConfiguration):
    """
    Configuration for creating a twisted interface between two nanoribbons with specified twist angle.

    Args:
        stack_components (List[SlabConfiguration]): List of two nanoribbons as slab configurations.
        angle (float): Twist angle in degrees for provenance.
    """

    angle: float = 0.0
