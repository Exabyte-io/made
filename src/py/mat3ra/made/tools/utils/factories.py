from mat3ra.utils.factory import BaseFactory


class PerturbationFunctionHelperFactory(BaseFactory):
    __class_registry__ = {
        "sine_wave": "mat3ra.made.tools.utils.helpers.SineWave",
    }
