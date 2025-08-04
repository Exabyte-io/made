import warnings

from pydantic.warnings import PydanticDeprecatedSince20

warnings.filterwarnings("ignore", category=DeprecationWarning, module="typing_extensions")
warnings.filterwarnings("ignore", category=UserWarning, module="pydantic")
warnings.filterwarnings("ignore", category=PydanticDeprecatedSince20)
