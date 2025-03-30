__path__ = __import__("pkgutil").extend_path(__path__, __name__)
# otherwise, `mat3ra.utils` path leads to an empty __init__.py file in the code.py package
