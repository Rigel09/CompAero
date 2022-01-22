from importlib_metadata import entry_points
from setuptools import setup, find_packages
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name="CompAero",
    version=get_version(os.path.join("CompAero", "__init__.py")),
    description="A python package for compressible aerodynamics",
    author="Rigel09",
    url="https://github.com/Rigel09/CompAero",
    packages=find_packages(),
    install_requires=["scipy", "numpy", "colorama", "matplotlib",],
    setup_requires=["pytest-runner", "black"],
    tests_require=["pytest"],
    license="MIT",
    long_description=open("README.md", "r", encoding="utf8").read(),
    entry_points={"console_scripts": ["compressible_calculator=CompAero.Calculator.main:main"]},
)
