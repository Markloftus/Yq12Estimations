from setuptools import setup
from Cython.Build import cythonize

setup(
    name="SRKmer",
    ext_modules=cythonize(
        "mainKmer.pyx",
        compiler_directives={"language_level": "3"},
    ),
)
