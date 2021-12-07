from distutils.core import setup, Extension

pyapproxmc_module = Extension(
    'pyapproxmc',
    sources=['src/pyapproxmc.cpp'],
    libraries=['approxmc', 'cryptominisat5'],
    language='C++', )

setup(
    name='pyapproxmc',
    ext_modules=[pyapproxmc_module]
)
