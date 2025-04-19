#
# ApproxMC
#
# Copyright (c) 2009-2017, Mate Soos. All rights reserved.
# Copyright (c) 2017, Pierre Vignet
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from setuptools import Extension, setup
import toml
import pathlib
from sys import platform

def _parse_toml(pyproject_path):
    pyproject_text = pyproject_path.read_text()
    pyproject_data = toml.loads(pyproject_text)
    return pyproject_data['project']['version']


def gen_modules(version):
    define_macros_val : list[tuple[str, str | None]] | None
    if platform == "win32" or platform == "cygwin":
        extra_compile_args_val = ['/std:c++17', "/DAPPMC_FULL_VERSION=\""+version+"\""]
        define_macros_val = [("TRACE", "")]

    else:
        extra_compile_args_val = ['-std=c++17']
        define_macros_val = [("TRACE", ""),("APPMC_FULL_VERSION", "\""+version+"\"")]

    modules = Extension(
        name = "pyapproxmc",
        sources = [
               "python/src/GitSHA1.cpp",
               "python/src/pyapproxmc.cpp",
               "src/approxmc.cpp",
               "src/appmc_constants.cpp",
               "src/counter.cpp",
           ],
        extra_compile_args = extra_compile_args_val,
        define_macros = define_macros_val,
        include_dirs = ["src/", "python/cryptominisat/src/", "python/arjun/src/", "python/sbva/src/", ],
        language = "c++",
    )
    return modules

if __name__ == '__main__':
    pyproject_path = pathlib.Path('pyproject.toml')
    version = _parse_toml(pyproject_path)
    modules = gen_modules(version)
    setup(
        ext_modules =  [modules],
    )
