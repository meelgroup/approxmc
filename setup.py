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


import sys
import os
from setuptools import Extension, setup
import sysconfig
import toml
import pathlib
from sys import platform

def _parse_toml(pyproject_path):
    pyproject_text = pyproject_path.read_text()
    pyproject_data = toml.loads(pyproject_text)
    return pyproject_data['project']['version']


picosatlib = ('picosatlib', {
    'sources': [
               "python/cryptominisat/src/picosat/picosat.c",
               "python/cryptominisat/src/picosat/version.c"],
    'language' : "c",
    'include_dirs' : ["python/cryptominisat/src/"]
    })


def gen_modules(version):
    if platform == "win32" or platform == "cygwin":
        extra_compile_args_val = ['/std:c++17', "/DCMS_LOCAL_BUILD=1", "/DAPPMC_FULL_VERSION=\""+version+"\""]
        define_macros_val = [("TRACE", "")]

    else:
        extra_compile_args_val = ['-std=c++17']
        define_macros_val = [('CMS_LOCAL_BUILD', 1),("TRACE", ""),("APPMC_FULL_VERSION", "\""+version+"\"")]

    modules = Extension(
        name = "pyapproxmc",
        sources = [
                   "python/src/GitSHA1.cpp",
                   "python/src/pyapproxmc.cpp",
                   "src/approxmc.cpp",
                   "src/appmc_constants.cpp",
                   "src/counter.cpp",
                   "python/cryptominisat/python/src/GitSHA1.cpp",
                   "python/cryptominisat/src/bva.cpp",
                   "python/cryptominisat/src/cardfinder.cpp",
                   "python/cryptominisat/src/ccnr_cms.cpp",
                   "python/cryptominisat/src/ccnr.cpp",
                   "python/cryptominisat/src/clauseallocator.cpp",
                   "python/cryptominisat/src/clausecleaner.cpp",
                   "python/cryptominisat/src/cnf.cpp",
                   "python/cryptominisat/src/completedetachreattacher.cpp",
                   "python/cryptominisat/src/cryptominisat_c.cpp",
                   "python/cryptominisat/src/cryptominisat.cpp",
                   "python/cryptominisat/src/datasync.cpp",
                   "python/cryptominisat/src/distillerbin.cpp",
                   "python/cryptominisat/src/distillerlitrem.cpp",
                   "python/cryptominisat/src/distillerlong.cpp",
                   "python/cryptominisat/src/distillerlongwithimpl.cpp",
                   "python/cryptominisat/src/frat.cpp",
                   "python/cryptominisat/src/gatefinder.cpp",
                   "python/cryptominisat/src/gaussian.cpp",
                   "python/cryptominisat/src/get_clause_query.cpp",
                   "python/cryptominisat/src/hyperengine.cpp",
                   "python/cryptominisat/src/intree.cpp",
                   "python/cryptominisat/src/lucky.cpp",
                   "python/cryptominisat/src/matrixfinder.cpp",
                   "python/cryptominisat/src/occsimplifier.cpp",
                   "python/cryptominisat/src/packedrow.cpp",
                   "python/cryptominisat/src/propengine.cpp",
                   "python/cryptominisat/src/reducedb.cpp",
                   "python/cryptominisat/src/sccfinder.cpp",
                   "python/cryptominisat/src/searcher.cpp",
                   "python/cryptominisat/src/searchstats.cpp",
                   "python/cryptominisat/src/sls.cpp",
                   "python/cryptominisat/src/solutionextender.cpp",
                   "python/cryptominisat/src/solverconf.cpp",
                   "python/cryptominisat/src/solver.cpp",
                   "python/cryptominisat/src/str_impl_w_impl.cpp",
                   "python/cryptominisat/src/subsumeimplicit.cpp",
                   "python/cryptominisat/src/subsumestrengthen.cpp",
                   "python/cryptominisat/src/varreplacer.cpp",
                   "python/cryptominisat/src/xorfinder.cpp",
                   "python/cryptominisat/src/oracle/oracle.cpp",
                   "python/arjun/src/arjun.cpp",
                   "python/arjun/src/backward.cpp",
                   "python/arjun/src/common.cpp",
                   "python/arjun/python/src/GitSHA1.cpp",
                   "python/arjun/src/simplify.cpp",
               ],
        extra_compile_args = extra_compile_args_val,
        define_macros = define_macros_val,
        include_dirs = ["src/", "python/cryptominisat/src/", "python/arjun/src/"],
        language = "c++",
    )
    return modules

if __name__ == '__main__':
    pyproject_path = pathlib.Path('pyproject.toml')
    version = _parse_toml(pyproject_path)
    modules = gen_modules(version)
    setup(
        ext_modules =  [modules],
        libraries = [picosatlib],
    )
