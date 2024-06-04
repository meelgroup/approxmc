/*
 ApproxMC

 Copyright (c) 2019, Mate Soos and Kuldeep S. Meel. All rights reserved
 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#pragma once

#include <cstdint>
#include <vector>
#include <string>

using std::vector;
using std::string;

#define verb_print(a, b) if (conf.verb >= a) cout << "c o " << b << endl
#define clear_toclear_seen() \
    do {\
      for(const auto& x: to_clear) seen[x] = 0;\
      to_clear.clear();} while (0)

namespace AppMCInt {

struct VarMap
{
    VarMap() {}
    VarMap(uint32_t _vars_to_inclusive, vector<uint32_t> _index_var_map) :
        vars_to_inclusive(_vars_to_inclusive),
        index_var_map(_index_var_map)
    {}

    uint32_t vars_to_inclusive = 0;
    vector<uint32_t> index_var_map;
};

class Constants
{
public:
    Constants();
    vector<double> probval;
    vector<VarMap> index_var_maps;

private:
    vector<string> sparseprobvalues;
    void readInSparseValues();
};

}
