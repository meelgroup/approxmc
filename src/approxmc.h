/*
 ApproxMC

 Copyright (c) 2019-2020, Mate Soos and Kuldeep S. Meel. All rights reserved
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


#ifndef APPROXMC_H__
#define APPROXMC_H__

#include "cryptominisat5/cryptominisat.h"


namespace ApproxMC {

struct AppMCPrivateData;
#ifdef _WIN32
class __declspec(dllexport) AppMC
#else
class AppMC
#endif
{
public:
    AppMC();
    ~AppMC();
    string get_version_info();
    void set_up_log(string log_file_name);
    void set_verbosity(uint32_t verb);
    void setup();
    void count();
    void set_sampling_set(const vector<uint32_t>& vars);

    uint32_t nVars();
    void new_vars(uint32_t num);
    void new_var();
    void add_clause(const vector<CMSat::Lit>& lits);
    void add_xor_clause(const vector<uint32_t>& vars, bool rhs);


private:
    ////////////////////////////
    // Do not bother with this, it's private
    ////////////////////////////
    AppMCPrivateData* data;
};

}

#endif
