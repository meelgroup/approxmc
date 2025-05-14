/*
 ApproxMC

 Copyright (c) 2019, Mate Soos and Kuldeep S. Meel. All rights reserved
 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi

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

#include <memory>
#include <vector>
#include <cstdint>
#include <string>
#include <gmpxx.h>
#include <cryptominisat5/solvertypesmini.h>

namespace AppMCInt {

class Config {
public:
    Config(const std::unique_ptr<CMSat::FieldGen>& _fg) : multiplier_weight(_fg->one()) {}
    uint32_t start_iter = 0;
    double epsilon = 0.80; //Tolerance.  CAV-2020 paper default
    double delta = 0.2;    //Confidence. CAV-2020 paper default
    int sparse = 0;
    unsigned verb = 0;
    unsigned verb_cls = 0;
    uint32_t seed = 1;
    int simplify = 1;
    double var_elim_ratio = 1.6;
    int reuse_models = 1;
    std::string logfilename = "";
    int dump_intermediary_cnf = 0;
    int debug = 0;
    int force_sol_extension = false;
    double appmc7_eps_cutoff = 5;

    std::vector<uint32_t> sampl_vars;
    bool sampl_vars_set = false;
    std::unique_ptr<CMSat::Field> multiplier_weight = nullptr;
};

}
