/*
 ApproxMC and gen_n_samples

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

#ifndef APPMCCONFIG
#define APPMCCONFIG

#include <vector>

struct AppMCConfig {
    uint32_t startiter = 0;
    double epsilon = 0.80;
    double delta = 0.2;
    uint32_t num_threads = 1;
    int sparse = 0;
    unsigned verb = 1;
    unsigned verb_appmc_cls = 0;
    uint32_t seed = 1;
    bool only_indep_samples = true;
    uint32_t multisample = 1;
    uint32_t samples = 0;
    int simplify = 1;
    double var_elim_ratio = 1.6;
    int reuse_models = 1;
    int force_sol_extension = 0;
    std::vector<uint32_t> sampling_set;
    double kappa = 0.638;      /* Corresponds to epsilon=16 */
    std::string sampleFilename;
    std::string logfilename = "";
    int cms_detach_xor = 1;
};

#endif //APPMCCONFIG
