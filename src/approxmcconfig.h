/*
 ApproxMC

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
    uint32_t start_iter = 0;
    uint32_t threshold = 72; //precision
    double epsilon = 0.80;
    uint32_t measurements = 9; //confidence of 0.81
    double delta = 0.2;
    uint32_t num_threads = 1;
    uint32_t startiter = 0;
    bool sparse = false;
    unsigned verb = 1;
    unsigned verb_appmc_cls = 0;
    uint32_t seed = 1;
    std::vector<uint32_t> sampling_set;
    double kappa = 0.638;      /* Corresponds to epsilon=16 */
    std::string logfilename = "";
};

#endif //APPMCCONFIG
