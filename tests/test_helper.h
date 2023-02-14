/******************************************
Copyright (c) 2016, Mate Soos

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
***********************************************/

#ifndef TEST_HELPER__H
#define TEST_HELPER__H

#include <vector>
#include <ostream>
#include <iostream>
#include <sstream>
#include <functional>
#include <cctype>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <cryptominisat5/cryptominisat.h>
#include <cryptominisat5/solvertypesmini.h>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::istringstream;
using std::stringstream;
using namespace CMSat;

class Xor
{
public:
    Xor()
    {}

    explicit Xor(const vector<uint32_t>& cl, const bool _rhs, const vector<uint32_t>& _clash_vars):
        rhs(_rhs)
        , clash_vars(_clash_vars)
    {
        for (uint32_t i = 0; i < cl.size(); i++) {
            vars.push_back(cl[i]);
        }
    }

    template<typename T>
    explicit Xor(const T& cl, const bool _rhs, const vector<uint32_t>& _clash_vars):
        rhs(_rhs)
        , clash_vars(_clash_vars)
    {
        for (uint32_t i = 0; i < cl.size(); i++) {
            vars.push_back(cl[i].var());
        }
    }

    explicit Xor(const vector<uint32_t>& cl, const bool _rhs, const uint32_t clash_var):
        rhs(_rhs)
    {
        clash_vars.push_back(clash_var);
        for (uint32_t i = 0; i < cl.size(); i++) {
            vars.push_back(cl[i]);
        }
    }

    vector<uint32_t>::const_iterator begin() const
    {
        return vars.begin();
    }

    vector<uint32_t>::const_iterator end() const
    {
        return vars.end();
    }

    vector<uint32_t>::iterator begin()
    {
        return vars.begin();
    }

    vector<uint32_t>::iterator end()
    {
        return vars.end();
    }

    bool operator<(const Xor& other) const
    {
        uint64_t i = 0;
        while(i < other.size() && i < size()) {
            if (other[i] != vars[i]) {
                return (vars[i] < other[i]);
            }
            i++;
        }

        if (other.size() != size()) {
            return size() < other.size();
        }
        return false;
    }

    const uint32_t& operator[](const uint32_t at) const
    {
        return vars[at];
    }

    uint32_t& operator[](const uint32_t at)
    {
        return vars[at];
    }

    void resize(const uint32_t newsize)
    {
        vars.resize(newsize);
    }

    vector<uint32_t>& get_vars()
    {
        return vars;
    }

    const vector<uint32_t>& get_vars() const
    {
        return vars;
    }

    size_t size() const
    {
        return vars.size();
    }

    bool empty() const
    {
        if (!vars.empty())
            return false;

        if (!clash_vars.empty())
            return false;

        if (rhs != false) {
            return false;
        }

        return true;
    }

    void merge_clash(const Xor& other, vector<uint16_t>& seen) {
        for(const auto& v: clash_vars) {
            seen[v] = 1;
        }

        for(const auto& v: other.clash_vars) {
            if (!seen[v]) {
                seen[v] = 1;
                clash_vars.push_back(v);
            }
        }

        for(const auto& v: clash_vars) {
            seen[v] = 0;
        }
    }


    bool rhs = false;
    vector<uint32_t> clash_vars;
    bool detached = false;
    vector<uint32_t> vars;
};

inline std::ostream& operator<<(std::ostream& os, const Xor& thisXor)
{
    for (uint32_t i = 0; i < thisXor.size(); i++) {
        os << Lit(thisXor[i], false);

        if (i+1 < thisXor.size())
            os << " + ";
    }
    os << " =  " << std::boolalpha << thisXor.rhs << std::noboolalpha;

    os << " -- clash: ";
    for(const auto& c: thisXor.clash_vars) {
        os << c+1 << ", ";
    }

    return os;
}

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

long int str_to_long_int(string& token)
{
    string trimmed = trim(token);
    size_t endptr;
    long i = std::stol(trimmed, &endptr);
    if (endptr != trimmed.size()) {
        cout << "Error, input token: '" << token << "' wasn't completely used up, wrong token!" << endl;
        exit(-1);
    }
    return i;
}

vector<Lit> str_to_cl(const string& data)
{
    vector<string> tokens;
    stringstream ss(data);
    string token;
    while (getline(ss,token, ','))
    {
        tokens.push_back(token);
    }

    vector<Lit> ret;
    for(string& token2: tokens) {
        long int i = str_to_long_int(token2);
        assert(i == (int)i);
        Lit lit(std::abs(i)-1, i < 0);
        ret.push_back(lit);
    }
    //cout << "input is: " << data << " LITs is: " << ret << endl;

    std::sort(ret.begin(), ret.end());
    return ret;
}

vector<uint32_t> str_to_vars(const string& data)
{
    vector<Lit> lits = str_to_cl(data);
    vector<uint32_t> vars;
    for(Lit lit: lits) {
        assert(lit.sign() == false);
        vars.push_back(lit.var());
    }
    return vars;
}

vector<Xor> str_to_xors(const string& data)
{
    vector<Xor> ret;
    stringstream ss(data);
    string token;
    while (getline(ss,token, ';'))
    {
        stringstream ss2(token);
        string token2;
        int at = 0;
        bool rhs = false;
        vector<uint32_t> vars;
        vector<uint32_t> clashes;
        while (getline(ss2,token2, '='))
        {
            //cout << "Token is: " << token2 << endl;
            if (at == 0) {
                vars = str_to_vars(token2);
            }
            if (at == 1) {
                uint32_t at2 = 0;
                stringstream ss3(token2);
                string token3;
                //cout << "parsing token2:" << token2 << endl;
                while (getline(ss3,token3, 'c')) {
                    if (at2 == 0) {
                        long r = str_to_long_int(token3);
                        assert(r >= 0 && r <= 1);
                        rhs = r;
                    } else if (at2 == 1) {
                        clashes = str_to_vars(token3);
                    }
                    assert(at2 < 2 && "We can have only at most one 'c' sign in an XOR");
                    at2++;
                }
            }
            assert(at < 2 && "We can only have one '=' sign in an XOR");
            at++;
        }

        assert(at == 2 && "You forgot the =0/1 from the XOR");
        ret.push_back(Xor(vars, rhs, clashes));
    }

    return ret;
}

vector<vector<Lit> > str_to_vecs(const string& data)
{
    vector<vector<Lit> > ret;
    stringstream ss(data);
    string token;
    while (getline(ss,token, ';'))
    {
        ret.push_back(str_to_cl(token));
    }

    return ret;
}

struct VecVecSorter
{
    bool operator()(const vector<Lit>&a, const vector<Lit>& b) const
    {
        if (a.size() != b.size()) {
            return a.size() < b.size();
        }

        for(size_t i = 0; i < a.size(); i++) {
            if (a[i] != b[i]) {
                return a[i] < b[i];
            }
        }
        return false;
    }
};

void check_fuzzy_equal(
    vector<vector<Lit> >& cls_expected,
    vector<vector<Lit> >& cls_actual)
{
    for(vector<Lit>& x: cls_actual) {
        std::sort(x.begin(), x.end());
    }
    for(vector<Lit>& x: cls_expected) {
        std::sort(x.begin(), x.end());
    }

    VecVecSorter sorter;
    std::sort(cls_actual.begin(), cls_actual.end(), sorter);
    std::sort(cls_expected.begin(), cls_expected.end(), sorter);

    EXPECT_EQ(cls_expected, cls_actual);
}

string print(const vector<vector<Lit> >& cls)
{
    std::stringstream ss;
    for(auto cl: cls) {
        ss << cl << endl;
    }
    return ss.str();
}

struct XorSorter
{
    bool operator()(const Xor& a, const Xor& b) const
    {
        if (a.size() != b.size())
            return a.size() < b.size();

        if (a.rhs != b.rhs) {
            return a.rhs < b.rhs;
        }

        for(size_t i = 0; i < a.size(); i++) {
            if (a[i] != b[i]) {
                return a[i] < b[i];
            }
        }

        return false;
    }
};


void sort_xor(Xor& x)
{
    std::sort(x.vars.begin(), x.vars.end());
    std::sort(x.clash_vars.begin(), x.clash_vars.end());
}

void check_xors_eq(const vector<Xor>& got_data, const std::string& expected)
{
    XorSorter xorsort;

    vector<Xor> expected_sorted = str_to_xors(expected);
    for(auto t: expected_sorted) {
        std::sort(t.begin(), t.end());
    }
    std::sort(expected_sorted.begin(), expected_sorted.end(), xorsort);

    vector<Xor> got_data_sorted = got_data;
    for(Xor& t: got_data_sorted) {
        sort_xor(t);
    }

    std::sort(got_data_sorted.begin(), got_data_sorted.end(), xorsort);
    EXPECT_EQ(expected_sorted.size(), got_data_sorted.size());
    for(size_t i = 0; i < expected_sorted.size(); i++) {
        EXPECT_EQ(expected_sorted[i].vars, got_data_sorted[i].vars);
        EXPECT_EQ(expected_sorted[i].rhs, got_data_sorted[i].rhs);
        EXPECT_EQ(expected_sorted[i].clash_vars, got_data_sorted[i].clash_vars);
    }
}

void check_xors_contains(const vector<Xor>& got_data, const std::string& expected)
{
    vector<Xor> expected_sorted = str_to_xors(expected);
    assert(expected_sorted.size() == 1);
    Xor expectedX = expected_sorted[0];
    std::sort(expectedX.begin(), expectedX.end());

    vector<Xor> got_data_sorted = got_data;
    for(auto& t: got_data_sorted) {
        sort_xor(t);
    }

    bool found = false;
    for(const Xor& x: got_data_sorted) {
        if (x.vars == expectedX.vars &&
            x.rhs == expectedX.rhs &&
            x.clash_vars == expectedX.clash_vars
        ) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
}

struct cnfdata {
    int64_t num_cls_per_header = -1;
    int64_t num_vars_per_header = -1;
    vector<vector<Lit>> cls;
    uint64_t num_vars = 0;
};

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

cnfdata cnf_file_read(std::string fname)
{
    cnfdata cnfdat;

    std::ifstream file(fname);
    std::string str;
    std::string file_contents;
    vector<Lit> cl;
    while (std::getline(file, str))
    {
        //cout << "CNF LINE: " << str << endl;
        if (str.find("cnf") != string::npos) {
            str.erase(0,5);
            vector<string> s = split(rtrim(ltrim(str)), ' ');
            assert(s.size() == 2);
            cnfdat.num_vars_per_header = std::stoi(s[0]);
            cnfdat.num_cls_per_header = std::stoi(s[1]);
            continue;
        }

        if (str.find("c ") == 0) {
            continue;
        }

        cl.clear();
        vector<string> s = split(rtrim(ltrim(str)), ' ');
        for(string& l: s) {
            if (l.length() == 0)
                continue;

            int x = std::stoi(l);
            if (x == 0) {
                break;
            }
            uint64_t var = std::abs(x)-1;
            cnfdat.num_vars = std::max(cnfdat.num_vars, var+1);
            bool sign = x < 0;
            cl.push_back(Lit(var, sign));
        }
        cnfdat.cls.push_back(cl);
    }
    return cnfdat;
}

bool cl_eq(const vector<Lit>& lits1, const vector<Lit>& lits2)
{
    if (lits1.size() != lits2.size())
        return false;



    vector<Lit> cl1_s = lits1;
    std::sort(cl1_s.begin(), cl1_s.end());

    vector<Lit> cl2_s = lits2;
    std::sort(cl2_s.begin(), cl2_s.end());
    for(size_t i = 0; i < cl1_s.size(); i++) {
        if (cl1_s[i] != cl2_s[i])
            return false;
    }
    return true;
}

bool cl_exists(const vector<vector<Lit> >& cls, const vector<Lit>& cl) {
    for(const vector<Lit>& cli: cls) {
        if (cl_eq(cli, cl)) {
            return true;
        }
    }
    return false;
}

// string print(const vector<Lit>& dat) {
//     std::stringstream m;
//     for(size_t i = 0; i < dat.size();) {
//         m << dat[i];
//         i++;
//         if (i < dat.size()) {
//             m << ", ";
//         }
//     }
//     return m.str();
// }

#endif //TEST_HELPER__H
