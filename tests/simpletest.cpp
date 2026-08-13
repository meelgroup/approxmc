/******************************************
Copyright (c) 2020, Mate Soos

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

#include "gtest/gtest.h"

#include "approxmc.h"
#include "cryptominisat5/solvertypesmini.h"
#include "test_helper.h"
#include <arjun/arjun.h>
#include <memory>
#include <cmath>
#include <string>
#include <vector>
using std::string;
using std::vector;

using namespace ApproxMC;
#include <vector>
using std::vector;


TEST(normal_interface, start)
{
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
    AppMC s(fg);
    s.set_sampl_vars({});
    SolCount c = s.count();
    EXPECT_EQ(1U, c.cellSolCount);
    EXPECT_EQ(0U, c.hashCount);
}

TEST(normal_interface, example1)
{
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
    AppMC s(fg);
    s.new_vars(2);
    s.add_clause(str_to_cl("-1, 2"));
    s.set_sampl_vars({0, 1});
    SolCount c = s.count();
    EXPECT_EQ(0U, c.hashCount);
    EXPECT_EQ(3U, c.cellSolCount);
}

TEST(normal_interface, example2)
{
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
    AppMC s(fg);
    s.new_vars(10);
    vector<uint32_t> sampl;
    for(uint32_t i = 0; i < 10; i++) sampl.push_back(i);
    s.set_sampl_vars(sampl);

    SolCount c = s.count();
    uint32_t x = std::pow(2, c.hashCount)*c.cellSolCount;
    EXPECT_EQ(std::pow(2, 10), x);
}

TEST(normal_interface, example3)
{
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
    AppMC s(fg);
    s.new_vars(10);
    vector<uint32_t> sampl;
    for(uint32_t i = 0; i < 10; i++) sampl.push_back(i);
    s.set_sampl_vars(sampl);

    s.add_clause(str_to_cl("-3"));
    SolCount c = s.count();
    uint32_t cnt = std::pow(2, c.hashCount)*c.cellSolCount;
    EXPECT_EQ(std::pow(2, 9), cnt);
}

TEST(normal_interface, example4)
{
    std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
    AppMC s(fg);
    s.new_vars(10);
    vector<uint32_t> sampl;
    for(uint32_t i = 0; i < 10; i++) sampl.push_back(i);
    s.set_sampl_vars(sampl);

    s.add_clause(str_to_cl("-3, 4"));
    s.add_clause(str_to_cl("3, -4"));
    SolCount c = s.count();
    uint32_t cnt = std::pow(2, c.hashCount)*c.cellSolCount;
    EXPECT_EQ(std::pow(2, 9), cnt);
}

// set_start_iter only seeds where the level search begins. Seeding it at or
// above the level the search would find must give the same answer as seeding it
// below, or not seeding it at all.
TEST(normal_interface, start_iter)
{
    auto run = [](uint32_t start_iter) {
        std::unique_ptr<FieldGen> fg = std::make_unique<ArjunNS::FGenMpq>();
        AppMC s(fg);
        s.set_seed(1);
        s.set_epsilon(0.8);
        s.set_delta(0.4); // one measurement, so a lost one leaves nothing
        s.set_start_iter(start_iter);
        s.new_vars(30);
        s.add_clause(str_to_cl("1, 2"));
        vector<uint32_t> sampl;
        for(uint32_t i = 0; i < 30; i++) sampl.push_back(i);
        s.set_sampl_vars(sampl);
        return s.count();
    };

    const SolCount ref = run(0);
    EXPECT_TRUE(ref.valid);
    EXPECT_GT(ref.hashCount, 2U);

    for(uint32_t offset = 0; offset <= 4; offset += 2) {
        const SolCount c = run(ref.hashCount + offset);
        EXPECT_TRUE(c.valid);
        EXPECT_EQ(ref.hashCount, c.hashCount);
        EXPECT_EQ(ref.cellSolCount, c.cellSolCount);
    }
    const SolCount below = run(ref.hashCount - 2);
    EXPECT_TRUE(below.valid);
    EXPECT_EQ(ref.hashCount, below.hashCount);
    EXPECT_EQ(ref.cellSolCount, below.cellSolCount);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
