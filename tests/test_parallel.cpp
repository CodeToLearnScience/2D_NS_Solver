// Unit tests for the j-slab decomposition logic (serial, no MPI needed).

#include <gtest/gtest.h>

#include <set>

#include "parallel/decomposition.hpp"

namespace {

using namespace ns::parallel;

TEST(Decomposition, RangesPartitionGlobalRows) {
    for (int ny : {1, 2, 7, 80, 81}) {
        for (int nranks = 1; nranks <= std::min(ny, 8); ++nranks) {
            auto ranges = SlabDecomposition::all_ranges(ny, nranks);
            ASSERT_EQ(static_cast<int>(ranges.size()), nranks);
            std::set<std::pair<int, int>> uniq;
            int covered = 0;
            for (auto& r : ranges) {
                EXPECT_LT(r.first, r.second) << "empty slab";
                uniq.insert(r);
                covered += r.second - r.first;
            }
            EXPECT_EQ(covered, ny) << "ny=" << ny << " np=" << nranks;
            EXPECT_EQ(static_cast<int>(uniq.size()), nranks) << "overlap";
            // contiguity & ordering
            int expect0 = 0;
            for (auto& r : ranges) {
                EXPECT_EQ(r.first, expect0);
                expect0 = r.second;
            }
            EXPECT_EQ(expect0, ny);
        }
    }
}

TEST(Decomposition, OwnershipFlags) {
    auto d0 = SlabDecomposition::make(80, 4, 0);
    EXPECT_TRUE(d0.owns_jmin());
    EXPECT_FALSE(d0.owns_jmax());
    EXPECT_FALSE(d0.owns_lower_neighbor());
    EXPECT_TRUE(d0.owns_upper_neighbor());

    auto d3 = SlabDecomposition::make(80, 4, 3);
    EXPECT_FALSE(d3.owns_jmin());
    EXPECT_TRUE(d3.owns_jmax());
    EXPECT_TRUE(d3.owns_lower_neighbor());
    EXPECT_FALSE(d3.owns_upper_neighbor());

    auto solo = SlabDecomposition::make(10, 1, 0);
    EXPECT_TRUE(solo.owns_jmin() && solo.owns_jmax());
}

TEST(Decomposition, RejectsInvalid) {
    EXPECT_THROW(SlabDecomposition::make(4, 8, 0), std::runtime_error);
    EXPECT_THROW(SlabDecomposition::make(10, 2, 5), std::runtime_error);
}

}  // namespace
