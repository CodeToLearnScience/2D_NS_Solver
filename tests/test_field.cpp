// Unit tests for Field / MultiField: ghost-aware indexing, memory layout,
// mdspan views, and the legacy ghost-fill convention.

#include <gtest/gtest.h>

#include <type_traits>

#include "fields/field.hpp"

namespace {

TEST(FieldLayout, GhostAwareIndexing) {
    ns::Field<double> f(4, 3, 2);
    EXPECT_EQ(f.ni(), 4);
    EXPECT_EQ(f.nj(), 3);
    EXPECT_EQ(f.ng(), 2);
    EXPECT_EQ(f.total_ni(), 8);
    EXPECT_EQ(f.total_nj(), 7);

    // interior + all ghost rings are addressable
    f(-2, -2) = 1.0;
    f(0, 0) = 2.0;
    f(3, 2) = 3.0;
    f(5, 4) = 4.0;  // ni+ng-1 = 5, nj+ng-1 = 4
    EXPECT_DOUBLE_EQ(f(-2, -2), 1.0);
    EXPECT_DOUBLE_EQ(f(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(f(3, 2), 3.0);
    EXPECT_DOUBLE_EQ(f(5, 4), 4.0);

    // operator() bounds are assert-guarded (active in debug builds); at()
    // must throw regardless of build type.
    EXPECT_THROW(f.at(-3, 0), std::out_of_range);
    EXPECT_THROW(f.at(6, 0), std::out_of_range);
    EXPECT_THROW(f.at(0, 5), std::out_of_range);  // j >= nj + ng
}

TEST(FieldLayout, ContiguousJFastestMatchingLegacy) {
    // (i,j) and (i,j+1) must be adjacent in memory; (i,j) to (i+1,j) strides
    // one full row -- identical to legacy double** layout x[i][j].
    ns::Field<double> f(5, 6, 2);
    const double* base = &f(1, 1);
    EXPECT_EQ(base + 1, &f(1, 2));
    EXPECT_EQ(base + f.stride_i(), &f(2, 1));
}

TEST(FieldLayout, FillAndGhostsCopy) {
    ns::Field<double> f(3, 3, 2);
    f.fill(0.0);
    // Distinct value per interior cell.
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) f(i, j) = 10 * i + j;

    f.fill_ghosts_copy();

    // Every ghost layer equals the nearest interior layer.
    for (int g = 1; g <= 2; ++g) {
        for (int j = -2; j < 5; ++j) {
            EXPECT_DOUBLE_EQ(f(-g, j), f(0, j)) << "i-ghost low";
            EXPECT_DOUBLE_EQ(f(2 + g, j), f(2, j)) << "i-ghost high";
        }
        for (int i = -2; i < 5; ++i) {
            EXPECT_DOUBLE_EQ(f(i, -g), f(i, 0)) << "j-ghost low";
            EXPECT_DOUBLE_EQ(f(i, 3 - 1 + g), f(i, 2)) << "j-ghost high";
        }
    }
    // Corners hold the nearest interior corner value (legacy convention).
    EXPECT_DOUBLE_EQ(f(-2, -2), f(0, 0));
    EXPECT_DOUBLE_EQ(f(4, 4), f(2, 2));
}

TEST(FieldLayout, MdspanViewCoversWholeCanvas) {
    ns::Field<double> f(2, 3, 1);
    f.fill(1.5);
    auto v = f.view();
    static_assert(std::is_same_v<decltype(v)::element_type, double>);
    EXPECT_EQ(v.extent(0), 4);
    EXPECT_EQ(v.extent(1), 5);
    // view index [0,0] corresponds to logical (-ng,-ng).
    v[0, 0] = 9.0;
    EXPECT_DOUBLE_EQ(f(-1, -1), 9.0);
    double sum = 0;
    for (int i = 0; i < v.extent(0); ++i)
        for (int j = 0; j < v.extent(1); ++j) sum += v[i, j];
    EXPECT_NEAR(sum, 20 * 1.5 - 1.5 + 9.0, 1e-12);
}

TEST(MultiFieldLayout, PlanesAreContiguousAndIndependent) {
    ns::MultiField<double> cv(4, 3, 4, 2);
    EXPECT_EQ(cv.nv(), 4);
    // Canvas origin (-ng,-ng) coincides with plane_data(v)[0].
    cv(0, -2, -2) = 1.0;
    cv(3, -2, -2) = 4.0;
    EXPECT_DOUBLE_EQ(cv.plane_data(0)[0], 1.0);
    EXPECT_DOUBLE_EQ(cv.plane_data(3)[0], 4.0);
    EXPECT_EQ(cv.plane_data(0) + cv.plane_size(), cv.plane_data(1));

    cv.fill(7.0);
    cv.fill_ghosts_copy();  // no-op after fill but must not crash
    EXPECT_DOUBLE_EQ(cv(2, -2, -2), 7.0);

    auto v1 = cv.view(1);
    v1[0, 0] = 3.25;
    EXPECT_DOUBLE_EQ(cv(1, -2, -2), 3.25);
}

}  // namespace
