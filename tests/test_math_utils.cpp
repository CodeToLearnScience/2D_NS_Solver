// Unit tests for the header-only math utilities in inc/basic_functions.h.
// The MinModLim tie-case test is a regression guard for the missing-return
// (undefined behavior) bug fixed during the Phase 0 audit.

#include <gtest/gtest.h>

#include <string>

#include "basic_functions.h"

namespace {

constexpr double kDx = 1.0;

TEST(MathUtils, MaxMinComparisons) {
    EXPECT_DOUBLE_EQ(Max2(1.0, 2.0), 2.0);
    EXPECT_DOUBLE_EQ(Max2(2.0, 1.0), 2.0);
    EXPECT_DOUBLE_EQ(Min2(1.0, -2.0), -2.0);
    EXPECT_DOUBLE_EQ(Max3(1.0, 5.0, 3.0), 5.0);
    EXPECT_DOUBLE_EQ(Max3(5.0, 3.0, 1.0), 5.0);
    EXPECT_DOUBLE_EQ(Min3(1.0, -5.0, 3.0), -5.0);
    EXPECT_DOUBLE_EQ(Max4(-4.0, 7.0, 7.0, 2.0), 7.0);
    EXPECT_DOUBLE_EQ(Min4(9.0, -9.0, 0.0, 3.0), -9.0);
}

TEST(MathUtils, SignIsThreeWay) {
    EXPECT_DOUBLE_EQ(Sign(0.3), 1.0);
    EXPECT_DOUBLE_EQ(Sign(-2.0e-9), -1.0);
    EXPECT_DOUBLE_EQ(Sign(0.0), 0.0);
}

TEST(MathUtils, ToStringWithPrecision) {
    EXPECT_EQ(ToStringWithPrecision(3.14159, 4), "3.142");
    EXPECT_EQ(ToStringWithPrecision(1.0, 16).substr(0, 1), "1");
}

TEST(MathUtils, GetStringBetweenDelimiters) {
    EXPECT_EQ(get_str_between_two_str(std::string("Iter_100_out.dat"), std::string("_"),
                                      std::string("_")),
              "100");
}

TEST(MinModLim, OppositeSlopesYieldZero) {
    // Ul > 0, Ur < 0 -> local extremum -> slope must be zeroed.
    EXPECT_DOUBLE_EQ(MinModLim(0.0, 1.0, 0.5, kDx), 0.0);
    // Ul == 0 -> product is zero -> zero.
    EXPECT_DOUBLE_EQ(MinModLim(1.0, 1.0, 2.0, kDx), 0.0);
}

TEST(MinModLim, SmallerMagnitudeSlopeWins) {
    // Ul = 1, Ur = 3 (same sign): minmod picks Ul = 1.
    EXPECT_DOUBLE_EQ(MinModLim(0.0, 1.0, 4.0, kDx), 1.0);
    // Ul = 3, Ur = 1: minmod picks Ur = 1.
    EXPECT_DOUBLE_EQ(MinModLim(-2.0, 1.0, 2.0, kDx), 1.0);
}

TEST(MinModLim, EqualSlopesReturnFiniteValue) {
    // Regression: |Ul| == |Ur| with same sign previously fell off the end of
    // the function (UB). Must now return a well-defined value equal to the slope.
    const double v = MinModLim(0.0, 1.0, 2.0, kDx);  // Ul = Ur = 1
    ASSERT_TRUE(v == v);                             // not NaN
    EXPECT_DOUBLE_EQ(v, 1.0);
}

TEST(ReconValue, VanishesForLinearProfile) {
    // Second-order reconstruction difference must be exact (zero up-slope
    // correction) for a linear profile u(i) = i*dx.
    const double u0 = 0.0, u1 = 1.0, u2 = 2.0, u3 = 3.0;
    EXPECT_NEAR(ReconValue(u0, u1, u2, u3, kDx), 0.0, 1e-14);
}

TEST(LogMean, ReturnsEndpointForEqualArguments) {
    EXPECT_DOUBLE_EQ(Log_Mean(2.0, 2.0), 2.0);
}

TEST(LogMean, NearEqualStatesMatchUntouchedApproximation) {
    // Small-jump branch kept from legacy code: 1.1 * F(u), u = f^2 with
    // xi = 1.2. Guards accidental changes to the tuned approximation.
    const double lm = Log_Mean(1.2, 1.0);
    EXPECT_NEAR(lm, 1.10305, 1e-4);
    // ... while still tracking the true logarithmic mean within ~1%.
    EXPECT_NEAR(lm, 0.2 / std::log(1.2), 0.01);
}

TEST(LogMean, LargeJumpIsExactAndPositive) {
    // Regression: for |xi - 1| > ~10% the legacy else-branch inverted f/ln(xi)
    // and returned NEGATIVE values (e.g. -1.56 instead of +1.4427 for (2,1)).
    const double true_lm = 1.0 / std::log(2.0);  // L(2,1) = 1/ln 2
    EXPECT_NEAR(Log_Mean(2.0, 1.0), true_lm, 1e-12);
    EXPECT_NEAR(Log_Mean(1.0, 2.0), true_lm, 1e-12);   // symmetric
    EXPECT_DOUBLE_EQ(Log_Mean(4.0, 2.0), 2.0 * Log_Mean(2.0, 1.0));  // homogeneous
}

}  // namespace
