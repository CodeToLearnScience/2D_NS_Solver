#pragma once

// 1-D slab domain decomposition along j (Phase 7).
//
// Global cell rows [0, ny_global) are split into nranks contiguous,
// non-overlapping ranges. Rank r owns [j0(r), j1(r)). i-direction is NOT
// decomposed, so imin/imax edges are fully local per rank.

#include <utility>
#include <vector>

namespace ns::parallel {

struct SlabDecomposition {
    int ny_global = 0;
    int nranks = 1;
    int rank = 0;

    int j0 = 0;         ///< first global cell row owned by this rank
    int j1 = 0;         ///< one past the last global cell row owned
    [[nodiscard]] int ny_local() const noexcept { return j1 - j0; }

    /// True if this rank owns a piece of the given edge (jmin/jmax) or,
    /// for i-edges, always true (i is not decomposed).
    [[nodiscard]] bool owns_jmin() const noexcept { return rank == 0; }
    [[nodiscard]] bool owns_jmax() const noexcept { return rank == nranks - 1; }
    [[nodiscard]] bool owns_lower_neighbor() const noexcept { return rank > 0; }
    [[nodiscard]] bool owns_upper_neighbor() const noexcept { return rank < nranks - 1; }

    [[nodiscard]] static SlabDecomposition make(int ny_global, int nranks, int rank);
    [[nodiscard]] static std::vector<std::pair<int, int>> all_ranges(int ny_global, int nranks);
};

}  // namespace ns::parallel
