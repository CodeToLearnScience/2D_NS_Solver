#include "parallel/decomposition.hpp"

#include <algorithm>
#include <stdexcept>

namespace ns::parallel {

SlabDecomposition SlabDecomposition::make(int ny_global, int nranks, int rank) {
    if (ny_global < 1 || nranks < 1 || rank < 0 || rank >= nranks)
        throw std::runtime_error("invalid decomposition arguments");
    if (nranks > ny_global)
        throw std::runtime_error("more ranks than cell rows is unsupported");

    SlabDecomposition d;
    d.ny_global = ny_global;
    d.nranks = nranks;
    d.rank = rank;

    const int base = ny_global / nranks;
    const int rem = ny_global % nranks;
    // first `rem` ranks get one extra row
    const int before = std::min(rank, rem);
    d.j0 = rank * base + before;
    d.j1 = d.j0 + base + (rank < rem ? 1 : 0);
    return d;
}

std::vector<std::pair<int, int>> SlabDecomposition::all_ranges(int ny_global, int nranks) {
    std::vector<std::pair<int, int>> r;
    r.reserve(nranks);
    for (int rk = 0; rk < nranks; ++rk) {
        auto d = make(ny_global, nranks, rk);
        r.emplace_back(d.j0, d.j1);
    }
    return r;
}

}  // namespace ns::parallel
