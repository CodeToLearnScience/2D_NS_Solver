// MPI halo-exchange test: runs under mpirun with 2 or 3 ranks.

#include <gtest/gtest.h>

#include <mpi.h>

#include "parallel/decomposition.hpp"
#include "parallel/halo.hpp"

namespace {

using namespace ns;

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    int rc = RUN_ALL_TESTS();
    MPI_Finalize();
    return rc;
}

namespace {

TEST(MpiHalo, GhostsMatchNeighbourRealRows) {
    int rank = 0, nranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    constexpr int NX = 6, NYG = 12, NG = 2;
    auto d = parallel::SlabDecomposition::make(NYG, nranks, rank);

    MultiField<double> f(3, NX, d.ny_local(), NG);
    for (int v = 0; v < 3; ++v)
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < d.ny_local(); ++j)
                f(v, i, j) = v * 100000 + (d.j0 + j) * 100 + i;

    parallel::HaloExchange hx(NX, d.ny_local(), NG,
                              d.owns_lower_neighbor() ? rank - 1 : -1,
                              d.owns_upper_neighbor() ? rank + 1 : -1);
    hx.exchange(f);

    for (int k = 1; k <= NG && d.j0 - k >= 0; ++k)
        for (int v = 0; v < 3; ++v)
            for (int i = 0; i < NX; ++i) {
                const double want = v * 100000 + (d.j0 - k) * 100 + i;
                ASSERT_NEAR(f(v, i, -k), want, 1e-12);
            }
    for (int k = 1; k <= NG && d.j1 - 1 + k < NYG; ++k)
        for (int v = 0; v < 3; ++v)
            for (int i = 0; i < NX; ++i) {
                const double want = v * 100000 + (d.j1 - 1 + k) * 100 + i;
                ASSERT_NEAR(f(v, i, d.ny_local() - 1 + k), want, 1e-12);
            }

    MPI_Barrier(MPI_COMM_WORLD);
}

}  // namespace
