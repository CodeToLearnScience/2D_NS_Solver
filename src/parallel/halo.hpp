#pragma once

// Halo exchange for j-slab decomposition: packs contiguous buffers, one
// matched MPI_Sendrecv per neighbour per field. Ghost ring k is filled from
// the neighbour's REAL rows (exact copies), which is what interior stencil
// faces on either side of a rank interface consume.
//
// Without NS_WITH_MPI every exchange compiles to a no-op, so the same solver
// code runs serially.

#include "fields/field.hpp"

#ifdef NS_WITH_MPI
#include <mpi.h>
#endif

namespace ns::parallel {

class HaloExchange {
public:
    HaloExchange() = default;
    HaloExchange(int nx, int ny_local, int ng, int lower_rank, int upper_rank,
                 int tag_base = 700);

    void exchange(MultiField<double>& f) const;
    void exchange(Field<double>& f) const;

private:
    void exchange_plane(double* plane, int stride_i, int ng_off, int plane_id) const;

    int nx_ = 0, ny_ = 0, ng_ = 0;
    int lower_ = -1, upper_ = -1, tag_base_ = 700;
    mutable std::vector<double> send_lo_, recv_lo_, send_hi_, recv_hi_;
};

}  // namespace ns::parallel
