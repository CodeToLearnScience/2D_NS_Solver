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

    /// True when this rank owns the global domain bottom (no lower neighbor).
    [[nodiscard]] bool at_global_bottom() const noexcept { return lower_ < 0; }
    /// True when this rank owns the global domain top (no upper neighbor).
    [[nodiscard]] bool at_global_top() const noexcept { return upper_ < 0; }

    /// Exchanges j-face metric rows (nf = ny+1 real face rows) so that ghost
    /// slots -1 and nf hold the TRUE neighbor face metrics instead of copies.
    /// xplane/yplane are strided arrays of nf + 2*ng doubles each.
    void exchange_face_metric_rows(double* xplane, double* yplane, int stride_i,
                                   int nf) const;

private:
    void exchange_plane(double* plane, int stride_i, int ng_off, int plane_id) const;

    int nx_ = 0, ny_ = 0, ng_ = 0;
    int lower_ = -1, upper_ = -1, tag_base_ = 700;
    mutable std::vector<double> send_lo_, recv_lo_, send_hi_, recv_hi_;
};

}  // namespace ns::parallel
