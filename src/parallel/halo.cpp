#include "parallel/halo.hpp"

#ifdef NS_WITH_MPI
#include <mpi.h>
#endif

#include <cstring>

namespace ns::parallel {
namespace {

// Pack ng real rows starting at row `j_first` (ascending) into dst[g*nx+i].
void pack_rows(const double* plane, int stride_i, int ng_off, int nx, int j_first,
               int ng, double* dst) {
    for (int g = 0; g < ng; ++g)
        for (int i = 0; i < nx; ++i)
            dst[static_cast<std::size_t>(g) * nx + i] =
                plane[static_cast<std::size_t>(i + ng_off) * stride_i +
                      static_cast<std::size_t>(j_first + g + ng_off)];
}

void unpack_upper(const double* src, double* plane, int stride_i, int ng_off, int nx,
                  int ny, int ng) {
    // upper ghost ring k (row ny-1+k) <- neighbour real row k-1 (ascending)
    for (int g = 0; g < ng; ++g)
        for (int i = 0; i < nx; ++i)
            plane[static_cast<std::size_t>(i + ng_off) * stride_i +
                  static_cast<std::size_t>(ny + g + ng_off)] =
                src[static_cast<std::size_t>(g) * nx + i];
}

void unpack_lower(const double* src, double* plane, int stride_i, int ng_off, int nx,
                  int ng) {
    // lower ghost ring k (row -k) <- neighbour's last real rows, delivered
    // ascending: ring 1 takes slot ng-1, ring 2 slot ng-2, ...
    for (int k = 1; k <= ng; ++k)
        for (int i = 0; i < nx; ++i)
            plane[static_cast<std::size_t>(i + ng_off) * stride_i +
                  static_cast<std::size_t>(ng_off - k)] =
                src[static_cast<std::size_t>(ng - k) * nx + i];
}

}  // namespace

HaloExchange::HaloExchange(int nx, int ny_local, int ng, int lower_rank,
                           int upper_rank, int tag_base)
    : nx_(nx), ny_(ny_local), ng_(ng), lower_(lower_rank), upper_(upper_rank),
      tag_base_(tag_base) {
    const std::size_t buf = static_cast<std::size_t>(ng_) * nx_;
    send_lo_.assign(buf, 0.0);
    recv_lo_.assign(buf, 0.0);
    send_hi_.assign(buf, 0.0);
    recv_hi_.assign(buf, 0.0);
}

void HaloExchange::exchange(MultiField<double>& f) const {
    for (int v = 0; v < f.nv(); ++v)
        exchange_plane(f.plane_data(v), f.stride_i(), f.ng(), v);
}

void HaloExchange::exchange(Field<double>& f) const {
    exchange_plane(f.data(), f.stride_i(), f.ng(), 0);
}

void HaloExchange::exchange_plane(double* plane, int stride_i, int ng_off,
                                  int plane_id) const {
#ifndef NS_WITH_MPI
    (void)plane;
    (void)stride_i;
    (void)ng_off;
    (void)plane_id;
#else
    pack_rows(plane, stride_i, ng_off, nx_, 0, ng_, send_lo_.data());          // bottom rows
    pack_rows(plane, stride_i, ng_off, nx_, ny_ - ng_, ng_, send_hi_.data());  // top rows

    const int t1 = tag_base_ + 2 * plane_id;      // tops travel upward
    const int t2 = tag_base_ + 2 * plane_id + 1;  // bottoms travel downward

    if (lower_ >= 0 && upper_ >= 0) {
        // Op1: my top rows -> upper; lower's top rows -> me
        MPI_Sendrecv(send_hi_.data(), static_cast<int>(send_hi_.size()), MPI_DOUBLE,
                     upper_, t1, recv_lo_.data(), static_cast<int>(recv_lo_.size()),
                     MPI_DOUBLE, lower_, t1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Op2: my bottom rows -> lower; upper's bottom rows -> me
        MPI_Sendrecv(send_lo_.data(), static_cast<int>(send_lo_.size()), MPI_DOUBLE,
                     lower_, t2, recv_hi_.data(), static_cast<int>(recv_hi_.size()),
                     MPI_DOUBLE, upper_, t2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (upper_ >= 0) {
        // first rank: send top rows up, receive upper's bottom rows down
        MPI_Sendrecv(send_hi_.data(), static_cast<int>(send_hi_.size()), MPI_DOUBLE,
                     upper_, t1, recv_hi_.data(), static_cast<int>(recv_hi_.size()),
                     MPI_DOUBLE, upper_, t2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (lower_ >= 0) {
        // last rank: receive lower's top rows, send bottom rows down
        MPI_Sendrecv(send_lo_.data(), static_cast<int>(send_lo_.size()), MPI_DOUBLE,
                     lower_, t2, recv_lo_.data(), static_cast<int>(recv_lo_.size()),
                     MPI_DOUBLE, lower_, t1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (upper_ >= 0) unpack_upper(recv_hi_.data(), plane, stride_i, ng_off, nx_, ny_, ng_);
    if (lower_ >= 0) unpack_lower(recv_lo_.data(), plane, stride_i, ng_off, nx_, ng_);
#endif
}

}  // namespace ns::parallel
