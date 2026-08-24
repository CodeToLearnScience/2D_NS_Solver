#pragma once

// Ghost-aware contiguous fields for the refactored solver.
//
// Index convention: interior cells span [0, ni) x [0, nj); ghost cells extend
// ng layers on every side, so operator() accepts i in [-ng, ni+ng). Storage is
// a single contiguous buffer, row-major with j fastest -- identical memory
// order to the legacy x[i][j] pointer-to-pointer layout, so parity porting and
// cache behaviour both line up.
//
// fill_ghosts_copy() reproduces the legacy ghost-filling convention exactly:
// every ghost layer copies the nearest interior layer (not staggered), and
// corner ghosts end up holding the nearest interior corner value.

#include <cassert>
#include <format>
#include <stdexcept>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace ns {

/// Minimal 2-D strided view with std::mdspan-compatible ergonomics.
/// Swap-in point: libstdc++ still lacks <mdspan> (checked on GCC 15); when it
/// lands, replace this alias set with the real thing -- the [i, j] subscript
/// syntax and extent() API below are deliberately identical.
template <typename T>
class View2D {
public:
    using element_type = T;

    View2D() noexcept = default;
    View2D(T* data, int ni, int nj) noexcept : data_(data), ni_(ni), nj_(nj) {}

    template <typename U>
        requires(std::is_const_v<T> && !std::is_const_v<U>)
    View2D(const View2D<U>& other) noexcept
        : data_(other.data()), ni_(other.extent(0)), nj_(other.extent(1)) {}

    [[nodiscard]] constexpr int extent(int rank) const noexcept {
        assert(rank == 0 || rank == 1);
        return rank == 0 ? ni_ : nj_;
    }
    [[nodiscard]] constexpr T* data() const noexcept { return data_; }

    T& operator[](int i, int j) const noexcept {
        assert(i >= 0 && i < ni_);
        assert(j >= 0 && j < nj_);
        return data_[static_cast<std::size_t>(i) * static_cast<std::size_t>(nj_) +
                     static_cast<std::size_t>(j)];
    }

private:
    T* data_ = nullptr;
    int ni_ = 0;
    int nj_ = 0;
};

template <typename T>
class Field {
public:
    Field() = default;
    Field(int ni, int nj, int ng) { allocate(ni, nj, ng); }

    void allocate(int ni, int nj, int ng) {
        ni_ = ni;
        nj_ = nj;
        ng_ = ng;
        data_.assign(canvas_size(), T{});
    }

    [[nodiscard]] constexpr int ni() const noexcept { return ni_; }
    [[nodiscard]] constexpr int nj() const noexcept { return nj_; }
    [[nodiscard]] constexpr int ng() const noexcept { return ng_; }
    [[nodiscard]] constexpr int total_ni() const noexcept { return ni_ + 2 * ng_; }
    [[nodiscard]] constexpr int total_nj() const noexcept { return nj_ + 2 * ng_; }
    /// Stride between consecutive i indices (elements), i.e. the row length.
    [[nodiscard]] constexpr int stride_i() const noexcept { return total_nj(); }
    [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }

    T& operator()(int i, int j) noexcept {
        assert(i >= -ng_ && i < ni_ + ng_);
        assert(j >= -ng_ && j < nj_ + ng_);
        return data_[static_cast<std::size_t>((i + ng_) * stride_i() + (j + ng_))];
    }
    const T& operator()(int i, int j) const noexcept {
        assert(i >= -ng_ && i < ni_ + ng_);
        assert(j >= -ng_ && j < nj_ + ng_);
        return data_[static_cast<std::size_t>((i + ng_) * stride_i() + (j + ng_))];
    }

    /// Bounds-checked access (always active, unlike operator()'s asserts);
    /// throws std::out_of_range on violation.
    T& at(int i, int j) {
        if (i < -ng_ || i >= ni_ + ng_ || j < -ng_ || j >= nj_ + ng_)
            throw std::out_of_range(std::format("Field index ({},{}) outside [{},{})x[{},{})",
                                                i, j, -ng_, ni_ + ng_, -ng_, nj_ + ng_));
        return (*this)(i, j);
    }
    const T& at(int i, int j) const {
        if (i < -ng_ || i >= ni_ + ng_ || j < -ng_ || j >= nj_ + ng_)
            throw std::out_of_range(std::format("Field index ({},{}) outside [{},{})x[{},{})",
                                                i, j, -ng_, ni_ + ng_, -ng_, nj_ + ng_));
        return (*this)(i, j);
    }

    [[nodiscard]] T* data() noexcept { return data_.data(); }
    [[nodiscard]] const T* data() const noexcept { return data_.data(); }

    /// View over the whole canvas; kernels index [0, total_ni) x [0, total_nj).
    [[nodiscard]] View2D<T> view() noexcept {
        return {data(), total_ni(), total_nj()};
    }
    [[nodiscard]] View2D<const T> view() const noexcept {
        return {data(), total_ni(), total_nj()};
    }

    void fill(const T& value) { data_.assign(data_.size(), value); }

    /// Legacy ghost convention: each ghost layer copies the nearest interior
    /// layer; corners end up with the nearest interior corner value.
    void fill_ghosts_copy() {
        for (int g = 1; g <= ng_; ++g) {
            for (int j = -ng_; j < nj_ + ng_; ++j) {
                (*this)(-g, j)      = (*this)(0, j);
                (*this)(ni_ - 1 + g, j) = (*this)(ni_ - 1, j);
            }
        }
        for (int g = 1; g <= ng_; ++g) {
            for (int i = -ng_; i < ni_ + ng_; ++i) {
                (*this)(i, -g)          = (*this)(i, 0);
                (*this)(i, nj_ - 1 + g) = (*this)(i, nj_ - 1);
            }
        }
    }

private:
    [[nodiscard]] std::size_t canvas_size() const noexcept {
        return static_cast<std::size_t>(total_ni()) * static_cast<std::size_t>(total_nj());
    }

    int ni_ = 0, nj_ = 0, ng_ = 0;
    std::vector<T> data_;
};

/// nv variable planes stacked contiguously; each plane has Field layout.
template <typename T>
class MultiField {
public:
    MultiField() = default;
    MultiField(int nv, int ni, int nj, int ng) { allocate(nv, ni, nj, ng); }

    void allocate(int nv, int ni, int nj, int ng) {
        nv_ = nv;
        ni_ = ni;
        nj_ = nj;
        ng_ = ng;
        data_.assign(static_cast<std::size_t>(nv_) * plane_size(), T{});
    }

    [[nodiscard]] constexpr int nv() const noexcept { return nv_; }
    [[nodiscard]] constexpr int ni() const noexcept { return ni_; }
    [[nodiscard]] constexpr int nj() const noexcept { return nj_; }
    [[nodiscard]] constexpr int ng() const noexcept { return ng_; }
    [[nodiscard]] constexpr int stride_i() const noexcept { return nj_ + 2 * ng_; }
    [[nodiscard]] constexpr std::size_t plane_size() const noexcept {
        return static_cast<std::size_t>(ni_ + 2 * ng_) *
               static_cast<std::size_t>(nj_ + 2 * ng_);
    }
    [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }

    T& operator()(int v, int i, int j) noexcept {
        assert(v >= 0 && v < nv_);
        assert(i >= -ng_ && i < ni_ + ng_);
        assert(j >= -ng_ && j < nj_ + ng_);
        return data_[static_cast<std::size_t>(v) * plane_size() +
                     static_cast<std::size_t>((i + ng_) * stride_i() + (j + ng_))];
    }
    const T& operator()(int v, int i, int j) const noexcept {
        assert(v >= 0 && v < nv_);
        assert(i >= -ng_ && i < ni_ + ng_);
        assert(j >= -ng_ && j < nj_ + ng_);
        return data_[static_cast<std::size_t>(v) * plane_size() +
                     static_cast<std::size_t>((i + ng_) * stride_i() + (j + ng_))];
    }

    [[nodiscard]] T* plane_data(int v) noexcept {
        assert(v >= 0 && v < nv_);
        return data_.data() + static_cast<std::size_t>(v) * plane_size();
    }
    [[nodiscard]] const T* plane_data(int v) const noexcept {
        assert(v >= 0 && v < nv_);
        return data_.data() + static_cast<std::size_t>(v) * plane_size();
    }

    [[nodiscard]] View2D<T> view(int v) noexcept {
        return {plane_data(v), ni_ + 2 * ng_, nj_ + 2 * ng_};
    }
    [[nodiscard]] View2D<const T> view(int v) const noexcept {
        return {plane_data(v), ni_ + 2 * ng_, nj_ + 2 * ng_};
    }

    void fill(const T& value) { data_.assign(data_.size(), value); }

    void fill_ghosts_copy() {
        for (int v = 0; v < nv_; ++v) {
            for (int g = 1; g <= ng_; ++g) {
                for (int j = -ng_; j < nj_ + ng_; ++j) {
                    (*this)(v, -g, j)          = (*this)(v, 0, j);
                    (*this)(v, ni_ - 1 + g, j) = (*this)(v, ni_ - 1, j);
                }
            }
            for (int g = 1; g <= ng_; ++g) {
                for (int i = -ng_; i < ni_ + ng_; ++i) {
                    (*this)(v, i, -g)          = (*this)(v, i, 0);
                    (*this)(v, i, nj_ - 1 + g) = (*this)(v, i, nj_ - 1);
                }
            }
        }
    }

private:
    int nv_ = 0, ni_ = 0, nj_ = 0, ng_ = 0;
    std::vector<T> data_;
};

}  // namespace ns
