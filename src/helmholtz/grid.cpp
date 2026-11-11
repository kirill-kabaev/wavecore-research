#include "wavecore/helmholtz/grid.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace wavecore::helmholtz {

Grid2D::Grid2D(size_type nx, size_type ny, double lx, double ly)
    : nx_{nx},
      ny_{ny},
      unknown_count_{0},
      lx_{lx},
      ly_{ly},
      hx_{0.0},
      hy_{0.0} {
    if (nx_ == 0 || ny_ == 0) {
        throw std::invalid_argument("Grid dimensions nx and ny must be positive");
    }

    if (!std::isfinite(lx_) || lx_ <= 0.0 || !std::isfinite(ly_) || ly_ <= 0.0) {
        throw std::invalid_argument("Domain lengths lx and ly must be finite and positive");
    }

    constexpr auto max_size = std::numeric_limits<size_type>::max();
    if (nx_ == max_size || ny_ == max_size) {
        throw std::overflow_error("Grid dimension is too large to include boundary intervals");
    }

    if (nx_ > max_size / ny_) {
        throw std::overflow_error("Grid unknown count nx * ny overflows size_type");
    }
    unknown_count_ = nx_ * ny_;

    hx_ = lx_ / (static_cast<double>(nx_) + 1.0);
    hy_ = ly_ / (static_cast<double>(ny_) + 1.0);
    if (!std::isfinite(hx_) || hx_ <= 0.0 || !std::isfinite(hy_) || hy_ <= 0.0) {
        throw std::invalid_argument("Computed grid spacing must be finite and positive");
    }
}

Grid2D::size_type Grid2D::nx() const noexcept {
    return nx_;
}

Grid2D::size_type Grid2D::ny() const noexcept {
    return ny_;
}

Grid2D::size_type Grid2D::unknown_count() const noexcept {
    return unknown_count_;
}

double Grid2D::lx() const noexcept {
    return lx_;
}

double Grid2D::ly() const noexcept {
    return ly_;
}

double Grid2D::hx() const noexcept {
    return hx_;
}

double Grid2D::hy() const noexcept {
    return hy_;
}

double Grid2D::x(size_type i) const {
    if (i >= nx_) {
        throw std::out_of_range("Grid x index is out of range");
    }
    return (static_cast<double>(i) + 1.0) * hx_;
}

double Grid2D::y(size_type j) const {
    if (j >= ny_) {
        throw std::out_of_range("Grid y index is out of range");
    }
    return (static_cast<double>(j) + 1.0) * hy_;
}

Grid2D::size_type Grid2D::linear_index(size_type i, size_type j) const {
    if (i >= nx_ || j >= ny_) {
        throw std::out_of_range("Grid indices are out of range");
    }
    return j * nx_ + i;
}

Grid2D::index_pair Grid2D::grid_indices(size_type p) const {
    if (p >= unknown_count_) {
        throw std::out_of_range("Linear grid index is out of range");
    }
    return {p % nx_, p / nx_};
}

}  // namespace wavecore::helmholtz
