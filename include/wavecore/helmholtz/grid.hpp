#pragma once

#include <cstddef>
#include <utility>

namespace wavecore::helmholtz {

// Uniform rectangular grid containing only interior unknowns.
//
// The public indices i and j are zero-based and correspond to the
// mathematical nodes x_{i+1} and y_{j+1}. Linear storage is x-fastest:
// p = j * nx + i. Therefore i varies in contiguous memory locations.
class Grid2D {
public:
    using size_type = std::size_t;
    using index_pair = std::pair<size_type, size_type>;

    Grid2D(size_type nx, size_type ny, double lx, double ly);

    [[nodiscard]] size_type nx() const noexcept;
    [[nodiscard]] size_type ny() const noexcept;
    [[nodiscard]] size_type unknown_count() const noexcept;

    [[nodiscard]] double lx() const noexcept;
    [[nodiscard]] double ly() const noexcept;
    [[nodiscard]] double hx() const noexcept;
    [[nodiscard]] double hy() const noexcept;

    [[nodiscard]] double x(size_type i) const;
    [[nodiscard]] double y(size_type j) const;

    [[nodiscard]] size_type linear_index(size_type i, size_type j) const;
    [[nodiscard]] index_pair grid_indices(size_type p) const;

private:
    size_type nx_;
    size_type ny_;
    size_type unknown_count_;
    double lx_;
    double ly_;
    double hx_;
    double hy_;
};

}  // namespace wavecore::helmholtz
