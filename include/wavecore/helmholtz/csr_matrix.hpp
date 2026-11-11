#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace wavecore::helmholtz {

// Immutable compressed sparse row matrix.
//
// Column indices must be strictly increasing within every row. This canonical
// representation excludes duplicate entries and makes the stored ordering
// deterministic.
class CsrMatrix {
public:
    using size_type = std::size_t;
    using value_type = std::complex<double>;
    using index_array = std::vector<size_type>;
    using value_array = std::vector<value_type>;

    CsrMatrix(size_type rows,
              size_type cols,
              index_array row_ptr,
              index_array col_idx,
              value_array values);

    [[nodiscard]] size_type rows() const noexcept;
    [[nodiscard]] size_type cols() const noexcept;
    [[nodiscard]] size_type nnz() const noexcept;

    [[nodiscard]] const index_array& row_ptr() const noexcept;
    [[nodiscard]] const index_array& col_idx() const noexcept;
    [[nodiscard]] const value_array& values() const noexcept;

    [[nodiscard]] value_array multiply(const value_array& vector) const;

private:
    void validate() const;

    size_type rows_;
    size_type cols_;
    index_array row_ptr_;
    index_array col_idx_;
    value_array values_;
};

}  // namespace wavecore::helmholtz
