#include "wavecore/helmholtz/csr_matrix.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace wavecore::helmholtz {

CsrMatrix::CsrMatrix(size_type rows,
                     size_type cols,
                     index_array row_ptr,
                     index_array col_idx,
                     value_array values)
    : rows_{rows},
      cols_{cols},
      row_ptr_{std::move(row_ptr)},
      col_idx_{std::move(col_idx)},
      values_{std::move(values)} {
    validate();
}

CsrMatrix::size_type CsrMatrix::rows() const noexcept {
    return rows_;
}

CsrMatrix::size_type CsrMatrix::cols() const noexcept {
    return cols_;
}

CsrMatrix::size_type CsrMatrix::nnz() const noexcept {
    return values_.size();
}

const CsrMatrix::index_array& CsrMatrix::row_ptr() const noexcept {
    return row_ptr_;
}

const CsrMatrix::index_array& CsrMatrix::col_idx() const noexcept {
    return col_idx_;
}

const CsrMatrix::value_array& CsrMatrix::values() const noexcept {
    return values_;
}

CsrMatrix::value_array CsrMatrix::multiply(const value_array& vector) const {
    if (vector.size() != cols_) {
        throw std::invalid_argument("Input vector size must equal the matrix column count");
    }

    value_array result(rows_, value_type{0.0, 0.0});
    for (size_type row = 0; row < rows_; ++row) {
        for (size_type entry = row_ptr_[row]; entry < row_ptr_[row + 1]; ++entry) {
            result[row] += values_[entry] * vector[col_idx_[entry]];
        }
    }
    return result;
}

void CsrMatrix::validate() const {
    constexpr auto max_size = std::numeric_limits<size_type>::max();
    if (rows_ == max_size) {
        throw std::overflow_error("Matrix row count overflows the required row_ptr size");
    }

    if (row_ptr_.size() != rows_ + 1) {
        throw std::invalid_argument("row_ptr size must equal rows + 1");
    }

    if (col_idx_.size() != values_.size()) {
        throw std::invalid_argument("col_idx and values sizes must match");
    }

    if (row_ptr_.front() != 0) {
        throw std::invalid_argument("row_ptr must start at zero");
    }

    const auto stored_nnz = values_.size();
    if (row_ptr_.back() != stored_nnz) {
        throw std::invalid_argument("row_ptr must end at nnz");
    }

    for (size_type position = 0; position < row_ptr_.size(); ++position) {
        if (row_ptr_[position] > stored_nnz) {
            throw std::invalid_argument("row_ptr entry exceeds nnz");
        }
        if (position > 0 && row_ptr_[position] < row_ptr_[position - 1]) {
            throw std::invalid_argument("row_ptr must be monotonically nondecreasing");
        }
    }

    for (size_type row = 0; row < rows_; ++row) {
        const auto row_begin = row_ptr_[row];
        const auto row_end = row_ptr_[row + 1];

        for (size_type entry = row_begin; entry < row_end; ++entry) {
            if (col_idx_[entry] >= cols_) {
                throw std::invalid_argument("Column index is outside the matrix bounds");
            }

            if (entry > row_begin && col_idx_[entry - 1] >= col_idx_[entry]) {
                throw std::invalid_argument(
                    "Column indices must be strictly increasing within each row");
            }

            if (!std::isfinite(values_[entry].real()) ||
                !std::isfinite(values_[entry].imag())) {
                throw std::invalid_argument("Matrix values must have finite real and imaginary parts");
            }
        }
    }
}

}  // namespace wavecore::helmholtz
