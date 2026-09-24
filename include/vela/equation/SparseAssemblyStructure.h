#pragma once

#include "vela/core/Types.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace vela {
/// Sweep-local structural data only. Serial access; matrices returned to callers
/// own their values. The scatter cursor is checked even when emission order changes.
struct SparseAssemblyStructure {
    struct Slot { Index row, column; Eigen::Index offset; };
    std::vector<std::uint64_t> key;
    SparseMatrixd pattern;
    std::vector<Slot> scatter;
    std::size_t builds=0, hits=0;

    bool matches(const std::vector<std::uint64_t>& candidate) const {
        return pattern.rows()!=0 && key==candidate;
    }
    void build(std::vector<std::uint64_t> identity, Index size,
               const std::vector<Eigen::Triplet<Real>>& entries) {
        SparseMatrixd next(size,size);
        next.setFromTriplets(entries.begin(),entries.end());
        next.makeCompressed();
        std::fill_n(next.valuePtr(),next.nonZeros(),0.);
        pattern=std::move(next);key=std::move(identity);scatter.clear();++builds;
    }
    SparseMatrixd zeroMatrix() const { return pattern; }
    void add(SparseMatrixd& matrix,std::size_t& cursor,Index row,Index column,Real value) {
        Eigen::Index offset;
        if(cursor<scatter.size() && scatter[cursor].row==row && scatter[cursor].column==column) {
            offset=scatter[cursor].offset;
        } else {
            const auto first=pattern.outerIndexPtr()[column],last=pattern.outerIndexPtr()[column+1];
            const auto* found=std::lower_bound(pattern.innerIndexPtr()+first,pattern.innerIndexPtr()+last,row);
            if(found==pattern.innerIndexPtr()+last || *found!=row)
                throw std::logic_error("Assembly contribution outside prepared sparse structure");
            offset=found-pattern.innerIndexPtr();
            const Slot slot{row,column,offset};
            if(cursor<scatter.size())scatter[cursor]=slot;else scatter.push_back(slot);
        }
        ++cursor;
        // Match the legacy omission of exact zero contributions (including -0).
        if(value!=0.)matrix.valuePtr()[offset]+=value;
    }
};
} // namespace vela
