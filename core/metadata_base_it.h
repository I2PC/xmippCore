#ifndef IT_BASE_METADATA_H
#define IT_BASE_METADATA_H

#include <memory>
#include "metadata_row_base.h"
#include "choose.h"

class MetaData;

/** Iterates over metadata rows */
template <bool IsConst>
struct MDBaseRowIterator {
    virtual std::unique_ptr<MDBaseRowIterator> clone() = 0;
    virtual void increment() = 0;
    virtual bool operator==(const MDBaseRowIterator& other) const = 0;
    virtual bool operator!=(const MDBaseRowIterator& other) const { return !(*this == other); }
    virtual typename TypeHelpers::choose<IsConst, MDRowConst&, MDRow&>::type operator*() = 0;
};

/** Iterates over metadata ids */
template <bool IsConst>
struct MDBaseIdIterator {
    virtual std::unique_ptr<MDBaseIdIterator> clone() = 0;
    virtual void increment() = 0;
    virtual bool operator==(const MDBaseIdIterator& other) const = 0;
    virtual bool operator!=(const MDBaseIdIterator& other) const { return !(*this == other); }
    virtual size_t operator*() = 0;
};

#endif
