#ifndef IT_BASE_METADATA_H
#define IT_BASE_METADATA_H

#include <memory>
#include "metadata_row_base.h"

class MetaData;

////////////////////////////// MetaData Iterator ////////////////////////////
/** Iterates over metadata rows indexes */
/*struct MDBaseIterator {
    virtual std::unique_ptr<MDBaseIterator> clone() = 0;
    virtual void increment() = 0;
    virtual bool operator==(const MDBaseIterator& other);
    virtual bool operator!=(const MDBaseIterator& other) { return !(*this == other); }
};*/

////////////////////////////// MetaData Row Iterator ////////////////////////////
/** Iterates over metadata rows */
struct MDBaseRowIterator {
    MDBaseRowIterator(const MetaData &md) {}
    virtual std::unique_ptr<MDBaseRowIterator> clone() = 0;
    virtual void increment() = 0;
    virtual bool operator==(const MDBaseRowIterator& other) = 0;
    virtual bool operator!=(const MDBaseRowIterator& other) { return !(*this == other); }
    virtual MDRow& operator*() = 0;
};

#endif
