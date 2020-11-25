#ifndef IT_VEC_METADATA_H
#define IT_VEC_METADATA_H

#include "metadata_base_it.h"

////////////////////////////// MetaData Iterator ////////////////////////////
/** Iterates over metadatas */
struct MDVecIterator : public MDBaseIterator {
    MDBaseIterator(const MetaData &md) {}
    MDBaseIterator& operator++() override;
    bool operator==(const MDBaseIterator& other) override;
};


////////////////////////////// MetaData Row Iterator ////////////////////////////
/** Iterates over metadata rows */
struct MDVecRowIterator : public MDBaseRowIterator {
    MDBaseRowIterator(const MetaData &md) {}
    MDBaseRowIterator operator++() override;
    bool operator==(const MDBaseRowIterator& other) override;
    MDRow& operator*() const override;
};

#endif
