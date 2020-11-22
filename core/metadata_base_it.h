#ifndef IT_VEC_METADATA_H
#define IT_VEC_METADATA_H

class MetaData;

////////////////////////////// MetaData Iterator ////////////////////////////
/** Iterates over metadatas */
struct MDBaseIterator {
    /** Empty constructor, creates an iterator from metadata */
    MDBaseIterator(const MetaData &md) {}
    virtual MDBaseIterator& operator++();
    virtual bool operator==(const MDBaseIterator& other);
    virtual bool operator!=(const MDBaseIterator& other) { return !(*this == other); }

    size_t objId; // for backwards-compatibility; TODO
};

////////////////////////////// MetaData Row Iterator ////////////////////////////
/** Iterates over metadata rows */
struct MDBaseRowIterator {
    MDBaseRowIterator(const MetaData &md) {}
    virtual MDBaseRowIterator operator++();
    virtual bool operator==(const MDBaseRowIterator& other);
    virtual bool operator!=(const MDBaseRowIterator& other) { return !(*this == other); }
    virtual MDRow& operator*() const;
};

#endif
