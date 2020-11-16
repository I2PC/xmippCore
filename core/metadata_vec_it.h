#ifndef IT_VEC_METADATA_H
#define IT_VEC_METADATA_H

////////////////////////////// MetaData Iterator ////////////////////////////
/** Iterates over metadatas */
class MDVecIterator
{
protected:
    size_t * objects;
    size_t size;

    /** Clear internal values to be used again*/
    void clear();
    /** Initialize internal values to NULL */
    void reset();
public:

    /** Internal function to initialize the iterator */
    void init(const MetaData &md, const MDQuery * pQuery=NULL);
    /** Empty constructor */
    MDIterator();
    /** Empty constructor, creates an iterator from metadata */
    MDIterator(const MetaData &md);
    /** Same as before but iterating over a query */
    MDIterator(const MetaData &md, const MDQuery &query);
    /** Destructor */
    ~MDIterator();

    /** This is the object ID in the metadata, usually starts at 1 */
    size_t objId;
    /** This is the index of the object, starts at 0 */
    size_t objIndex;
    /** Function to move to next element.
     * return false if there aren't more elements to iterate.
     */
    bool moveNext();
    /** Function to check if exist next element
     */
    bool hasNext();
}
;//class MDIterator


////////////////////////////// MetaData Row Iterator ////////////////////////////
/** Iterates over metadata rows */
class MDVecRowIterator
{
protected:

    // Current row data.
    MDRow currentRow;

    // Flag set to true if a new row has been retrieved.
    bool rowReturned;

public:

    /** Empty constructor, creates an iterator from metadata */
    MDRowIterator(MetaData &md);

    // Get current row data.
    MDRow *getRow(void);

    /** Function to move to next element.
     * return false if there aren't more elements to iterate.
     */
    bool moveNext(MetaData &md);

    /** Function to check if exist next element */
    bool hasNext();
}
;//class MDRowIterator



#endif
