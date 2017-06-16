#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>

using namespace seqan;

int main(int argc, char** argv)
{
    FragmentStore<> fragStore;
    if(!loadContigs(fragStore, argv[1]))
        return 1;
    if(!loadReads(fragStore, argv[2]))
        return 1;
    
    typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
    typedef GetValue<TReadSeqStore>::Type TReadSeq;
    typedef Index<TReadSeqStore, IndexQGram<SimpleShape, OpenAddressing > > TIndex;
    typedef Fibre<TIndex, FibreSA>::Type const      TSA;
    typedef Infix<TSA>::Type                        TOccurrences;
    typedef Iterator<String<Dna5Q> >::Type TIter;
    typedef Shape<Dna5Q, SimpleShape > Shape;
    typedef Value<Shape>::Type hValueType;
    
    TIndex index(fragStore.readSeqStore);
    Shape t_shape;
    TOccurrences occ;
    hValueType sum=0;
    unsigned kmerLength = atof(argv[3]);
    resize(indexShape(index), kmerLength);
    resize(t_shape, kmerLength);
    
    double start=sysTime();
    indexCreate(index, FibreSADir()); 
    std::cout << "Index creat time: " << sysTime() - start << std::endl;
    unsigned count=0;
    
    start=sysTime();   
    for (unsigned k = 0; k < length(fragStore.contigStore); ++k){
        TIter it = begin(fragStore.contigStore[k].seq);
        hashInit(t_shape, it);
        unsigned itLength = end(fragStore.contigStore[k].seq) - begin(fragStore.contigStore[k].seq);
        for (unsigned j = 0; j < itLength - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j); 
            count += countOccurrences(index, t_shape);  
        }
    }
    std::cout << count << std::endl;
    std::cout << "Time: " << sysTime() - start <<std::endl;
}

