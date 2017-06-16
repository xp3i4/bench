#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>
#include <bitset>
#include <vector>
// Evaluate performance of hash functions


using namespace seqan;

const unsigned shapelength = 5;
const unsigned shapeweight = 3;
const unsigned blocklimit = 32;

typedef Iterator<String<Dna> >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength, shapeweight> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;

typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength, shapeweight>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
typedef Index<StringSet<DnaString>, IndexQGram<SimpleMShape, OpenAddressing > > TMIndex;

void uTest(StringSet<DnaString> &genome, StringSet<DnaString> &read )
{
    TIndex_u index(genome);   
    indexCreate(index, FibreSADir());
    TShape shape;
    for (uint64_t k = 0; k < length(read); k++)
    {
        hashInit(shape, begin(read[k]));
        for (uint64_t j = 0; j < length(read[k]); j++)
        {
            hashNext(shape, begin(read[k]) + j);
            for (uint64_t n = index.dir[getOccurrenced(index, shape.hValue)]; n < index.dir[getOccurrenced(index, shape.hValue) + 1]; n++)  
            {
                sum ^= getSA(index.dir[n]);
            }
        }
    }
}

int main(int argc, char** argv)
{
    if (argc < 3)
        return 1;
    //time = sysTime();
    SeqFileIn rFile(toCString(argv[1]));
    SeqFileIn gFile(toCString(argv[2]));
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    StringSet<DnaString> genome;
    readRecords(ids, reads, rFile);
    readRecords(ids, genome, gFile);
    unsigned step = atof(toCString(argv[3]));
    unsigned l = atof(toCString(argv[4]));
    std::cout << "read done sysTime " << sysTime() - time << std::endl;
    uTest(reads, genome);
    return 0;
}
