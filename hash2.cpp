#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>
#include <bitset>

using namespace seqan;

const unsigned shapelength = 30;
const unsigned shapeweight = 23;

typedef Iterator<DnaString >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength, shapeweight> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;

int hTest(StringSet<DnaString>  &  genome)
{
    std::cout << "hTest() " << std::endl;
    TShape shape;
    double time = sysTime();
    for (uint64_t j = 0; j < length(genome); j++)
    {
        TIter it = begin(genome[j]);
        hashInit(shape, it);
        for (uint64_t k = 0; k < length(genome[j]); k++)
        {
            hashNext(shape, it + k);           
            if (shape.hValue ^ xy2h(shape, shape.XValue, shape.YValue))
                std::cout << "    " << shape.hValue << " " << xy2h(shape, shape.XValue, shape.YValue) << std::endl;
            //std::cout << "    " << std::bitset<64>(shape.hValue) << " " << std::bitset<64>(xy2h(shape, shape.XValue, shape.YValue)) << std::endl;
        }
    }
    std::cout << "    End hTest() " << sysTime() - time << " " << shape.XValue << std::endl;
    return 0;
}

int main(int argc, char** argv)
{
   
    SeqFileIn rFile(toCString(argv[2]));
    SeqFileIn gFile(toCString(argv[1]));
    
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    StringSet<DnaString> genome;
    StringSet<String<uint64_t> > hs;

    readRecords(ids, reads, rFile);
    readRecords(ids, genome, gFile);

    hTest(genome);
}

