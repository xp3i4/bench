//![header]
/*==========================================================================
  SeqAn - The Library for Sequence Analysis
  http://www.seqan.de
  ===========================================================================
  Copyright (C) 2010
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  ===========================================================================
  Minimizer benchmark 
  ===========================================================================*/
//![header]
//![includes]
#include <cstdio>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>

#include <time.h>
#include <math.h>

using namespace seqan;
//![includes]

const unsigned SHAPE_LENGTH = 23;
const unsigned SHAPE_WEIGHT = 5;

//![typedefs]
typedef FragmentStore<>::TReadSeqStore TReadSeqStore;
typedef Value<TReadSeqStore>::Type TReadSeq;
typedef FragmentStore<>::TContigStore TContigStore;
typedef Value<TContigStore>::Type TContigStoreElement;
typedef TContigStoreElement::TContigSeq TContigSeq;
typedef Index<TReadSeqStore, IndexQGram<UngappedShape<SHAPE_LENGTH>, OpenAddressing > > TIndex2;
typedef Index<TReadSeqStore, IndexQGram< MinimizerShape<SHAPE_LENGTH, SHAPE_WEIGHT> > > TIndex;

typedef Shape<Dna5Q, UngappedShape<SHAPE_LENGTH> > Ungapped_L_Shape;
typedef Shape<Dna5Q, MinimizerShape<SHAPE_LENGTH, SHAPE_WEIGHT> > MiniShape;
typedef Value<Ungapped_L_Shape>::Type    THashValue;
typedef Value<MiniShape>::Type TMiniHashValue;

typedef Iterator<String<Dna5Q> >::Type TIter;
typedef Infix<String<Dna5Q> >::Type  TTextInfix;
String<Dna5Q> text;
typedef Fibre<TIndex, FibreSA>::Type const      TSA;
typedef typename Size< typename Fibre< TIndex, FibreSA>::Type const >::Type TOccurrences;

//![main-input]
int main(int argc, char** argv)
{
    std::cout << SHAPE_LENGTH <<std::endl;
    // 0) Handle command line arguments.
    if (argc < 2)
    {
        std::cerr << "Invalid number of arguments." << std::endl
                  << "USAGE: minimapper GENOME.fasta READS.fasta " << std::endl;
        return 1;
    }

    // 1) Load contigs and reads.
    FragmentStore<> fragStore;
    if (!loadContigs(fragStore, argv[1]))
        return 1;

    if (!loadReads(fragStore, argv[2]))
        return 1;

    TIndex index(fragStore.readSeqStore);
    TIndex2 index2(fragStore.readSeqStore);

    for (unsigned i = 0; i < length(fragStore.contigStore); ++i)
    { 
        double timeStart; 
        TIter it = begin(fragStore.contigStore[i].seq);
        unsigned itLength= end(fragStore.contigStore[i].seq) - begin(fragStore.contigStore[i].seq);
        timeStart = sysTime();
        indexCreate(index, FibreSADir()); 
        std::cout << "indexCreate time = " << (sysTime() - timeStart) * 1000 << std::endl;
 
        #define N 1024000
        //TOccurrences occ;
        TOccurrences count[N];
        for (int k = 0; k < N; k++)
            count[k] = 0;
        Ungapped_L_Shape u_shape;
        MiniShape t_shape; 
        
        timeStart = sysTime();
        hashInit(t_shape, it);
        for (unsigned j = 0; j < itLength - SHAPE_LENGTH + 1; j++)
        {
            hashNext(t_shape, it + j); 
            //std::cout << t_shape.hValue << std::endl;
            count[t_shape.hValue]++;
        }
        for (int k = 0; k < N; k++)
            if (count[k] != 0)
            std::cout << k << " " << count[k] << std::endl;
        std::cout << " hash time = " << (sysTime() - timeStart) * 1000 << std::endl;
    

/////////////////////////////ungappd hash//////////////////////
//        timeStart = sysTime();
//        indexCreate(index2, FibreSADir());
//        std::cout << "index2 Create time = " << (sysTime() - timeStart) * 1000 << std::endl;
//        
//        timeStart = sysTime();
//        hashInit(u_shape, it);
//        n = 0;
//        for (unsigned j = 0; j < itLength - SHAPE_LENGTH + 1; j++)
//        {
//            hashNext(u_shape, it + j);    
//            if (u_shape.hValue > 0)
//                n++;
//        }
//        std::cout << n << " ungapped hash time = " << (sysTime() - timeStart) * 1000 << std::endl;
//
//        occ = 0;    
//        hashInit(u_shape, it); 
//        timeStart = sysTime();
//        for (unsigned j = 0; j < itLength - SHAPE_LENGTH + 1; j++)
//        {
//            hashNext(u_shape, it + j);
//            occ += countOccurrences(index, u_shape);
//        }
//        std::cout << occ << " ungapped hash + countOccurrences time = " << (sysTime() - timeStart) * 1000 << std::endl;
//
//
    } 
    return 0;
}
//![main-output]
