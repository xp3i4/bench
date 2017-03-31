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
const unsigned blocklimit = 32;


typedef Iterator<String<Dna> >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength, shapeweight> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;

typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength, shapeweight>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
typedef Value<TShape>::Type HValue;


template <unsigned TSpan, unsigned TWeight>
void _qgramClearDir(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing> > & index)
{
    typedef Shape<Dna, MinimizerShape<TSpan, TWeight> > TM_Shape;
    typedef typename Value<TM_Shape>::Type HValue;
    resize (indexDir(index), _fullDirLength(index) + lengthSum(indexText(index)) + 2);
    index.start = _fullDirLength(index);
    index._Empty_Dir_ = 0;
    for (HValue k = 0; k < length(index.dir); k++) 
    {
        index.dir[k] = _bitEmpty;
    }
    std::cout << "        _qgramClearDir():" << std::endl;
    std::cout << "            _fullDirLength(index) = " << _fullDirLength(index) << std::endl;
    std::cout << "            lengh(index.dir) = " << length(index.dir) << std::endl;
    std::cout << "            End _qgramClearDir()" << std::endl;
}

template <unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<Dna, MinimizerShape<TSpan, TWeight> > TM_Shape;
    typedef typename Value<TM_Shape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue, HValue, HValue> HTuple;
    typedef String<HTuple> Stringtuple;

    TM_Shape shape;
    Stringtuple hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();
    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    std::cout << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    std::cout << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;
    for(HValue k = 0; k < length(indexText(index)); k++)
    {
        TIter it = begin(indexText(index)[k]);
        hashInit(shape, it);
        for (HValue j = 0; j < length(indexText(index)[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            hs[m++] = std::make_tuple(shape.XValue, shape.hValue, shape.YValue, k, j);
        }
    }
    std::cout << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1, (HValue)1);
    std::sort(begin(hs), end(hs),
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) > std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    hs[length(hs) - 1] = std::make_tuple((HValue)1, (HValue)1, (HValue)0, (HValue)1, (HValue)1);
    std::cout << "            sort sysTime(): " << sysTime() - time << std::endl;
    HValue countx = 1, counth = 1, countdh = 1, countb = 1, hk = 0;
    for (HValue k = 1; k < length(hs); k++)
    {
        //std::cout << k << std::endl;
        if (std::get<1>(hs[k]) != std::get<1>(hs[k - 1]))
        { _setBodyNode(index.dir[index.start + hk], std::get<2>(hs[k-1]), _BodyType_code, counth);
            //std::cout << counth << std::endl;
            hk++;
            countb++;
            countdh++;
        }
        //else
        counth++;

        if (std::get<0>(hs[k]) != std::get<0>(hs[k - 1]))
        {
            if (countb < blocklimit)
            {
                requestDir(index.dir, index.start, _makeHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(index.start + hk - countb));
                for (HValue j = 0; j < countb; j++)
                    _setBodyType_Begin(index.dir[index.start + hk - countb]);
            }
            else
            {
                hk -= countb;
                requestDir(index.dir, index.start, _makeVtlHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(index.start + hk));
                for (HValue j = k - countx; j < k; j++)
                    if (std::get<1>(hs[j]) != std::get<1>(hs[j + 1]))
                    {
                        requestDir(index.dir, index.start, _makeHVlHeadNode(std::get<1>(hs[j])), _makeEmptyNode(index.start+hk));
                        _setBodyType_Begin(index.dir[index.start + hk]);
                        hk++;
                    }
            }
            countb = 0;
            countx = 1;
        }
        else
        {
            countx++;
        }
    }
    
    resize(index.dir, index.start + countdh + 10);
    index._Empty_Dir_ = index.start + countdh + 1;
    std::cout << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}

template <unsigned TSpan, unsigned TWeight>
void createQGramIndex(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > >& index)
{
    double time = sysTime(); 
    std::cout << "    createQGramIndex() sysTime(): " << sysTime() - time << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams(index);
    std::cout << "        index.dir[index._Empty_Dir_] = " << index.dir[index._Empty_Dir_] << std::endl << "        length(index.dir) = " << length(index.dir) << std::endl;
    std::cout << "        End createQGramIndex() sysTime(): " << sysTime() - time << std::endl;
}

int mTest1(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape;
    TIndex index(reads);
    uint64_t sum = 0, dirh = 0;
    double time = sysTime();
    std::cout << "mTest1(): " << std::endl;
    createQGramIndex(index);
    std::cout << "    lenght Dir = " << length(index.dir) << " length Text = " << lengthSum(indexText(index)) << std::endl;
    std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
    for(uint64_t k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            sum = index.dir[getDir(index, shape)];
        }
    }
    sum=_getBodyCounth(sum);
    std::cout << "    sum = " << sum << std::endl;
    std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    End mTest1()" << std::endl;

    return 0;
}

int uTest(StringSet<DnaString> & reads, StringSet<DnaString> & genome, float alpha)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t sum=0, s=0;
    double time = sysTime();
    std::cout << "uTest():\n";
    std::cout << "    fullDirLength " << _fullDirLength(index) << std::endl; 
    resize(indexDir(index), _fullDirLength(index)   , 0); 
    resize(index.bucketMap.qgramCode, _fullDirLength(index)  , BucketMap<uint64_t>::EMPTY); 
    for (unsigned k = 0; k < length(index.bucketMap.qgramCode); k++) {
        index.bucketMap.qgramCode[k] = -1;
    }

    for(unsigned k = 0; k < length(reads); k++)
    {
        TIter it = begin(reads[k]);
        s+=length(reads[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(reads[k]) - t_shape.span + 1; j++)
        {
            hashNext(t_shape, it + j);
            index.dir[requestBucket(index.bucketMap, t_shape.hValue)]++;
        }
    }
    std::cout << "    getBucket start sysTime(): " << sysTime() - time<< std::endl;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j);
            sum+=index.dir[getBucket(index.bucketMap, t_shape.hValue)];
        }
    }
    std::cout << "    sum = " << sum << std::endl;
    std::cout << "    getBucket end sysTime(): "<< sysTime() - time<< std::endl;
    std::cout << "    End uTest()" << std::endl;
    return 0;
}

int umTest(StringSet<DnaString> & reads,  StringSet<DnaString> & genome)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t s=0, sum=0;
    std::cout << "umTest()" << std::endl;
    resize(indexDir(index), _fullDirLength(index)   , 0);
    resize(index.bucketMap.qgramCode, _fullDirLength(index)  , BucketMap<uint64_t>::EMPTY);
    for (unsigned k = 0; k < length(index.bucketMap.qgramCode); k++)
    {
        index.bucketMap.qgramCode[k] = BucketMap<uint64_t>::EMPTY;
    }

    TShape shape;
    TIndex index1(reads);
    createQGramIndex(index1); 

    double time=sysTime();
    for(unsigned k = 0; k < length(reads); k++)
    {
        TIter it = begin(reads[k]);
        sum+=length(reads[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(reads[k]) - t_shape.span + 1; j++)
        {
            hashNext(t_shape, it + j);

            index.dir[requestBucket(index.bucketMap, t_shape.hValue)]++;
        }
    }
    std::cout << "    h2y(shape, BucketMap<uint64_t>::EMPTY) = " << h2y(shape, BucketMap<uint64_t>::EMPTY) << std::endl;
    uint64_t v1, v2;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j);
            hashNext(shape, it + j);
            if ((uint64_t)t_shape.hValue - (uint64_t)shape.hValue)
                std::cout << "    hValue unequal" << j << " " <<  t_shape.hValue << " " <<  shape.hValue << std::endl;
            v1=h2y(shape, index.bucketMap.qgramCode[getBucket(index.bucketMap, t_shape.hValue)]);
            v2=_getBodyValue(index1.dir[getDir(index1, shape)]);
            if (v1 != v2)
                if (v1 ^ h2y(shape, BucketMap<uint64_t>::EMPTY) || v2 ^ _bitEmpty)
                    std::cout << "    YValue unequal " << t_shape.hValue << " " << getDir(index1, shape) << " " << shape.hValue << " " << v2 << std::endl;
        }
    }
    std::cout << "    s = " << s << std::endl;
    std::cout << "    End umTest() " << std::endl;
    return 0;
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
    //std::cout << "read done sysTime " << sysTime() - time << std::endl;
    //umTest(reads, genome);
    uTest(reads, genome, 1.8);
    mTest1(reads, genome);
    return 0;
}
