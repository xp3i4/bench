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
typedef Value<TShape>::Type HvalueType;



int _f(String<uint64_t> c, uint64_t empty=_dirEmpty)
{
    unsigned j = 0;
    uint64_t s = 0, s1=0;
    for (unsigned k = 0; k < length(c); k++)
        if (c[k]!=empty)
        s++;
    for (unsigned k = 0; k < length(c); k++)
    {
        if(c[k]!=empty)
            s1++;
        if ((float)k/length(c) > j * 0.01)
        {
            std::cout <<  (float)s1/s << " ";
            j++;
        }
    }
    return 1;
}

/*
//void _f1(string<pair<uint64_t, uint64_t> > hs, empty=_dirempty )
//{
//
//}
template < typename tdir,  typename tspec, typename tstring, typename tstepsize>
inline void
_qgramcountqgrams(tdir &dir, StringSet<tstring, tspec> const &text, tshape shape, tstepsize stepsize){
    seqan_checkpoint
        typedef typename Iterator<tstring const, standard>::type    tIterator;
        typedef typename value<tdir>::type                        tsize;
        typedef typename std::tuple<tsize, tsize, tsize>          thashtuple;

        //string<pair<tsize, tsize> > hs,hs1;
        string<thashtuple> hs, hs1;
        tsize num_qgrams = 0;

        for (unsigned seqno = 0; seqno < length(text); seqno++)
        {
            if (length(text[seqno]) < length(shape) || empty(shape))
                continue;
            num_qgrams += (length(text[seqno]) - length(shape)) / stepsize + 1;
        }

        resize(hs, num_qgrams);

        for (unsigned seqno = 0; seqno < length(text); seqno++)
        {
            tIterator it = begin(text[seqno], standard());
            hashInit(shape, it);
            if (stepsize == 1){
                for (tsize k = 0; k < num_qgrams; k++){
                    hashNext(shape, it + k);
                    hs[k] = std::make_tuple(shape.XValue, shape.hValue, shape.YValue);
                }
                        }

        //else{
        //}
        }
        std::sort(begin(hs), end(hs),
                [](const thashtuple &a, const thashtuple &b)
                {return (std::get<0>(a) > std::get<0>(b))||((std::get<0>(a) == std::get<0>(b)) && std::get<1>(a) > std::get<1>(b));});
        uint64_t count = 0, hsk = 0;
        for (unsigned k = 1; k < length(hs); k++)
        {
            if(std::get<1>(hs[k-1]) != std::get<1>(hs[k]))
                count++;
        }
        resize(hs1, count);
        for (unsigned k = 0; k < length(hs); k++)
        {
            if(std::get<1>(hs[k-1]) != std::get<1>(hs[k]))
            {
                hs1[hsk++] = std::make_tuple((count << 39)+(std::get<0>(hs[k])<<3), shape.hValue, shape.YValue);
                count = 1;
                std::cout << std::get<0>(hs1[k]) << std::endl;
            }
            else
            {
                count++;
            }
        }

    //setvalue(hs, hashNext());
    //sortvalue(hs)
    //count(hs);
    //request(hs1);
}
*/

//template <typename tdir>
//inline void _qgramcleardir(tdir &dir, typename tparalleltag)
//{
//
//}

/*
struct XHYValue{
    uint64_t i1,i2,i3,i4,i5;    
    inline int assignValue(uint64_t v1, uint64_t v2){//, uint64_t v3, uint64_t v4, uint64_t v5){
        i1=v1;i2=v2;
        //i3=v3;i4=v4;
        //i5=v5; 
    }
};
template <typename TDir>
int _qgramCountQGrams(TDir &dir, StringSet<DnaString>  &StringSet, TShape shape)
{

    std::cout << "mtest() _qgramcountqgrams " << length(StringSet) << " length dir" << length(dir)<< std::endl;
    typedef std::tuple<HvalueType, HvalueType, HvalueType, HvalueType, HvalueType> tuple;
    typedef String<std::tuple<HvalueType, HvalueType, HvalueType, HvalueType, HvalueType> > StringTuple;
    typedef Iterator<StringTuple>::Type TIter2;

    typedef Pair<HvalueType, HvalueType> pairh;
    typedef String<pairh> StringPair;
    StringTuple hs, hs1;
    StringPair hs2;
    String<XHYValue> hs3;
    HvalueType  m=0;
    uint64_t sum=0, count=0, pre;

    double start=sysTime();
    resize(hs, lengthSum(StringSet) - length(StringSet) * (shape.span - 1) + 1);
    resize(hs2, lengthSum(StringSet) - length(StringSet) * (shape.span - 1) + 1);
    resize(hs3, lengthSum(StringSet) - length(StringSet) * (shape.span - 1) + 1);

    std::cout << lengthSum(StringSet) << "length\nhash " << sysTime() - start <<"\n";
    for(HvalueType k = 0; k < length(StringSet); k++)
    {
        TIter it = begin(StringSet[k]);
        hashInit(shape, it);
        for (HvalueType j = 0; j < length(StringSet[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
        //    std::cout << std::bitset<64>(shape.XValue) << " " << std::bitset<64>(shape.hValue) << std::endl;
            hs[m++] = std::make_tuple(shape.XValue, shape.hValue, shape.YValue, k, j);
            assignValueI1(hs2[m-1], shape.XValue);
            assignValueI2(hs2[m-1], shape.hValue);
            //hs3[m-1].assignValue(shape.XValue, shape.hValue);//, shape.YValue, k, j);

            
            //std::cout << std::get<0>(hs[m -1]) << " " << std::get<1>(hs[m-1]) << std::endl;

        }
    }
    std::cout << std::get<0>(hs[length(hs)-2]) << " " << std::get<1>(hs[length(hs)-2]) << std::endl;
    std::cout << "hashtime " << sysTime() - start << std::endl;


    hs[length(hs) - 1] = std::make_tuple((HvalueType)0, (HvalueType)0, (HvalueType)0, (HvalueType)1, (HvalueType)1);
    std::cout << " sort1 " << sysTime() - start << std::endl;
    std::sort(begin(hs2), end(hs2),
        [](const pairh &a, const pairh &b)
        {return a.i1 > b.i1; });
    std::cout << "sort " << sysTime() -start << std::endl;
    std::sort(begin(hs), end(hs),
        [](const tuple &a, const tuple &b)
        {return std::get<0>(a) > std::get<0>(b);});//||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b));});
    std::cout << "sort time " << sysTime() - start << std::endl;

    hs[length(hs) - 1] = std::make_tuple((HvalueType)1, (HvalueType)1, (HvalueType)0, (HvalueType)1, (HvalueType)1);
    for (unsigned k = 1 ; k < length(hs); k++)
    {
        //if (k > length(hs) - 100)
        //std::cout << std::get<0>(hs[k]) << " " << std::bitset<64>(std::get<1>(hs[k])) << std::endl;
        if(std::get<0>(hs[k]) != std::get<0>(hs[k - 1]) )

        //if((std::get<0>(hs[k]) != std::get<0>(hs[k - 1])) && (std::get<1>(hs[k]) == std::get<1>(hs[k - 1] )))

            count++;
    }
    std::cout << "count " << count << "  " << lengthSum(StringSet) <<  " " << length(hs1) << sysTime() - start << std::endl;

    resize(hs1, count - 1);
    unsigned hsk=0;
    count = 1;
    sum = 0;
    for (uint64_t k = 1; k < length(hs); k++)
    {
        //std::cout << std::get<0>(hs[k]) << std::endl;
        if (std::get<0>(hs[k]) != std::get<0>(hs[k - 1]))
        {
            for (uint64_t j = k - count; j < k; j++)
                if (count < blocklimit)
                {
                    hs1[hsk++] = std::make_tuple((std::get<0>(hs[j]) << _bitvalue_end) + ((k - j) << 52) + _bitheadcode, (std::get<1>(hs[j]) << _bitvalue_end) + _bitheadcode,
                                (std::get<2>(hs[j]) << _bitvalue_end) + ((k - j) << _bitlength_end) + _bitbodycode, std::get<3>(hs[j]), std::get<4>(hs[j]));
                }
                else
                {
                    hs1[hsk++] = std::make_tuple((std::get<0>(hs[j]) << _bitvalue_end) +((k - j) << 52) + _bitvtlheadcode, (std::get<1>(hs[j]) << _bitvalue_end) + _bitheadcode,
                                (std::get<2>(hs[j]) << _bitvalue_end) + (1 << _bitlength_end) + _bitbodycode, std::get<3>(hs[j]), std::get<4>(hs[j]));
                }
            sum += count;
            count = 1;
        }
        else
            if (std::get<1>(hs[k]) != std::get<1>(hs[k - 1]))
                count++;
    }
//    std::cout << length(hs1) << " " << hsk << " " << sum << " sort " << sysTime() - start << std::endl;
//    std::sort(begin(hs1), end(hs1),
//        [](const tuple &a, const tuple &b)
//        //{return (std::get<0>(a) > std::get<0>(b) || (std::get<0>(a) == std::get<0>(b) &&
//        //        std::get<3>(a) == std::get<3>(b) && std::get<4>(a) > std::get<4>(b)));});
//        {return (std::get<0>(a) > std::get<0>(b));} );
//
//    std::cout << sysTime() - start << std::endl;
//    uint64_t dirh = 0, hk = 0, blocklength;
//    bool flag = true;
//    uint64_t key = ((uint64_t)1 << 52) - 1;
//    sum = 0;
//    unsigned u=0;
//
//    while (hk < length(hs1))
//    {
//        blocklength = std::get<0>(hs1[hk]) >> 52;
//        if(blocklength < blocklimit)
//        {
//            dirh = requestbucket(dir, std::get<0>(hs1[hk]) & key, blocklength + 1);
//            dir[dirh] = (std::get<0>(hs[hk] ) & key);
//            for (unsigned j = 1; j <= blocklength; j++)
//            {
//                dir[dirh + j] = std::get<2>(hs1[hk++]);
//            }
//        }
//        else
//        {
//            dirh = requestbucket(dir, std::get<0>(hs1[hk]) & key, 1);
//            dir[dirh] = (std::get<0>(hs1[hk]) & key);
//            for (unsigned j = 0; j < blocklength; j++ )
//            {
//                dirh = requestbucket(dir, std::get<1>(hs1[hk]), 2);
//                dir[dirh] = std::get<1>(hs1[hk]);
//                dir[dirh + 1] = std::get<2>(hs1[hk++]);
//            }
//        }
//        sum += dir[dirh];
//    }
//
//    while (hk < length(hs1))
//    {
//        blocklength = std::get<0>(hs1[hk]) >> 52;
//        //requestbucket(dir, );
//        _fillblock(dir, it2);
//    }
//
//    std::cout << dir[hk] << sum  << length(hs1)  << std::endl << sysTime() - start << std::endl;
    return 0;
}
*/


template <typename TDir>
int _qgramCountQGrams(TDir& dir, StringSet<DnaString>& StringSet, TShape shape, uint64_t dirStart)
{

    std::cout << "mtest() _qgramcountqgrams " << length(StringSet) << " length dir" << length(dir)<< std::endl;
    typedef std::tuple<HvalueType, HvalueType, HvalueType, HvalueType, HvalueType> HTuple;
    typedef String<HTuple> Stringtuple;
    Stringtuple hs, hs1;
    String<TShape> hs3;
    HvalueType  m=0;
    uint64_t sum=0;

    double start=sysTime();
    resize(hs, lengthSum(StringSet) - length(StringSet) * (shape.span - 1) + 1);
    resize(hs3, lengthSum(StringSet) - length(StringSet) * (shape.span - 1) + 1);

    std::cout << lengthSum(StringSet) << "length\n";
    for(HvalueType k = 0; k < length(StringSet); k++)
    {
        TIter it = begin(StringSet[k]);
        hashInit(shape, it);
        for (HvalueType j = 0; j < length(StringSet[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            hs[m++] = std::make_tuple(shape.XValue, shape.hValue, shape.YValue, k, j);
        }
    }
    std::cout << "hashtime " << sysTime() - start << std::endl;
    hs[length(hs) - 1] = std::make_tuple((HvalueType)0, (HvalueType)0, (HvalueType)0, (HvalueType)1, (HvalueType)1);
    std::sort(begin(hs), end(hs),
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) > std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    std::cout << sysTime() - start << std::endl;
    hs[length(hs) - 1] = std::make_tuple((HvalueType)1, (HvalueType)1, (HvalueType)0, (HvalueType)1, (HvalueType)1);

    unsigned hsk=0;
    uint64_t countx = 1, counth = 1, countb = 1, hk = 0;
    std::cout << "done" << std::endl;
    for (uint64_t k = 1; k < length(hs); k++)
    {
        if (std::get<1>(hs[k]) != std::get<1>(hs[k - 1]))
        {
            _setBodyNode(dir[dirStart + hk], std::get<2>(hs[k-1]), _BodyType_code, counth);
            //std::cout << counth << std::endl;
            hk++;
            countb++;
        }
        else
            counth++;

        if (std::get<0>(hs[k]) != std::get<0>(hs[k - 1]))
        {
            if (countb < blocklimit)
            {
                requestDir(dir, dirStart, _makeHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(dirStart + hk - countb));
                for (unsigned j = 0; j < countb; j++)
                    _setBodyType_Begin(dir[dirStart + hk - countb]);
            }
            else
            {
                hk -= countb;
                requestDir(dir, dirStart, _makeVtlHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(dirStart + hk));
                for (uint64_t j = k - countx; j < k; j++)
                    if (std::get<1>(hs[j]) != std::get<1>(hs[j + 1]))
                    {
                        requestDir(dir, dirStart, _makeHVlHeadNode(std::get<1>(hs[j])), _makeEmptyNode(dirStart+hk));
                        _setBodyType_Begin(dir[dirStart + hk]);
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
    std::cout << "done" << std::endl;
    resize(dir, dirStart + counth + 10);
    return 0;
}


/*
    template<typename TDir, typename THashValue>
    inline THashValue
    _fillBlock(TDir & dir, TIter2 it)//, Tag<TParallelTag> parallelTag)
    {
        typedef unsigned long TSize;
        TSize hlen = length(dir);
        if (hlen == 0ul) return code;

        TSize h1 = _hashFunction(dir, code >> _bitValue_END);
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif
        TSize delta = 0;
        (void)delta;
        unsigned count = 0;
        for (auto k = 0; k < len; k++)
        {
            switch (dir[h1+k] & _bitCode){
                case 1:
                    h1 = (0 + h1 + (dir[h1 + k] >> _bitLength_END) + delta) &  hlen ;
                    ++delta; k=0;
                    count++;
                    break;
                case 2:
                    h1 = (0 + h1 + (dir[h1 + k + 1] >> _bitLength_END) + delta) &  hlen ;
                    ++delta; k=0;
                    count++ ;
                    break;
                case 3:
                    std::cout << k << std::endl;
                    h1 = (h1 + 1 + delta) &  hlen ;
                    ++delta; k=0;
                    count++;
                    break;
            }
        }
        return h1;
    }
*/

int mTest1(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape;
    TIndex index(reads);
    //unsigned count = 0;
    uint64_t sum = 0, dirh = 0;
    //resize(indexDir(index), _fullDirLength(index));
    resize (indexDir(index), _fullDirLength(index) + lengthSum(reads));
    index.start = _fullDirLength(index);
    index._Empty_Dir_ = index.start - 2;

    std::cout << "mTest(): index dir length " << _fullDirLength(index)<< " length(index.dir)= " << length(index.dir) << std::endl;

    for (unsigned k = 0; k < length(indexDir(index)); k++)
    {
        indexDir(index)[k] = _bitEmpty;
    }
    double start = sysTime();
    _qgramCountQGrams(indexDir(index), indexText(index), shape, index.start);
    std::cout << " _qgramCountQGrams() length(index.dir) = " << length(index.dir) << std::endl;
    //for (uint64_t k = 0; k < length(indexDir(index)); k++)
    //{
    //    std::cout << indexDir(index)[k] << std::endl;
    //}

//    std::cout << sysTime() - start << std::endl;
//
    std::cout << sysTime() - start << std::endl;
    for(uint64_t k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            //sum = _getBodyCounth(index.dir[getBucket(indexDir(index), shape.XValue, shape.YValue, shape.hValue)+1] - index.dir[getBucket(indexDir(index), shape.XValue, shape.YValue, shape.hValue)]);
             sum = index.dir[getDir(index, shape)];
        }
    }
    sum=_getBodyCounth(sum);
//    for (uint64_t k = 0; k < _fullDirLength(index); k++)
//    {
//        if(index.dir[k] == 4637156639017092940)
//        std::cout << k;
//    }

    std::cout << sum << " " << sysTime() - start << std::endl;

    return 0;
}

/*
int mTest(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    StringSet<String<uint64_t> > h1s;
    String<Pair<uint64_t, uint64_t> > hs,hs1;
    uint64_t sum=0, count=0, pre;
    TShape t_shape;
    TIndex index(reads);
    unsigned kmerLength = t_shape.span;
    float alpha = 1.8;
    resize(indexDir(index), (_fullDirLength(index)<<1) - 1);

    std::cout << "mTest() " << length(indexText(index)[0]) << " " << length(indexDir(index))<< std::endl;

    //_qgramCountQgrams(indexDir(index), indexText(index), t_shape, 1);

    for (unsigned k = 0; k < length(indexDir(index)); k++)
    {
        indexDir(index)[k] = _dirEmpty;
    }

    String<unsigned> cc;
    resize(cc, lent64_t j = k - count; j < k; j++)
    for (unsigned k =0; k < length(cc); k++)
        cc[k]=0;

    double start=sysTime();
    std::cout << start << std::endl;
    for(unsigned k = 0; k < length(reads); k++)
    {
        TIter it = begin(reads[k]);
        hashInit(t_shape, it);
        resize(hs, length(reads[k]));
        for (uint64_t j = 0; j < length(reads[k]) - kmerLength + 1; j++)
        {
            sum+=hashNext(t_shape, it + j);
            //std::cout << std::bitset<64>(t_shape.hValue) << " " << std::bitset<64>(t_shape.XValue) << std::endl;
            setValueI1(hs[j], t_shape.XValue);
            setValueI2(hs[j], t_shape.hValue);
//            std::cout << std::bitset<64>(t_shape.hValue) << " " << std::bitset<64>(t_shape.XValue) << std::endl;
        }
    }
    std::cout << "hashTime " << sysTime() - start << std::endl;
    pre = hs[0].i1;
    count = 0;
    unsigned cp_hs[10000] = {0};
    for (unsigned k = 0; k < length(hs); k++)
    {
        //std::cout << hs[k] << std::endl;
        if(pre != hs[k].i1)
        {
            cp_hs[count]++;
            pre=hs[k].i1;
            count=1;
        }
        else
            count++;
    }
    unsigned sum_1 = 0;
    for (unsigned k = 0; k < 10000; k++)
    {
        if (cp_hs[k] != 0)
            sum_1 += cp_hs[k] * k;

    }
    unsigned tmp = 0;
    for (unsigned k = 0; k < 10000; k++)
    {
        if (cp_hs[k]!=0)
        {
            tmp += cp_hs[k] * k ;
            std::cout << k << " " << (float)tmp / sum_1 << std::endl;
        }
    }
    std::cout << std::endl;
    std::sort(begin(hs), end(hs),
        [](const Pair<uint64_t, uint64_t> &a, const Pair<uint64_t, uint64_t> &b)
        {return (a.i1 > b.i1)||((a.i1 == b.i1) && a.i2 > b.i2);});
    for (unsigned k = 1 ; k < length(hs); k++)
    {
        if(hs[k].i1 != hs[k-1].i1 || hs[k].i2 != hs[k-1].i2)
            count++;
    }
    resize(hs1, count);
    unsigned hsk=0, cp[10000]={0};
    pre = hs[0].i1;
    count = 0;
    for (unsigned k = 1 ; k < length(hs); k++)
    {
        if(hs[k].i1 != hs[k-1].i1 || hs[k].i2 != hs[k-1].i2)
        {
            hs1[hsk++] = hs[k];
        }

    }
    for(unsigned k = 0; k < length(hs1); k++)
    {
        if (hs1[k].i1 == pre)
        {
            count++;
        }
        else
        {
            cp[count]++;
            for (unsigned j = k - count; j < k; j++)
            {
                hs1[j].i1 += (count << 56);
            }
            if (count > 800)
                std::cout << count << " " << std::bitset<64>(hs1[k-1].i1) << std::endl;//" " << std::bitset<64>(hs1[k-1].i2) << std::endl;

            //std::cout << std::endl;
            count=1;
            pre = hs1[k].i1;
        }

        //std::cout << hs1[k] << " " << hs1[k].i1 << std::endl;
    }
    for (unsigned j = length(hs1) - count; j < length(hs1); j++)
    {
        hs1[j].i1 += (count << 56);
    }

//    unsigned sum_cp=0;
//    tmp=0;
//    for (unsigned j = 0; j < 10000; j++)
//    {
//        sum_cp+=cp[j] * j;
//    }
//    for (unsigned j = 0; j < 10000; j++)
//    {
//        if (cp[j] != 0)
//        {
//            tmp += cp[j] * j;
//            std::cout << j << " " << (float)tmp / sum_cp << std::endl;
//        }
//    }
    std::sort(begin(hs1), end(hs1),
        [](const Pair<uint64_t, uint64_t> &a, const Pair<uint64_t, uint64_t> &b)
        {return (a.i1 > b.i1)||((a.i1 == b.i1) && a.i2 > b.i2);});

    std::cout << sysTime() - start << std::endl;
    unsigned l=0,cm=0;
    uint64_t s=0;

    for (unsigned k = 0; k < length(hs1); k++)
    {
        //std::cout << (hs1[k].i1 >> 56)  <<" " << k << std::endl;
        if((hs1[k].i1 >> 56) < blockLimit)
        {
            cc[k]=requestBucket(indexDir(index), hs1[k].i1 << 8 >> 8, hs1[k].i1 >> 56);
            s+=cc[k];
        }
        else
        {
            cm++;
        }
    }
    std::cout << "length " << length(hs1) << " " << cm << std::endl;
    for (unsigned k = 0; k < cm; k+=(hs1[k].i1 >> 56))
    {
        s+=requestBucket(indexDir(index), hs1[k].i1, (uint64_t)1);
        for (unsigned j = k; j < k + (hs1[k].i1 >> 56); j++)
        {
            s+=requestBucket(indexDir(index), hs1[k].i2, (uint64_t)1);
            //std::cout << j << " " << k << std::endl;
        }
    }
//    for (unsigned k = 0; k < length(hs1); k++)
//    {
//        if ((hs1 >> 56) > blockLimit)
//
//    }
    //requestBucket(indexDir(index), hs1);
    std::sort(begin(cc), end(cc),[](const unsigned &a, const unsigned &b){return a > b;});
    std::cout << s << " " << length(hs1)  << std::endl << sysTime() - start << std::endl;
    //_f(indexDir(index));
    return 0;
}
*/
int uTest(StringSet<DnaString> & reads, StringSet<DnaString> & genome, float alpha)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t s=0, sum=0;
    //float alpha = 1.6;
    std::cout << "uTest()\nfullDirLength " << _fullDirLength(index) << " alpha" << alpha << std::endl;
    resize(indexDir(index), _fullDirLength(index)   , 0);
    resize(index.bucketMap.qgramCode, _fullDirLength(index)  , BucketMap<uint64_t>::EMPTY);
    //resize(indexDir(index), length(reads[0]) * alpha , 0);
    //resize(index.bucketMap.qgramCode, length(reads[0]) * alpha, BucketMap<uint64_t>::EMPTY);
    for (unsigned k = 0; k < length(index.bucketMap.qgramCode); k++)
    {
        index.bucketMap.qgramCode[k] = -1;
        //std::cout << index.bucketMap.qgramCode[k] << std::endl;
    }

    double start=sysTime();
    for(unsigned k = 0; k < length(reads); k++)
    {
        TIter it = begin(reads[k]);
        sum+=length(reads[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(reads[k]) - t_shape.span + 1; j++)
        {
            hashNext(t_shape, it + j);
            //s+=_hashFunction(index.bucketMap, t_shape.hValue);
            //s+=((t_shape.hValue * 43) ^ (t_shape.hValue >> 20)) + t_shape.hValue;

            //std::cout << s << std::endl;
            index.dir[requestBucket(index.bucketMap, t_shape.hValue)]++;
            //s+=cc[j];
        }
    }
    std::cout << sysTime() - start << std::endl;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j);
            //if(index.bucketMap.qgramCode[getBucket(index.bucketMap, t_shape.hValue)] != t_shape.hValue)
            //    std::cout <<index.bucketMap.qgramCode[getBucket(index.bucketMap, t_shape.hValue)] << " " << t_shape.hValue << std::endl;
            //else
            //    std::cout << t_shape.hValue<< std::endl;
            s+=index.dir[getBucket(index.bucketMap, t_shape.hValue)];
        }
    }
    //std::sort(begin(cc), end(cc),[](const unsigned &a, const unsigned &b){return a > b;});
    std::cout << s << std::endl << sysTime() - start << std::endl;
    //_f(index.bucketMap.qgramCode, BucketMap<uint64_t>::EMPTY);
  //  start = sysTime();
  //  for(unsigned k = 0; k < length(genome); k++)
  //  {
  //      TIter it = begin(genome[k]);
  //      hashInit(t_shape, it);
  //      for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
  //      {
  //          hashNext(t_shape, it + j);
  //          s += getBucket(index.bucketMap, t_shape.hValue);
  //          //s+=cc[j];
  //      }
  //  }

  //  //std::sort(begin(cc), end(cc),[](const unsigned &a, const unsigned &b){return a > b;});
  //  std::cout << s << std::endl << sysTime() - start << std::endl;

    return 0;
}

int umTest(StringSet<DnaString> & reads,  StringSet<DnaString> & genome)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t s=0, sum=0;
    std::cout << "umTest()\nfullDirLength " << std::endl;
    resize(indexDir(index), _fullDirLength(index)   , 0);
    resize(index.bucketMap.qgramCode, _fullDirLength(index)  , BucketMap<uint64_t>::EMPTY);
    for (unsigned k = 0; k < length(index.bucketMap.qgramCode); k++)
    {
        index.bucketMap.qgramCode[k] = -1;
    }

    TShape shape;
    TIndex index1(reads);
    uint64_t dirh = 0;
    resize (indexDir(index1), _fullDirLength(index1) + lengthSum(reads)+2);
    index1.start = _fullDirLength(index1);
    index1._Empty_Dir_ = length(index1) - 2;

    for (unsigned k = 0; k < length(indexDir(index1)); k++)
    {
        indexDir(index1)[k] = _bitEmpty;
    }
    _qgramCountQGrams(indexDir(index1), indexText(index1), shape, index1.start);


    double start=sysTime();
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
    std::cout << sysTime() - start << std::endl;
    uint64_t v1, v2;
    std::cout << index1.start << std::endl;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        hashInit(shape, it);
        //hash(t_shape, it);
        std::cout << " INIT " << t_shape.hValue << " " << shape.hValue  << std::endl;
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j);
            hashNext(shape, it + j);
            if ((uint64_t)t_shape.hValue - (uint64_t)shape.hValue)
                std::cout << j << "hValue " << t_shape.hValue << " " <<  shape.hValue << std::endl;
            v1=h2y(shape, index.bucketMap.qgramCode[getBucket(index.bucketMap, t_shape.hValue)]);
            v2=_getBodyValue(index1.dir[getDir(index1, shape)]);
            //v2=_getBodyValue(index1.dir[getBucket(index1.dir, index1.start, shape.XValue, shape.YValue, shape.hValue)]);
            if (v1 != v2)
                if (v1 != 327679 || v2 != 0)
                    std::cout << "unequal " << t_shape.hValue << " " << getDir(index1, shape) << " " << shape.hValue << " " << v2 << std::endl;
        }
    }
    std::cout << s << std::endl << sysTime() - start << std::endl;
    return 0;
}

int main(int argc, char** argv)
{
    if (argc < 3)
        return 1;
    SeqFileIn rFile(toCString(argv[1]));
    SeqFileIn gFile(toCString(argv[2]));
    StringSet<CharString> ids;
    StringSet<DnaString> reads;
    StringSet<DnaString> genome;
    readRecords(ids, reads, rFile);
    readRecords(ids, genome, gFile);
    std::cout << "read done" << std::endl;
    //float alpha = 1.8;

    //mTest(reads, genome)    ;
    //for (unsigned k = 0; k < 20; k+=2)
    //{
    //    alpha = 2 + (float)k / 10;
    //    std::cout << alpha << std::endl;
    uTest(reads, genome, 1.8);
    //    std::cout << std::endl;
    //}
    //mTest1(reads, genome);
    //umTest(reads, genome);
    TShape shape;
    std::cout << h2y(shape, (uint64_t)-1) << std::endl;
}
