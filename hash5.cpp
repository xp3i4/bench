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

const unsigned shapelength = 30;
const unsigned shapeweight = 23;
const unsigned blocklimit = 32;

typedef Iterator<String<Dna> >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength, shapeweight> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;

typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength, shapeweight>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
typedef Index<StringSet<DnaString>, IndexQGram<SimpleMShape, OpenAddressing > > TMIndex;

typedef Value<TShape>::Type HValue;

void _setSANode(uint64_t & node, uint64_t const & k, uint64_t const & j)
{
    node = k << _bitSANum + j;
}

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

template <typename HValue>
inline void _compress(HValue & k, HValue &j, HValue &YValue)
{
   YValue += (k << 36) + (j << 18);
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

    std::vector<uint64_t> hs2;

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
            //hs[m++] = std::Pair(shape.XValue, _compress(k, j, shape.YValue));
        }
    }
    std::cout << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1, (HValue)1);
    std::sort(begin(hs), end(hs),
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) > std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    hs[length(hs) - 1] = std::make_tuple((HValue)1, (HValue)1, (HValue)0, (HValue)1, (HValue)1);
    
    std::cout << "            sort sysTime(): " << sysTime() - time << std::endl;
    HValue countx = 1, counth = 1, tmp = 0, countdh = 0, countb = 1, hk = 0;
    resize(index.sa, length(hs));
    for (HValue k = 1; k < length(hs); k++)
    {
        //std::cout << k << std::endl;
        index.sa[k] = (std::get<0>(hs[k]) << 40)  + std::get<4>(hs[k]);
        if (std::get<1>(hs[k]) != std::get<1>(hs[k - 1]))
        { _setBodyNode(index.dir[index.start + hk], std::get<2>(hs[k-1]), _BodyType_code, tmp);
            //std::cout << counth << std::endl;
            hk++;
            countb++;
            countdh++;
            tmp = counth;
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

    std::cout << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    std::cout << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}

template <unsigned TSpan, unsigned TWeight>
void createQGramIndexDirOnly(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > >& index)
{
    double time = sysTime(); 
    std::cout << "    createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams(index);
    std::cout << "        index.dir[index._Empty_Dir_] = " << index.dir[index._Empty_Dir_] << std::endl << "        length(index.dir) = " << length(index.dir) << std::endl;
    std::cout << "        End createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
}

int mTest1(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape;
    TIndex index(reads);
    uint64_t sum = 0, p = 0;
    double time = sysTime();
    std::cout << "mTest1(): " << "    start " << time << std::endl;
    createQGramIndexDirOnly(index);
    std::cout << "    lenght Dir = " << length(index.dir) << " length Text = " << lengthSum(indexText(index)) << std::endl;
    std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
    for(uint64_t k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
        //for (uint64_t j = 0; j < 100; j++) 
        {
            hashNext(shape, it + j);
            p = getDir(index, shape);
            sum += index.dir[p + 1] - index.dir[p];
        }
    }
    sum=_getBodyCounth(sum);
    std::cout << "    sum = " << sum << std::endl;
    std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    End mTest1()" << std::endl;

    return 0;
}

int mTest2(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape shape;
    TIndex index(reads);
    
    uint64_t sum = 0;
    double time = sysTime();
    std::cout << "mTest2(): " << std::endl;
    for (unsigned j = 24; j < 31; j++)
    for (unsigned k = 22; k < 23; k++)
    {
        std::cout << j << " " << k << std::endl;
        //resize(shape, 30, k);
        //resize(index.shape, 30, k);
        shape.span=j;
        shape.weight=k;
        index.shape.weight=k;
        index.shape.span=j;
        createQGramIndexDirOnly(index);
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
    } 

    return 0;
}

int uTest(StringSet<DnaString> & reads, StringSet<DnaString> & genome)
{
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t sum=0, count=0, p = 0;
    double time = sysTime();
    std::cout << "uTest():\n" << "    start " << time << std::endl;
    std::cout << "    fullDirLength " << _fullDirLength(index) << std::endl; 

    indexCreate(index, FibreDir());
    std::cout << "    getBucket start sysTime(): " << sysTime() - time<< std::endl;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        //for (uint64_t j = 0; j < 100; j++)
        {
            hashNext(t_shape, it + j);
            p = getBucket(index.bucketMap, t_shape.hValue);
            sum += index.dir[p + 1] - index.dir[p];
        }
    }
    std::cout << "    sum = " << sum << " count = " << count << std::endl;
    std::cout << "    getBucket end sysTime(): "<< sysTime() - time<< std::endl;
    std::cout << "    End uTest()" << std::endl;
    return 0;
}

int umTest(StringSet<DnaString> & reads,  StringSet<DnaString> & genome)
{
    
    std::cout << "umTest()" << std::endl;
    TShape_u t_shape;
    TIndex_u index(reads);
    unsigned kmerLength = t_shape.span;
    uint64_t s=0, sum=0;
   
    TShape shape;
    TIndex index1(reads);
    
    indexCreate(index, FibreDir());
    createQGramIndexDirOnly(index1);

    std::cout << "    h2y(shape, BucketMap<uint64_t>::EMPTY) = " << h2y(shape, BucketMap<uint64_t>::EMPTY) << std::endl;
    uint64_t v1, v2, p1, p2;
    for(unsigned k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(t_shape, it);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - kmerLength + 1; j++)
        {
            hashNext(t_shape, it + j);
            hashNext(shape, it + j);
            p1 = getBucket(index.bucketMap, t_shape.hValue);
            p2 = getDir(index1, shape);
            if (index.dir[p1+1] - index.dir[p1] != _getBodyCounth(index1.dir[p2 + 1]) - _getBodyCounth(index1.dir[p2]))
                std::cout << index.dir[p1+1] - index.dir[p1] << " " << _getBodyCounth(index1.dir[p2 + 1]) - _getBodyCounth(index1.dir[p2]) << " " << shape.hValue<< std::endl;
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


uint64_t hashBit(uint64_t key) {
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}
/*
inline uint64_t hashInit1(TShape & shape, TIter it)
{
    hashInit(shape, it);
    return shape.hValue;
}

inline uint64_t hashNext(TShape & shape, TIter it)
{
    sha
}
*/

uint64_t _min(TShape me, uint64_t hValue)
{
    uint64_t v1;
    me.XValue = 1<<63;
    for (unsigned k = 64-(me.span << 1) ; k <= 64 - (me.weight << 1); k+=2)
    {   
        v1 = hValue << k >> (64-(me.weight<<1));
        if(me.XValue > v1) 
        {   
            me.XValue=v1;
        }   
    }   
    return me.XValue;
}

int hTest(StringSet<DnaString> & reads)
{
    std::cout << "hTest() " << std::endl;
    TShape shape;
    uint64_t sum = 0;
    double time = sysTime();
    for (uint64_t j = 0; j < length(reads); j++)   
    { 
        TIter it = begin(reads[j]);
        hashInit(shape, it);
        for (uint64_t k = 0; k < length(reads[j]); k++)
        {
            hashNext(shape, it + k);
            //sum^=_min(shape, hashBit(shape.hValue));
            std::cout << std::bitset<64> (_min(shape, hashBit(shape.hValue))) << " " << std::bitset<64> (hashBit(shape.hValue)) << std::endl;
        }
    }
    
    std::cout << "    hashNext() " << sysTime() - time << "\n    " << sum << std::endl;
}

void _sort1(String<Pair<uint64_t, uint64_t> > & arr, String<Pair<uint64_t, uint64_t> > &hs, uint64_t const & p_bit, uint64_t const & l)
{
    uint64_t l_move = 64, r_move = 64 - p_bit, count[1<<p_bit];
    //String<Pair<uint64_t, uint64_t>* > output, tmp;
    String<uint32_t> out, tmp;
    //String<Pair<uint64_t, uint64_t> > tmp1 = arr;
    resize(out, length(arr));
    resize(tmp, length(arr));
    for (uint64_t k = 0; k < length(tmp); k++)
        tmp[k] = k;
    for (uint64_t j = 0; j < l; j++) 
    {  
        l_move -= p_bit; 
        for (uint64_t k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (uint64_t k = 0; k < length(tmp); k++) 
        {
            count[arr[tmp[k]].i1 << l_move >> r_move]++;
        }
        for (uint64_t k = 1; k < (1 << p_bit); k++)
        {
            count[k] += count[k - 1];  
        } 
        for (int64_t k = length(out)-1; k >=0; k-- )
        {
            out[--count[arr[tmp[k]].i1 << l_move >> r_move]] = tmp[k];
            //count[arr[k] << l_move >> r_move]--;
        }
        for(uint64_t k = 0; k < length(out); k++)
        tmp[k] = out[k];
    }
    for (uint64_t k = 0; k < length(tmp); k++)
        hs[k] = arr[out[k]];
}

void _sort(String<Pair<uint64_t, uint64_t> > & arr, uint64_t const & p_bit, uint64_t const & l)
{
    uint64_t l_move = 64, r_move = 64 - p_bit, count[1<<p_bit];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, length(arr));
    for (uint64_t j = 0; j < l; j++) 
    {  
        l_move -= p_bit; 
        for (uint64_t k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (uint64_t k = 0; k < length(arr); k++) 
            count[arr[k].i1 << l_move >> r_move]++;
        for (uint64_t k = 1; k < (1 << p_bit); k++)
        {
            count[k] += count[k - 1];  
        } 
        for (int64_t k = length(arr)-1; k >=0; k-- )
        {
            output[--count[arr[k].i1 << l_move >> r_move]] = arr[k];
            //count[arr[k] << l_move >> r_move]--;
        }
        
        arr = output;
    }
}

void _sort2(String<Pair<uint64_t, uint64_t> > & arr, uint64_t const & p_bit, uint64_t const & l)
{
    uint64_t l_move = 64, r_move = 64 - p_bit, count[1<<p_bit];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, length(arr));
    for (uint64_t j = 0; j < l; j++) 
    {  
        l_move -= p_bit; 
        for (uint64_t k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (uint64_t k = 0; k < length(arr); k++) 
            count[arr[k].i2 << l_move >> r_move]++;
        for (uint64_t k = 1; k < (1 << p_bit); k++)
        {
            count[k] += count[k - 1];  
        } 
        for (int64_t k = length(arr)-1; k >=0; k-- )
        {
            output[--count[arr[k].i2 << l_move >> r_move]] = arr[k];
            //count[arr[k] << l_move >> r_move]--;
        }
        
        arr = output;
    }
}


void _sort(String<uint64_t> & arr, uint64_t const & p_bit, uint64_t const & l)
{
    uint64_t l_move = 64, r_move = 64 - p_bit, count[1<<p_bit];
    String<uint64_t> output;
    resize(output, length(arr));
    for (uint64_t j = 0; j < l; j++) 
    {  
        l_move -= p_bit; 
        for (uint64_t k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (uint64_t k = 0; k < length(arr); k++) 
            count[arr[k] << l_move >> r_move]++;
        for (uint64_t k = 1; k < (1 << p_bit); k++)
        {
            count[k] += count[k - 1];  
        } 
        for (int64_t k = length(arr)-1; k >=0; k-- )
        {
            output[--count[arr[k] << l_move >> r_move]] = arr[k];
            //count[arr[k] << l_move >> r_move]--;
        }
        
        arr = output;
    }
}

void radix_sort(uint64_t *begin, uint64_t *end)
{
    uint64_t *begin1 = new uint64_t[end - begin];
    uint64_t *end1 = begin1 + (end - begin);
    for (uint64_t shift = 0; shift < 32; shift += 8) {
        size_t count[0x100] = {};
        for (uint64_t *p = begin; p != end; p++)
            count[(*p >> shift) & 0xFF]++;
        uint64_t *bucket[0x100], *q = begin1;
        for (uint64_t i = 0; i < 0x100; q += count[i++])
            bucket[i] = q;
        for (uint64_t *p = begin; p != end; p++)
            *bucket[(*p >> shift) & 0xFF]++ = *p;
        std::swap(begin, begin1);
        std::swap(end, end1);
    }
    delete[] begin1;
}

int sTest(StringSet<DnaString> & reads)
{
    typedef typename Value<TShape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue, HValue, HValue> HTuple;
    typedef String<HTuple> Stringtuple;

    TShape shape;
    TIndex index(reads);
    Stringtuple hs, hs1;
    HValue  m = 0, sum = 0;

    String<uint64_t> hs2;
    resize(hs2, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);
    
    String<Pair<uint64_t, uint64_t> > hs3;
    resize(hs3, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);


    double time = sysTime();
    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    std::cout << " sTest() " << sysTime() - time << std::endl;
    std::cout << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;
    for(HValue k = 0; k < length(indexText(index)); k++)
    {
        TIter it = begin(indexText(index)[k]);
        hashInit(shape, it);
        for (HValue j = 0; j < length(indexText(index)[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            hs[m++] = std::make_tuple(shape.XValue, shape.hValue, shape.YValue, k, j);
            hs2[m - 1] = shape.XValue;
            hs3[m - 1].i1 = shape.XValue;
            hs3[m - 1].i2 = shape.YValue;
        }
    }

    // another little test  
    std::sort(begin(hs3), end(hs3), [](Pair<uint64_t, uint64_t> const & a,  Pair<uint64_t, uint64_t> const & b){return a.i1<b.i1;});

    //std::sort(begin(hs2), end(hs2), [](uint64_t & a, uint64_t&b){return a < b;});
    time = sysTime();
    uint64_t step = 0, first = hs3[0].i1;
    for (uint64_t k = 0 ; k < length(hs3); k++)
    {
        if (first != hs3[k].i1)
        {
            if (step > 2)
                std::sort(begin(hs3)+k - step , begin(hs3)+k, [](Pair<uint64_t, uint64_t> const & a, Pair<uint64_t, uint64_t> const & b){return a.i2<b.i2;});
            step=1;
            first = hs3[k].i1;
        }
        else 
            step++;
        //for (uint64_t  n=0; n < 32; n++)
        //    std::cout << hs2[k + n] << " " << n << std::endl;
    }
    //for (uint64_t k = 0; k < length(hs2); k++)
    //{
    //    if (hs2[k] > hs2[k+1] && k % step != step-1)
    //        std::cout << hs2[k] << " " << hs2[k+1] << " " << k % step << std::endl;
    //}
    std::cout << "hs2" << hs3[10] << " " << hs3[11] << " " << hs3[12] << std::endl;
    std::cout << sysTime() - time << std::endl;
    // little test end
    
    unsigned count = 0;
    for(uint64_t j = 0; j < 100000; j++ )
    {
        unsigned flag = 0;
        for (uint64_t k = j + 1; k < 100000; k++)
        {
            if (hs2[k] != hs2[j])
            {
                flag = 1;
                continue;
            }
            if (flag && hs2[k] == hs2[j])
                count++;
        }
    } 
    String<Pair<uint64_t, uint64_t> > hs22;
    resize(hs22, length(hs2));
    std::cout << sysTime() - time << std::endl;
    unsigned j = 0;
    for (uint64_t k = 0; k < length(hs2)-1; k++)
    {
        if(hs2[k] != hs2[k+1])
        {
            hs22[j++].i1 = hs2[k];
            hs22[j++].i2 = k;
        }
    }
    resize(hs22, j);
    _sort(hs22, 10, 5);
    std::cout << j << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1, (HValue)1);
    String<uint64_t> tmp = hs2;
    String<Pair<uint64_t, uint64_t> > tmp1=hs3, hs4=tmp1;
    time = sysTime();
    std::cout << "std::sort() start " << sysTime() - time  << std::endl;
    std::sort (begin(hs3), end(hs3), [](Pair<uint64_t, uint64_t> const & a, Pair<uint64_t, uint64_t> const & b){return a.i1 > b.i1;});
    //radix_sort(&hs2[0], &hs2[length(hs2)-1]);
    //for (uint64_t k = 0; k < length(hs2); k++)
    //std::cout << hs2[k] << " std::sort() end " << sysTime() - time<< std::endl;
    //resize(tmp, length(hs2));
    //for (uint64_t k = 4; k <= shape.weight; k++)
    for (uint64_t k = 2; k <= 20; k++)
    {
        hs3 = tmp1;
        uint64_t l = 2*shape.weight/k;
        if (l * k != 2*shape.weight)
            l++;
        uint64_t l1 = 20/k;
        if (l1*k != 20)  
            l1++;
        std::cout << k << " " << l1 << " \n"; 
        time = sysTime();
        //for (uint64_t k = 0; k < length(hs2); k++)
        //    sum += hs2[k]; 
        //String<uint64_t> t;
        //resize(t, 10) ;
        //t[0] = 3;
        //t[1] = 5;
        //t[2] = 2;
        //t[3] = 4;
        //t[4] = 8;
        //t[5] = 7;
        //t[6] = 6;
        //t[7] = 9;
        //t[8] = 0;
        //t[9] = 1;
        _sort2(hs3, 10, 2); 
        _sort(hs3, 10, 5);
        //for (uint64_t k = 0; k<length(t); k++)
        //std::cout << t[k] << " ";
            //std::cout << hs2[k] << std::endl;
        std::cout << "            sort sysTime(): " << hs3[2] << " " << sysTime() - time << std::endl;
    } 

    return 0;
}

inline void _countSort(String<uint64_t> & arr, String<uint64_t> & output, uint64_t const & l_move, uint64_t const & r_move, uint64_t const & p_bit, uint64_t const & l)
{
    uint64_t count[1 << p_bit] = {0};
    for (uint64_t k = 0; k < length(arr); k++) 
        count[arr[k] << l_move >> r_move]++;
    for (uint64_t k = 1; k < (1 << p_bit); k++)
    {
        count[k] += count[k - 1];  
    } 
    for (int64_t k = length(arr)-1; k >=0; k-- )
    {
        output[count[arr[k] << l_move >> r_move] - 1] = arr[k];
        count[arr[k] << l_move >> r_move]--;
    }

    arr = output;
}

//void _sort(String<uint64_t> & arr, uint64_t const & p_bit, uint64_t const & l)
//{
//    uint64_t l_move = 64, r_move = 64 - p_bit, count[1<<p_bit] = {0}, sum[1<<p_bit + 1] = {0};
//    String<uint64_t> output;
//    resize(output, length(arr));
//    for (uint64_t j = 0; j < l; j++) 
//    {  
//        l_move -= p_bit; 
//        for (uint64_t k = 0; end < 1<<p_bit; k++) 
//        {
//            for (uint64_t n = begin; n < end; n++ )
//                count[arr[n] << l_move >> r_move]++;
//            for (uint64_t n = 1; n < (1 << p_bit); n++)
//            {
//                sum[n] += count[n - 1];  
//            } 
//            for (int64_t n = begin; n < end; n++ )
//            {
//                output[count[arr[n] << l_move >> r_move] + sum[n]] = arr[n];
//                count[arr[n] << l_move >> r_move]--;
//            }
//            if (sum[k] ^ sum[k+1])
//            { 
//                begin = sum[k];
//                end = sum[k + 1];
//            }
//            else 
//                continue;
//        } 
//        arr = output;
//    }
//}



//void _sort(String<uint64_t> & arr, uint64_t const & p_bit, uint64_t const & l)
//{
//    uint64_t l_move = 64, r_move = 64 - p_bit;
//    String<uint64_t> output;
//    std::cout << "radixSort() " << std::endl;
//    resize(output, length(arr));
//    for (uint64_t k = 0; k < l; k++) 
//    {  
//        l_move -= p_bit; 
//        _countSort(arr, output, l_move, r_move, p_bit, l);
//    }
//}
 
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
    //uTest(reads, genome);
    mTest1(reads, genome);
    //mTest2(reads, genome);
    //hTest(reads);
    //sTest(reads);
    return 0;
}
