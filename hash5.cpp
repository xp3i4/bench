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

const unsigned shapelength = 25;
const unsigned shapeweight = 18;
const unsigned blocklimit = 32;

typedef Iterator<String<Dna> >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength, shapeweight> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;

typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength, shapeweight>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
typedef Index<StringSet<DnaString>, IndexQGram<SimpleMShape, OpenAddressing > > TMIndex;

typedef Value<TShape>::Type HValue;
/*
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

void _sort2(String<Pair<uint64_t, uint64_t> > & arr, unsigned const & p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    uint64_t count[1<<p_bit];
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
            count[k] += count[k - 1];  
        for (int64_t k = length(arr)-1; k >=0; k-- )
            output[--count[arr[k].i1 << l_move >> r_move]] = arr[k];
            //count[arr[k] << l_move >> r_move]--;
        
        arr = output;
    }
}
*/

template <typename TIt>
inline void _sort3(TIt const & begin, const TIt & end, unsigned const & p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    //uint64_t count[1<<p_bit];
    uint64_t count[1024];
    //int count[1024];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, end - begin);
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i1 << l_move >> r_move]++;
        for (int k = 1; k < (1 << p_bit); k++)
            count[k] += count[k - 1];
        for (int64_t k = end - begin - 1; k >=0; k-- )
            output[--count[(begin + k)->i1 << l_move >> r_move]] = *(begin + k);
        for (int64_t k = 0; k < end - begin; k++)
            *(begin + k)  = output[k];
    }
}

/*
template <typename TIt>
inline void _sort3I2(TIt const & begin, const TIt & end, unsigned const & p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    uint64_t count[1<<p_bit];
    //String<Pair<uint64_t, uint64_t> > output;
    Pair<uint64_t, uint64_t> output[end - begin];
    //resize(output, end - begin);
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (uint64_t k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (uint64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i2 << l_move >> r_move]++;
        for (uint64_t k = 1; k < (1 << p_bit); k++)
            count[k] += count[k - 1];
        for (int64_t k = end - begin -1; k >=0; k-- )
            output[--count[(begin + k)->i2 << l_move >> r_move]] = *(begin + k);
        for (uint64_t k = 0; k < end - begin; k++)
            *(begin + k)  = output[k];
    }
}

template <typename TIt>
inline void _sort3I22(TIt const & begin, const TIt & end, unsigned const & p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    uint64_t count[1<<p_bit];
    //String<Pair<uint64_t, uint64_t> > output;
    Pair<uint64_t, uint64_t> output[end - begin];
    //resize(output, end - begin);
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (unsigned k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (uint64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i2 << l_move >> r_move]++;
        for (unsigned k = (1 << p_bit) - 2; k >=0;  k--)
            count[k] += count[k + 1];
        for (int64_t k = 0; k < end - begin; k++)
            output[--count[(begin + k)->i2 << l_move >> r_move]] = *(begin + k);
        for (uint64_t k = 0; k < end - begin; k++)
            *(begin + k)  = output[k];
    }
}
*/

template <typename TIt>
inline void _insertSort(TIt const & begin, TIt const & end )
{
    Pair<uint64_t, uint64_t> key;
    for (int j = 1; j < end - begin; j++)
    {
        key = *(begin + j);
        int k = j - 1;
        while (k >= 0)
        {
            if (((begin + k)->i2 < key.i2))
                *(begin+k+1) = *(begin+k);
            else
            {
                break;
            }
            k--;
        }    
        *(begin+k+1) = key;
    }
}

template <typename TIt>
inline void _insertSort2(TIt const & begin, TIt const & end )
{
    Pair<uint64_t, uint64_t> key;
    for (int j = 1; j < end - begin; j++)
    {
        key = *(begin + j);
        int k = j - 1;
        while (k >= 0)
        {
            if (((begin + k)->i2 < key.i2))
                *(begin+k+1) = *(begin+k);
            else
            {
                break;
            }
            k--;
        }    
        *(begin+k+1) = key;
    }
}

/*
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
*/

//-----------------------------------------------------------------------------
// test if sort functions are correct
//-----------------------------------------------------------------------------
int sTest4(StringSet<DnaString> & reads)
{
    typedef String<Pair<uint64_t, uint64_t>> StringPair;

    TShape shape;    
    StringPair hs, hs_std, tmp;
    uint64_t m = 0;
    std::cout << " sTest4()" << std::endl;

    resize(hs, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1));
    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            hs[m].i1 = shape.XValue;
            hs[m++].i2 = shape.YValue;
        }
    }
    tmp = hs;

    hs = tmp; 
    _insertSort(begin(hs), end(hs));
    
    hs_std = tmp;
    std::sort(begin(hs_std), end(hs_std), [](Pair<uint64_t, uint64_t> & a, Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
    
    std::cout << "done " << std::endl;
    for (uint64_t n = 0; n < length(hs) - 1; n++)
    {
        if (hs[n].i2 != hs_std[n].i2)
        {
            std::cout << "    Error " << hs[n] << " " << hs_std[n]<< std::endl;
        }
    }
    return 0;
}
/*
int sTest1(StringSet<DnaString> & reads)
{
    typedef String<Pair<uint64_t, uint64_t>> StringPair;

    TShape shape;    
    StringPair hs, tmp;
    uint64_t m = 0;
    std::cout << " sTest1()" << std::endl;

    resize(hs, lengthSum(reads) - length(reads) * (shape.span - 1) + 1);
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1) + 1);
    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            hs[m].i1 = shape.XValue;
            hs[m++].i2 = shape.YValue;
        }
    }
    tmp = hs;
    double time = sysTime();
    for (uint64_t k = 1; k <= shape.weight * 2; k++)
    {
        uint64_t l = 2*shape.weight/k;
        if (l * k != 2*shape.weight)
            l++;
        hs = tmp; 
        time = sysTime();
        _sort2(hs, k, l);
        std::cout << "    End sort sysTime(): " << k << " " << l << " " << hs[2] << " " << sysTime() - time << std::endl;
        for (uint64_t n = 0; n < length(hs) - 1; n++)
            if (hs[n+1].i1 < hs [n].i1)
            {
                std::cout << "    Error " << hs[n] << std::endl;
            }
    }
}
*/

void _createValueArray(StringSet<DnaString> & reads, String<Pair<uint64_t, uint64_t> > & hs, unsigned & step, unsigned & l)
{
    typedef String<Pair<uint64_t, uint64_t>> StringPair;

    TShape shape;    
    StringPair tmp;
    String<uint64_t> tmp3;
    uint64_t p = 0, q=0, c = 0, n = -1, pre = ~0, count = 0, mask = ((uint64_t)1 << 63);
    std::cout << "    _createValueArray() " << std::endl;
    double time = sysTime();
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(tmp3, lengthSum(reads) - length(reads) * (shape.span - 1));
    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            tmp3[p] = shape.YValue;
            if (pre ^ shape.XValue)
            {
                tmp[q].i1 = shape.XValue;
                tmp[q].i2 = p;
                pre = shape.XValue;
                tmp3[p] |= mask;
                q++;
            }
            p++;            
        }
    }
    resize(tmp, q);
    *(end(tmp3)) |= (mask);

    std::cout << "        loading time " << sysTime() - time << std::endl;
    c = tmp[0].i2;  
    p = q = count = 0;
    n = -1;
    _sort3(begin(tmp), end(tmp), step, l);         // sort parameters 1
    std::cout << "        sort xvalue time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
    for (uint64_t q = 0;  q < length(hs) - 1; q++)
    {
        if (tmp3[c] & mask)
        {
            c = tmp[++n].i2;
        }
        hs[q].i1 = tmp[n].i1;
        hs[q].i2 = tmp3[c] & (~mask);
        c++;
    }
    std::cout << "        xvalue expand " << sysTime() - time << std::endl;
    hs[length(hs)-1].i2 |= mask;
    pre = hs[0].i1;
    for (uint64_t k = 0; k < length(hs); k++)
    {
        if (pre ^ hs[k].i1)
        {   
            pre = hs[k].i1;
            if (count < 20)                   // sort parameters
                _insertSort(begin(hs) + k - count, begin(hs) + k);
            else 
                std::sort(begin(hs) + k -count, begin(hs) + k, [](Pair<uint64_t, uint64_t> & a, 
                Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
            count = 0;
        }
        count++;
   } 


   std::cout << "        End sort sysTime(): " <<  sysTime() - time << std::endl;
}

//--------------------------------------------------------------------
// create sa array at the same time
//--------------------------------------------------------------------
void _createValueArray2(StringSet<DnaString> & reads, String<Pair<uint64_t, uint64_t> > & hs, unsigned & step, unsigned & l)
{
    typedef String<Pair<uint64_t, uint64_t>> StringPair;

    TShape shape;    
    StringPair tmp;
    String<uint64_t> tmp3;
    uint64_t p = 0, q=0, c = 0, n = -1, pre = ~0, count = 0, mask = ((uint64_t)1 << 63);
    std::cout << "    _createValueArray() " << std::endl;
    double time = sysTime();
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(tmp3, lengthSum(reads) - length(reads) * (shape.span - 1));
    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            _setBodyNode(tmp3[p], shape.YValue, (uint64_t)1, _createSANode(j, k));
            if (pre ^ shape.XValue)
            {
                tmp[q].i1 = shape.XValue;
                tmp[q].i2 = p;
                pre = shape.XValue;
                tmp3[p] |= mask;
                q++;
            }
            p++;            
        }
    }
    resize(tmp, q);
    *(end(tmp3)) |= (mask);

    std::cout << "        loading time " << sysTime() - time << std::endl;
    c = tmp[0].i2;  
    p = q = count = 0;
    n = -1;
    _sort3(begin(tmp), end(tmp), step, l);         // sort parameters 1
    std::cout << "        sort xvalue time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
    for (uint64_t q = 0;  q < length(hs) - 1; q++)
    {
        if (tmp3[c] & mask)
        {
            c = tmp[++n].i2;
        }
        hs[q].i1 = tmp[n].i1;
        hs[q].i2 = tmp3[c] & (~mask);
        c++;
    }
    std::cout << "        xvalue expand " << sysTime() - time << std::endl;
    hs[length(hs)-1].i2 |= mask;
    pre = hs[0].i1;
    for (uint64_t k = 0; k < length(hs); k++)
    {
        if (pre ^ hs[k].i1)
        {   
            pre = hs[k].i1;
            if (count < 20)                   // sort parameters
                _insertSort(begin(hs) + k - count, begin(hs) + k);
            else 
                std::sort(begin(hs) + k -count, begin(hs) + k, [](Pair<uint64_t, uint64_t> & a, 
                Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
            count = 0;
        }
        count++;
   } 


   std::cout << "        End sort sysTime(): " <<  sysTime() - time << std::endl;
}

//-----------------------------------------------------------------------------
// test if _creatValueArray() is correct
//-----------------------------------------------------------------------------
void sTest3(StringSet<DnaString> & reads, unsigned & step, unsigned l)
{
    std::cout << "sTest3() " << std::endl; 
    String<Pair<uint64_t, uint64_t> > hs, hs_std;
    TShape shape;
    resize(hs, lengthSum(reads) - length(reads) * (shape.span - 1) + 1);
    resize(hs_std, lengthSum(reads) - length(reads) * (shape.span - 1));
    _createValueArray(reads, hs, step, l);
    for (uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            hs_std[k].i1 = shape.XValue;
            hs_std[k].i2 = shape.YValue;
        }
    }
    std::sort(begin(hs_std), end(hs_std), [](Pair<uint64_t, uint64_t> & a, 
    Pair<uint64_t, uint64_t> & b){return a.i1 < b.i1 || (a.i1 == b.i1 && a.i2 > b.i2);}); 

    for (uint64_t k = 0; k < length(hs) - 1; k++)
    {
        if (hs[k] != hs_std[k])
        {
            std::cout << hs[k] << " " << hs_std[k] << std::endl;
        }
    }
    std::cout << "    End sTest3 " << std::endl;

}

void sTest3_(StringSet<DnaString> & reads, unsigned & step, unsigned l)
{
    std::cout << "sTest3_() " << std::endl; 
    String<Pair<uint64_t, uint64_t> > hs, hs_std;
    TShape shape;
    resize(hs, lengthSum(reads) - length(reads) * (shape.span - 1) + 1);
    resize(hs_std, lengthSum(reads) - length(reads) * (shape.span - 1));
    _createValueArray2(reads, hs, step, l);
    for (uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            hs_std[k].i1 = shape.XValue;
            hs_std[k].i2 = shape.YValue;
        }
    }
    std::sort(begin(hs_std), end(hs_std), [](Pair<uint64_t, uint64_t> & a, 
    Pair<uint64_t, uint64_t> & b){return a.i1 < b.i1 || (a.i1 == b.i1 && a.i2 > b.i2);}); 

    for (uint64_t k = 0; k < length(hs) - 1; k++)
    {
        if (hs[k].i1 != hs_std[k].i1 || _getBodyValue(hs[k].i2) != hs_std[k].i2)
        {
            std::cout << hs[k] << " " << hs_std[k] << std::endl;
        }
    }
    std::cout << "    End sTest3_ " << std::endl;

}

//-----------------------------------------------------------------------------
// prototype of _createValueArray()
//-----------------------------------------------------------------------------

int sTest2(StringSet<DnaString> & reads)
{
    typedef String<Pair<uint64_t, uint64_t> > StringPair;

    TShape shape;    
    StringPair hs, tmp, tmp1, hs_std;
    String<uint64_t> tmp3;
    uint64_t p = 0, q=0, c = 0, n = 0, pre = ~0, count = 0, mask = (uint64_t)1 << 63;
    std::cout << " sTest2()" << std::endl;

    resize(hs, lengthSum(reads) - length(reads) * (shape.span - 1)+1);
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(tmp3, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(hs_std, lengthSum(reads) - length(reads) * (shape.span - 1));

    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            tmp3[p] = shape.YValue;
            //std::cout << shape.XValue << " " << shape.YValue << std::endl;
            hs_std[p].i1 = shape.XValue;
            hs_std[p].i2 = shape.YValue;
            if (pre ^ shape.XValue)
            {
                tmp[q].i1 = shape.XValue;
                tmp[q].i2 = p;
                pre = shape.XValue;
                tmp3[p] |= mask;
                q++;
            }
            p++;            
        }
    }
    *end(tmp3) |= mask;
    //tmp3 = hs;
    //hs = hs_std;
    double time = sysTime();
    std::sort(begin(hs_std), end(hs_std), [](Pair<uint64_t, uint64_t> & a, Pair<uint64_t, uint64_t> & b){return a.i1 < b.i1 ||(a.i1==b.i1 && a.i2> b.i2);});
    std::cout << sysTime() - time << std::endl;
    resize(tmp, q);
    resize(tmp1, q);
    tmp1 = tmp;
    time = sysTime();
    for (uint64_t k = 8; k < 9; k++)
    //for (uint64_t k = 5; k <= shape.weight * 2; k++)
    {
       // uint64_t l = 2*shape.weight/k;
        //if (l * k != 2*shape.weight)
        //    l++;
        //tmp = tmp1; 
        p = 0;
        q = 0;
        n = -1;
        count = 0;
        time = sysTime();
        //_sort2(tmp, k, l);
        //_sort2(hs, k, l);
        _sort3(begin(tmp), end(tmp), 8, 6);
        //std::sort(begin(tmp), end(tmp), [](Pair<uint64_t, uint64_t> & a, Pair<uint64_t, uint64_t> & b){return a.i1 > b.i1;});
        c = tmp[0].i2;
        //std::cout << "& " << (tmp3[c] & mask) << " " << std::bitset<64>(tmp3[c]) << " " << std::bitset<64>(~mask)<< std::endl;

        for (uint64_t q = 0;  q < length(hs) - 1; q++)
        {
            if (tmp3[c] & mask)
            {
                c = tmp[++n].i2;
            }
            hs[q].i1 = tmp[n].i1;
            hs[q].i2 = tmp3[c] & (~mask);
            //std::cout << n << " " << tmp3[c] << " " << c << std::endl;
            c++;
        }
        hs[length(hs)-1].i2 |= mask;
        pre = hs[0].i1;
        for (uint64_t k = 0; k < length(hs); k++)
        {
            if (k < 100)
                std::cout <<pre << " " << hs[k].i1 << count << std::endl;
            if (pre ^ hs[k].i1)
            {   
                pre = hs[k].i1;
                if (count < 20)
                    _insertSort(begin(hs) + k - count, begin(hs) + k);
                else 
                    //_sort3I22(begin(hs) + k - count, begin(hs) + k, 6, 3);
                //if (count > 2)
                    std::sort(begin(hs) + k -count, begin(hs) + k, [](Pair<uint64_t, uint64_t> & a, Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
                count = 0;
            }
            count++;
        }
        


        std::cout << "    End sort sysTime(): " <<  sysTime() - time << std::endl;

        for (uint64_t n = 0; n < length(hs) - 1; n++)
        { 
            if (hs[n] != hs_std[n])
            {
                std::cout << n << " " <<  hs[n] << " " << hs_std[n] << std::endl;
            }
            //std::cout << hs[n]<< " " << hs_std[n]<< std::endl;
            //if (hs[n+1].i1 < hs[n].i1)
            //{
            //    
            //    std::cout << "    Error" << hs[n - 1] << " " << hs[n] << " " << hs[n + 1] << std::endl;
            //    break;
            //    //return 1;
            //}
        }
    }
    return 0;
}
/*
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
*/

/*
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
*/

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
    unsigned step = atof(toCString(argv[3]));
    unsigned l = atof(toCString(argv[4]));
    //std::cout << "read done sysTime " << sysTime() - time << std::endl;
    //umTest(reads, genome);
    //uTest(reads, genome);
    //mTest1(reads, genome);
    //mTest2(reads, genome);
    //hTest(reads);
    //sTest(reads);
    //sTest1(reads);
    //sTest2(reads);
    //while (true)
    //{
        std::cout << " step = " << step << " l = " << l << std::endl;
        //sTest3(reads, step, l);
        sTest3_(reads, step, l);
    //}
    //sTest4(reads);
    return 0;
}
