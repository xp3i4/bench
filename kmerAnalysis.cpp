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

const unsigned shapelength =30;
const unsigned shapeweight = 22;
const unsigned blocklimit = 32;


typedef Iterator<String<Dna> >::Type TIter;
typedef Shape<Dna, MinimizerShape<shapelength, shapeweight> > TShape;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;
typedef Shape<Dna, SimpleShape> TSpShape;

typedef Index<StringSet<DnaString>, IndexQGram<MinimizerShape<shapelength, shapeweight>, OpenAddressing > > TIndex;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
typedef Index<StringSet<DnaString>, IndexQGram<SimpleMShape, OpenAddressing > > TMIndex;
//typedef Index<StringSet<DnaString>, IndexQGram<SimpleSpShape, OpenAddressing > > TSpIndex;


typedef Value<TShape>::Type HValue;

static uint64_t const _Anchor_Empty_ = 1 << 63;
static uint64_t const _AnchorBit_i1_ = 10;
static uint64_t const _AnchorMask_i1_ = (1 <<  _AnchorBit_i1_) - 1;

inline void _initAnchor(String<uint64_t> & anchor)
{
    for (uint64_t k = 0; k < length(anchor); k++)
        anchor[k] = _Anchor_Empty_;
}

inline int _initAnchor(String<uint64_t> & anchor, uint64_t const & mask)
{
    for (uint64_t k = 0; k < mask + 1; k++)
        anchor[k] = _Anchor_Empty_;
    return 0;
}

inline int _initAnchor(uint64_t *anchor, uint64_t const & mask)
{
    for (uint64_t k = 0; k < mask + 1; k++)
        anchor[k] = _Anchor_Empty_;
    return 0;
}

inline uint64_t _getAnchor_i1 (uint64_t & node)
{
    return node >> _AnchorBit_i1_;
}

inline uint64_t _createAnchorNode (uint64_t & code)
{
    return code << _AnchorBit_i1_;
}

inline uint64_t _getAnchor_i2 (uint64_t & node)
{
    return node & _AnchorMask_i1_;
}

inline void _updateValue(uint64_t & value, uint64_t code)
{
    
}

  inline uint64_t   _h1(uint64_t val)
     {
     uint64_t key = val;
           key = (~key) + (key << 21); // key = (key << 21) - key - 1;
   key = key ^ (key >> 24);
   key = (key + (key << 3)) + (key << 8); // key * 265
   key = key ^ (key >> 14);
   key = (key + (key << 2)) + (key << 4); // key * 22
   key = key ^ (key >> 28);
   key = key + (key << 31);
   return key;
     }

inline uint64_t _requestAnchor(String<uint64_t> & anchor, uint64_t const & code, uint64_t const & mask /*length of anchor - 1*/)
{
    uint64_t k = 0, x = code & mask; 
    while((_getAnchor_i1(anchor[x]) ^code) & (_getAnchor_i1(anchor[x]) ^ _Anchor_Empty_))
    {
        x = (x + ++k) & mask;
    }
    anchor[x] = (code << _AnchorBit_i1_) + ((anchor[x] + 1 ) & _AnchorMask_i1_); 
    return x;
}

inline uint64_t _requestAnchor(uint64_t * anchor, uint64_t const & code, uint64_t const & mask /*length of anchor - 1*/)
{
    uint64_t k = 0, x = code & mask; 
    while((_getAnchor_i1(anchor[x]) ^code) & (_getAnchor_i1(anchor[x]) ^ _Anchor_Empty_))
    {
        x = (x + ++k) & mask;
    }
    anchor[x] = (code << _AnchorBit_i1_) + ((anchor[x] + 1 ) & _AnchorMask_i1_); 
    return x;
//    while (true)
//    {
//        if (_getAnchor_i1(anchor[x]) ==  code)
//        {
//            anchor[x]++;
//            return x;
//        }
//        else if (_getAnchor_i1(anchor[x]) == _Anchor_Empty_)
//        {
//            anchor[x] = (code << _AnchorBit_i1_) + 1; 
//            return x;
//        }
//        x = (x + ++k) & mask;
//    }
}

inline uint64_t _requestAnchor2(uint64_t * anchor, uint64_t const & code, uint64_t const & mask /*length of anchor - 1*/)
{
    uint64_t k = 0, x = code & mask; 
    while((_getAnchor_i1(anchor[x]) ^code) & (_getAnchor_i1(anchor[x]) ^ _Anchor_Empty_))
    {
        x = (x + ++k) & mask;
    }
    anchor[x] = (code << _AnchorBit_i1_) + ((anchor[x] + 1 ) & _AnchorMask_i1_); 
    return x;
}


inline uint64_t _getAnchor(String<uint64_t> & anchor, uint64_t const & code, uint64_t const & mask /*length of anchor - 1*/)
{
    uint64_t k = 0, x = code & mask;
    while ((_getAnchor_i1(anchor[x]) ^ code) & (_getAnchor_i1(anchor[x]) ^ _Anchor_Empty_))
    {
        x = (x + ++k) & mask;
    }
    return x;
}

inline uint64_t _getAnchor(uint64_t * anchor, uint64_t const & code, uint64_t const & mask /*length of anchor - 1*/)
{
    uint64_t k = 0, x = code & mask;
    while ((_getAnchor_i1(anchor[x]) ^ code) & (_getAnchor_i1(anchor[x]) ^ _Anchor_Empty_))
    {
        x = (x + ++k) & mask;
    }
    return x;
}



void opKmerAnalysis2(StringSet<DnaString> & genome, StringSet<DnaString> & reads)
{
    
    TShape shape; 
    TIndex_u index(genome);
    uint64_t mask  = 131072 - 1, dn, count=0, sum = 0;
    uint64_t anchor[mask + 1];
    double time = sysTime();
    std::cout << "mnKmerAnalysis() \n";

    //createQGramIndexDirOnly(index);
    indexCreate(index, FibreSADir());
    _initAnchor(anchor, mask);
    for (uint64_t j = 0; j < length(reads); j++)  
    //for (uint64_t j = 0; j < 3000; j++)
    //for (uint64_t j = 27; j< 28; j++)
    {
        TIter it = begin(reads[j]);
        hashInit(shape, it);
        _initAnchor(anchor, mask); 
        for (uint64_t k = 0; k < (length(reads[j]) - shape.span + 1); k++) 
        {
            hashNext(shape, it + k);
            dn = getBucket(index.bucketMap, shape.hValue); 
            uint64_t pre = ~0; 
            for (uint64_t n = index.dir[dn]; n < index.dir[dn + 1]; n++)
            {
                if (index.sa[n].i2 - pre > 100)
                {
                    uint64_t x=_requestAnchor(anchor, (index.sa[n].i2 - k) >> 9, mask);
                    pre = index.sa[n].i2;
                //sum ^= (_getSA_i2(index.sa[n]) - k) >> 9;
                //if ((index.sa[n].i2 - k) >>9==298768 && j==13)
                //if(j==2 && index.dir[dn+1]-index.dir[dn] < 32)
                if(j== 9)
                    std::cout << shape.hValue << " " << k << " " << index.sa[n].i2<< std::endl;
                
                }
            }
        }
        uint64_t max1 = 0, max2 = 0, k1 = 0, k2 =0, countk = 0;
        //for (uint64_t k = 0; k < length(anchor); k++)
        for (uint64_t k = 0; k < mask + 1; k++) 
        {
            if (anchor[k] != _Anchor_Empty_)
            {
                if (max1 < (anchor[k] & _AnchorMask_i1_))
                {
                    max1 = _getAnchor_i2(anchor[k]);
                    k1 = _getAnchor_i1(anchor[k]);
                }
                else if(max2 < (anchor[k] & _AnchorMask_i1_))
                {
                    max2 = _getAnchor_i2(anchor[k]);
                    k2 = _getAnchor_i1(anchor[k]);
                }
                //if (j==44)
                //std::cout << (_getAnchor_i1(anchor[k]) << 9) << " " << _getAnchor_i2(anchor[k]) << std::endl;

                //anchor[k] = _Anchor_Empty_;
            }
            //else 
            //    countk++;
        }
        //if (j < 100)
        //std::cout << " " << j << " " << length(reads[j]) << " " << (k1 << 9) << " " << max1 << " " << k2 << " " << max2 << " " << countk / (float)length(anchor) << std::endl;
        //std:: cout << k << " " << (anchor[k] & _AnchorMask_i1_) << std::endl; //<< " " << std::bitset<64>(anchor[k]) << std::endl;
        //std::cout << j << " " << std::endl;
        sum ^= max1;
        if (max1 == 0)
            count++;
        
    }
    std::cout << "    count % = " << (float)count  << " " << length(reads) << " " << (float)count / length(reads) << std::endl;
    std::cout << "    anchor[0] = " << anchor[0] << std::endl;
    std::cout << "    End mnKmerAnalysis() systime() = " << sysTime() - time << std::endl;
}

int opKmerAnalysis(TIndex_u & index, StringSet<DnaString> & genome)
{
    TShape_u shape;
    unsigned step = 6, N=1; // actual step = 2 << step;
    String<float> count;
    String<uint64_t> anchor;
    String<Pair<uint64_t, uint64_t> > anchor1;
    resize(count, 100000000);
    uint64_t mask = ((1<<32)-1);
    uint64_t flag = 0;
    //resize(count, length(index) >> step);
    resize(anchor, lengthSum(index.sa));
    resize(anchor1, length(index.sa));
    for (uint64_t j = 0; j < length(genome); j++) 
    {
        
        std::cout << length(genome[j]) << " ";
         
        for(uint64_t k = 0; k < length(count); k++)
        {
            count[k] = 0;
            //anchor[k] = 0;
        }
        //std::cout << "> " << j + 1 << length(genome[j]);
        
        for (uint64_t k = 0; k < length(anchor); k++)
            anchor[k] = 0;
        for (uint64_t k = 0; k < length(anchor1); k++)
        {
            anchor1[k].i1 = 0;
            anchor1[k].i2 = 0;
        }
        TIter it=begin(genome[j]);
        hashInit(shape, it);
        for (uint64_t k = 0; k < length(genome[j])/N - shape.span + 1; k++)
        {
            flag = 0;
            hashNext(shape, it+k);
            float occ = index.dir[getBucket(index.bucketMap, shape.hValue) + 1] - index.dir[getBucket(index.bucketMap, shape.hValue)];
            for (uint64_t n = index.dir[getBucket(index.bucketMap, shape.hValue)]; n < index.dir[getBucket(index.bucketMap, shape.hValue) + 1] ; n++)
            {
                count[(uint64_t)index.sa[n].i2 /length(genome[j])]+=1/occ*N;//<< std::endl;
                //std::cout << "index.sa " << index.sa[n].i2 - k << " " << length(anchor) << std::endl;
                //if (anchor[index.sa[n].i2 - k] == 0)
                //    anchor[index.sa[n].i2 - k] = k<<32;
                //anchor[index.sa[n].i2 - k] =  (anchor[index.sa[n].i2 - k] & (~mask)) + k;
                ////std::cout << "done " << std::endl;
                //uint64_t l = (anchor[index.sa[n].i2 - k] & mask) - (anchor[index.sa[n].i2 - k] >>32);
                //if ((anchor[index.sa[n].i2 - k] & mask) - (anchor[index.sa[n].i2 - k] >>32) > 500)
                //    std::cout << k << " anchor " << index.sa[n].i2 -k << " " << l << std::endl;
                if (anchor1[(index.sa[n].i2 - k)/100].i1 == 0)
                    anchor1[(index.sa[n].i2 - k)/100].i1 = k;
                anchor1[(index.sa[n].i2 - k)/100].i2 = k;
                uint64_t l = anchor1[(index.sa[n].i2 - k)/100].i2 - anchor1[(index.sa[n].i2 - k)/100].i1;
                if (l > 100)
                { 
                    std::cout << j << " anchor1 " << k << " " << (index.sa[n].i2 - k)/100 << " " << l << std::endl;
                    flag = 1;
                    break;
                }
            }
            if (flag)
                break;
        }
        //uint64_t max = 0;
        //unsigned K = 0;
        //for (uint64_t k = 0; k < length(count); k++)
        //    if(count[k] != 0) 
        //    {
        //        //std::cout << k << " " << count[k] << std::endl;
        //        if(max < count[k])
        //        {
        //            max = count[k];
        //            K = k;
        //        }
        //    }
        //
        //std::cout << " K " << " " << K << std::endl;//k << " " << count[k] << std::endl;//<< " " << (double)length(genome[j]) / (count+1) << std::endl;
    }
}

/*
void _detectAnchor(TIndex_u & index, StringSet<DnaString> genome)
{
    TShape_u shape;
    String<Pair<uint64_t, uint64_t> > anchor;
    uint64_t dir=0;
    for (uint64_t j = 0; j < length(genome); j++) 
    {
        resize(count, length(index.text)/length(genome[j])+1, 0);
        TIter it=begin(genome[j]);
        hashInit(shape, it);
        for (uint64_t k = 0; k < length(genome[j]) - shape.span + 1; k++)
        {
            hashNext(shape, it+k);
            dir = getBucket(index.bucketMap, shape.hValue);
            for (uint64_t n = index.dir[dir]; n < index.dir[dir + 1] ; n++)
            {
                anchor[index.sa[n]] 
            }
        }
    }
}
*/

void mnKmerAnalysis(StringSet<DnaString> & genome, StringSet<DnaString> & reads)
{
    
    TShape shape; 
    TIndex index(genome);
    String<uint64_t> anchor;
    uint64_t mask  = 131072 - 1, dn, count=0, sum = 0;
    double time = sysTime();
    std::cout << "mnKmerAnalysis() \n";

    resize(anchor, mask + 1); //4096=2^12
    //resize(anchor, lengthSum(genome));
    createQGramIndexDirOnly(index);
    _initAnchor(anchor, mask);
    //for (uint64_t j = 0; j < length(reads); j++)  
    for (uint64_t j = 0; j < 10000; j++)
    //for (uint64_t j = 27; j< 28; j++)
    {
        TIter it = begin(reads[j]);
        hashInit(shape, it);
        
        //_initAnchor(anchor, mask);
        //_initAnchor(anchor);
        for (uint64_t k = 0; k < (length(reads[j]) - shape.span + 1); k++) 
        {
            hashNext(shape, it + k);
            dn = getDir(index, shape); 
            
            for (uint64_t n = _getBodyCounth(index.dir[dn]); n < _getBodyCounth(index.dir[dn + 1]); n++)
            {
                _requestAnchor(anchor, (_getSA_i2(index.sa[n]) - k) >> 9, mask);
            //sum ^= (_getSA_i2(index.sa[n]) - k)>>9;
            }
        }
        //std::sort(begin(anchor), end(anchor), [](uint64_t & a, uint64_t & b){return (a & _AnchorMask_i1_ ) > (b & _AnchorMask_i1_);});
        uint64_t max1 = 0, max2 = 0, k1 = 0, k2 =0, countk = 0;
        for (uint64_t k = 0; k < length(anchor); k++)
        {
            if (anchor[k] != _Anchor_Empty_)
            {
                if (max1 < (anchor[k] & _AnchorMask_i1_))
                {
                    //max1 = anchor[k] & _AnchorMask_i1_;
                    max1 = _getAnchor_i2(anchor[k]);
                    k1 = _getAnchor_i1(anchor[k]);
         //           std::cout << " k1 " << k1 << std::endl;
                }
                anchor[k] = _Anchor_Empty_;
            }
            //else 
            //    countk++;
        }
            //else if(max2 < (anchor[k] & _AnchorMask_i1_))
            //{

            //    max2 = _getAnchor_i2(anchor[k]);
            //    k2 = _getAnchor_i1(anchor[k]);
            //}
        if (j < 100)
        std::cout << " " << j << " " << length(reads[j]) << " " << (k1 << 9) << " " << max1 << " " << k2 << " " << max2 << " " << countk / (float)length(anchor) << std::endl;
        //std:: cout << k << " " << (anchor[k] & _AnchorMask_i1_) << std::endl; //<< " " << std::bitset<64>(anchor[k]) << std::endl;
        //std::cout << j << " " << std::endl;
        sum ^= max1;
        if (max1 == 0)
            count++;
        
    }
    std::cout << "    count % = " << (float)count  << " " << length(reads) << " " << (float)count / length(reads) << std::endl;
    std::cout << "    anchor[0] = " << anchor[0] << std::endl;
    std::cout << "    End mnKmerAnalysis() systime() = " << sysTime() - time << std::endl;
}

void mnKmerAnalysis2(StringSet<DnaString> & genome, StringSet<DnaString> & reads)
{
    
    TShape shape; 
    TIndex index(genome);
    uint64_t mask  = 131072 - 1, dn, count=0, sum = 0;
    uint64_t anchor[mask + 1];
    double time = sysTime();
    std::cout << "mnKmerAnalysis() \n";

    createQGramIndexDirOnly(index);
    _initAnchor(anchor, mask);
    for (uint64_t j = 0; j < length(reads); j++)  
    //for (uint64_t j = 0; j < 3000; j++)
    //for (uint64_t j = 27; j< 28; j++)
    {
        TIter it = begin(reads[j]);
        hashInit(shape, it);
        _initAnchor(anchor, mask); 
        for (uint64_t k = 0; k < (length(reads[j]) - shape.span + 1); k++) 
        {
            hashNext(shape, it + k);
            dn = getDir(index, shape); 
            
            for (uint64_t n = _getBodyCounth(index.dir[dn]); n < _getBodyCounth(index.dir[dn + 1]); n++)
            {
                uint64_t x=_requestAnchor(anchor, (_getSA_i2(index.sa[n]) - k) >> 9, mask);
                //sum ^= (_getSA_i2(index.sa[n]) - k) >> 9;
                //if ((_getSA_i2(index.sa[n]) - k) >>9<<9==122704896 && j==91)
                //    std::cout << shape.hValue << " " << k << " " << _getSA_i2(index.sa[n]) << std::endl;
            }
        }
        uint64_t max1 = 0, max2 = 0, k1 = 0, k2 =0, countk = 0;
        //for (uint64_t k = 0; k < length(anchor); k++)
        for (uint64_t k = 0; k < mask + 1; k++) 
        {
            if (anchor[k] != _Anchor_Empty_)
            {
                if (max1 < (anchor[k] & _AnchorMask_i1_))
                {
                    max1 = _getAnchor_i2(anchor[k]);
                    k1 = _getAnchor_i1(anchor[k]);
                }
                else if(max2 < (anchor[k] & _AnchorMask_i1_))
                {
                    max2 = _getAnchor_i2(anchor[k]);
                    k2 = _getAnchor_i1(anchor[k]);
                }
                if (j==2)
                std::cout << (_getAnchor_i1(anchor[k]) << 9) << " " << _getAnchor_i2(anchor[k]) << std::endl;

                //anchor[k] = _Anchor_Empty_;
            }
            //else 
            //    countk++;
        }
        if (j < 100)
        std::cout << " " << j << " " << length(reads[j]) << " " << (k1 << 9) << " " << max1 << " " << k2 << " " << max2 << " " << countk / (float)length(anchor) << std::endl;
        //std:: cout << k << " " << (anchor[k] & _AnchorMask_i1_) << std::endl; //<< " " << std::bitset<64>(anchor[k]) << std::endl;
        //std::cout << j << " " << std::endl;
        sum ^= max1;
        if (max1 == 0)
            count++;
        
    }
    std::cout << "    count % = " << (float)count  << " " << length(reads) << " " << (float)count / length(reads) << std::endl;
    std::cout << "    anchor[0] = " << anchor[0] << std::endl;
    std::cout << "    End mnKmerAnalysis() systime() = " << sysTime() - time << std::endl;
}


int main(int argc, char** argv)
{
    if (argc < 3)
        return 1;
    double time;
    SeqFileIn gFile(toCString(argv[1]));
    SeqFileIn rFile(toCString(argv[2]));
    StringSet<CharString> ids;
    StringSet<DnaString> genome;
    StringSet<DnaString> reads;
    readRecords(ids, genome, gFile);
    readRecords(ids, reads, rFile);
    TIndex_u index_u(reads);
    //time = sysTime();
    //indexCreate(index_u, FibreSADir());
    //std::cout << sysTime() - time << std::endl;
    //opKmerAnalysis(index_u, genome);
    opKmerAnalysis2(genome, reads);
    //mnKmerAnalysis(genome, reads);
    //mnKmerAnalysis2(genome, reads);
    return 0;
} 
