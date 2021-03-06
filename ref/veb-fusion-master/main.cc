#include "veb.h"
#include "fusion.h"
#include <set>
#include <vector>
#include "util.h"

int main(int argc, char* argv[])
{
    FusionTree<unsigned long> test;
    
    //veb<unsigned long long> test;
    /*
    srand(time(NULL));
    std::unordered_set<int> inserts;
    std::cout << "Start Inserts" << std::endl;
    for(int i = 0; i < 500000; ++i)
    {
        //  test.insert(((unsigned long long)rand()<<32)|rand());
        int t = rand();
        
        if(inserts.count(t) > 0)
        {
            //std::cout << "duplicate value " << t << std::endl;
            //  exit(1);
        }
        else {
            test.insert(t);
            inserts.insert(t);
        }
             
    }
    std::cout << "Initialize" << std::endl;
    test.initialize();
    unsigned long next = rand();
    std::cout << "Start successor" << std::endl;
    std::cout << "suc " << next << ":" << test.successor(next) << std::endl;
    */

    //for (int i = 1; i < 26; ++i) test.insert(i);
    //std::cout << "Initialize" << std::endl;
    //test.initialize();
    //test.traverse();
    //std::cout << "suc " << 13 << ":" << test.successor(13) << std::endl;
    //std::cout << "suc " << 7 << ":" << test.successor(7) << std::endl;
    //std::cout << "suc " << 20 << ":" << test.successor(20) << std::endl;

    std::vector<unsigned long> keys;
    keys.push_back(4);
    keys.push_back(8);
    keys.push_back(13);
    std::vector<int> bs = get_impor_bits(keys);
    //for(int i = 0; i < bs.size(); ++i) {
    //    std::cout << bs[i] << std::endl;
    //}
    unsigned long m;
    std::vector<int> ms = get_m(bs, m);
    //for(int i = 0; i < ms.size(); ++i) {
    //    std::cout << ms[i] << std::endl;
    //}
    //std::cout << m << std::endl;
    std::cout << bs.back() << std::endl << ms.back() << std::endl << bs.front() << std::endl << ms.front() << std::endl;

    // std::set<unsigned long long> rbTree;

    // std::cout << "Start RB Tree" << std::endl;
    // for(int i = 0; i < 500000; ++i)
    // {
    //     rbTree.insert(rand());
    // }
    // unsigned  long long n = rand();
    // std::cout << "Start successor" << std::endl;
    // unsigned long long suc = *(rbTree.upper_bound(n));
    // std::cout << "suc " << n << ":" << suc << std::endl;


    // std::cout << "Test important bits" << std::endl;
    // std::vector<unsigned long long > keys;
    // keys.push_back(16);
    // keys.push_back(20);
    // keys.push_back(21);
    // keys.push_back(29);
    // std::vector<int> msbs = get_impor_bits(keys);
    // for(int i = 0; i < msbs.size(); ++i)
    // {
    //     std::cout << msbs[i] << std::endl;
    // }

    // unsigned long long mask = get_mask<unsigned long long>(msbs);
    // std::cout << "mask:" << mask << std::endl;

    // unsigned long  long m;

    // std::vector<int> m_bits = get_m(msbs, m);

    // for(int i = 0; i < msbs.size(); ++i)
    //     std::cout << msbs[i]+m_bits[i] << std::endl;

    // std::cout << "gap:" << msbs.back()+m_bits.back()-msbs.front()-m_bits.front() << std::endl;
    // std::cout << m << std::endl;

    // unsigned long long mh = mul_high(2305843009213693952LL, 2305843009213693952LL);

    // std::cout << mh << std::endl;
    
}
