package neotango;

import java.lang.Math;
import java.lang.Integer;
import java.lang.Long;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Collections;

class Utils {

    static int high(long x) {
        return (int) (x >> 32 & 0xffffffff);
    }

    static int low(long x) {
        return (int) (x & 0xffffffff); 
    }

    // 0 indexed msb -- returns the index
    static int get_msb(long arg) {
        if (arg == 0) return -1;
        return 63-Long.numberOfLeadingZeros(arg);
    }

    // finding the sketch bits
    static ArrayList<Integer> get_impor_bits(ArrayList<Long> keys) {
        ArrayList<Integer> ret = new ArrayList<Integer>(keys.size());
        HashSet<Integer> check = new HashSet<Integer>(keys.size());

        for(int i = 0; i < keys.size(); ++i) { 
            for(int j = 0; j < keys.size(); ++j) {
                long xored = keys.get(i).longValue() ^ keys.get(j).longValue();
                if(xored!=0) {
                    Integer val = new Integer(get_msb(xored));
                    if(!check.contains(val)) {
                        ret.add(val);
                        check.add(val);
                    }
                }
            }
        }
        Collections.sort(ret/*, Collections.reverseOrder()*/);
        return ret;
    }

    // create a bitmask from bit indices
    static long get_mask(ArrayList<Integer> mask_bits) {
        long ret = 0;
        for(Integer i: mask_bits) ret |= 1L << i;
        return ret;
    }
    
    // create bitmask from pairwise sums of bit indices
    static long get_combo_mask(ArrayList<Integer> a, ArrayList<Integer> b) {
        ArrayList<Integer> temp = new ArrayList<Integer>(a.size());
        for(int i = 0; i < a.size(); ++i) 
            temp.add(new Integer(a.get(i).intValue() + b.get(i).intValue()));
        return get_mask(temp);
    }

    // tests if a bit can be an m_j bit
    private static boolean m_bit_works(ArrayList<Integer> m_primes, ArrayList<Integer> b_bits, int mt) {
        int i = m_primes.size();
        int r = b_bits.size();
        for(int j = 0; j < i; ++j) {
            for(int k = 0; k < r; ++k) {
                for(int l = 0; l < r; ++l) {
                    if(mt == m_primes.get(j).intValue() + b_bits.get(k).intValue() - b_bits.get(l).intValue()) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // finds the m_j bits to make up the multiplication constant m for approximate sketching
    static ArrayList<Integer> get_m_bits(ArrayList<Integer> b_bits) {
        int r = b_bits.size();
        int r3 = r*r*r;
        ArrayList<Integer> m_primes = new ArrayList<Integer>(r);
        ArrayList<Integer> m_bits = new ArrayList<Integer>(r);

        int mt = 0;
        boolean okay = false;
        for(int i = 0; i < r; ++i) {
            boolean found = false;
            while(!found) {
                if(m_bit_works(m_primes, b_bits, mt)) {
                    found = true;
                    m_primes.add(new Integer(mt));
                    int mt_final = 64 - b_bits.get(i).intValue() + i*r3;
                    mt_final = mt_final / r3 * r3;
                    mt_final += mt;
                    m_bits.add(new Integer(mt_final));
                } else {
                    ++mt;
                }
            }
        }
        //for(Integer i : m_primes) System.out.println(i);
        return m_bits;
    }

    // uses the m_bits to compute the multiplication cnstant m for approximate sketching
    static long get_m(ArrayList<Integer> m_bits) {
        long m = 0;
        for(Integer mj: m_bits) 
            m |= 1L << mj;
        return m;
    }

    // result is 32 high order bits of a multiplication
    static long mul_high(long a, long b) {
        long a1 = a >> 32;
        long a0 = a & 0xffffffff;
        long b1 = b >> 32;
        long b0 = b & 0xffffffff;
   
        long res = a1*b1;
        res += (a1*b0) >> 32;
        res += (a0*b1) >> 32;

        return res;
    }

    // approximate sketch in fusion tree
    static long approx_sketch(long m, long x, long b_mask, long bm_mask, int shift_dis) {
        long x_prime = x & b_mask;
        long as = mul_high(x_prime, m);
        return as & bm_mask;
    }

    // performs parallel comparison on the fusion tree node
    static int par_comp(long sketch_node, long sketch_k, long sketch_maskh, long sketch_maskl, int k, int gap) {
        long diff = sketch_node - sketch_k;
        diff &= sketch_maskh;
        int msb = get_msb(diff);

        if(msb < 0) return k; //bigger than all
        int ind = msb/(gap+1);
        return k-ind-1;
    }   

    public static void main(String[] args) {
        /*
        long x = 0x0000000100000001L;
        System.out.println(x);
        System.out.println(high(x));
        System.out.println(low(x));
        System.out.println(get_msb(x));
        */
        ArrayList<Long> keys = new ArrayList<Long>();
        keys.add(new Long(4)); //0100
        keys.add(new Long(8)); //1000
        keys.add(new Long(13)); //1101
        ArrayList<Integer> b_bits = get_impor_bits(keys); // 2 3
        for(Integer i: b_bits) System.out.println(i);
        ArrayList<Integer> m_bits = get_m_bits(b_bits); //56 66
        long m = get_m(m_bits);
        for(Integer i: m_bits) System.out.println(i);
        System.out.println(m); //72057594037927940 = 2^56 + 2^2
    }
}
