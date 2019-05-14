package neotango;

import java.util.ArrayList;
import java.lang.Math;
import java.lang.Integer;
import java.lang.Long;
//import neotango.Utils;

class FusionTree {
    
    private class Node {
        // should be sorted
        ArrayList<Long> keys;
        ArrayList<Node> children;

        boolean isLeaf;
        int n; 

        ArrayList<Integer> bs; // important bits
        ArrayList<Integer> ms; // m_j bits

        long m; // multiplicative contant for approximate sketching
        long b_mask; // mask out the important bits
        long bm_mask; // mask out the important bits and m_j bits
        int sketch_gap; // space a sketch takes up
        long sketches; // word containing all sketches
        
        long sketch_maskl; // mask with 1 at start of sketch
        long sketch_maskh; // mask with 1 at end of sketch

        Node(int max_keys) {
            this.keys = new ArrayList<Long>(max_keys);
            this.children = new ArrayList<Node>(max_keys+1);
            this.n = 0;
        }

        Node() {this.n = 0;}

        void precompute() {
            if (this.n == 0) return;
            this.bs = Utils.get_impor_bits(keys);
            this.ms = Utils.get_m_bits(this.bs);
            this.m = Utils.get_m(this.ms);
            
            this.b_mask = Utils.get_mask(this.bs);
            this.bm_mask = Utils.get_combo_mask(this.bs, this.ms);

            int r = this.bs.size();
            int bb = this.bs.get(r-1).intValue();
            int mb = this.ms.get(r-1).intValue();
            int bf = this.bs.get(0).intValue();
            int mf = this.ms.get(0).intValue();
            this.sketch_gap = bb + mb - bf - mf;
            this.sketch_gap = (this.sketch_gap == 0) ? 1 : this.sketch_gap;

            for(int j = 0; j < this.n; ++j) {
                long cur_key = this.keys.get(this.n - j - 1).longValue();
                long sketch = Utils.approx_sketch(this.m, cur_key, this.b_mask, this.bm_mask);
                this.sketches |= (sketch | (1L << this.sketch_gap)) << j*(this.sketch_gap + 1);
                this.sketch_maskl |= (1L << j*(this.sketch_gap + 1));
                this.sketch_maskh |= (1L << this.sketch_gap) << j*(this.sketch_gap + 1);
            }
        }
    }

    Node root;
    int max_keys;
    int min_keys;

    public FusionTree() {
         this.max_keys = 2;
         this.min_keys = max_keys/2;
         this.root = new Node(max_keys);
         this.root.isLeaf = true;
    }

    // x is a non-full parent, i is point to split on
    private void split_child(Node x, int i) {
        Node z = new Node(max_keys);
        Node y = x.children.get(i);

        int pivot = max_keys/2;
        z.n = max_keys - pivot - 1;
        
        for(int j = 0; j < z.n; ++j) {
            z.keys.set(j, y.keys.get(pivot+j));
        }

        if(!y.isLeaf) {
            for(int j = 0; j < z.n+1; ++j) {
                z.children.set(j, y.children.get(pivot+j+1));
            }
        }

        x.children.add(i+1, z);
        x.children.remove(x.children.size()-1);

        x.keys.add(i, y.keys.get(pivot));
        x.keys.remove(x.keys.size()-1);
    }

    // 
    private void insert_nonfull(Node x, long k) {
        int i = x.n;

        if(x.isLeaf) {
            while(i >= 1 && k < x.keys.get(i-1)) 
                i--;
            x.keys.add(i, new Long(k));
            x.n++;
        } else {
            while(i >= 1 && k < x.keys.get(i-1)) {
                i--;
            }
            i++;

            if(x.children.get(i-1).n == max_keys)
            {
                split_child(x, i-1);
                if( k > x.keys.get(i-1))
                    i++;
            }
            insert_nonfull(x.children.get(i-1), k);
        }
    }
   
    // finds a query successor in the tree 
    private long fusion_successor(Node x, long k) {
        if(x.n == 0) {
            if(x.isLeaf)
                return Long.MIN_VALUE;
            return fusion_successor(x.children.get(0), k);
        }

        long app_sketch = Utils.approx_sketch(x.m,
                                              k,
                                              x.b_mask,
                                              x.bm_mask);

        long repeat_sketch = app_sketch * x.sketch_maskl;
        
        int i1 = Utils.par_comp(x.sketches,
                                repeat_sketch,
                                x.sketch_maskh,
                                x.sketch_maskl,
                                x.n,
                                x.sketch_gap);

        long y = 0;
        
        if(i1 < x.n) {
            long lcp1 = k ^ x.keys.get(i1);
            long msb1 = Long.numberOfLeadingZeros(lcp1);
            y = msb1;
        }
        
        if(i1 != 0) {
            long lcp2 = k ^ x.keys.get(i1-1);
            long msb2 = Long.numberOfLeadingZeros(lcp2);
            y = (y > msb2) ? y : msb2; // max
        }

        // set y to next least significant bit
        // ie. does it go left or right
        y = 63 - y;

        // this puts a 1 at first bit which is different
        long e_mask = 1L << y;
        // This makes all bits after that 1s
        e_mask -= 1;

        // this puts a 1 at the first bit after
        long new_mask = 1L << y;
        new_mask ^= -1L; // hopefully this has one 0
        
        long e = k | e_mask; 
        if(((1L << y) & k) != 0) {
            e &= new_mask;
        } else {
            e &= -1L ^ e_mask;
            e |= 1L<<y;
        }

        app_sketch = Utils.approx_sketch(x.m,
                                         e,
                                         x.b_mask,
                                         x.bm_mask);
        
        int i = Utils.par_comp(x.sketches,
                               app_sketch*x.sketch_maskl,
                               x.sketch_maskh,
                               x.sketch_maskl,
                               x.n,
                               x.sketch_gap);

        // this check is necessary due to singlteon nodes being weird
        // This check could be at the top but I want to give the fusion
        // buddies a chance.
        if(x.n != i && k > x.keys.get(i)) {
            i++;
        }

        //std::cout << "index after is:  " << i << std::endl;
        if(x.n != i) {
            System.out.println(k + " is less than " + x.keys.get(i) + '\n');
            //std::cout << "test:" << (k < x->keys[i]) << std::endl; 
        } else {
            System.out.println(k + " is bigger than all " + x.keys.get(x.n-1) + '\n');
            //std::cout << "test:" << (k > x->keys[x->n-1]) << std::endl;
        }
        
        if(x.isLeaf) {
            if(i == x.n)
                return Long.MIN_VALUE;
            
            return x.keys.get(i).longValue();
        }

        long succ = fusion_successor(x.children.get(i), k);
        if(succ == Long.MIN_VALUE) {
            if(i == x.n)
                return succ;
            return x.keys.get(i);
        }
        else return succ;
    }

    private void initialize(Node x) {
        x.precompute();
        if(!x.isLeaf) {
            for(int i = 0; i < x.n+1; ++i) {
                initialize(x.children.get(i));
            }
        }
    }

    public void insert(long k) {
        Node r = root;
        //if(r.isLeaf) System.out.println("root is leaf");
        if(root.n == max_keys) {
            Node s = new Node(max_keys);
            root = s;
            s.isLeaf = false;
            s.n = 0;
            s.children.add(r);
            split_child(s, 0);
            insert_nonfull(s, k);
        } else {
            insert_nonfull(r, k);
        }
    }

    long successor(long k) {
        return fusion_successor(root, k); //successor(root, k);
    }

    public void initialize() {
        initialize(root);
    }
}
