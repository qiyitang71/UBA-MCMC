package jltl2dstar;

import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2ba.MyBitSet;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public class GFG extends NBA{
    private HashMap<Integer, HashMap<APElement, MyBitSet>> accEdges;
    /**
     * Constructor.
     * @param apset The underlying APSet
     */
    public GFG (APSet apset)
    {
        super(apset);
        accEdges = new HashMap<>();
    }
    public void addAccEdge(int from,APElement ap,int to){

        if(accEdges.get(from) == null){
            HashMap<APElement, MyBitSet> hm = new HashMap<>();
            MyBitSet bs = new MyBitSet(this.getAPSize());
            bs.set(to);
            hm.put(ap,bs);
            accEdges.put(from, hm);
        }else{
            HashMap<APElement, MyBitSet> hm = accEdges.get(from);
            if(hm.get(ap) == null){
                MyBitSet bs = new MyBitSet(this.getAPSize());
                bs.set(to);
                hm.put(ap,bs);
            }else{
                MyBitSet bs = hm.get(ap);
                bs.set(to);
            }
        }
    }
    public void addAccEdge(int from,APMonom apm,int to){
        APSet ap_set = this.getAPSet();

        for (Iterator<APElement> it = apm.APElementIterator(ap_set); it.hasNext(); ) {
            APElement cur = it.next();
            addAccEdge(from, cur, to);
            // System.out.println("State " + _state.getName() + " added edge to " + state.getName() + " through " + cur);
        }
    }


}
