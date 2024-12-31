package jltl2dstar;

import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2ba.MyBitSet;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Vector;
// for transition-based acceptance
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

    public boolean isAcc(int from,APElement ap,int to){
        if(accEdges.get(from) == null){
            return false;
        }else{
            HashMap<APElement, MyBitSet> hm = accEdges.get(from);
            if(hm.get(ap) == null){
                return false;
            }else{
                MyBitSet bs = hm.get(ap);
                return bs.get(to);
            }
        }
    }

    public HashMap<Integer, HashMap<APElement, MyBitSet>> getAccEdges(){
        return accEdges;
    }

    /** Print the NBA as a HOA automaton to out */
    public void print_hoa(PrintStream out) {
        APSet _apset = this.getAPSet();
        out.println("HOA: v1");
        out.println("States: "+size());
        _apset.print_hoa(out);
        out.println("Start: "+getStartState().getName());
        out.println("Acceptance: 1 Inf(0)");
        out.println("acc-name: Buchi");
        out.println("properties: trans-labels explicit-labels trans-acc");
        out.println("--BODY--");
        for (NBA_State state : this.getIndexSet()) {
            out.print("State: "+state.getName());  // id
            out.println((state.isFinal() ? " {0}" : ""));

            for (Map.Entry<APElement, MyBitSet> edge : state) {
                APElement label = edge.getKey();
                String labelString = "["+label.toStringHOA(_apset.size())+"]";
                MyBitSet to_states = edge.getValue();
                for (Integer to : to_states) {
                    out.print(labelString);
                    out.print(" ");
                    out.println((isAcc(state.getName(),label,to) ? to + " {0}" : to));
                }
            }
        }
        out.println("--END--");
    }

}
