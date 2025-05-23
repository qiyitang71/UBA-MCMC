package explicit.gfg;

import java.awt.Point;
import java.util.*;
import java.util.Map.Entry;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.RCDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import common.IterableBitSet;
import common.IterableStateSet;
import explicit.*;
import explicit.uba.DTMCUBAProduct;
import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2ba.MyBitSet;
import jltl2dstar.GFG;
import jltl2dstar.NBA;
import jltl2dstar.NBA_State;
import parser.State;
import parser.VarList;
import prism.Pair;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismSettings;

public class DTMCGFGProduct extends MDPSimple {

    protected HashMap<Integer, Set<APElement>> accEdges;
    protected int ubaSize;
    protected int invMap[];
    protected int dtmcSize;

    private PrismLog mainLog;
    private PrismSettings settings;
    private int verbosity = 0;

    private enum MatrixType {MATRIX_DENSE, MATRIX_SPARSE, MATRIX_RC};
    private explicit.gfg.DTMCGFGProduct.MatrixType matrixType = explicit.gfg.DTMCGFGProduct.MatrixType.MATRIX_RC;

    private APSet actions;
    //public DTMCUBAProduct dtmcProduct;

    HashMap<Integer, Integer> stateMap = new HashMap<>();
    ArrayList<Integer>  splitterIds = new ArrayList<>();

    HashMap<Integer, Integer>  equivDTMC = new HashMap<>();;
    HashMap<Integer, BitSet> equivGFG;


    private DoubleMatrix2D newMatrix(int rows, int columns) {
        switch (matrixType) {
            case MATRIX_DENSE:
                return new DenseDoubleMatrix2D(rows, columns);
            case MATRIX_RC:
                return new RCDoubleMatrix2D(rows, columns);
            case MATRIX_SPARSE:
                return new SparseDoubleMatrix2D(rows, columns);
            default:
                throw new IllegalStateException();
        }
    }

    public int getDTMCState(final int prodState) {
        return invMap[prodState]/ubaSize;
    }

    public int getUBAState(final int prodState) {
        return invMap[prodState]  % ubaSize;
    }

    public HashMap<Integer, Integer>  getEquivDTMCStates(){
        return equivDTMC;
    }
    public HashMap<Integer, BitSet> getEquivGFG(){
        return equivGFG;
    }

    public DTMCGFGProduct(PrismComponent parent, DTMC dtmc, GFG uba, Vector<BitSet> labelBS, BitSet statesOfInterest) throws PrismException {
        this.mainLog = parent.getLog();
        this.settings = parent.getSettings();
        this.dtmcSize = dtmc.getNumStates();

        verbosity = settings.getInteger(PrismSettings.PRISM_UBA_VERBOSITY);

        if (verbosity >= 2) {
            mainLog.println("Start constructing product");
        }

        ubaSize = uba.getStateCount();
        int numAPs = uba.getAPSize();
        int modelNumStates = dtmc.getNumStates() + 1;
        int prodNumStates = modelNumStates * ubaSize;
        int s_1, s_2, q_1, q_2;
        MyBitSet s_labels = new MyBitSet(numAPs);
        List<State> prodStatesList = null;

        VarList newVarList = null;

        if (dtmc.getVarList() != null) {
            newVarList = (VarList) dtmc.getVarList().clone();
        }

        // Create a (simple, mutable) model of the appropriate type
        setVarList(newVarList);

        // Encoding:
        // each probabiblistic state s' = <s, q> = 2*(s * ubaSize + q)
        // each nondeterministic state s' = <s, q> = 2*(s * ubaSize + q) + 1
        // s(s') = s'/ (ubaSize)
        // q(s') = (s') % (ubaSize)

        LinkedList<Point> queue = new LinkedList<Point>();
        int map[] = new int[prodNumStates];
        Arrays.fill(map, -1);

        if (dtmc.getStatesList() != null) {
            prodStatesList = new ArrayList<State>();
        }

        //check which APs are projected to the same original letter
        Map<String, Set<Integer>> labelStringMap = new HashMap<>();
        for (int k = 0; k < numAPs; k++) {
            String key = uba.getAPSet().getAP(k).split("_")[0];
            if (labelStringMap.containsKey(key)) {
                labelStringMap.get(key).add(k);
            } else {
                Set<Integer> set = new HashSet<>();
                set.add(k);
                labelStringMap.put(key, set);
            }
        }

        Map<Integer, Set<Integer>> labelMap = new HashMap<>();
        for (int k = 0; k < numAPs; k++) {
            String key = uba.getAPSet().getAP(k).split("_")[0];
            labelMap.put(k, labelStringMap.get(key));
        }



        // We need results for all states of the original model in statesOfInterest
        // We thus explore states of the product starting from these states.
        // These are designated as initial states of the product model
        // (a) to ensure reachability is done for these states; and
        // (b) to later identify the corresponding product state for the original states
        //     of interest
        if (verbosity >= 2) {
            mainLog.println("labelBS: " + labelBS);
            mainLog.println("labelMap: " + labelMap);
        }


        int newDtmcInit = dtmc.getNumStates();
        for (int s_0 : new IterableStateSet(statesOfInterest, dtmc.getNumStates())) {
            // Get BitSet representing APs (labels) satisfied by state s_0
            if (verbosity >= 2) {
                mainLog.println("statesOfInterest: " + statesOfInterest);
                mainLog.println("s_labels: " + s_labels);
            }

            // Find corresponding initial state in DA
            NBA_State ubaStartState = uba.getStartState();
            MyBitSet destinations = new MyBitSet(uba.size());

            for (int k = 0; k < numAPs; k++) {
                if (!labelBS.get(k).isEmpty()) {
                    for (int k2 : labelMap.get(k)) {
                        MyBitSet s_labels2 = new MyBitSet(numAPs);
                        s_labels2.set(k2, labelBS.get(k).get(s_0));
                        destinations.or(ubaStartState.getEdge(new APElement(s_labels2)));
                    }
                }
            }

            if (destinations.isEmpty()) {
                // Language is empty starting from this DTMC start state -> just skip
                continue;
            }
            if (verbosity >= 2) {
                mainLog.println("LMC start state = " + s_0 + ", destinations =" + destinations);
            }

            for(Iterator<Integer> initialStatesIterator = destinations.iterator(); initialStatesIterator.hasNext();) {
                // Add (initial) state to product
                Integer q_0 = initialStatesIterator.next();
                queue.add(new Point(s_0,q_0.intValue()));
                addState();
                addInitialState(getNumStates() - 1);
                if (verbosity >= 2) {
                    mainLog.println("INITIAL: " + (getNumStates()-1) + "->(" + s_0 + "," + q_0.intValue() + ")" + " in constructing GFG DTMC product ");
                }
                map[s_0 * ubaSize + q_0] = getNumStates() - 1;
            }
        }

        // Product states
        //System.out.println("Map = " + Arrays.toString(map));
        BitSet visited = new BitSet(prodNumStates);
        while (!queue.isEmpty()) {
            Point p = queue.pop();
            s_1 = p.x;
            q_1 = p.y;
            visited.set(s_1 * ubaSize + q_1);


            // Go through transitions from state s_1 in original model
            for (APElement ap : uba.nba_i_getAPSet().elements()) {

                // Find corresponding successor in UBA
                NBA_State ubaState = uba.get(q_1);
                MyBitSet destinations = ubaState.getEdge(ap);
                if (destinations.isEmpty() || ap.cardinality() != 1) continue;

                Iterator<Map.Entry<Integer, Double>> iter = ((DTMC) dtmc).getTransitionsIterator(s_1);

                while (iter.hasNext()) {
                    Map.Entry<Integer, Double> e = iter.next();
                    s_2 = e.getKey();
                    double prob = e.getValue();

                    //check can continue
                    int idx = ap.nextSetBit(0);
                    //mainLog.println("idx: " + idx);

                    boolean isSet = false;
                    for (int i : labelMap.get(idx)) {
                        if (labelBS.get(i).get(s_2)){
                            isSet = true;
                            break;
                            }
                    }
                    if (!isSet) continue;


                    //MyBitSet destinations = ubaState.getEdge(new APElement(s_labels));
                    Distribution dist = new Distribution();
                    for (Integer destination : destinations) {
                        q_2 = destination;
                        if (q_2 < 0) {
                            throw new PrismException("The deterministic automaton is not complete (state " + q_1 + ")");
                        }
                        // Add state/transition to model
                        if (!visited.get(s_2 * ubaSize + q_2) && map[s_2 * ubaSize + q_2] == -1) {
                            queue.add(new Point(s_2, q_2));
                            //Create probabilistic state
                            addState();
                            map[s_2 * ubaSize + q_2] = getNumStates() - 1;
                            if (prodStatesList != null) {
                                // Store state information for the product
                                prodStatesList.add(dtmc.getStatesList().get(s_2));
                            }
                        }
                        dist.add(map[s_2 * ubaSize + q_2], (1.0 * prob) / (destinations.cardinality()));
                    }
                    addActionLabelledChoice(map[s_1 * ubaSize + q_1], dist, ap);
                }
            }
        }


        // Build a mapping from state indices to states (s,q), encoded as (s * daSize + q)
        //System.out.println("Map = " + Arrays.toString(map));
        invMap = new int[getNumStates()];
        for (int i = 0; i < map.length; i++) {
            if (map[i] != -1) {
                invMap[map[i]] = i;
            }
        }

        findDeadlocks(false);

        if (prodStatesList != null) {
            setStatesList(prodStatesList);
        }

        //final MDPSparse productModel;
        //mainLog.println("Converting product model to MDPSparse");

        //Lift acceptance condition
        accEdges = new HashMap<>();

        for (int i = 0; i < getNumStates(); i++) {
            int q = getUBAState(i);
            //APSet apset = new APSet();
            if (uba.getAccEdges().get(q) != null) {
                Set<APElement> apset = uba.getAccEdges().get(q).keySet();
                accEdges.put(i, apset);
            }
        }

        actions = uba.getAPSet();
        if (verbosity >= 3) {
            //mainLog.println("UBA accepting transitions: " + uba.getAccEdges());
            //mainLog.println("Accepting transitions: " + accEdges);
            //uba.print_hoa(System.out);
            //mainLog.println("Product: " + this.toString());
            mainLog.println("Labels: " + uba.getAPSet());
            mainLog.print("invMap: ");
            for (int i = 0; i < this.numStates; i++) {
                mainLog.print(invMap[i] + ", ");
            }
            mainLog.print("\n");
        }

        //// generate equivalent GFG states
        equivGFG = new HashMap<>();
        Iterator<NBA_State> iter = uba.iterator();
        while (iter.hasNext()) {
            NBA_State s = iter.next();
            for (APElement ap : uba.nba_i_getAPSet().elements()) {
                if (ap.cardinality() != 1) continue;
                MyBitSet succs = s.getEdge(ap);
                if (succs.cardinality() > 1) {
                    for (int t : IterableBitSet.getSetBits(succs)) {
                        //if(bisMap.containsKey(t)) break;
                        equivGFG.put(t, succs);
                    }
                }

            }
        }

        if(verbosity >= 3) {
            mainLog.println("equivGFG = " + equivGFG);
        }

        /// compute equivalent LMC states

        /// compute equivalent LMC states
        // Step 0: build dtmc successor map

        Map<Integer, Set<Integer>> symbolMap = new HashMap<>();
        int index = 0;
        for(String key: labelStringMap.keySet()){
            symbolMap.put(index, labelStringMap.get(key));
            //FIX
        }

        HashMap<Integer, HashMap<Integer,Integer>> succMap = new HashMap<>(); //successors map for dtmc: S x Symbols -> 2^S
        for(int s = 0; s < dtmc.getNumStates(); s++ ){
            HashMap<Integer,Integer> succLetterMap = new HashMap<>();
            Iterator<Map.Entry<Integer, Double>> iter2 = ((DTMC) dtmc).getTransitionsIterator(s);
            while(iter2.hasNext()){
                int t = iter2.next().getKey();
                for(int letter = 0; letter < symbolMap.size(); letter++){
                    for(int ap: symbolMap.get(letter)){
                        if(!labelBS.get(ap).isEmpty()){
                            if(labelBS.get(ap).get(t)){
                                succLetterMap.put(letter, t);
                            }
                        }
                    }
                }
            }
            for(int letter = 0; letter < symbolMap.size(); letter++){
                if(!succLetterMap.containsKey(letter)){
                    succLetterMap.put(letter, dtmc.getNumStates());
                }
            }
            succMap.put(s,succLetterMap);
        }
        {
            HashMap<Integer, Integer> succLetterMap = new HashMap<>();
            for(int letter = 0; letter < symbolMap.size(); letter++){
                succLetterMap.put(letter, dtmc.getNumStates());
            }
            succMap.put(dtmc.getNumStates(),succLetterMap);
        }
        if(verbosity >= 2) {
            mainLog.println("succMap = " + succMap.size());
        }

        // Step 1: Partition states into final and non-final sets
        Set<Integer> nonFinalStates = new HashSet<>();
        for (int state = 0; state < dtmc.getNumStates(); state++) {
            nonFinalStates.add(state);
        }

        Set<Integer> finalStates = new HashSet<>();
        finalStates.add(dtmc.getNumStates());

        List<Set<Integer>> partitions = new ArrayList<>();
        partitions.add(finalStates);
        partitions.add(nonFinalStates);

//        if(verbosity >= 2) {
//            mainLog.println("partitions = " + partitions);
//        }
        // Step 2: Refine partitions using transition-based distinctions
        Queue<Set<Integer>> workQueue = new LinkedList<>(partitions);
        while (!workQueue.isEmpty()) {
            Set<Integer> currentSplitter = workQueue.poll();

            for (int letter = 0; letter < symbolMap.size(); letter++) {
                Map<Integer, Set<Integer>> reverseTransitions = new HashMap<>();
                for (int state = 0; state < dtmc.getNumStates()+1; state++) {
                    int targetState = succMap.get(state).get(letter);
                    if (currentSplitter.contains(targetState)) {
                        reverseTransitions.putIfAbsent(state, new HashSet<>());
                        reverseTransitions.get(state).add(state);
                    }
                }

                List<Set<Integer>> newPartitions = new ArrayList<>();
                for (Set<Integer> partition : partitions) {
                    Set<Integer> intersect = new HashSet<>();
                    Set<Integer> difference = new HashSet<>();

                    for (int state : partition) {
                        if (reverseTransitions.containsKey(state)) {
                            intersect.add(state);
                        } else {
                            difference.add(state);
                        }
                    }

                    if (!intersect.isEmpty() && !difference.isEmpty()) {
                        newPartitions.add(intersect);
                        newPartitions.add(difference);
                        workQueue.add(intersect);
                        workQueue.add(difference);
                    } else {
                        newPartitions.add(partition);
                    }
                }
                partitions = newPartitions;
            }
        }

        if(verbosity >= 3) {
            mainLog.println("partitions = " + partitions);
        }

        equivDTMC = new HashMap<>();
        for(int idx = 0; idx < partitions.size(); idx++){
            for(int s: partitions.get(idx)){
                equivDTMC.put(s, idx);
            }
        }

        /*
        equivDTMC = new HashMap<>();
        for(int s = 0; s < dtmc.getNumStates(); s++){
            BitSet bs = new BitSet(dtmc.getNumStates());
            bs.set(0,dtmc.getNumStates());
            equivDTMC.put(s, bs);
        }
        */

        if(verbosity >= 2) {
            mainLog.println("equivDTMC size = " + partitions.size());
        }
    }


    public int getNumDtmcStates(){
        return this.dtmcSize;
    }
    public HashMap<Integer, Set<APElement>> getAccEdges() {
        return accEdges;
    }

    DoubleMatrix2D positivityMatrixForMCC(BitSet mcc, Map<Integer,Integer> map, boolean subtractIdentity) throws PrismException {

        int size = 0;
        if (map == null) {
            map = new LinkedHashMap<Integer,Integer>();
        }
        for (int i : IterableBitSet.getSetBits(mcc)) {
            map.put(i, size);
            size++;
        }

        if (verbosity >= 3) {
            mainLog.println("MCC mapping = "+map);
        }

        DoubleMatrix2D matrix = newMatrix(size,size);

        for (Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < getNumChoices(from); c++) {
                for (Iterator<Entry<Integer, Double>> it = this.getTransitionsIterator(from,c); it.hasNext(); ) {
                    Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();

                    if (mcc.get(to)) {
                        matrix.setQuick(map.get(from), map.get(to), matrix.get(map.get(from),map.get(to))+ probMove.getValue());
                    }
                }
            }
        }

        if (verbosity >= 2) {
            //mainLog.println("A = \n"+matrix.toString());
        }

        if (subtractIdentity) {
            for (int i = 0; i < size; i++) {
                matrix.setQuick(i,i,  matrix.getQuick(i, i) - 1.0);
            }

            if (verbosity >= 2) {
               // mainLog.println("A - I = \n"+matrix.toString());
            }
        }

        return matrix;
    }

    DoubleMatrix2D valueMatrixForSCC(BitSet scc, Set<Integer> cut, LinkedHashMap<Integer,Integer> map) throws PrismException {

        int size = 0;
        for (int i : IterableBitSet.getSetBits(scc)) {
            map.put(i, size);
            size++;
        }

        // sanity check:
        for (int i : cut) {
            if (!map.containsKey(i)) {
                throw new PrismException("Implementation error: Cut is not contained in probStates(SCC)");
            }
        }

        if (verbosity >= 3) {
            mainLog.println("MCC mapping = "+map);
        }

        // one row additionally for cut
        DoubleMatrix2D matrix = newMatrix(size+1, size);

        for (Map.Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < getNumChoices(from); c++) {
                for(Iterator<Map.Entry<Integer, Double>> it = this.getTransitionsIterator(from, c); it.hasNext();) {
                    Map.Entry<Integer,Double> probMove = it.next();
                    int to = probMove.getKey();

                    if (scc.get(to)) {
                        matrix.setQuick(map.get(from), map.get(to), matrix.get(map.get(from), map.get(to)) + probMove.getValue());
                    }
                }
            }
        }

        if (verbosity >= 2) {
            //mainLog.println("A = \n"+matrix.toString());
        }

        for (int i = 0; i < size; i++) {
            matrix.setQuick(i,i,  matrix.getQuick(i, i) - 1.0);
        }

        if (verbosity >= 2) {
            //mainLog.println("A - I = \n"+matrix.toString());
        }

        // Add cut
        for (int i : cut) {
            matrix.setQuick(size, map.get(i), 1.0);
            mainLog.println("cut i = "+i + ", index = " + map.get(i));
        }

        if (verbosity >= 2) {
            //mainLog.println("(A - I) + cut condition = \n"+matrix.toString());
        }

        return matrix;
    }

    Pair<DoubleMatrix2D,DoubleMatrix2D> reachabilityEquationSystem(BitSet target, BitSet unknown, StateValues knownValues, LinkedHashMap<Integer,Integer> map) throws PrismException {

        int size = 0;
        for (int i : IterableBitSet.getSetBits(unknown)) {
            map.put(i, size);
            size++;
        }

        if (verbosity >= 3) {
            mainLog.println("unknown states mapping = "+map);
        }

        DoubleMatrix2D matrix = newMatrix(size, size);
        DoubleMatrix2D B = newMatrix(size, 1);

        for (Map.Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < getNumChoices(from); c++) {
                for (Iterator<Map.Entry<Integer, Double>> it = this.getTransitionsIterator(from, c); it.hasNext(); ) {
                    Map.Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();

                    if (unknown.get(to)) {
                        matrix.setQuick(map.get(from), map.get(to), matrix.getQuick(map.get(from), map.get(to)) + probMove.getValue());
                    } else {
                        double value = B.getQuick(map.get(from), 0);
                        value -= probMove.getValue() * (Double) knownValues.getValue(to);
                        B.setQuick(map.get(from), 0, value);
                    }

                }
            }
        }

        if (verbosity >= 2) {
            //mainLog.println("A = \n"+matrix.toString());
        }

        for (int i = 0; i < size; i++) {
            matrix.setQuick(i,i,  matrix.getQuick(i, i) - 1.0);
        }

        if (verbosity >= 2) {
            //mainLog.println("A - I = \n" + matrix.toString());
            //mainLog.println("B = \n" + B.toString());
        }

        return new Pair<DoubleMatrix2D, DoubleMatrix2D>(matrix, B);
    }

        public Set<Integer> getProbStatesSuccessors(final int probState, final APElement symbol) {
            Set<Integer> result = new HashSet<Integer>();

            for(int c = 0; c < getNumChoices(probState); c++){
                APElement ap = (APElement) getAction(probState, c);
                //mainLog.println("ap = " + ap.toString() + ", ap.IntegerList()" + ap.IntegerList());
                String a1 = actions.getAP(ap.IntegerList().get(0));
                String a2 = actions.getAP(symbol.IntegerList().get(0));
                //mainLog.println("a1 = " + a1 + ", a2 = " + a2 + ", a1.split()[0] = " + a1.split("_")[0]);

                if(a1.split("_")[0].equals( a2.split("_")[0])){
                    result.addAll(getProbStatesSuccessors(probState,c));
                    //mainLog.println("result = " + result);
                }
            }
            return result;
        }

//        public Set<Integer> getProbStatesSuccessors(final int probState) {
//            Set<Integer> result = new HashSet<Integer>();
//            for(int c = 0; c < getNumChoices(probState); c++) {
//                APElement symbol = (APElement) getAction(probState, c);
//                result.addAll(getProbStatesSuccessors(probState, symbol));
//            }
//            return result;
//        }

        public Set<Integer> getProbStatesSuccessors(final int probState) {
            Set<Integer> result = new HashSet<Integer>();

            for(Iterator<Integer> probSucc = getSuccessorsIterator(probState); probSucc.hasNext();) {
                result.add(probSucc.next());
            }

            return result;
        }

    public Set<Integer> getProbStatesLetterSuccessors(final int probState, final int letter) {
        Set<Integer> result = new HashSet<Integer>();

        for(Iterator<Integer> probSucc = getSuccessorsIterator(probState); probSucc.hasNext();) {
            int succ = probSucc.next();
            if(getDTMCState(succ) == letter){
                result.add(succ);
            }
        }
        return result;
    }
        public Set<Integer> getProbStatesSuccessors(final int probState, final int symbol) {
            Set<Integer> result = new HashSet<Integer>();

            for (Iterator<Integer> probSucc = getSuccessorsIterator(probState, symbol); probSucc.hasNext();) {
                Integer succState = probSucc.next();
                result.add(succState);
            }

            return result;
        }

    public Set<Integer> getProbStatesSuccessors(final int probState, List<APElement> word) {
        Set<Integer> current = new LinkedHashSet<Integer>();
        current.add(probState);

        for (APElement symbol : word) {
            Set<Integer> next = new LinkedHashSet<Integer>();
            for (int state : current) {
                next.addAll(getProbStatesSuccessors(state, symbol));
            }
            current = next;
        }

        return current;
    }


    public Set<Integer> getProbStatesSuccessors(Set<Integer> probStates, List<APElement> word) {
        Set<Integer> result = new HashSet<Integer>();
        for (Integer probState : probStates) {
            result.addAll(getProbStatesSuccessors(probState, word));
        }
        return result;
    }


    public int[] getInvMap() {
        return invMap;
    }

    public DTMCSimple projectToDTMC(){
        DTMCSimple dtmc = new DTMCSimple(this.getNumStates());
        for (int from = 0; from < this.getNumStates(); from++) {
            for (int c = 0; c < getNumChoices(from); c++) {
                for (Iterator<Map.Entry<Integer, Double>> it = this.getTransitionsIterator(from, c); it.hasNext(); ) {
                    Map.Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();
                    double value = dtmc.getTransitions(from).get(to);
                    dtmc.setProbability(from, to, value + probMove.getValue());
                }
            }
        }
        return dtmc;
    }
}



