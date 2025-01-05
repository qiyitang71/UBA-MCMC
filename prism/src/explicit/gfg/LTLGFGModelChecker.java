package explicit.gfg;

import java.io.PrintStream;
import java.util.*;
import java.util.Map.Entry;

import acceptance.AcceptanceOmega;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import automata.DA;
import automata.DASimplifyAcceptance;
import common.IterableBitSet;
import common.IterableStateSet;
import common.PathUtil;
import common.StopWatch;
import cern.colt.function.IntIntDoubleFunction;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.QRDecomposition;
import explicit.*;
import explicit.uba.LTL2UBA;
import explicit.uba.SharedWord;
import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2ba.MyBitSet;
import jltl2dstar.GFG;
import jltl2dstar.NBA;
import jltl2dstar.NBA_State;
import parser.ast.Expression;
import parser.type.TypeDouble;
import prism.ModelType;
import prism.Pair;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import prism.PrismSettings;
import prism.PrismUtils;

public class LTLGFGModelChecker extends PrismComponent
{
    // TODO: stub
    public class LTLProduct {
        private DTMCGFGProduct prod;

        public LTLProduct(DTMCGFGProduct prod) {
            this.prod = prod;
        }

        public DTMCGFGProduct getProduct() {
            return prod;
        }

        public HashMap<Integer, Set<APElement>> getAccEdges() {
            return prod.getAccEdges();
        }

        public StateValues projectToOriginalModel(DTMC model, StateValues prodValues) throws PrismLangException {
            StateValues result = new StateValues(TypeDouble.getInstance(), 1);

            for (int i : prod.getInitialStates()) {
                int dtmcState = prod.getDTMCState(i);

                Double value = (Double) result.getValue(dtmcState);
                if (value == null)
                    value = 0.0;

                value += (Double) prodValues.getValue(i);
                result.setDoubleValue(dtmcState, value);
            }

            return result;
        }
    };

    private boolean isUFA = false;

    /** The underlying model checker */
    private ProbModelChecker mc;

    private int verbosity = 0;

    private boolean sanityCheck = false;

    private HashMap<Integer, BitSet> bisMap = new HashMap<>();

    public LTLGFGModelChecker(ProbModelChecker mc)
    {
        super(mc);
        this.mc = mc;
        verbosity = getSettings().getInteger(PrismSettings.PRISM_UBA_VERBOSITY);
        sanityCheck = getSettings().getBoolean(PrismSettings.PRISM_UBA_SANITY);
    }

    /**
     * Construct an unambiguous Buchi automaton (UBA) for an LTL formula, having first extracted maximal state formulas
     * and model checked them with the passed in model checker. The maximal state formulas are assigned labels
     * (L0, L1, etc.) which become the atomic propositions in the resulting DA. BitSets giving the states which
     * satisfy each label are put into the vector {@code labelBS}, which should be empty when this function is called.
     * <br>
     *
     * @param mc a ProbModelChecker, used for checking maximal state formulas
     * @param model the model
     * @param expr a path expression, i.e. the LTL formula
     * @param labelBS empty vector to be filled with BitSets for subformulas
     * @return the unambiguous NBA
     */
    public GFG constructUBAForLTLFormula(ProbModelChecker mc, DTMC model, Expression expr, Vector<BitSet> labelBS) throws PrismException
    {
        Expression ltl;
        GFG uba;
        long time;

        if (Expression.containsTemporalRewardBounds(expr)) {
            throw new PrismException("Can not handle reward bounds via deterministic automata.");
        }

        if (Expression.containsTemporalTimeBounds(expr)) {
            if (model.getModelType().continuousTime()) {
                throw new PrismException("Automaton construction for time-bounded operators not supported for " + model.getModelType()+".");
            }
        }

        if (Expression.isHOA(expr)) {
            LTL2GFG ltl2uba = new LTL2GFG(this);
            time = System.currentTimeMillis();
            mainLog.println("Parsing and constructing HOA automaton for "+expr);
            Vector<Expression> apExpressions = new Vector<Expression>();
            uba = ltl2uba.fromExpressionHOA(expr, PathUtil.getDirectoryForRelativePropertyResource(mc.getModulesFile(), mc.getPropertiesFile()), apExpressions);

            mainLog.println("Determining states satisfying atomic proposition expressions of the automaton...");
            for (int i=0; i<uba.getAPSet().size(); i++) {
                Expression label = apExpressions.get(i);
                label.typeCheck();
                if(mc.checkExpression(model, label, null) != null) {
                    BitSet labelStates = mc.checkExpression(model, label, null).getBitSet();
                    labelBS.add(labelStates);
                }else{
                    labelBS.add(new BitSet());
                }
                //uba.getAPSet().renameAP(i, "L"+i);//do not rename
            }
        } else {
            mainLog.println("\nLTL to GFG not supported");

            // Model check maximal state formulas
            LTLModelChecker ltlMC = new LTLModelChecker(this);
            ltl = ltlMC.checkMaximalStateFormulas(mc, model, expr.deepCopy(), labelBS);

            // Convert LTL formula to UBA
            mainLog.println("\nBuilding unambiguous Buchi automaton (for " + ltl + ")...");
            time = System.currentTimeMillis();
            LTL2GFG ltl2uba = new LTL2GFG(this);
            uba = ltl2uba.convertLTLFormulaToUBA(ltl, mc.getConstantValues());
        }
        mainLog.println("GFG has " + uba.size() + " states.");
        //checkForCanonicalAPs(uba.getAPSet(), labelBS.size());
        time = System.currentTimeMillis() - time;
        mainLog.println("Time for GFG translation: " + time / 1000.0 + " seconds.");
        // If required, export UBA
        if (settings.getExportPropAut()) {
            mainLog.println("Exporting UBA to file \"" + settings.getExportPropAutFilename() + "\"...");
            PrintStream out = PrismUtils.newPrintStream(settings.getExportPropAutFilename());
            uba.print(out, settings.getExportPropAutType());
            out.close();
        }

        /*
        isUFA = checkUpwardClosedness(uba) && !settings.getBoolean(PrismSettings.PRISM_UBA_PURE);
        if(isUFA) {
            mainLog.println("GFG is actually an (upward-closed) UFA");
        }*/

        return uba;
    }

    /**
     *
     *
     */

    public Integer findExtensionDFS(int probState, int cur, List<Integer> word, Map<Integer, Set<Integer>> succForZ) {
        return null;
    }

    /**
     * Returns true, if every final state in an NBA has an edge to a final state for every symbol
     * @param uba the UBA, that should be checked
     * @return true, if every final state an edge to a final state for every symbol
     */
    public boolean checkUpwardClosedness(GFG uba)
    {
        for(int finalState : uba.getFinalStates()) {
            //Check whether all final states have an edge to a final state for every symbol
            NBA_State st = uba.get(finalState);
            for(Iterator<APElement> it = uba.getAPSet().elementIterator(); it.hasNext();) {
                //Check, whether there exists an edge to a final state
                if(!st.getEdge(it.next()).intersects(uba.getFinalStates())) {
                    return false;
                }
            }
        }
        return true;
    }


    public LTLProduct constructProduct(DTMC model, GFG uba, Vector<BitSet> labelBS, BitSet statesOfInterest) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("computing GFG-DTMC product");
        DTMCGFGProduct prod = new DTMCGFGProduct(this, model, uba, labelBS, statesOfInterest);

        if (settings.isExportDTMCUBAProduct()) {
            mainLog.println("Exporting DTMC GFG product to file \"" + settings.getExportDTMCUBAProductFilename() + "\" …");
            prod.exportToDotFile(settings.getExportDTMCUBAProductFilename());
            for(int i = 0; i < prod.getNumStates(); i++) {
                System.out.println(i + " -> (" + prod.getDTMCState(i) + "," + prod.getUBAState(i) + ")");
            }
        }

        LTLProduct ltlProd = new LTLProduct(prod);

        timer.stop(" (product has "+prod.getNumStates()+" states)");
        return ltlProd;
    }

    public StateValues computeValues(LTLProduct product) throws PrismException
    {
        StateValues result = new StateValues(TypeDouble.getInstance(), product.getProduct());
        // the states in positive MCCs
        BitSet S_MCC_pos = new BitSet();

        StopWatch timer = new StopWatch(mainLog);

        BitSet unknownStates = new BitSet();
        if (!isUFA) {
            timer.start("computing MCC in GFG-DTMC product");
            List<BitSet> mccs = computeMCCs(product.getProduct());
            int mccCount = mccs.size();
            timer.stop(" (found "+mccCount+" non-trivial MCCs)");

            int posMccCount = 0;
            int mccIndex = 0;
            timer.start("computing MCC probabilities for positive MCCs");
            for (BitSet mcc : mccs) {
                mccIndex++;

                if (verbosity >= 2) {
                    mainLog.println("\n MCC " + mccIndex + ": " + mcc);
                } else if (verbosity >= 1) {
                    mainLog.println("\nMCC " + mccIndex + " has " + mcc.cardinality() + " states");
                }

                boolean isPositive = handleMCC(product, mccIndex, mcc, result, S_MCC_pos);;
                /*
                if (getSettings().getBoolean(PrismSettings.PRISM_UBA_POWER)) {
                    isPositive = handleMCCPower(product, mccIndex, mcc, result, S_MCC_pos);
                } else {
                    isPositive = handleMCC(product, mccIndex, mcc, result, S_MCC_pos);
                }*/
                if (isPositive) posMccCount++;
            }
            timer.stop(" ("+posMccCount+" positive MCCs, known probabilities for "+S_MCC_pos.cardinality()+" states)");

            if (verbosity >= 2) {
                mainLog.println("Partial result (all MCCs):");
                result.printFiltered(mainLog, S_MCC_pos, true, false, false, true);
            }
            //*******************

            timer.start("determining states with probability zero");
            PredecessorRelation pre = product.getProduct().getPredecessorRelation(this, true);
            BitSet S_prePositiveMCC = pre.calculatePreStar(null, S_MCC_pos, new BitSet());
            BitSet S_zero = new BitSet();
            S_zero.flip(0, product.getProduct().getNumStates());
            S_zero.andNot(S_prePositiveMCC);

            unknownStates.flip(0, product.getProduct().getNumStates());
            unknownStates.andNot(S_zero);
            unknownStates.andNot(S_MCC_pos);

            // we are only interested in the prob states

            timer.stop(" ("+S_zero.cardinality()+" zero prob. states, "+unknownStates.cardinality()+" remaining unknown)");
        }

        /*if (isUFA) {
            for (int state : IterableBitSet.getSetBits(product.getFinalStates())) {
                // we set the value to 1 for all the known states ( = final states )
                result.setDoubleValue(state, 1.0);
            }

            // we calculate all states that can reach a final states, as those have P>0
            PredecessorRelation pre = product.getProduct().getPredecessorRelation(this, true);
            BitSet preFinal = pre.calculatePreStar(null, product.getFinalStates(), null);
            preFinal.andNot(product.getFinalStates());

            // we determine the reachability probability for reaching the final states
            // for those states that can reach them
            if (getSettings().getString(PrismSettings.PRISM_UBA_REACH_METHOD).equals("GS")) {
                computeReachabilityGS(product, result, product.getFinalStates(), preFinal);
            } else {
                computeReachability(product, result, product.getFinalStates(), preFinal);
            }
        } else {

        }*/

        if (!unknownStates.isEmpty()) {
            if (getSettings().getString(PrismSettings.PRISM_UBA_REACH_METHOD).equals("GS")) {
                computeReachabilityGS(product, result, S_MCC_pos, unknownStates);
            } else {
                computeReachability(product, result, S_MCC_pos, unknownStates);
            }
        }
        if (verbosity >= 2) {
            mainLog.println("Result:");
            result.printFiltered(mainLog, null, true, false, false, true);
        }

        return result;
    }

    private boolean handleMCC(LTLProduct product, int mccIndex, BitSet mcc, StateValues result, BitSet knownValues) throws PrismException
    {
        boolean isPositive = checkIsMCCPositive(product, mccIndex, mcc);
        if (!isPositive)
            return false;

        int probState = mcc.nextSetBit(0);
        Set<Integer> cut = generateCut(product.getProduct(), probState, mcc);
        if (verbosity >= 2) {
            mainLog.println("Cut: " + cut);

            for (Integer cutState : cut) {
                DTMCGFGProduct prod = product.getProduct();
                mainLog.print("(" + prod.getDTMCState(cutState) + "," + prod.getUBAState(cutState) +  ")" );
            }
            mainLog.println();
        }

        computeMCCProbs(product, mccIndex, mcc, cut, result, knownValues);

        return isPositive;
    }

    //TODO: check it has positive solution with x(c) = x(d) for all c ~ d
    private boolean checkIsMCCPositive(LTLProduct product, int mccIndex, BitSet mcc) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("checking whether MCC " + mccIndex + " is positive");
        DTMCGFGProduct prod = product.getProduct();
        HashMap<Integer, Set<APElement>>  accEdges = product.getAccEdges();
        // first, check that MCC intersects accepting edges
        boolean acc = false;
        for (int state : new IterableStateSet(mcc, prod.getNumStates())) {
            Set<APElement> set = accEdges.get(state);
            if(set == null) continue;
            for(APElement ap: set){
                for(int i = 0; i < prod.getNumChoices(state);i++) {
                    if(!prod.getAction(state, i).equals(ap)) continue;
                    if(prod.allSuccessorsInSet(state, i, mcc)){
                        acc = true;
                        break;
                    }
                }
            }
        }

        if(!acc){
            if (verbosity >= 1) timer.stop(" (MCC is zero, no accepting edges)");
            return false;
        }

        StopWatch timerMatrix = new StopWatch(mainLog);
        timerMatrix.start("building positivity matrix");
        DoubleMatrix2D mccMatrix = prod.positivityMatrixForMCC(mcc, null, true);
        assert(mccMatrix.rows() == mccMatrix.columns());
        assert(mccMatrix.rows() == mcc.cardinality());
        if (verbosity >= 1) {
            timerMatrix.stop();
        }

        int rows = mccMatrix.rows();

        boolean positive;
        String posMethod = getSettings().getString(PrismSettings.PRISM_UBA_POS_METHOD);
        switch (posMethod) {
            case "SVD": {
                Algebra algebra = new Algebra();
                int rank = algebra.rank(mccMatrix);
                if (rank < rows-1 || rank > rows) {
                    throw new PrismException("Strange things are going on (rank = " + rank + ", matrix has size " + rows + "x" + rows +")..");
                }

                if (verbosity >= 1) mainLog.println("Rank of MCC " + mccIndex + " = " + rank + ", full rank is " + rows);
                positive = rank < rows;
                break;
            }
            case "QR": {
                QRDecomposition qr = new QRDecomposition(mccMatrix);
                positive = !hasFullRankWithTolerance(qr);
                break;
            }
            default:
                throw new PrismException("Unknown UBA method for checking MCC positivity: "+posMethod);
        }

        if (verbosity >= 1) timer.stop(" (MCC is " + (positive ? "positive" : "zero")+", checked via " + posMethod + ")");

        return positive;
    }

    private void computeMCCProbs(LTLProduct product, int mccIndex, BitSet mcc, Set<Integer> cut, StateValues result, BitSet knownValues) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("building matrix for positive MCC probabilities, MCC " + mccIndex);

        LinkedHashMap<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
        DoubleMatrix2D matrix = product.getProduct().valueMatrixForSCC(mcc, cut, map);
        assert(matrix.rows() + 1 == matrix.columns());
        assert(matrix.rows() < mcc.cardinality() + 1);

        Algebra algebra = new Algebra();
        DoubleMatrix2D B = new DenseDoubleMatrix2D(matrix.rows(), 1);
        B.setQuick(matrix.rows()-1, 0, 1);
        //if (verbosity >= 2) mainLog.println("B = \n" + B);

        if (verbosity >= 1) timer.stop();

        if (verbosity >= 1) mainLog.println("Solving MCC " + mccIndex + " probability values...");
        timer.start("solving equation system for positive MCC " + mccIndex);
        DoubleMatrix2D solution = algebra.solve(matrix, B);
        assert (solution.rows() == matrix.rows());
        assert (solution.columns() == 1);

        double cutSum = 0.0;
        double minValue = 0.0;
        double maxValue = 0.0;
        boolean first = true;
        for (Entry<Integer, Integer> entry : map.entrySet()) {
            int productIndex = entry.getKey();
            int solutionIndex = entry.getValue();

            double value = solution.getQuick(solutionIndex, 0);
            if (sanityCheck && value == 0.0) {
                throw new PrismException("Something strange going on (probability in positive MCC is zero for state "+productIndex+")");
            }
            result.setDoubleValue(productIndex, value);
            assert(!knownValues.get(productIndex));
            knownValues.set(productIndex);

            if (cut.contains(productIndex)) {
                cutSum += value;
            }
            if (first) {
                minValue = maxValue = value;
                first = false;
            } else {
                minValue = Math.min(minValue, value);
                maxValue = Math.max(maxValue, value);
            }
        }
        if (verbosity >= 1) timer.stop();

        if (verbosity >= 1) {
            mainLog.println("Sum of probabilities for the cut C = "+cutSum + " for MCC "+mccIndex);
            mainLog.println("Probabilities in MCC " +mccIndex + " are in the range ["+minValue+","+maxValue+"]");
        }
    }

    private boolean hasFullRankWithTolerance(QRDecomposition qr)
    {
        DoubleMatrix2D Rdiag = qr.getR();
        //mainLog.println(Rdiag);

        for (int j = 0; j < Rdiag.columns(); j++) {
            double value = Rdiag.getQuick(j,j);
            // mainLog.println(" Rdiag["+j+"]  = " + value);
            if (value < 1E-10) return false;
        }
        return true;
    }

    private void computeReachability(LTLProduct product, StateValues result, BitSet s_MCC_pos, BitSet unknownStates) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("solving linear equation system for remaining unknown states");
        LinkedHashMap<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
        Pair<DoubleMatrix2D, DoubleMatrix2D> lgs = product.getProduct().reachabilityEquationSystem(s_MCC_pos, unknownStates, result, map);

        Algebra algebra = new Algebra();
        DoubleMatrix2D solution = algebra.solve(lgs.first, lgs.second);

        if (verbosity >= 2) mainLog.println("Solution = " + solution);

        for (Entry<Integer, Integer> entry : map.entrySet()) {
            int productIndex = entry.getKey();
            int solutionIndex = entry.getValue();

            double value = solution.getQuick(solutionIndex, 0);
            if (value == 0.0) {
                throw new PrismException("Something strange going on (probability during reachability is zero for state "+productIndex+")");
            }
            result.setDoubleValue(productIndex, value);

        }
        timer.stop();
    }

    private void computeReachabilityGS(LTLProduct product, StateValues result, BitSet s_MCC_pos, BitSet unknownStates) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("solving linear equation system for remaining unknown states (GS)");

        BitSet yes = new BitSet(); // empty
        BitSet no  = new BitSet(); // empty
        BitSet known = new BitSet();
        known.set(0, product.getProduct().getNumStates());
        known.andNot(unknownStates);

        // while the product is not actually a DTMC (probabilities can sum to > 1),
        // the Gauss-Seidel methods of DTMCModelChecker does not care
        // we just feed the results for the known states
        double[] values = result.getDoubleArray();
//        MDPModelChecker productMdp = new MDPModelChecker(mc);
//      ModelCheckerResult r = productMdp.computeReachProbs(product.getProduct(),);
        DTMCModelChecker productMc = new DTMCModelChecker(mc);
        ModelCheckerResult r = productMc.computeReachProbsGaussSeidel(product.getProduct().projectToDTMC(), no, yes, values, known);

        if (verbosity >= 2){
            mainLog.print("Solution = ");
            for(int  i= 0; i < values.length; i++){
                mainLog.print(values[i] + ",");
            }
            mainLog.println();
        }

        for (int i : new IterableBitSet(unknownStates)) {
            double value = r.soln[i];
            if (value == 0.0) {
                throw new PrismException("Something strange going on (probability during reachability is zero for state "+i+")");
            }
            result.setDoubleValue(i, value);
        }
        timer.stop();
    }

    private void generateBisimilarMap (DTMCGFGProduct product, final BitSet mcc){
        bisMap = new HashMap<>();        //int s = mcc.nextSetBit(0);

        for (int s : IterableBitSet.getSetBits(mcc)) {
            for(int c = 0; c < product.getNumChoices(s); c++){
                if(product.allSuccessorsInSet(s,c,mcc) && product.getProbStatesSuccessors(s,c).size() > 1){
                    BitSet bset = new BitSet(product.getNumStates());
                    for(int t: product.getProbStatesSuccessors(s,c)){
                        bset.set(t);
                    }
                    for(int t: product.getProbStatesSuccessors(s,c)){
                        if(bisMap.containsKey(t)) break;
                        bisMap.put(t, bset);
                    }
                }
            }
        }
        if(verbosity >= 2) {
            mainLog.println("bisMap = " + bisMap);
            for(int i: bisMap.keySet()){
                mainLog.print(i + " == " );
                for(Integer j: IterableBitSet.getSetBits(bisMap.get(i))){
                    mainLog.print( j + "(" +product.getDTMCState(j)+ "," + product.getUBAState(j) + ")" );
                }
            }
        }
    }

        private Set<Integer> generateCut(DTMCGFGProduct product, final int probState, final BitSet mcc) throws PrismException  {
        StopWatch timer = new StopWatch(mainLog);

        if (!mcc.get(probState)) {
            throw new PrismException("Initial state for cut generation is not in MCC");
        }

        timer.start("generating cut");
        // tracks, for every state q the set delta(q, z)
        Map<Integer,Set<Integer>> succsForZ = new HashMap<Integer,Set<Integer>>();

        // partition the mcc
        generateBisimilarMap(product, mcc);

        // for every state, the state can be reached via z=epsilon
        for(Integer probState_: IterableBitSet.getSetBits(mcc)) {
            Set<Integer> reach = new HashSet<Integer>();
            reach.add(probState_);
            succsForZ.put(probState_, reach);
        }

        boolean found = true;
        int iterations = 0;

        long countVisited = 0;

        if (verbosity >= 2) mainLog.println("START WITH: (" + probState + "," + probState + ")");
        while (found) {
            //Search for an extension
            found = false;
            ++iterations;

            // the set of expanded states
            Set<ProductState> visited = new HashSet<ProductState>();

            // the BFS queue, consisting of a state in the self-product and the word used for reaching it
            Queue<Pair<ProductState, SharedWord<APElement>>> queue = new LinkedList<Pair<ProductState, SharedWord<APElement>>>();
            queue.add(new Pair<ProductState, SharedWord<APElement>>(new ProductState(probState, probState),
                    new SharedWord<APElement>()));
            while(!queue.isEmpty()) {
                Pair<ProductState, SharedWord<APElement>> entry = queue.poll();
                ProductState current = entry.first;
                if (!visited.contains(current)) {
                    SharedWord<APElement> word = entry.second;
                    Integer left = current.getFirstState();
                    Integer right = current.getSecondState();
                    if (!mcc.get(left) || !mcc.get(right)) {
                        // we don't want to leave the MCC
                        continue;
                    }

                    //mainLog.println("product.getNumChoices(left) = " + product.getNumChoices(left));

                    for(int c = 0; c < product.getNumChoices(left); c++) {
                        APElement ap = (APElement) product.getAction(left, c);
                        Set<Integer> leftSuccs = product.getProbStatesSuccessors(left,ap);
                        Set<Integer> rightSuccs;
                        if (left == right) {
                            rightSuccs = leftSuccs;
                        } else {
                            rightSuccs = product.getProbStatesSuccessors(right,ap);
                        }
                        //mainLog.println("leftSuccs = " + leftSuccs + "; rightSuccs = " + rightSuccs);

                        visited.add(current);
                        countVisited++;
                        for (Integer leftSucc : leftSuccs) {
                            for (Integer rightSucc : rightSuccs) {
                                if (!mcc.get(leftSucc) || !mcc.get(rightSucc)) {
                                    // we don't want to leave the SCC
                                    //mainLog.println("leftSucc = " + leftSucc + "; rightSucc = " + rightSucc);
                                    continue;
                                }

                                // if the move has a different symbol...
//                                if (product.getDTMCState(leftSucc) != product.getDTMCState(rightSucc)) {
//                                    continue;
//                                }

                                // ... we extend the word and add to the BFS queue
                                SharedWord<APElement> curWord = word.append(ap);

                                boolean isLeftBis = leftSucc == probState || ( (bisMap.get(probState) == null)? false: bisMap.get(probState).get(leftSucc));
                                boolean isRightBis = rightSucc == probState || ((bisMap.get(probState) == null)? false: bisMap.get(probState).get(rightSucc));


                                if ((!isLeftBis  &&(rightSucc == probState && !succsForZ.get(leftSucc).isEmpty()) ||
                                                (!isRightBis  &&leftSucc == probState && !succsForZ.get(rightSucc).isEmpty()))) {
                                    //Found extension
                                    found = true;
                                    if (verbosity >= 2) {
                                        String wordInfo = curWord.size() > 3 ? "length(y) = " + curWord.size() : "y = " + curWord;
                                        mainLog.println("FOUND EXTENSION: (" + leftSucc + "," + rightSucc + ") with word " + wordInfo + ";");
                                        String cutInfo;
                                        if (leftSucc == probState) {
                                            if (succsForZ.get(rightSucc).size() > 3) {
                                                cutInfo = succsForZ.get(rightSucc).size() + " states";
                                            } else {
                                                cutInfo = succsForZ.get(rightSucc).toString();
                                            }
                                        } else {
                                            if (succsForZ.get(leftSucc).size() > 3) {
                                                cutInfo = succsForZ.get(leftSucc).size() + " states";
                                            } else {
                                                cutInfo = succsForZ.get(leftSucc).toString();
                                            }
                                        }
                                        mainLog.println("add " + cutInfo + " to cut");
                                    }

                                    if (verbosity >= 2) {
                                        mainLog.println("Visited " + visited.size() + " states for finding an extension.");
                                    }

                                    visited.clear();
                                    queue.clear();
                                    //Update succsforZ
                                    Map<Integer, Set<Integer>> succsForZ_ = new HashMap<Integer, Set<Integer>>(succsForZ.size());
                                    List<APElement> y = curWord.getWord();
                                    for (Integer probState_ : IterableBitSet.getSetBits(mcc)) {
                                        Set<Integer> newReach = new HashSet<Integer>();
                                        Set<Integer> succs = product.getProbStatesSuccessors(probState_, y);
                                        for (Integer succState : succs) {
                                            if(succsForZ.get(succState) != null) {
                                                newReach.addAll(succsForZ.get(succState));
                                            }
                                        }
                                        succsForZ_.put(probState_, newReach);
                                    }

                                    if (sanityCheck) {
                                        if (!succsForZ_.get(probState).containsAll(succsForZ.get(probState))) {
                                            throw new PrismException("The potential cut has got smaller.");
                                        }
                                    }

                                    succsForZ = succsForZ_;
                                } else {
                                    // we have to continue, add new product state and curWord to BFS queue
                                    queue.add(new Pair<ProductState, SharedWord<APElement>>(new ProductState(leftSucc, rightSucc),
                                            curWord));
                                }
                            }
                        }
                    }
                }
            }
            if(!found) {
                if (verbosity >= 2) {
                    mainLog.println("DID NOT FIND EXTENSION, but visited " + visited.size() + " states");
                }
            }
        }

        Set<Integer> tempC = succsForZ.get(probState);

        if (verbosity >= 1) {
            timer.stop(" ("+iterations+" iterations, " + countVisited + " extension checks, cut tempC has " + tempC.size() + " states)");
        }

        //get rid of bisimilar states
            Set<Integer> C = new HashSet<Integer>();
            for (int c: tempC){
                Set<Integer> setC = new HashSet<>();
                boolean isContained = false;
                if(bisMap.get(c) != null){
                    for(int s: IterableBitSet.getSetBits(bisMap.get(c)) ) {
                        if (C.contains(s)) {
                            isContained = true;
                            break;
                        }
                    }
                }
                if(!isContained)
                    C.add(c);
            }

            if (verbosity >= 1) {
                timer.stop(" ("+iterations+" iterations, " + countVisited + " extension checks, cut C has " + C.size() + " states)");
            }

            return C;
    }


    /**
     * Validates that the atomic propositions
     * conform to the standard values that PRISM expects:
     *   L0, ..., Ln-1 (in arbitrary order)
     * if there are {@code n} expected atomic propositions.
     * <br/>
     * The automaton may actually have less atomic propositions than expected,
     * e.g., if the given atomic proposition does not influence the acceptance
     * of a run in the automaton.
     * <br/>
     * If there is an error, throws a {@code PrismException} detailing the problem.
     * @param expectedNumberOfAPs the expected number of atomic propositions
     */
    private void checkForCanonicalAPs(APSet aps, int expectedNumberOfAPs) throws PrismException {
        BitSet seen = new BitSet();
        for (String ap : aps) {
            if (!ap.substring(0,1).equals("L")) {
                throw new PrismException("In UBA, unexpected atomic proposition "+ap+", expected L0, L1, ...");
            }
            try {
                int index = Integer.parseInt(ap.substring(1));
                if (seen.get(index)) {
                    throw new PrismException("In UBA, duplicate atomic proposition "+ap);
                }
                if (index < 0) {
                    throw new PrismException("In UBA, unexpected atomic proposition "+ap+", expected L0, L1, ...");
                }
                if (index >= expectedNumberOfAPs) {
                    throw new PrismException("In UBA, unexpected atomic proposition "+ap+", expected highest index to be "+(expectedNumberOfAPs-1));
                }
                seen.set(index);
            } catch (NumberFormatException e) {
                throw new PrismException("In UBA, unexpected atomic proposition "+ap+", expected L0, L1, ...");
            }
        }
        // We are fine with an empty apList or an apList that lacks some of the expected Li.
    }

    private List<BitSet> computeSCCs(DTMCGFGProduct prod) throws PrismException {
        SCCConsumerStore sccs = new SCCConsumerStore();
        SCCComputer sccc = SCCComputer.createSCCComputer(this, prod, sccs);
        sccc.computeSCCs();
        return sccs.getSCCs();
    }

    private List<BitSet> computeMCCs(DTMCGFGProduct prod) throws PrismException {
//        ECComputer ecs = ECComputer.createECComputer(this, prod);
//        ecs.computeMECStates();

        ECComputerDefault ecs = new ECComputerDefault(this, prod);
        ecs.computeMECStates();
        return ecs.getMECStates();
    }

    public BitSet computeAccBSCCs(Model model, final BitSet finalStates) throws PrismException
    {
        final BitSet result = new BitSet();

        // Compute bottom strongly connected components (BSCCs)
        // and check using the following SCCConsumerBSCCs:
        SCCConsumerBSCCs sccConsumer = new SCCConsumerBSCCs() {

            @Override
            public void notifyNextBSCC(BitSet bscc) {
                if (bscc.intersects(finalStates)) {
                    result.or(bscc);
                }
            }
        };

        SCCComputer sccComputer = SCCComputer.createSCCComputer(this, model, sccConsumer);
        sccComputer.computeSCCs();

        // now, the result is ready
        return result;
    }

    public StateValues computeWithDA(DTMC model, NBA uba, Vector<BitSet> labelBS, BitSet statesOfInterest) throws PrismException {
        LTLModelChecker mcLtl = new LTLModelChecker(this);
        DA<BitSet, ? extends AcceptanceOmega> da = uba.createPrismDAFromDeterministicNBA();

        da = DASimplifyAcceptance.simplifyAcceptance(this, da, AcceptanceType.RABIN, AcceptanceType.REACH);
        mainLog.println("The UBA is actually deterministic, we continue with standard model checking with a "+da.getAutomataType()+" with "+da.size()+" states.");

        explicit.LTLModelChecker.LTLProduct<DTMC> product = mcLtl.constructProductModel(da, model, labelBS, statesOfInterest);

        // Find accepting states + compute reachability probabilities
        BitSet acc;
        if (product.getAcceptance() instanceof AcceptanceReach) {
            mainLog.println("\nSkipping BSCC computation since acceptance is defined via goal states...");
            acc = ((AcceptanceReach)product.getAcceptance()).getGoalStates();
        } else {
            mainLog.println("\nFinding accepting BSCCs...");
            acc = mcLtl.findAcceptingBSCCs(product.getProductModel(), product.getAcceptance());
        }
        mainLog.println("\nComputing reachability probabilities...");
        DTMCModelChecker mcProduct = new DTMCModelChecker(this);
        mcProduct.inheritSettings(mc);
        StateValues probsProduct = StateValues.createFromDoubleArray(mcProduct.computeReachProbs(product.getProductModel(), acc).soln, product.getProductModel());

        // Mapping probabilities in the original model
        StateValues probs = product.projectToOriginalModel(probsProduct);
        probsProduct.clear();

        return probs;
    }

    public StateValues checkProbPathFormulaLTL(DTMC model, Expression expr, boolean qual, BitSet statesOfInterest) throws PrismException
    {
        StateValues probs;
        LTLProduct product;

        Vector<BitSet> labelBS = new Vector<BitSet>();

        GFG uba = constructUBAForLTLFormula(mc, model, expr, labelBS);

        /*
        if (!getSettings().getBoolean(PrismSettings.PRISM_UBA_PURE) &&
                uba.isDeterministic()) {
            return computeWithDA(model, uba, labelBS, statesOfInterest);
        }
        */

        product = constructProduct(model, uba, labelBS, statesOfInterest);
        probs = computeValues(product);

        StateValues result = product.projectToOriginalModel(model, probs);

        return result;
    }
}
