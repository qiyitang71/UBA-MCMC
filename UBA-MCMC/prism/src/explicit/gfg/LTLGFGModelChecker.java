package explicit.gfg;

import java.io.PrintStream;
import java.util.*;
import java.util.Map.Entry;

import acceptance.AcceptanceOmega;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import automata.DA;
import automata.DASimplifyAcceptance;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
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
import explicit.uba.SharedWord;
import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2dstar.GFG;
import jltl2dstar.NBA;
import jltl2dstar.NBA_State;
import parser.ast.Expression;
import parser.type.TypeDouble;
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

    private HashMap<Integer, Integer> equivDTMC = new HashMap<>();
    private HashMap<Integer, BitSet> equivGFG = new HashMap<>();

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


                /// /////////////////CHEAT!

//                BitSet labelStates = new BitSet(model.getNumStates());
//
//                for(int s = 0; s < model.getNumStates(); s++){
//                    if (s%4 == i/2){
//                        labelStates.set(s);
//                    }
//                }
//                labelBS.add(labelStates);

                /// /////////////////

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
        equivGFG = prod.getEquivGFG();
        equivDTMC = prod.getEquivDTMCStates();
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

                boolean isPositive;//handleMCC(product, mccIndex, mcc, result, S_MCC_pos);;

                if (getSettings().getBoolean(PrismSettings.PRISM_GFG_POWER)) {
                    isPositive = handleMCCPower(product, mccIndex, mcc, result, S_MCC_pos);
                } else {
                    isPositive = handleMCCQuotient(product, mccIndex, mcc, result, S_MCC_pos);//handleMCCQuotient(product, mccIndex, mcc, result, S_MCC_pos);
                }
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

        boolean isPositive = checkIsMCCPositiveQuotient(product, mccIndex, mcc);//checkIsMCCPositive(product, mccIndex, mcc);
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

    private boolean handleMCCQuotient(LTLProduct product, int mccIndex, BitSet mcc, StateValues result, BitSet knownValues) throws PrismException
    {
//        if(!checkIsMCCPositiveQuotient(product, mccIndex, mcc)){
//            return false;
//        }

        StopWatch timer = new StopWatch(mainLog);
        timer.start("checking whether MCC " + mccIndex + " is positive");
        DTMCGFGProduct prod = product.getProduct();
        HashMap<Integer, Set<APElement>> accEdges = product.getAccEdges();

        if (verbosity >= 2){
            mainLog.println("MCC states:");

            for (int i : IterableBitSet.getSetBits(mcc)) {
                int ubaState = prod.getUBAState(i);
                int lmcState = prod.getDTMCState(i);
                mainLog.print("("+ubaState+", "+ lmcState + ") ");
            }
            mainLog.println();
        }


        // first, check that MCC intersects accepting edges
        boolean acc = false;
        for (int state : new IterableStateSet(mcc, prod.getNumStates())) {
            Set<APElement> set = accEdges.get(state);
            if (set == null) continue;
            for (APElement ap : set) {
                for (int i = 0; i < prod.getNumChoices(state); i++) {
                    if (!prod.getAction(state, i).equals(ap)) continue;
                    if (prod.allSuccessorsInSet(state, i, mcc)) {
                        acc = true;
                        break;
                    }
                }
            }
        }

        if (!acc) {
            if (verbosity >= 1) timer.stop(" (MCC is zero, no accepting edges)");
            return false;
        }



        //FIX only handle LMCs with the same languages
        for(int idx = 0; idx < product.getProduct().getNumDtmcStates(); idx++){
            equivDTMC = new HashMap<>();
            equivDTMC.put(idx, 0);
        }

        StopWatch timerMatrix = new StopWatch(mainLog);
        timerMatrix.start("building positivity matrix");
        Map<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
        //DoubleMatrix2D sccMatrix = product.getProduct().positivityMatrixForMCC(mcc, map, false);


        Map<Integer, Integer> equivMap = new LinkedHashMap<Integer, Integer>();

        int size = 0;
        for (int i : IterableBitSet.getSetBits(mcc)) {
            int ubaState = prod.getUBAState(i);
            int lmcState = prod.getDTMCState(i);
            boolean isSet = false;
            for (Entry<Integer, Integer> states : map.entrySet()) {
                int s = states.getKey();
                if(((equivGFG.containsKey(ubaState) && equivGFG.get(ubaState).get(prod.getUBAState(s))) || prod.getUBAState(s) == ubaState) && (Objects.equals(equivDTMC.get(lmcState), equivDTMC.get(prod.getDTMCState(s))))){
                    mainLog.println("s = " + prod.getUBAState(s) + ", " + prod.getDTMCState(s) + ", i = " + ubaState + ", " + lmcState);

                    equivMap.put(i, map.get(s));
                    isSet = true;
                    break;
                }
            }
            if(!isSet){
                mainLog.println("new: i = " + ubaState + ", " + lmcState);
                map.put(i, size);
                equivMap.put(i, size);
                size++;
            }
        }

        if (verbosity >= 3) {
            mainLog.println("MCC mapping = " + map);
            for (Entry<Integer, Integer> states : map.entrySet()) {
                int prodstate = states.getKey();
                int ubaState = prod.getUBAState(prodstate);
                int lmcState = prod.getDTMCState(prodstate);
                mainLog.print("("+ ubaState + ", " + lmcState + ") ");
            }
            mainLog.println();
        }

        //DoubleMatrix2D matrix = new SparseDoubleMatrix2D(size,size);
        DoubleMatrix2D sccMatrix = new SparseDoubleMatrix2D(size,size);
        for (Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < prod.getNumChoices(from); c++) {
                for (Iterator<Entry<Integer, Double>> it = prod.getTransitionsIterator(from,c); it.hasNext(); ) {
                    Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();
                    if (mcc.get(to)) {
                        int matrixTo = equivMap.get(to);
                        sccMatrix.setQuick(map.get(from), matrixTo, sccMatrix.get(map.get(from),matrixTo)+ probMove.getValue());
                    }
                }
            }
        }

        if (verbosity >= 2) {

            mainLog.println("sccMatrix: ");
            for (int i = 0; i< sccMatrix.rows(); i++) {
                for (int j = 0; j < sccMatrix.columns(); j++) {
                    mainLog.print(sccMatrix.getQuick(i, j) + ", ");
                }
                mainLog.println();
            }
            mainLog.println();
        }

        assert(sccMatrix.rows() == sccMatrix.columns());
        assert(sccMatrix.rows() == mcc.cardinality());
        if (verbosity >= 1) {
            timerMatrix.stop();
        }

        //check positivity

        if (isEigenValueBig(sccMatrix)) {
            mainLog.println("MCC has eigenvalue bigger than 1, MCC is zero...");
            return false;
        }

        if (isEigenValueOne(sccMatrix)) {
            mainLog.println("has eigenvalue 1...");
        }else {
            mainLog.println("MCC does not have eigenvalue 1, MCC is zero...");
            return false;
        }

        int probState = mcc.nextSetBit(0);
        Set<Integer> cut = generateCut(product.getProduct(), probState, mcc);
        if (verbosity >= 2) {
            mainLog.println("Cut: " + cut);

            for (Integer cutState : cut) {
                mainLog.print("(" + prod.getDTMCState(cutState) + "," + prod.getUBAState(cutState) +  ")" );
            }
            mainLog.println();
        }


        //new matrix solve equation system

        // one row additionally for cut
        DoubleMatrix2D matrix = new SparseDoubleMatrix2D(size+1, size);

        for (int i=0; i < size +1; i++) {
            for(int j = 0; j < size; j++){
                matrix.setQuick(i, j, sccMatrix.getQuick(i, j));
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
            matrix.setQuick(size, equivMap.get(i), 1.0);
        }

        if (verbosity >= 2) {
            //mainLog.println("(A - I) + cut condition = \n"+matrix.toString());
        }

        //solve equation
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
        for (Entry<Integer, Integer> entry : equivMap.entrySet()) {
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

        //computeMCCProbs(product, mccIndex, mcc, cut, result, knownValues);

        return true;
    }


    private boolean handleMCCPower(LTLProduct product, int mccIndex, BitSet mcc, StateValues result, BitSet knownValues) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("checking whether SCC " + mccIndex + " is positive and computing eigenvector");
        // first, check that SCC intersects accepting edges
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
        Map<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
        //DoubleMatrix2D sccMatrix = product.getProduct().positivityMatrixForMCC(mcc, map, false);


        Map<Integer, Integer> equivMap = new LinkedHashMap<Integer, Integer>();

        int size = 0;
        for (int i : IterableBitSet.getSetBits(mcc)) {
            int ubaState = prod.getUBAState(i);
            int lmcState = prod.getDTMCState(i);
            boolean isSet = false;
            for (Entry<Integer, Integer> states : map.entrySet()) {
                int s = states.getKey();
                if(equivGFG.containsKey(ubaState) && equivGFG.get(ubaState).get(prod.getUBAState(s)) && (equivDTMC.get(lmcState) == equivDTMC.get(prod.getDTMCState(s)))){
                    equivMap.put(i, map.get(s));
                    isSet = true;
                    break;
                }
            }
            if(!isSet){
                map.put(i, size);
                equivMap.put(i, size);
                size++;
            }
        }

        if (verbosity >= 3) {
            mainLog.println("MCC mapping = "+map);
        }

        //DoubleMatrix2D matrix = new SparseDoubleMatrix2D(size,size);
        DoubleMatrix2D sccMatrix = new SparseDoubleMatrix2D(size,size);
        for (Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < prod.getNumChoices(from); c++) {
                for (Iterator<Entry<Integer, Double>> it = prod.getTransitionsIterator(from,c); it.hasNext(); ) {
                    Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();
                    if (mcc.get(to)) {
                        int matrixTo = equivMap.get(to);
                        sccMatrix.setQuick(map.get(from), matrixTo, sccMatrix.get(map.get(from),matrixTo)+ probMove.getValue());
                    }
                }
            }
        }

        assert(sccMatrix.rows() == sccMatrix.columns());
        assert(sccMatrix.rows() == mcc.cardinality());
        if (verbosity >= 1) {
            timerMatrix.stop();
        }

        // M' = (M+I) / 2 to avoid periodicity problem
        // do it using forEachNonZero to exploit sparseness of the matrix
        sccMatrix.forEachNonZero(new IntIntDoubleFunction() {
            @Override
            public double apply(int row, int column, double value)
            {
                return value / 2;
            }
        });
        // add I/2 to diagonal
        for (int i = 0; i < sccMatrix.rows(); i++) {
            sccMatrix.setQuick(i, i, sccMatrix.getQuick(i, i) + 0.5);
        }

        if (verbosity >= 2) {
            //mainLog.println("(M+I)/2 =" + sccMatrix);
        }

        // do iterations
        //int n = mcc.cardinality();
        int n = size;
        int maxIters = Integer.MAX_VALUE;//getSettings().getInteger(PrismSettings.PRISM_MAX_ITERS);
        boolean absolute = (getSettings().getString(PrismSettings.PRISM_TERM_CRIT).equals("Absolute"));
        double epsilon = getSettings().getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);

        DenseDoubleMatrix1D oldX = new DenseDoubleMatrix1D(n);
        // initialize with 1
        for (int i=0; i < n; i++) {
            oldX.setQuick(i, 1.0);
        }
        boolean positive = false;
        for (int iter = 0; iter < maxIters; iter++) {
            DenseDoubleMatrix1D newX = new DenseDoubleMatrix1D(n);
            sccMatrix.zMult(oldX, newX);

            if (verbosity >= 2) {
                //mainLog.println("Iteration " + iter+ ": "+newX.toString());
            }

            // check whether the we have a "drop" in all elements of the vector
            boolean allSmaller = true;
            boolean converged = true;
            for (int i = 0; i < n; i++) {
                double oldValue = oldX.getQuick(i);
                double newValue = newX.getQuick(i);

                if (PrismUtils.doublesAreClose(oldValue, newValue, epsilon, absolute)) {
                    allSmaller = false;
                } else {
                    converged = false;
                    if (!(newValue < oldValue)) {
                        allSmaller = false;
                    }
                }

                if ((!allSmaller) && (!converged)) {
                    break;
                }
            }

            if (allSmaller) {
                // we know that the vector will converge against zero, we are done
                if (verbosity >= 1)
                    timer.stop(" (SCC is zero, "+iter+" iterations)");
                return false;
            }

            oldX = newX;

            if (converged) {
                positive = true;
                if (verbosity >= 1)
                    timer.stop(" (SCC is positive, " + iter + " iterations)");
                break;
            }
        }

        if (!positive) {
            throw new PrismException("SCC analysis (power method) did not converge within " + maxIters + " iterations");
        }

        if (verbosity >= 2) {
            //mainLog.println("Eigenvector = " +oldX);
        }

        // compute cut
        int probState = mcc.nextSetBit(0);
        Set<Integer> cut = generateCut(product.getProduct(), probState, mcc);
        if (verbosity >= 2) {
            mainLog.println("Cut: " + cut);

            for (Integer cutState : cut) {
                mainLog.print("(" + prod.getDTMCState(cutState) + "," + prod.getUBAState(cutState) +  ")" );
            }
            mainLog.println();
        }

        // we have to weight the values...
        double sumCut = 0.0;
        for (int productIndex : cut) {
            int solutionIndex = map.get(equivMap.get(productIndex));
            sumCut += oldX.getQuick(solutionIndex);
        }
        double alpha = 1 / sumCut;
        if (verbosity >= 1) {
            mainLog.println("Weighing the eigenvector with alpha = " + alpha + " to obtain probabilities");
        }

        double cutSum = 0.0;
        double minValue = 0.0;
        double maxValue = 0.0;
        boolean first = true;
        for (Entry<Integer, Integer> entry : map.entrySet()) {
            int productIndex = entry.getKey();
            int solutionIndex = entry.getValue();

            double value = oldX.getQuick(solutionIndex);
            if (sanityCheck && value == 0.0) {
                throw new PrismException("Something strange going on (probability in positive SCC is zero for state "+productIndex+")");
            }
            // scale
            value = value * alpha;
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
//		if (verbosity >= 1) timer.stop();

        if (verbosity >= 1) {
            mainLog.println("Sum of probabilities for the cut C = "+cutSum + " for SCC "+mccIndex);
            mainLog.println("Probabilities in SCC " +mccIndex + " are in the range ["+minValue+","+maxValue+"]");
        }

        return true;
    }


    private boolean handleMCCPower2(LTLProduct product, int mccIndex, BitSet mcc, StateValues result, BitSet knownValues) throws PrismException
    {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("checking whether SCC " + mccIndex + " is positive and computing eigenvector");
        // first, check that SCC intersects accepting edges
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

        //FIX only handle LMCs with the same languages
        for(int idx = 0; idx < product.getProduct().getNumDtmcStates(); idx++){
            equivDTMC = new HashMap<>();
            equivDTMC.put(idx, 0);
        }

        StopWatch timerMatrix = new StopWatch(mainLog);
        timerMatrix.start("building positivity matrix");
        Map<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
        //FIX
        //DoubleMatrix2D sccMatrix = product.getProduct().positivityMatrixForMCC(mcc, map, false);


        Map<Integer, Integer> equivMap = new LinkedHashMap<Integer, Integer>();

        int size = 0;
        for (int i : IterableBitSet.getSetBits(mcc)) {
            int ubaState = prod.getUBAState(i);
            int lmcState = prod.getDTMCState(i);
            boolean isSet = false;
            for (Entry<Integer, Integer> states : map.entrySet()) {
                int s = states.getKey();
                if(((equivGFG.containsKey(ubaState) && equivGFG.get(ubaState).get(prod.getUBAState(s))) || prod.getUBAState(s) == ubaState) && (equivDTMC.get(lmcState) == equivDTMC.get(prod.getDTMCState(s)))){
                    equivMap.put(i, map.get(s));
                    isSet = true;
                    break;
                }
            }
            if(!isSet){
                map.put(i, size);
                equivMap.put(i, size);
                size++;
            }
        }

        if (verbosity >= 3) {
            mainLog.println("MCC mapping = "+map);
        }

        DoubleMatrix2D sccMatrix = new SparseDoubleMatrix2D(size,size);
        for (Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < prod.getNumChoices(from); c++) {
                for (Iterator<Entry<Integer, Double>> it = prod.getTransitionsIterator(from,c); it.hasNext(); ) {
                    Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();
                    if (mcc.get(to)) {
                        int matrixTo = equivMap.get(to);
                        sccMatrix.setQuick(map.get(from), matrixTo, sccMatrix.get(map.get(from),matrixTo)+ probMove.getValue());
                    }
                }
            }
        }

        if (verbosity >= 3) {
            mainLog.println("MCC mapping = " + map);
        }

        assert(sccMatrix.rows() == sccMatrix.columns());
        assert(sccMatrix.rows() == mcc.cardinality());
        if (verbosity >= 1) {
            timerMatrix.stop();
        }

        /*
        if (isEigenValueBig(sccMatrix)) {
            mainLog.println("MCC has eigenvalue bigger than 1, MCC is zero...");
            return false;
        }

        if (isEigenValueOne(sccMatrix)) {
            mainLog.println("has eigenvalue 1...");
        }else {
            mainLog.println("MCC does not have eigenvalue 1, MCC is zero...");
            return false;
        }
         */


        // M' = (M+I) / 2 to avoid periodicity problem
        // do it using forEachNonZero to exploit sparseness of the matrix
        sccMatrix.forEachNonZero(new IntIntDoubleFunction() {
            @Override
            public double apply(int row, int column, double value)
            {
                return value / 2;
            }
        });
        // add I/2 to diagonal
        for (int i = 0; i < sccMatrix.rows(); i++) {
            sccMatrix.setQuick(i, i, sccMatrix.getQuick(i, i) + 0.5);
        }

        if (verbosity >= 2) {
            //mainLog.println("(M+I)/2 =" + sccMatrix);
        }

        // do iterations
        //FIX
        //int n = mcc.cardinality();
        int n = size;
        int maxIters = getSettings().getInteger(PrismSettings.PRISM_MAX_ITERS);
        boolean absolute = (getSettings().getString(PrismSettings.PRISM_TERM_CRIT).equals("Absolute"));
        double epsilon = getSettings().getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);

        DenseDoubleMatrix1D oldX = new DenseDoubleMatrix1D(n);
        // initialize with 1
        for (int i=0; i < n; i++) {
            oldX.setQuick(i, 1.0);
        }
        boolean positive = false;
        for (int iter = 0; iter < maxIters; iter++) {
            DenseDoubleMatrix1D newX = new DenseDoubleMatrix1D(n);
            sccMatrix.zMult(oldX, newX);

            if (verbosity >= 2) {
                //mainLog.println("Iteration " + iter+ ": "+newX.toString());
            }

            // check whether the we have a "drop" in all elements of the vector
            boolean allSmaller = true;
            boolean converged = true;
            for (int i = 0; i < n; i++) {
                double oldValue = oldX.getQuick(i);
                double newValue = newX.getQuick(i);

                if (PrismUtils.doublesAreClose(oldValue, newValue, epsilon, absolute)) {
                    allSmaller = false;
                } else {
                    converged = false;
                    if (!(newValue < oldValue)) {
                        allSmaller = false;
                    }
                }

                if ((!allSmaller) && (!converged)) {
                    break;
                }
            }

            if (allSmaller) {
                // we know that the vector will converge against zero, we are done
                if (verbosity >= 1)
                    timer.stop(" (SCC is zero, "+iter+" iterations)");
                return false;
            }

            oldX = newX;

            if (converged) {
                positive = true;
                if (verbosity >= 1)
                    timer.stop(" (SCC is positive, " + iter + " iterations)");
                break;
            }
        }

        if (!positive) {
            mainLog.println("MCC does not converge, MCC is zero...");
            return false;
            //throw new PrismException("SCC analysis (power method) did not converge within " + maxIters + " iterations");
        }

//        if (verbosity >= 2) {
//            mainLog.println("Eigenvector = " +oldX.getQuick(0) + oldX.getQuick(1));
//        }

        // compute cut
        int probState = mcc.nextSetBit(0);
        Set<Integer> cut = generateCut(product.getProduct(), probState, mcc);
        if (verbosity >= 2) {
            mainLog.println("Cut: " + cut);

            for (Integer cutState : cut) {
                mainLog.print("(" + prod.getDTMCState(cutState) + "," + prod.getUBAState(cutState) +  ")" );
            }
            mainLog.println();
        }
        if (verbosity >= 2) {
            mainLog.println("Eigenvector = " +oldX.getQuick(0) + oldX.getQuick(1));
        }
        // we have to weight the values...
        double sumCut = 0.0;
        for (int productIndex : cut) {
            int tmp = equivMap.get(productIndex);
            if (verbosity >= 2) {
                mainLog.println("map.size = " + map.size() + ", tmp=" + tmp + "map.getindex" + map.get(productIndex));
            }
            //int solutionIndex = map.get(tmp);
            sumCut += oldX.getQuick(tmp);//productIndex
            //FIX
            //sumCut += oldX.getQuick(productIndex);//tmp
        }
        mainLog.println("sumCut = " + sumCut);

        double alpha = 1 / sumCut;
        if(alpha <= 0 || Double.isInfinite(alpha) || Double.isInfinite(sumCut) ){
            mainLog.println("MCC does not converge, MCC is zero...");
            return false;
        }
        if (verbosity >= 1) {
            mainLog.println("Weighing the eigenvector with alpha = " + alpha + " to obtain probabilities");
        }

        double cutSum = 0.0;
        double minValue = 0.0;
        double maxValue = 0.0;
        boolean first = true;
        for (Entry<Integer, Integer> entry : map.entrySet()) {
            int productIndex = entry.getKey();
            int solutionIndex = entry.getValue();

            double value = oldX.getQuick(solutionIndex);
            if (sanityCheck && value == 0.0) {
                throw new PrismException("Something strange going on (probability in positive SCC is zero for state "+productIndex+")");
            }
            // scale
            value = value * alpha;
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
//		if (verbosity >= 1) timer.stop();

        if (verbosity >= 1) {
            mainLog.println("Sum of probabilities for the cut C = "+cutSum + " for SCC "+mccIndex);
            mainLog.println("Probabilities in SCC " +mccIndex + " are in the range ["+minValue+","+maxValue+"]");
        }

        return true;
    }


    private boolean checkIsMCCPositiveQuotient(LTLProduct product, int mccIndex, BitSet mcc) throws PrismException {
        StopWatch timer = new StopWatch(mainLog);
        timer.start("checking whether MCC " + mccIndex + " is positive");
        DTMCGFGProduct prod = product.getProduct();
        HashMap<Integer, Set<APElement>> accEdges = product.getAccEdges();
        // first, check that MCC intersects accepting edges
        boolean acc = false;
        for (int state : new IterableStateSet(mcc, prod.getNumStates())) {
            Set<APElement> set = accEdges.get(state);
            if (set == null) continue;
            for (APElement ap : set) {
                for (int i = 0; i < prod.getNumChoices(state); i++) {
                    if (!prod.getAction(state, i).equals(ap)) continue;
                    if (prod.allSuccessorsInSet(state, i, mcc)) {
                        acc = true;
                        break;
                    }
                }
            }
        }

        if (!acc) {
            if (verbosity >= 1) timer.stop(" (MCC is zero, no accepting edges)");
            return false;
        }

        StopWatch timerMatrix = new StopWatch(mainLog);
        timerMatrix.start("building positivity matrix");
        Map<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
        //DoubleMatrix2D sccMatrix = product.getProduct().positivityMatrixForMCC(mcc, map, false);


        /////////////////////compute language equivalent states
        /// compute equivalent LMC states
        // Step 0: build dtmc successor map
        Map<Integer, Map<Integer, Set<Integer>>> succMap = new HashMap<>(); //successors map for dtmc: S x Symbols -> 2^S
        Set<Integer> mccSet =new HashSet<>();
        for(int s: IterableBitSet.getSetBits(mcc)) {
            mccSet.add(s);
        }
        
        for (int s: IterableBitSet.getSetBits(mcc)) {
            Map<Integer, Set<Integer>> succLetterMap = new HashMap<>();
            for(int letter = 0; letter < prod.getNumDtmcStates(); letter++){
                Set<Integer> possible = prod.getProbStatesLetterSuccessors(s, letter);
                mainLog.println(s + ", possible = " + possible);

                Set<Integer> t = new HashSet<>();
                for(int tmp: possible){
                    if(mcc.get(tmp)){
                        mainLog.println(s + ", letter = " +letter + ", " + tmp);
                        t.add(tmp);
                    }
                }
                succLetterMap.put(letter, t);
            }
            succMap.put(s, succLetterMap);
        }
        {
            HashMap<Integer, Set<Integer>> succLetterMap = new HashMap<>();
            for (int letter = 0; letter < prod.getNumDtmcStates(); letter++) {
                Set<Integer> tmp = new HashSet<>();
                tmp.add(-1);
                succLetterMap.put(letter, tmp);
            }
            succMap.put(-1, succLetterMap);
        }
        if(verbosity >= 2) {
            mainLog.println("succMap = " + succMap.size());
            for(int i: succMap.keySet()){
                mainLog.print( i + " -> " + succMap.get(i) + ", ");
            }
            mainLog.println();

        }


        // Initial partition: accepting vs non-accepting
        Set<Integer> accepting = new HashSet<>(mccSet);
        Set<Integer> nonAccepting = new HashSet<>();
        nonAccepting.add(-1);

        List<Set<Integer>> partitions = new ArrayList<>();
        if (!accepting.isEmpty()) partitions.add(accepting);
        if (!nonAccepting.isEmpty()) partitions.add(nonAccepting);

        boolean refined;
        do {
            refined = false;
            List<Set<Integer>> newPartitions = new ArrayList<>();

            for (Set<Integer> group : partitions) {
                Map<Map<Integer, Set<Set<Integer>>>, Set<Integer>> splitter = new HashMap<>();

                for (int state : group) {
                    Map<Integer, Set<Integer>> succLetterMap = succMap.getOrDefault(state, Collections.emptyMap());
                    Map<Integer, Set<Set<Integer>>> image = new HashMap<>();

                    for(int letter = 0; letter < prod.getNumDtmcStates(); letter++){
                        Set<Integer> succs = succLetterMap.getOrDefault(letter, Collections.emptySet());
                        Set<Set<Integer>> targetBlocks = new HashSet<>();
                        for (Set<Integer> block : partitions) {
                            for (int s : succs) {
                                if (block.contains(s)) {
                                    targetBlocks.add(block);
                                    break;
                                }
                            }
                        }
                        image.put(letter, targetBlocks);
                    }

                    splitter.computeIfAbsent(image, k -> new HashSet<>()).add(state);
                }

                if (splitter.size() == 1) {
                    newPartitions.add(group);
                } else {
                    refined = true;
                    newPartitions.addAll(splitter.values());
                }
            }

            partitions = newPartitions;
        } while (refined);

        if (verbosity >= 3) {
            mainLog.println("partitions = " + partitions);
        }

        Map<Integer, Integer> equivMap = new LinkedHashMap<Integer, Integer>();
        for(int idx = 0; idx < partitions.size(); idx++){
            boolean first = true;
            for(int s: partitions.get(idx)){
                if(s == -1) continue;
                if(first){
                    map.put(s, idx);
                    first = false;
                }
                equivMap.put(s, idx);
            }
        }
        int size = map.size();

        /////////////////////end language equivalence check


        /*
        Map<Integer, Integer> equivMap = new LinkedHashMap<Integer, Integer>();

        int size = 0;
        for (int i : IterableBitSet.getSetBits(mcc)) {
            int ubaState = prod.getUBAState(i);
            int lmcState = prod.getDTMCState(i);
            boolean isSet = false;
            mainLog.println("(ubaState, lmcState) = " + "(" + ubaState + ", " + lmcState + ")");

            for (Entry<Integer, Integer> states : map.entrySet()) {
                int s = states.getKey();
                if(((equivGFG.containsKey(ubaState) && equivGFG.get(ubaState).get(prod.getUBAState(s))) || prod.getUBAState(s) == ubaState) && (Objects.equals(equivDTMC.get(lmcState), equivDTMC.get(prod.getDTMCState(s))))){
                    mainLog.println("equivalent to: "+ "(" + prod.getUBAState(s) + ", " + prod.getDTMCState(s) + ")");

                    equivMap.put(i, map.get(s));
                    isSet = true;
                    break;
                }
            }
            if(!isSet){
                mainLog.println("new");
                map.put(i, size);
                equivMap.put(i, size);
                size++;
            }
        }
        */


        if (verbosity >= 3) {
            mainLog.println("MCC mapping = "+map);
        }

        //DoubleMatrix2D matrix = new SparseDoubleMatrix2D(size,size);
        DoubleMatrix2D sccMatrix = new SparseDoubleMatrix2D(size,size);
        for (Entry<Integer, Integer> states : map.entrySet()) {
            int from = states.getKey();
            for(int c = 0; c < prod.getNumChoices(from); c++) {
                for (Iterator<Entry<Integer, Double>> it = prod.getTransitionsIterator(from,c); it.hasNext(); ) {
                    Entry<Integer, Double> probMove = it.next();
                    int to = probMove.getKey();
                    if (mcc.get(to)) {
                        int matrixTo = equivMap.get(to);
                        sccMatrix.setQuick(map.get(from), matrixTo, sccMatrix.get(map.get(from),matrixTo)+ probMove.getValue());
                    }
                }
            }
        }
        if (verbosity >= 3) {
            mainLog.println("MCC mapping = " + map);
            for (Entry<Integer, Integer> states : map.entrySet()) {
                int prodstate = states.getKey();
                int ubaState = prod.getUBAState(prodstate);
                int lmcState = prod.getDTMCState(prodstate);
                mainLog.print("("+ubaState+", "+ lmcState + ") ");
            }
            mainLog.println();
        }

        if (verbosity >= 2) {

            mainLog.println("sccMatrix: ");
            for (int i = 0; i< sccMatrix.rows(); i++) {
                for (int j = 0; j < sccMatrix.columns(); j++) {
                    mainLog.print(sccMatrix.getQuick(i, j) + ", ");
                }
                mainLog.println();
            }
            mainLog.println();
        }
        assert(sccMatrix.rows() == sccMatrix.columns());
        assert(sccMatrix.rows() == mcc.cardinality());
        if (verbosity >= 1) {
            timerMatrix.stop();
        }


        if (isEigenValueBig(sccMatrix)) {
            mainLog.println("MCC has eigenvalue bigger than 1, MCC is zero...");
            return false;
        }

        if (isEigenValueOne(sccMatrix)) {
            mainLog.println("has eigenvalue 1...");
            return true;
        }else {
            mainLog.println("MCC does not have eigenvalue 1, MCC is zero...");
            return false;
        }

         /*

        int rows = sccMatrix.rows();

        boolean positive;
        String posMethod = getSettings().getString(PrismSettings.PRISM_UBA_POS_METHOD);
        switch (posMethod) {
            case "SVD": {
                Algebra algebra = new Algebra();
                int rank = algebra.rank(sccMatrix);
                if (rank < rows-1 || rank > rows) {
                    throw new PrismException("Strange things are going on (rank = " + rank + ", matrix has size " + rows + "x" + rows +")..");
                }

                if (verbosity >= 1) mainLog.println("Rank of SCC " + mccIndex + " = " + rank + ", full rank is " + rows);
                positive = rank < rows;
                break;
            }
            case "QR": {
                QRDecomposition qr = new QRDecomposition(sccMatrix);
                positive = !hasFullRankWithTolerance(qr);
                break;
            }
            default:
                throw new PrismException("Unknown UBA method for checking SCC positivity: "+posMethod);
        }

        if (verbosity >= 1) timer.stop(" (MCC is " + (positive ? "positive" : "zero")+", checked via " + posMethod + ")");

        return positive;
          */
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

        if(isEigenValueOne(mccMatrix)){
            mainLog.println("has eigenvalue 1...");
            return true;
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
    //check whether eigenvalue is greater than 1
    private boolean isEigenValueBig(DoubleMatrix2D matrix){
        EigenvalueDecomposition eig = new EigenvalueDecomposition(matrix);
        DoubleMatrix1D realParts = eig.getRealEigenvalues();
        DoubleMatrix1D imagParts = eig.getImagEigenvalues();
        double epsilon = getSettings().getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);

        boolean hasEigenvalueGreaterThanOne = false;

        for (int i = 0; i < realParts.size(); i++) {
            double real = realParts.get(i);
            double imag = imagParts.get(i);
            double abs = Math.hypot(real, imag); // sqrt(real^2 + imag^2)

            if (abs-epsilon > 1.0) {
                hasEigenvalueGreaterThanOne = true;
                break;
            }
        }
        return hasEigenvalueGreaterThanOne;
    }

    //check whether eigenvalue is greater than 1
    private boolean isEigenValueOne(DoubleMatrix2D matrix){
        EigenvalueDecomposition eig = new EigenvalueDecomposition(matrix);
        DoubleMatrix1D realParts = eig.getRealEigenvalues();
        DoubleMatrix1D imagParts = eig.getImagEigenvalues();
        boolean absolute = (getSettings().getString(PrismSettings.PRISM_TERM_CRIT).equals("Absolute"));
        double epsilon = getSettings().getDouble(PrismSettings.PRISM_TERM_CRIT_PARAM);
        boolean hasEigenvalueOne = false;

        for (int i = 0; i < realParts.size(); i++) {
            double real = realParts.get(i);
            double imag = imagParts.get(i);
            double abs = Math.hypot(real, imag); // sqrt(real^2 + imag^2)
            if (PrismUtils.doublesAreClose(abs, 1, 1e-5, absolute)) {
                hasEigenvalueOne = true;
                break;
            }
        }
        return hasEigenvalueOne;
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
        if (verbosity >= 2) {

            mainLog.println("Matrix: ");
            for (int i = 0; i< matrix.rows(); i++) {
                for (int j = 0; j < matrix.columns(); j++) {
                    mainLog.print(matrix.getQuick(i, j) + ", ");
                }
                mainLog.println();
            }
            mainLog.println();
            mainLog.print("B: ");
            for (int i = 0; i< matrix.rows(); i++){
                mainLog.print(B.getQuick(i,0) + ", ");
            }
            mainLog.println();
        }

        if (verbosity >= 1) timer.stop();

        if (verbosity >= 1) mainLog.println("Solving MCC " + mccIndex + " probability values...");
        timer.start("solving equation system for positive MCC " + mccIndex);

        //********* new solver
        DoubleMatrix2D At = algebra.transpose(matrix);
        DoubleMatrix2D AtA = algebra.mult(At, matrix);
        DoubleMatrix2D AtB = algebra.mult(At, B);

        // Now solve AtA * x = AtB
        DoubleMatrix2D solution = algebra.solve(AtA, AtB);
        //*********

        //DoubleMatrix2D solution = algebra.solve(matrix, B);
        assert (solution.rows() == matrix.rows());
        assert (solution.columns() == 1);

        if (verbosity >= 2) {
            mainLog.println("Solution: " + solution.rows() + ", " + solution.columns());
            for (int i = 0; i < solution.rows(); i++) {
                mainLog.print(solution.getQuick(i, 0) + ", ");
            }
            mainLog.println();
        }

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

        //if (verbosity >= 2) mainLog.println("Solution = " + solution);

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

    private Set<Integer> generateCut(DTMCGFGProduct product, final int probState, final BitSet mcc) throws PrismException  {
        StopWatch timer = new StopWatch(mainLog);

        if (!mcc.get(probState)) {
            throw new PrismException("Initial state for cut generation is not in MCC");
        }

        timer.start("generating cut");
        // tracks, for every state q the set delta(q, z)
        Map<Integer,Set<Integer>> succsForZ = new HashMap<Integer,Set<Integer>>();

        // partition the mcc
        //generateBisimilarMap(product, mcc);

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
                if (verbosity >= 4) mainLog.println("queue entry: (" + current.getFirstState() + "," + current.getSecondState() + ")");

                if (!visited.contains(current)) {
                    SharedWord<APElement> word = entry.second;
                    Integer left = current.getFirstState();
                    Integer right = current.getSecondState();

                    visited.add(current);
                    countVisited++;

                    if (!mcc.get(left) || !mcc.get(right)) {
                        // we don't want to leave the MCC
                        continue;
                    }

                    //mainLog.println("product.getNumChoices(left) = " + product.getNumChoices(left));

                    for(int c = 0; c < product.getNumChoices(left); c++) {
                        APElement ap = (APElement) product.getAction(left, c);
                        Set<Integer> leftSuccs = product.getProbStatesSuccessors(left,ap);
                        Set<Integer> rightSuccs;
                        if (left.equals(right)) {
                            rightSuccs = leftSuccs;
                        } else {
                            rightSuccs = product.getProbStatesSuccessors(right,ap);
                        }
                        //mainLog.println("leftSuccs = " + leftSuccs + "; rightSuccs = " + rightSuccs);



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


                                boolean isLeftBis = (leftSucc == probState) || (equivGFG.get(product.getUBAState(probState)) != null && equivGFG.get(product.getUBAState(probState)).get(product.getUBAState(leftSucc)) && (Objects.equals(equivDTMC.get(product.getDTMCState(probState)), equivDTMC.get(product.getDTMCState(leftSucc)))));
                                boolean isRightBis = (rightSucc == probState) || (equivGFG.get(product.getUBAState(probState)) != null && equivGFG.get(product.getUBAState(probState)).get(product.getUBAState(rightSucc))&& (Objects.equals(equivDTMC.get(product.getDTMCState(probState)), equivDTMC.get(product.getDTMCState(rightSucc)))));


                                if ((!isLeftBis  && (rightSucc == probState) && !succsForZ.get(leftSucc).isEmpty()) ||
                                                (!isRightBis  && (leftSucc == probState) && !succsForZ.get(rightSucc).isEmpty())) {
                                    //Found extension
                                    found = true;
                                    if (verbosity >= 3) {
                                        String wordInfo = curWord.size() > 3 ? "length(y) = " + curWord.size() : "y = " + curWord;
                                        mainLog.println("FOUND EXTENSION: (" + leftSucc + "," + rightSucc + ") with word " + wordInfo + ";");
                                        String cutInfo;
                                        if (leftSucc == probState) {
                                            if (succsForZ.get(rightSucc).size() > 4) {
                                                cutInfo = succsForZ.get(rightSucc).size() + " states";
                                            } else {
                                                cutInfo = succsForZ.get(rightSucc).toString();
                                            }
                                        } else {
                                            if (succsForZ.get(leftSucc).size() > 4) {
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
                                   ////////////////////////stop here
                                    if(succsForZ_.get(probState).size() == succsForZ.get(probState).size()){
                                        Set<Integer> tempC = succsForZ.get(probState);

                                        if (verbosity >= 1) {
                                            timer.stop(" ("+iterations+" iterations, " + countVisited + " extension checks, cut tempC has " + tempC.size() + " states)");
                                        }

                                        //get rid of bisimilar states
                                        Set<Integer> C = new HashSet<Integer>();
                                        for (int cc: tempC){
                                            boolean isContained = false;
                                            if(equivGFG.get(product.getUBAState(cc)) != null){
                                                for(int t: C){
                                                    if(equivGFG.get(product.getUBAState(cc)).get(product.getUBAState(t)) && (equivDTMC.get(product.getDTMCState(cc))==equivDTMC.get(product.getDTMCState(t)))){
                                                        isContained = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            if(!isContained)
                                                C.add(cc);
                                        }

                                        if (verbosity >= 1) {
                                            timer.stop(" ("+iterations+" iterations, " + countVisited + " extension checks, cut C has " + C.size() + " states)");
                                        }

                                        return C;
                                    }
                                    /// /////////////////////////////////////////
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
                boolean isContained = false;
                if(equivGFG.get(product.getUBAState(c)) != null){
                    for(int t: C){
                        if(equivGFG.get(product.getUBAState(c)).get(product.getUBAState(t)) && (equivDTMC.get(product.getDTMCState(c))==equivDTMC.get(product.getDTMCState(t)))){
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

    private List<BitSet> computeMCCs(DTMCGFGProduct prod) throws PrismException {
//        ECComputer ecs = ECComputer.createECComputer(this, prod);
//        ecs.computeMECStates();

        ECComputerDefault ecs = new ECComputerDefault(this, prod);
        ecs.computeMECStates();
        return ecs.getMECStates();
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
