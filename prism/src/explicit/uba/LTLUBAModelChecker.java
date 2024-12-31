package explicit.uba;

import java.io.PrintStream;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;
import java.util.Vector;

import acceptance.AcceptanceOmega;
import acceptance.AcceptanceReach;
import acceptance.AcceptanceType;
import automata.DA;
import automata.DASimplifyAcceptance;
import common.IterableBitSet;
import common.PathUtil;
import common.StopWatch;
import cern.colt.function.IntIntDoubleFunction;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.QRDecomposition;
import explicit.DTMC;
import explicit.DTMCModelChecker;
import explicit.LTLModelChecker;
import explicit.Model;
import explicit.ModelCheckerResult;
import explicit.PredecessorRelation;
import explicit.ProbModelChecker;
import explicit.ProductState;
import explicit.SCCComputer;
import explicit.SCCConsumerBSCCs;
import explicit.SCCConsumerStore;
import explicit.StateModelChecker;
import explicit.StateValues;
import jltl2ba.APElement;
import jltl2ba.APSet;
import jltl2ba.MyBitSet;
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

public class LTLUBAModelChecker extends PrismComponent
{
	// TODO: stub
	public class LTLProduct {
		private DTMCUBAProduct prod;
		
		public LTLProduct(DTMCUBAProduct prod) {
			this.prod = prod;
		}
		
		public DTMCUBAProduct getProduct() {
			 return prod;
		}
		
		public MyBitSet getFinalStates() {
			return prod.getFinalStates();
		}
		
		public StateValues projectToOriginalModel(DTMC model, StateValues prodValues) throws PrismLangException {
			StateValues result = new StateValues(TypeDouble.getInstance(), model);

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

	public LTLUBAModelChecker(ProbModelChecker mc)
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
	public NBA constructUBAForLTLFormula(ProbModelChecker mc, DTMC model, Expression expr, Vector<BitSet> labelBS) throws PrismException
	{
		Expression ltl;
		NBA uba;
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
			LTL2UBA ltl2uba = new LTL2UBA(this);
			time = System.currentTimeMillis();
			mainLog.println("Parsing and constructing HOA automaton for "+expr);
			Vector<Expression> apExpressions = new Vector<Expression>();
			uba = ltl2uba.fromExpressionHOA(expr, PathUtil.getDirectoryForRelativePropertyResource(mc.getModulesFile(), mc.getPropertiesFile()), apExpressions);

			mainLog.println("Determining states satisfying atomic proposition expressions of the automaton...");
			for (int i=0; i<uba.getAPSet().size(); i++) {
				Expression label = apExpressions.get(i);
				label.typeCheck();
				BitSet labelStates = mc.checkExpression(model, label, null).getBitSet();
				labelBS.add(labelStates);
				uba.getAPSet().renameAP(i, "L"+i);
			}
			mainLog.println("labelBS: " + labelBS);

		} else {
			// Model check maximal state formulas
			LTLModelChecker ltlMC = new LTLModelChecker(this);
			ltl = ltlMC.checkMaximalStateFormulas(mc, model, expr.deepCopy(), labelBS);

			// Convert LTL formula to UBA
			mainLog.println("\nBuilding unambiguous Buchi automaton (for " + ltl + ")...");
			time = System.currentTimeMillis();
			LTL2UBA ltl2uba = new LTL2UBA(this);
			uba = ltl2uba.convertLTLFormulaToUBA(ltl, mc.getConstantValues());
		}
		mainLog.println("UBA has " + uba.size() + " states.");
		checkForCanonicalAPs(uba.getAPSet(), labelBS.size());
		time = System.currentTimeMillis() - time;
		mainLog.println("Time for UBA translation: " + time / 1000.0 + " seconds.");
		// If required, export UBA
		if (settings.getExportPropAut()) {
			mainLog.println("Exporting UBA to file \"" + settings.getExportPropAutFilename() + "\"...");
			PrintStream out = PrismUtils.newPrintStream(settings.getExportPropAutFilename());
			uba.print(out, settings.getExportPropAutType());
			out.close();
		}
		
		isUFA = checkUpwardClosedness(uba) && !settings.getBoolean(PrismSettings.PRISM_UBA_PURE);
		if(isUFA) {
			mainLog.println("UBA is actually an (upward-closed) UFA");
		}

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
	public boolean checkUpwardClosedness(NBA uba)
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

	
	public LTLProduct constructProduct(DTMC model, NBA uba, Vector<BitSet> labelBS, BitSet statesOfInterest) throws PrismException
	{
		StopWatch timer = new StopWatch(mainLog);
		timer.start("computing UBA-DTMC product");
		DTMCUBAProduct prod = new DTMCUBAProduct(this, model, uba, labelBS, statesOfInterest);

		if (settings.isExportDTMCUBAProduct()) {
			mainLog.println("Exporting DTMC UBA product to file \"" + settings.getExportDTMCUBAProductFilename() + "\" …");
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
		// the states in positive SCCs
		BitSet S_SCC_pos = new BitSet();

		StopWatch timer = new StopWatch(mainLog);

		BitSet unknownStates = new BitSet();
		if (!isUFA) {
			timer.start("computing SCC in UBA-DTMC product");
			List<BitSet> sccs = computeSCCs(product.getProduct());
			int sccCount = sccs.size();
			timer.stop(" (found "+sccCount+" non-trivial SCCs)");

			int posSccCount = 0;
			int sccIndex = 0;
			timer.start("computing SCC probabilities for positive SCCs");
			for (BitSet scc : sccs) {
				sccIndex++;

				if (verbosity >= 2) {
					mainLog.println("\nSCC " + sccIndex + ": " + scc);
				} else if (verbosity >= 1) {
					mainLog.println("\nSCC " + sccIndex + " has " + scc.cardinality() + " states");
				}

				boolean isPositive;
				if (getSettings().getBoolean(PrismSettings.PRISM_UBA_POWER)) {
					isPositive = handleSCCPower(product, sccIndex, scc, result, S_SCC_pos);
				} else {
					isPositive = handleSCC(product, sccIndex, scc, result, S_SCC_pos);
				}
				if (isPositive) posSccCount++;
			}
			timer.stop(" ("+posSccCount+" positive SCCs, known probabilities for "+S_SCC_pos.cardinality()+" states)");

			if (verbosity >= 2) {
				mainLog.println("Partial result (all SCCs):");
				result.printFiltered(mainLog, S_SCC_pos, true, false, false, true);
			}
		
		
			timer.start("determining states with probability zero");
			PredecessorRelation pre = product.getProduct().getPredecessorRelation(this, true);
			BitSet S_prePositiveSCC = pre.calculatePreStar(null, S_SCC_pos, new BitSet());
			BitSet S_zero = new BitSet();
			S_zero.flip(0, product.getProduct().getNumStates());
			S_zero.andNot(S_prePositiveSCC);

			unknownStates.flip(0, product.getProduct().getNumStates());
			unknownStates.andNot(S_zero);
			unknownStates.andNot(S_SCC_pos);

			// we are only interested in the prob states

			timer.stop(" ("+S_zero.cardinality()+" zero prob. states, "+unknownStates.cardinality()+" remaining unknown)");
		}
		
		if (isUFA) {
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
			if (!unknownStates.isEmpty()) {
				if (getSettings().getString(PrismSettings.PRISM_UBA_REACH_METHOD).equals("GS")) {
					computeReachabilityGS(product, result, S_SCC_pos, unknownStates);
				} else {
					computeReachability(product, result, S_SCC_pos, unknownStates);
				}
			}
		}

		if (verbosity >= 2) {
			mainLog.println("Result:");
			result.printFiltered(mainLog, null, true, false, false, true);
		}

		return result;
	}
	
	private boolean handleSCC(LTLProduct product, int sccIndex, BitSet scc, StateValues result, BitSet knownValues) throws PrismException
	{
		boolean isPositive = checkIsSCCPositive(product, sccIndex, scc);
		if (!isPositive)
			return false;

		int probState = scc.nextSetBit(0);
		Set<Integer> cut = generateCut(product.getProduct(), probState, scc);
		if (verbosity >= 2) {
			mainLog.println("Cut: " + cut);

			for (Integer cutState : cut) {
				DTMCUBAProduct prod = product.getProduct();
				mainLog.print("(" + prod.getDTMCState(cutState) + "," + prod.getUBAState(cutState) +  ")" );
			}
			mainLog.println();
		}

		computeSCCProbs(product, sccIndex, scc, cut, result, knownValues);

		return isPositive;
	}

	private boolean handleSCCPower(LTLProduct product, int sccIndex, BitSet scc, StateValues result, BitSet knownValues) throws PrismException
	{
		StopWatch timer = new StopWatch(mainLog);
		timer.start("checking whether SCC " + sccIndex + " is positive and computing eigenvector");
		// first, check that SCC intersects F
		if (!scc.intersects(product.getFinalStates())) {
			if (verbosity >= 1) timer.stop(" (SCC is zero, no final state)");
			return false;
		}

		StopWatch timerMatrix = new StopWatch(mainLog);
		timerMatrix.start("building positivity matrix");
		Map<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
		DoubleMatrix2D sccMatrix = product.getProduct().positivityMatrixForSCC(scc, map, false);
		assert(sccMatrix.rows() == sccMatrix.columns());
		assert(sccMatrix.rows() == scc.cardinality());
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
			mainLog.println("(M+I)/2 =" + sccMatrix);
		}

		// do iterations
		int n = scc.cardinality();

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
			throw new PrismException("SCC analysis (power method) did not converge within " + maxIters + " iterations");
		}

		if (verbosity >= 2) {
			mainLog.println("Eigenvector = " +oldX);
		}

		// compute cut
		int probState = scc.nextSetBit(0);
		Set<Integer> cut = generateCut(product.getProduct(), probState, scc);
		if (verbosity >= 2) {
			mainLog.println("Cut: " + cut);

			for (Integer cutState : cut) {
				DTMCUBAProduct prod = product.getProduct();
				mainLog.print("(" + prod.getDTMCState(cutState) + "," + prod.getUBAState(cutState) +  ")" );
			}
			mainLog.println();
		}

		// we have to weight the values...
		double sumCut = 0.0;
		for (int productIndex : cut) {
			int solutionIndex = map.get(productIndex);
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
			mainLog.println("Sum of probabilities for the cut C = "+cutSum + " for SCC "+sccIndex);
			mainLog.println("Probabilities in SCC " +sccIndex + " are in the range ["+minValue+","+maxValue+"]");
		}

		return true;
	}

	private boolean checkIsSCCPositive(LTLProduct product, int sccIndex, BitSet scc) throws PrismException
	{
		StopWatch timer = new StopWatch(mainLog);
		timer.start("checking whether SCC " + sccIndex + " is positive");
		// first, check that SCC intersects F
		if (!scc.intersects(product.getFinalStates())) {
			if (verbosity >= 1) timer.stop(" (SCC is zero, no final state)");
			return false;
		}

		StopWatch timerMatrix = new StopWatch(mainLog);
		timerMatrix.start("building positivity matrix");
		DoubleMatrix2D sccMatrix = product.getProduct().positivityMatrixForSCC(scc, null, true);
		assert(sccMatrix.rows() == sccMatrix.columns());
		assert(sccMatrix.rows() == scc.cardinality());
		if (verbosity >= 1) {
			timerMatrix.stop();
		}

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
			
			if (verbosity >= 1) mainLog.println("Rank of SCC " + sccIndex + " = " + rank + ", full rank is " + rows);
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

		if (verbosity >= 1) timer.stop(" (SCC is " + (positive ? "positive" : "zero")+", checked via " + posMethod + ")");

		return positive;
	}

	private void computeSCCProbs(LTLProduct product, int sccIndex, BitSet scc, Set<Integer> cut, StateValues result, BitSet knownValues) throws PrismException
	{
		StopWatch timer = new StopWatch(mainLog);
		timer.start("building matrix for positive SCC probabilities, SCC " + sccIndex);

		LinkedHashMap<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
		DoubleMatrix2D matrix = product.getProduct().valueMatrixForSCC(scc, cut, map);
		assert(matrix.rows() + 1 == matrix.columns());
		assert(matrix.rows() < scc.cardinality() + 1);

		Algebra algebra = new Algebra();
		DoubleMatrix2D B = new DenseDoubleMatrix2D(matrix.rows(), 1);
		B.setQuick(matrix.rows()-1, 0, 1);
		//if (verbosity >= 2) mainLog.println("B = \n" + B);
		
		if (verbosity >= 1) timer.stop();

		if (verbosity >= 1) mainLog.println("Solving SCC " + sccIndex + " probability values...");
		timer.start("solving equation system for positive SCC " + sccIndex);
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
				throw new PrismException("Something strange going on (probability in positive SCC is zero for state "+productIndex+")");
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
			mainLog.println("Sum of probabilities for the cut C = "+cutSum + " for SCC "+sccIndex);
			mainLog.println("Probabilities in SCC " +sccIndex + " are in the range ["+minValue+","+maxValue+"]");
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

	private void computeReachability(LTLProduct product, StateValues result, BitSet s_SCC_pos, BitSet unknownStates) throws PrismException
	{
		StopWatch timer = new StopWatch(mainLog);
		timer.start("solving linear equation system for remaining unknown states");
		LinkedHashMap<Integer, Integer> map = new LinkedHashMap<Integer, Integer>();
		Pair<DoubleMatrix2D, DoubleMatrix2D> lgs = product.getProduct().reachabilityEquationSystem(s_SCC_pos, unknownStates, result, map);

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

	private void computeReachabilityGS(LTLProduct product, StateValues result, BitSet s_SCC_pos, BitSet unknownStates) throws PrismException
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
		DTMCModelChecker productMc = new DTMCModelChecker(mc);
		ModelCheckerResult r = productMc.computeReachProbsGaussSeidel(product.getProduct(), no, yes, values, known);

		if (verbosity >= 2) mainLog.println("Solution = " + values);

		for (int i : new IterableBitSet(unknownStates)) {
			double value = r.soln[i];
			if (value == 0.0) {
				throw new PrismException("Something strange going on (probability during reachability is zero for state "+i+")");
			}
			result.setDoubleValue(i, value);
		}
		timer.stop();
	}

	private Set<Integer> generateCut(DTMCUBAProduct product, final int probState, final BitSet scc) throws PrismException  {
		StopWatch timer = new StopWatch(mainLog);

		if (!scc.get(probState)) {
			throw new PrismException("Initial state for cut generation is not in SCC");
		}

		timer.start("generating cut");
		// tracks, for every state q the set delta(q, z)
		Map<Integer,Set<Integer>> succsForZ = new HashMap<Integer,Set<Integer>>();


		// for every state, the state can be reached via z=epsilon
		for(Integer probState_: IterableBitSet.getSetBits(scc)) {
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
			Queue<Pair<ProductState,SharedWord<Integer>>> queue = new LinkedList<Pair<ProductState, SharedWord<Integer>>>();
			queue.add(new Pair<ProductState, SharedWord<Integer>>(new ProductState(probState, probState),
			                                                      new SharedWord<Integer>()));
			while(!queue.isEmpty()) {
				Pair<ProductState, SharedWord<Integer>> entry = queue.poll();
				ProductState current = entry.first;
				if (!visited.contains(current)) {
					SharedWord<Integer> word = entry.second;
					Integer left = current.getFirstState();
					Integer right = current.getSecondState();
					if (!scc.get(left) || !scc.get(right)) {
						// we don't want to leave the SCC
						continue;
					}

					Set<Integer> leftSuccs = product.getProbStatesSuccessors(left);
					Set<Integer> rightSuccs;
					if (left == right) {
						rightSuccs = leftSuccs;
					} else {
						rightSuccs = product.getProbStatesSuccessors(right);
					}
					visited.add(current);
					countVisited++;
					for (Integer leftSucc : leftSuccs) {
						for(Integer rightSucc : rightSuccs) {
							if(!scc.get(leftSucc) || !scc.get(rightSucc)) {
								// we don't want to leave the SCC
								continue;
							}
							// if the move has a different symbol...
							if (product.getDTMCState(leftSucc) != product.getDTMCState(rightSucc)) {
								continue;
							}

							// ... we extend the word and add to the BFS queue
							SharedWord<Integer> curWord = word.append(Integer.valueOf(product.getDTMCState(leftSucc)));
							
							if (leftSucc != rightSucc &&
							    ((rightSucc == probState && !succsForZ.get(leftSucc).isEmpty()) ||
							     (leftSucc == probState && !succsForZ.get(rightSucc).isEmpty()))) {
								//Found extension 
								found = true;
								if (verbosity >= 2) {
									String wordInfo = curWord.size() > 3 ? "length(y) = " + curWord.size() : "y = " + curWord;
									mainLog.println("FOUND EXTENSION: (" + leftSucc + "," + rightSucc + ") with word " + wordInfo + ";");
									String cutInfo;
									if(leftSucc == probState) {
										if(succsForZ.get(rightSucc).size() > 3) {
											cutInfo = succsForZ.get(rightSucc).size() + " states";
										} else {
											cutInfo = succsForZ.get(rightSucc).toString();
										}
									} else {
										if(succsForZ.get(leftSucc).size() > 3) {
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
								Map<Integer,Set<Integer>> succsForZ_ = new HashMap<Integer,Set<Integer>>(succsForZ.size());
								List<Integer> y = curWord.getWord();
								for(Integer probState_ : IterableBitSet.getSetBits(scc)) {
									Set<Integer> newReach = new HashSet<Integer>();
									Set<Integer> succs = product.getProbStatesSuccessors(probState_, y);
									for(Integer succState : succs) {
										newReach.addAll(succsForZ.get(succState));
									}
									succsForZ_.put(probState_, newReach);
								}
								
								if (sanityCheck) {
									if(!succsForZ_.get(probState).containsAll(succsForZ.get(probState))) {
										throw new PrismException("The potential cut has got smaller.");
									}
								}
								
								succsForZ = succsForZ_;
							} else {
								// we have to continue, add new product state and curWord to BFS queue
								queue.add(new Pair<ProductState, SharedWord<Integer>>(new ProductState(leftSucc, rightSucc),
								                                                      curWord));
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
		
		Set<Integer> C = succsForZ.get(probState);
		
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

	private List<BitSet> computeSCCs(DTMCUBAProduct prod) throws PrismException {
		SCCConsumerStore sccs = new SCCConsumerStore();
		SCCComputer sccc = SCCComputer.createSCCComputer(this, prod, sccs);
		sccc.computeSCCs();
		return sccs.getSCCs();
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

		NBA uba = constructUBAForLTLFormula(mc, model, expr, labelBS);

		if (!getSettings().getBoolean(PrismSettings.PRISM_UBA_PURE) &&
		    uba.isDeterministic()) {
			return computeWithDA(model, uba, labelBS, statesOfInterest);
		}

			product = constructProduct(model, uba, labelBS, statesOfInterest);
		probs = computeValues(product);
		
		StateValues result = product.projectToOriginalModel(model, probs);

		return result;
	}
}
