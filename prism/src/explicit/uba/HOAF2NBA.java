

package explicit.uba;

import java.io.FileInputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import prism.PrismException;
import jhoafparser.ast.AtomAcceptance;
import jhoafparser.ast.AtomLabel;
import jhoafparser.ast.BooleanExpression;
import jhoafparser.consumer.HOAConsumer;
import jhoafparser.consumer.HOAConsumerException;
import jhoafparser.parser.HOAFParser;
import jhoafparser.parser.generated.ParseException;
import jhoafparser.util.ImplicitEdgeHelper;
import jltl2ba.APSet;
import jltl2dstar.APMonom;
import jltl2dstar.NBA;


/**
 * A HOAConsumer for jhoafparser that constructs a jltl2dstar.NBA from the parsed automaton.
 * <br>
 * The automaton has to be deterministic and complete, with state-based acceptance and
 * labels (explicit/implicit) on the edges.
 * <br>
 * There are (currently) more restrictions on the automaton:
 * <ul>
 * <li>The Start and States headers have to be present</li>
 * <li>At least one state in the automaton.
 * <li>All explicit edge labels have to be in disjunctive normal form (disjunction of conjunctive clauses)</li>
 * <li>At most 30 atomic propositions</li>
 * </ul>
 */
public class HOAF2NBA implements HOAConsumer {

	/** The resulting nondeterministic Buechi automaton */
	private NBA nba;
	/** The set of atomic propositions of the automaton (in APSet form) */
	private APSet aps = new APSet();

	/** Size, i.e. number of states */
	private int size;
	/** Do we know the number of states? Is provided by the optional HOA States-header */
	private boolean knowSize = false;

	/** Start state (index) */
	private int startState;
	/** Do we know the start state? Is provided by the HOA Start-header */
	private boolean knowStartState = false;

	/** The condition name from the acc-name header (optional) */
	private String accName;

	/** The helper for handling implicit edges */
	private ImplicitEdgeHelper implicitEdgeHelper = null;

	/** The properties of the automaton */
	Set<String> properties = new TreeSet<String>();

	/** Clear the various state information */
	public void clear() {
		aps = new APSet();

		implicitEdgeHelper = null;

		size = 0;
		knowSize = false;

		startState = 0;
		knowStartState = false;

		accName = null;
	}

	/** Constructor */
	public HOAF2NBA() {
	}

	@Override
	public boolean parserResolvesAliases() {
		return true;
	}

	@Override
	public void notifyHeaderStart(String version) throws HOAConsumerException {
		// NOP
	}

	@Override
	public void setNumberOfStates(int numberOfStates)
			throws HOAConsumerException {
		size = numberOfStates;
		knowSize = true;
		if (numberOfStates == 0) {
			throw new HOAConsumerException("Automaton with zero states, need at least one state");
		}
	}

	@Override
	public void addStartStates(List<Integer> stateConjunction)
			throws HOAConsumerException {
		if(stateConjunction.size() > 1 || knowStartState) {
			throw new HOAConsumerException("Currently, only a single initial state is supported");
		}
		startState = stateConjunction.get(0).intValue();
		knowStartState = true;
	}

	@Override
	public void addAlias(String name, BooleanExpression<AtomLabel> labelExpr)
			throws HOAConsumerException {
		// NOP, aliases are already resolved
	}

	@Override
	public void setAPs(List<String> aps) throws HOAConsumerException {
		if (aps.size() > 30) {
			throw new HOAConsumerException("Automaton has "+aps.size()+" atomic propositions, at most 30 are supported");
		}

		for (String ap : aps) {
			this.aps.addAP(ap);
		}
		
		this.aps.powersetSize();
	}

	@Override
	public void setAcceptanceCondition(int numberOfSets,
			BooleanExpression<AtomAcceptance> accExpr)
			throws HOAConsumerException {
		// ignore
	}

	@Override
	public void provideAcceptanceName(String name, List<Object> extraInfo)
			throws HOAConsumerException {
		accName = name;
		if (!name.equals("Buchi")) {
			throw new HOAConsumerException("Can only parse Buchi automata where acc-name is provided");
		}
	}

	@Override
	public void setName(String name) throws HOAConsumerException {
		// NOP
	}

	@Override
	public void setTool(String name, String version) throws HOAConsumerException {
		// NOP
	}

	@Override
	public void addProperties(List<String> properties)
			throws HOAConsumerException {
		if(properties.contains("univ-branch")) {
			throw new HOAConsumerException("A HOAF with universal branching is not deterministic");
		}
		
		if(properties.contains("state-labels")) {
			throw new HOAConsumerException("Can't handle state labelling");
		}

		// store properties
		this.properties.addAll(properties);
	}

	@Override
	public void addMiscHeader(String name, List<Object> content)
			throws HOAConsumerException {
		if (name.substring(0,1).toUpperCase().equals(name.substring(0,1))) {
			throw new HOAConsumerException("Unknown header "+name+" potentially containing semantic information, can not handle");
		}
	}

	@Override
	public void notifyBodyStart() throws HOAConsumerException {
		if (!knowSize) {
			throw new HOAConsumerException("Can currently only parse automata where the number of states is specified in the header");
		}
		if (!knowStartState) {
			throw new HOAConsumerException("Can currently only parse automata where the initial state specified (Start header)");
		}
		if (startState >= size) {
			throw new HOAConsumerException("Initial state "+startState+" is out of range");
		}

		if (accName == null) {
			throw new HOAConsumerException("Can only parse Buchi automaton where acc-name is provided");
		}

		nba = new NBA(aps);
		for (int i = 0; i< size; i++) {
			nba.nba_i_newState();
		}
		nba.nba_i_setStartState(startState);

		implicitEdgeHelper = new ImplicitEdgeHelper(aps.size());
	}


	@Override
	public void addState(int id, String info,
			BooleanExpression<AtomLabel> labelExpr, List<Integer> accSignature)
			throws HOAConsumerException {
		implicitEdgeHelper.startOfState(id);

		if(labelExpr != null) {
			throw new HOAConsumerException("State "+id+" has a state label, currently only supports labels on transitions");
		}
		
		if (id >= size) {
			throw new HOAConsumerException("Illegal state index "+id+", out of range");
		}

		nba.nba_i_setFinal(id, false);
		if (accSignature != null) {
			if (accSignature.contains(0)) {
				nba.nba_i_setFinal(id, true);
			}
		}
	}

	@Override
	public void addEdgeImplicit(int stateId, List<Integer> conjSuccessors,
			List<Integer> accSignature) throws HOAConsumerException {
		if (conjSuccessors.size() != 1) {
			throw new HOAConsumerException("Not an NBA, state "+stateId+" has transition with conjunctive target");
		}

		if (accSignature != null) {
			throw new HOAConsumerException("NBA has transition-based acceptance (state "+stateId+", currently only state-labeled acceptance is supported");
		}

		int to = conjSuccessors.get(0);

		APMonom m = new APMonom();
		long tmp = implicitEdgeHelper.nextImplicitEdge();
		int index = 0;
		while (tmp != 0) {
			if (tmp % 2 == 1) {
				m.andAP(index, true);
			} else {
				m.andAP(index, false);
			}
			tmp = tmp >> 1L;
			index++;
		}
		nba.nba_i_addEdge(stateId, m, to);
	}
	
	/**
	 * Returns a list of APMonoms for the expression. The expression currently has to be in
	 * disjunctive normal form. Returns one APMonom for each clause of the DNF.
	 */
	private List<APMonom> labelExpressionToAPMonom(BooleanExpression<AtomLabel> expr) throws HOAConsumerException {
		List<APMonom> result = new ArrayList<APMonom>();
		
		switch (expr.getType()) {
		case EXP_AND:
		case EXP_ATOM:
		case EXP_NOT: {
			APMonom monom = new APMonom();
			labelExpressionToAPMonom(expr, monom);
			result.add(monom);
			return result;
		}
		case EXP_TRUE:
			result.add(new APMonom(true));
			return result;
		case EXP_FALSE:
			result.add(new APMonom(false));
		case EXP_OR:
			result.addAll(labelExpressionToAPMonom(expr.getLeft()));
			result.addAll(labelExpressionToAPMonom(expr.getRight()));
			return result;
		}
		throw new UnsupportedOperationException("Unsupported operator in label expression: "+expr);
	}

	
	/**
	 * Returns a single APMonom for a single clause of the overall DNF formula.
	 * Modifies APMonom result such that in the end it is correct.
	 */
	private void labelExpressionToAPMonom(BooleanExpression<AtomLabel> expr, APMonom result) throws HOAConsumerException {
		try {
			switch (expr.getType()) {
			case EXP_TRUE:
			case EXP_FALSE:
			case EXP_OR:
				throw new HOAConsumerException("Complex transition labels are not yet supported, only disjunctive normal form: "+expr);

			case EXP_AND:
				labelExpressionToAPMonom(expr.getLeft(), result);
				labelExpressionToAPMonom(expr.getRight(), result);
				return;
			case EXP_ATOM: {
				int apIndex = expr.getAtom().getAPIndex();
				if (result.isSet(apIndex) && result.getValue(apIndex)!=true) {
					throw new HOAConsumerException("Complex transition labels are not yet supported, transition label evaluates to false");
				}
				result.setValue(apIndex, true);
				return;
			}
			case EXP_NOT: {
				if (!expr.getLeft().isAtom()) {
					throw new HOAConsumerException("Complex transition labels are not yet supported, only conjunction of (negated) labels");
				}
				int apIndex = expr.getLeft().getAtom().getAPIndex();
				if (result.isSet(apIndex) && result.getValue(apIndex)!=false) {
					throw new HOAConsumerException("Complex transition labels are not yet supported, transition label evaluates to false");
				}
				result.setValue(apIndex, false);
				return;
			}
			}
		} catch (PrismException e) {
			throw new HOAConsumerException("While parsing, APMonom exception: "+e.getMessage());
		}
	}
	
	@Override
	public void addEdgeWithLabel(int stateId,
			BooleanExpression<AtomLabel> labelExpr,
			List<Integer> conjSuccessors, List<Integer> accSignature)
			throws HOAConsumerException {

		if (conjSuccessors.size() != 1) {
			throw new HOAConsumerException("Not an NBA, state "+stateId+" has transition with conjunctive target");
		}

		if (accSignature != null) {
			throw new HOAConsumerException("NBA has transition-based acceptance (state "+stateId+", currently only state-labeled acceptance is supported");
		}

		if (labelExpr == null) {
			throw new HOAConsumerException("Missing label on transition");
		}

		int to = conjSuccessors.get(0);

		for (APMonom monom : labelExpressionToAPMonom(labelExpr)) {
			nba.nba_i_addEdge(stateId, monom, to);
		}
	}

	@Override
	public void notifyEndOfState(int stateId) throws HOAConsumerException
	{
		implicitEdgeHelper.endOfState();
	}

	@Override
	public void notifyEnd() throws HOAConsumerException {
		clear();
	}

	@Override
	public void notifyAbort() {
		clear();
		
	}
	
	public NBA getNBA() {
		return nba;
	}
	
	public Set<String> getProperties() {
		return properties;
	}

	@Override
	public void notifyWarning(String warning) throws HOAConsumerException
	{
		// warnings are fatal
		throw new HOAConsumerException(warning);
	}

	/**
	 * Parses an NBA HOA automaton from the input stream, returns the automaton
	 * and stores all attached properties in the propertyStore 
	 * @param input the input stream
	 * @param propertyStore a Set for receiving the properties of the automaton (may be {@code null})
	 * @return the automaton
	 */
	public static NBA parseNBA(InputStream input, Set<String> propertyStore) throws ParseException {
		HOAF2NBA consumerNBA = new HOAF2NBA();
		HOAFParser.parseHOA(input, consumerNBA);
		
		if (propertyStore != null) {
			propertyStore.addAll(consumerNBA.getProperties());
		}
		return consumerNBA.getNBA();
	}

	/**
	 * Parses an UBA HOA automaton from the input stream, returns the automaton.
	 * If the automaton does not have the 'unambiguous' property, throws an exception.
	 * @param input the input stream
	 * @return the automaton
	 */
	public static NBA parseUBA(InputStream input) throws ParseException {
		Set<String> propertyStore = new TreeSet<String>();

		NBA nba = parseNBA(input, propertyStore);
		if (!propertyStore.contains("unambiguous") && !propertyStore.contains("deterministic")) {
			throw new ParseException("Automaton is required to be marked as 'unambiguous' or 'deterministic' HERE IN UBA");
		}
	
		return nba;
	}

	/** Command-line interface for reading, parsing and printing a HOA automaton (for testing) */
	public static void main(String args[])
	{
		int rv = 0;
		InputStream input = null;
		try {
			if (args.length != 2) {
				System.err.println("Usage: input-file output-file\n\n Filename can be '-' for standard input/output");
				System.exit(1);
			}
			if (args[0].equals("-")) {
				input = System.in;
			} else {
				input = new FileInputStream(args[0]);
			}

			PrintStream output;
			String outfile = args[1];
			if (outfile.equals("-")) {
				output = System.out;
			} else {
				output = new PrintStream(outfile);
			}

			Set<String> properties = new TreeSet<String>();
			NBA result = parseNBA(input, properties);
			if (result == null) {
				throw new PrismException("Could not construct NBA");
			}

			// should we try and simplify?
			// result = DASimplifyAcceptance.simplifyAcceptance(result, acceptance.AcceptanceType.REACH);

			result.print_hoa(output);

			output.println("\n----- Properties: ---------\n");
			for (String p : properties) {
				output.println(p);
			}
		} catch (Exception e) {
			System.err.println(e.toString());
			rv = 1;
		}
		
		if (rv != 0) {
			System.exit(rv);
		}
	}
}
