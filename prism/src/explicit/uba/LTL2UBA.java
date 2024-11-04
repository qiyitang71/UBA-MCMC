
package explicit.uba;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import common.PathUtil;

import jhoafparser.parser.generated.ParseException;
import jltl2ba.APSet;
import jltl2ba.SimpleLTL;
import jltl2dstar.NBA;
import parser.Values;
import parser.ast.Expression;
import parser.visitor.ExpandStepBoundsSyntactically;
import parser.ast.ExpressionHOA;
import parser.ast.ExpressionLabel;
import parser.ast.ExpressionUnaryOp;
import prism.Pair;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismNotSupportedException;
import prism.PrismSettings;

/**
 * Infrastructure for constructing unambiguous Buchi automata (UBA) for LTL formulas.
 */
public class LTL2UBA extends PrismComponent
{

	public LTL2UBA(PrismComponent parent) throws PrismException
	{
		super(parent);
	}

	/**
	 * Convert an LTL formula into an unambiguous Buchi automaton.
	 * The LTL formula is represented as a PRISM Expression,
	 * in which atomic propositions are represented by ExpressionLabel objects.
	 * @param ltl the formula
	 * @param constantValues the values of constants, may be {@code null}
	 */
	public NBA convertLTLFormulaToUBA(Expression ltl, Values constantValues) throws PrismException
	{
		NBA result = null;

		boolean useExternal = useExternal();
		if (!useExternal) {
			throw new PrismNotSupportedException("Internal UBA generation not supported, use external tool");
		}

		boolean containsTemporalBounds = Expression.containsTemporalTimeBounds(ltl);

		if (containsTemporalBounds) {
			// remove time bounds syntactically
			ExpandStepBoundsSyntactically visitor = new ExpandStepBoundsSyntactically(constantValues);
			ltl = (Expression)ltl.deepCopy().accept(visitor);

			mainLog.println("LTL formula with syntactically expanded bounds: "+ltl);
		}
		result = convertLTLFormulaToUBAWithExternalTool(ltl, constantValues);

		if (result == null) {
			throw new PrismNotSupportedException("Could not convert LTL formula to UBA");
		}

		return result;
	}

	public NBA convertLTLFormulaToUBAWithExternalTool(Expression ltl, Values constants)
			throws PrismException
	{
		// reuse LTL2DA tool setting
		// TODO: extra LTL2UBA setting?
		String ltl2ubaTool = getSettings().getString(PrismSettings.PRISM_LTL2DA_TOOL);

		SimpleLTL ltlFormula = ltl.convertForJltl2ba();

		// switch from the L0, L1, ... APs of PRISM to the
		// safer p0, p1, ... APs for the external tool
		SimpleLTL ltlFormulaSafeAP = ltlFormula.clone();
		ltlFormulaSafeAP.renameAP("L", "p");

		NBA result = null;

		try {

			String syntax = getSettings().getString(PrismSettings.PRISM_LTL2DA_SYNTAX);
			if (syntax == null || syntax.isEmpty()) {
				throw new PrismException("No LTL syntax option provided");
			}
			String ltlOutput;
			switch (syntax) {
			case "LBT":
				ltlOutput = ltlFormulaSafeAP.toStringLBT();
				break;
			case "Spin":
				ltlOutput = ltlFormulaSafeAP.toStringSpin();
				break;
			case "Spot":
				ltlOutput = ltlFormulaSafeAP.toStringSpot();
				break;
			case "Rabinizer":
				ltlFormulaSafeAP = ltlFormulaSafeAP.toBasicOperators();
				ltlOutput = ltlFormulaSafeAP.toStringSpot();
				break;
			default:
				throw new PrismException("Unknown LTL syntax option \"" + syntax + "\"");
			}

			File ltl_file = File.createTempFile("prism-ltl-external-", ".ltl", null);
			File uba_file = File.createTempFile("prism-ltl-external-", ".hoa", null);
			File tool_output = File.createTempFile("prism-ltl-external-", ".output", null);

			FileWriter ltlWriter = new FileWriter(ltl_file);
			ltlWriter.write(ltlOutput);
			ltlWriter.close();

			List<String> arguments = new ArrayList<String>();
			arguments.add(ltl2ubaTool);

			getLog().print("Calling external LTL->UBA tool: ");
			for (String s : arguments) {
				getLog().print(" " + s);
			}
			getLog().println();

			getLog().print("LTL formula (in " + syntax + " syntax):  ");
			getLog().println(ltlOutput);
			getLog().println();

			arguments.add(ltl_file.getAbsolutePath());
			arguments.add(uba_file.getAbsolutePath());

			ProcessBuilder builder = new ProcessBuilder(arguments);
			builder.redirectOutput(tool_output);
			builder.redirectErrorStream(true);

			// if we are running under the Nailgun environment, setup the
			// environment to include the environment variables of the Nailgun client
			prism.PrismNG.setupChildProcessEnvironment(builder);

			Process p = builder.start();
			p.getInputStream().close();

			int rv;
			while (true) {
				try {
					rv = p.waitFor();
					break;
				} catch (InterruptedException e) {
				}
			}
			if (rv != 0) {
				throw new PrismException("Call to external LTL->UBA tool failed, return value = " + rv + ".\n"
						+ "To investigate, please consult the following files:" + "\n LTL formula:                     " + ltl_file.getAbsolutePath()
						+ "\n Automaton output:                " + uba_file.getAbsolutePath() + "\n Tool output (stdout and stderr): "
						+ tool_output.getAbsolutePath() + "\n");
			}

			tool_output.delete();

			try {
				result = HOAF2NBA.parseUBA(new FileInputStream(uba_file));

				if (result == null) {
					throw new PrismException("Could not construct UBA");
				}
				checkAPs(ltlFormulaSafeAP, result.getAPSet());

				// rename back from safe APs, i.e., p0, p1, ... to L0, L1, ...
				APSet automatonAPs = result.getAPSet();
				APSet newAPs = new APSet();
				for (int i = 0; i < automatonAPs.size(); i++) {
					if (automatonAPs.getAP(i).startsWith("p")) {
						String renamed = "L" + automatonAPs.getAP(i).substring("p".length());
						newAPs.addAP(renamed);
					}
				}
				result.switchAPSet(newAPs);
			} catch (ParseException e) {
				throw new PrismException("Parse error: " + e.getMessage() + ".\n" + "To investigate, please consult the following files:\n"
						+ " LTL formula:        " + ltl_file.getAbsolutePath() + "\n Automaton output: " + uba_file.getAbsolutePath() + "\n");
			} catch (PrismException e) {
				throw new PrismException(e.getMessage() + ".\n" + "To investigate, please consult the following files:" + "\n LTL formula: "
						+ ltl_file.getAbsolutePath() + "\n Automaton output: " + uba_file.getAbsolutePath() + "\n");
			}

			uba_file.delete();
			ltl_file.delete();
		} catch (IOException e) {
			throw new PrismException(e.getMessage());
		}


		return result;
	}


	/**
	 * Reads a HOA automaton from the given file and returns the parsed automaton.
	 * Throws a PrismException if the automaton does not have one of the allowed acceptance conditions.
	 */
	public NBA readHOA(String filename) throws PrismException
	{
		NBA result = null;
		try {
			mainLog.println("Reading HOA automaton from "+filename+"...");
			InputStream input = new FileInputStream(filename);
			result = HOAF2NBA.parseUBA(input);

			if (result == null) {
				throw new PrismException("Could not construct UBA");
			}
		} catch (ParseException e) {
			throw new PrismException("Reading from " + filename + ", parse error: "+e.getMessage());
		} catch (PrismException e) {
			throw new PrismException("Reading from " + filename + ", " + e.getMessage());
		} catch (IOException e) {
			throw new PrismException("I/O error: " + e.getMessage());
		}

		return result;
	}

	/**
	 * Handles a ExpressionHOA expression, possibly negated.
	 * Reads the HOA automaton from the specified file and returns the parsed automaton.
	 *
	 * If the expression is of the form "! HOA: {....}" will attempt to negate the automation
	 * via an appropriate switch of the acceptance condition.
	 *
	 * Throws a PrismException if the automaton does not have one of the allowed acceptance conditions.
	 * @param expr the expression
	 * @param resourceDirectory path to the directory that is used for loading automata with relative paths ({@code null} = current working dir)
	 */
	public NBA fromExpressionHOA(Expression expr, Path resourceDirectory, Vector<Expression> apExpressions) throws PrismException
	{
		boolean negate = false;
		while (Expression.isNot(expr) || Expression.isParenth(expr)) {
			if (Expression.isNot(expr)) {
				negate = !negate;
				expr = ((ExpressionUnaryOp)expr).getOperand();
			} else if (Expression.isParenth(expr)) {
				expr = ((ExpressionUnaryOp)expr).getOperand();
			} else {
				throw new PrismException("implementation error");
			}
		}
		if (!(expr instanceof ExpressionHOA)) {
			throw new PrismException("Expected HOA path automaton expression, found "+expr);
		}

		if (negate) {
			throw new PrismNotSupportedException("Negated HOA path formulas not supported via UBA");
		}

		ExpressionHOA exprHOA = (ExpressionHOA) expr;

		File automatonFile = new File(exprHOA.getAutomatonFile().getText());
		if (!automatonFile.isAbsolute()) {
			if (resourceDirectory == null) {
				mainLog.printWarning("Loading automaton file with relative path relative to current working directory...");
			}
			automatonFile = PathUtil.resolvePath(resourceDirectory, automatonFile.toPath()).toFile();
		}

	    NBA nba = readHOA(automatonFile.getPath());

		assert(apExpressions.size() == 0);
		for (String ap : nba.getAPSet()) {
			// default: ap corresponds to label expression
			apExpressions.add(new ExpressionLabel(ap));
		}

		// handle renames
		for (Pair<String, Expression> rename : exprHOA.getRenames()) {
			int i = nba.getAPSet().indexOf(rename.getKey());

			if (i == -1) {
				throw new PrismException("Can not rename label \"" + rename.getKey() +"\", not an atomic proposition of the automaton");
			} else {
				apExpressions.set(i, rename.getValue());
			}
		}

		return nba;
	}

	/** Check whether we should use an external LTL->DA/UBA tool */
	private boolean useExternal()
	{
		String ltl2da_tool = getSettings().getString(PrismSettings.PRISM_LTL2DA_TOOL);
		if (ltl2da_tool != null && !ltl2da_tool.isEmpty()) {
			return true;
		}
		return false;
	}

	/** Check the atomic propositions of the (externally generated) automaton */
	private void checkAPs(SimpleLTL ltl, APSet automatonAPs) throws PrismException
	{
		APSet ltlAPs = ltl.getAPs();
		for (String ap : automatonAPs) {
			if (!ltlAPs.hasAP(ap)) {
				throw new PrismException("Generated automaton has extra atomic proposition \"" + ap + "\"");
			}
		}
		// It's fine for the automaton to not have APs that occur in the formula, e.g., for
		// p0 | !p0, the external tool could simplify to 'true' and omit all APs
	}

}
