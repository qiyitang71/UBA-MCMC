package explicit.gfg;

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
import jltl2dstar.GFG;
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
public class LTL2GFG extends PrismComponent
{

    public LTL2GFG(PrismComponent parent) throws PrismException
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
    //not supported
    public GFG convertLTLFormulaToUBA(Expression ltl, Values constantValues) throws PrismException
    {
        GFG result = null;
        return result;
    }

    public GFG convertLTLFormulaToUBAWithExternalTool(Expression ltl, Values constants)
            throws PrismException
    {
        GFG result = null;
        return result;
    }


    /**
     * Reads a HOA automaton from the given file and returns the parsed automaton.
     * Throws a PrismException if the automaton does not have one of the allowed acceptance conditions.
     */
    public GFG readHOA(String filename) throws PrismException
    {
        GFG result = null;
        try {
            mainLog.println("Reading HOA automaton from "+filename+"...");
            InputStream input = new FileInputStream(filename);
            result = HOAF2GFG.parseUBA(input);

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
    public GFG fromExpressionHOA(Expression expr, Path resourceDirectory, Vector<Expression> apExpressions) throws PrismException
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

        GFG nba = readHOA(automatonFile.getPath());

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
