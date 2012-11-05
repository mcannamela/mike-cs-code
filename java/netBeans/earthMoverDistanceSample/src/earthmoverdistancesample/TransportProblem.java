package earthmoverdistancesample;

import java.util.ArrayList;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.linear.LinearConstraint;
import org.apache.commons.math3.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optimization.linear.Relationship;
import org.apache.commons.math3.optimization.linear.SimplexSolver;

/**
 *
 * @author mcannamela
 */
public class TransportProblem implements Runnable{
    private Simple2DArrayInterface costMatrix;
    private Signature supplySignature;
    private Signature demandSignature;
    
    ArrayList<LinearConstraint> constraints = new ArrayList<>();
    SimplexSolver solver = new SimplexSolver();
    LinearObjectiveFunction objectiveFunction;
    private static final boolean restrictNonNegative = true;
    
    private PointValuePair solution;

    
    public TransportProblem(Signature supplySignature, Signature demandSignature) {    
        this.supplySignature = supplySignature;
        this.demandSignature = demandSignature;
        
        setCostMatrix();
    }
    
    @Override
    public void run() {
        makeLinearConstraints();
        makeObjectiveFunction();
        solution = solver.optimize(objectiveFunction, constraints, GoalType.MINIMIZE, restrictNonNegative);
    }

    public PointValuePair getSolution() {
        return solution;
    }
    
    private void setCostMatrix(){
        costMatrix = supplySignature.distanceMatrix(demandSignature);
    }
    private void validateSupply(){
        if (supplySignature.getWeightsNorm()<demandSignature.getWeightsNorm()){
            supplySignature.normalizeWeights(demandSignature.getWeightsNorm());
        }
    }
    
    private double[] getObjectiveFunctionCoeffs(){
        ArrayList<Double> coeffList = costMatrix.flatten();
        double[] coeffs;
        coeffs = Simple2DDoubleArray.ArrayListToDoubleArray(coeffList);
        return coeffs;
    } 
    
    public void setDemandSignature(Signature demandSignature) {
        this.demandSignature = demandSignature;
        setCostMatrix();
        validateSupply();
    }

    public void setSupplySignature(Signature supplySignature) {
        this.supplySignature = supplySignature;
        setCostMatrix();
        validateSupply();
    }
    
    private void makeObjectiveFunction(){
        double[] coeffs = getObjectiveFunctionCoeffs();
        double constantTerm = 0.0;
        objectiveFunction = new LinearObjectiveFunction(coeffs,constantTerm);
    }
    
    private void makeLinearConstraints(){
        constraints.clear();
        int nRows,nColumns;
        
        LinearConstraint currentConstraint;
        double[] lhsCoeffs;
        
        nRows = costMatrix.shape()[0];
        for (int i=0;i<nRows;i++){
            lhsCoeffs = Simple2DDoubleArray.ArrayListToDoubleArray(
                    costMatrix.onesInRow(i).flatten()  );
            currentConstraint = new LinearConstraint(lhsCoeffs, Relationship.LEQ,supplySignature.weightAt(i));
            constraints.add(currentConstraint);
        }
        
        nColumns = costMatrix.shape()[1];
        for (int j=0;j<nColumns;j++){
            lhsCoeffs = Simple2DDoubleArray.ArrayListToDoubleArray(
                    costMatrix.onesInColumn(j).flatten()  );
            currentConstraint = new LinearConstraint(lhsCoeffs, Relationship.EQ,demandSignature.weightAt(j));
            constraints.add(currentConstraint);
        }
    }
    
}
