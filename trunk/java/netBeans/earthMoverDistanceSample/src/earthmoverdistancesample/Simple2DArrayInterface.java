/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package earthmoverdistancesample;

/**
 *
 * @author wichtelwesen
 */
public interface Simple2DArrayInterface {

    Simple2DArrayInterface elementWiseOperate(Simple2DArrayInterface other, 
                                            doubleOperatorInterface operator);
    
    Simple2DArrayInterface elementWiseMultiply(Simple2DArrayInterface other);
    
    Simple2DArrayInterface elementWiseAdd(Simple2DArrayInterface other);
   
    double[] flatten();

    void setValueAt(int row, int column, double value);

    int[] shape();

    double sum();

    double[] sumOfColumns();

    double[] sumOfRows();

    double valueAt(int row, int column);
    
}
