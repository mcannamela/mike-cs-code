/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package earthmoverdistancesample;

import java.util.ArrayList;

/**
 *
 * @author wichtelwesen
 */
public interface AbstractSimple2DArray<T> {

    AbstractSimple2DArray elementWiseOperate(AbstractSimple2DArray<T> other, 
                                            BinaryOperatorInterface<T> operator);
    
    AbstractSimple2DArray elementWiseMultiply(AbstractSimple2DArray other);
    
    AbstractSimple2DArray elementWiseAdd(AbstractSimple2DArray other);
   
    ArrayList<T> flatten();
    
    AbstractSimple2DArray ravel(ArrayList<T> flatArray, int nRows, int nColumns);
    
    AbstractSimple2DArray ravel(ArrayList<T> flatArray);

    void setValueAt(int row, int column, T value);

    int[] shape();

    T sum();

    ArrayList<T> sumOfColumns();

    ArrayList<T> sumOfRows();

    T valueAt(int row, int column);
    
}
