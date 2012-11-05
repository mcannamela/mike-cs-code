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
public interface Simple2DArrayInterface<T> {

    Simple2DArrayInterface elementWiseOperate(Simple2DArrayInterface<T> other, 
                                            BinaryOperatorInterface<T> operator);
    
    Simple2DArrayInterface elementWiseMultiply(Simple2DArrayInterface other);
    
    Simple2DArrayInterface elementWiseAdd(Simple2DArrayInterface other);
    
    Simple2DDoubleArray onesInRow(int row);
    Simple2DDoubleArray onesInColumn(int column);
   
    ArrayList<T> flatten();
    
    Simple2DArrayInterface ravel(ArrayList<T> flatArray, int nRows, int nColumns);
    
    Simple2DArrayInterface ravel(ArrayList<T> flatArray);

    void setValueAt(int row, int column, T value);

    int[] shape();

    T sum();

    ArrayList<T> sumOfColumns();

    ArrayList<T> sumOfRows();

    T valueAt(int row, int column);
    
}
