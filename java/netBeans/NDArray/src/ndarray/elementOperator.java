/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

/**
 *
 * @author Michael
 */
abstract public class elementOperator{
    NDDoubleArray own;
    NDDoubleArray other;
    NDDoubleArray result;
    
    public elementOperator(){
    }

    public elementOperator(NDDoubleArray own, NDDoubleArray other, NDDoubleArray result) {
        this.own = own;
        this.other = other;
        this.result = result;
    }
    
    abstract public void operate(int[] ownIndex, int[] otherIndex, int[] resultIndex);
    
}