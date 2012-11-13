/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

/**
 *
 * @author Michael
 */
abstract public class elementOperator<ARRAY_TYPE>{
    ARRAY_TYPE own;
    ARRAY_TYPE other;
    ARRAY_TYPE result;
    
    public elementOperator(){
    }

    public elementOperator(ARRAY_TYPE own, ARRAY_TYPE other, ARRAY_TYPE result) {
        this.own = own;
        this.other = other;
        this.result = result;
    }
    
    /*
     * carry out the operation on one element of the arrays
     * ownIndex, otherIndex, resultIndex - pre-flattened indices into the respective
     *                                      arrays
     */
    abstract public void operate(int ownIndex, int otherIndex, int resultIndex);
    
}
