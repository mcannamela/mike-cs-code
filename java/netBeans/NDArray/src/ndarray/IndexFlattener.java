/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

/**
 *
 * @author Michael
 */
abstract public class IndexFlattener {
    int[] flattenerStrides;
    
    int flatten(int[] nDIndex){
        int idx = nDIndex[0]*flattenerStrides[0];
            for (int i = 1; i < flattenerStrides.length; i++) {
                idx += flattenerStrides[i] * nDIndex[i];
            }
            return idx;
    }
    
    int[] getFlattenerStrides(){
        return NDEntity.arrayCopy(flattenerStrides);
    }
}
