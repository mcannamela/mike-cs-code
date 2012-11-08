/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.ArrayList;

/**
 *
 * @author mcannamela
 */
public class SliceCounter extends NDCounter{
    protected ArrayList<Slice> slices;
    

    SliceCounter(ArrayList<Slice> slices) {
        this.slices = slices;     
        shape = new int[this.slices.size()];
        for(int i=0;i<nDimensions();i++){
            shape[i]=this.slices.get(i).nIndices();
        }
        initNElements();
        initStrides();
        initNDIndex();
    }

    @Override
    protected int[] getCurrentIndex(){
        int[] sliceIndex = newIndex();
        for (int i=0;i<nDimensions();i++){
            sliceIndex[i] = slices.get(i).get(nDIndex[i]);
        }
        return sliceIndex;
    }
    
    
    

    
}
