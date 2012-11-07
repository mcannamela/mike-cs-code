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
class SliceList {
    int nDimensions;
    int[] shape;
    ArrayList<Slice> slices;

    SliceList(ArrayList<Slice> slices) {
        this.slices = slices;
        nDimensions = this.slices.size();
        shape = new int[nDimensions];
        for(int i=0;i<nDimensions;i++){
            shape[i]=this.slices.get(i).nIndices();
        }
        
    }


    public NDIterator iterator() {
        return new NDIterator(shape);
    }
    public int[] shape(){
        return NDEntity.idxCopy(shape);
    }
    public int[] getSlicePosition(int[] nDIndex){
        int[] slicePosition = new int[nDimensions];
        for(int i=0;i<nDimensions;i++){
            slicePosition[i] = slices.get(i).get(nDIndex[i]);
        }
        return slicePosition;
    }

    
}
