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
    protected Slice[] slices;
    
    SliceCounter(){
        super();
        this.slices = new Slice[]{new Slice(0,1)};
    }
    
    SliceCounter(int[] shape){
        super(shape);
        this.slices =new Slice[nDimensions()];
        for (int i=0;i<nDimensions();i++){
            this.slices[i]=new Slice(0,shape[i]);
        }
    }

    SliceCounter(Slice[] slices) {
        this.slices = slices;     
        shape = new int[this.slices.length];
        for(int i=0;i<nDimensions();i++){
            shape[i]=this.slices[i].nIndices();
        }
        initNElements();
        initStrides();
        initNDIndex();
    }

    @Override
    protected int[] getCurrentIndex(){
        int[] sliceIndex = newIndex();
        for (int i=0;i<nDimensions();i++){
            sliceIndex[i] = slices[i].get(nDIndex[i]);
        }
        return sliceIndex;
    }
}
