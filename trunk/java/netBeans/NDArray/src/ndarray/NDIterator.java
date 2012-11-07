/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.NoSuchElementException;

/**
 *
 * @author mcannamela
 */
public class NDIterator extends NDEntity {
    private int[] nDIndex;
    private int   index;

    public NDIterator(int[] shape) {
        this.shape = shape;
        initNDIndex();
        initStrides();
        initNElements();
    }
    
    private void initNDIndex(){
        nDIndex = new int[shape.length];
        for (int i=0;i<nDIndex.length;i++){
            nDIndex[i]=0;
        }
    }

    public void setStrides(int[] strides) {
        this.strides = strides;
    }
    public void setStrideAt(int index, int stride){
        strides[index] = stride;
    }
    
    public boolean hasNext(){
        return index<(nElements-1);
    }
    private boolean dimensionHasNext(int dimension){
        
        return nDIndex[dimension]<(shape[dimension]-1);
    }
    
    public int[] next() throws NoSuchElementException{
        int dimCounter = shape.length;
        int[] idx = idxCopy(nDIndex);
        
        if (hasNext()){
            recursiveIncrement(dimCounter-1);
            index++;
        }
        else{
            throw new NoSuchElementException("this iterator is exhausted");
        }
        
        return idx;
    }
    
    private void recursiveIncrement(int activeDimension){
        if (dimensionHasNext(activeDimension)){
            nDIndex[activeDimension]+=1;
        }
        else{
            nDIndex[activeDimension] = 0;
            if (activeDimension>0){
                recursiveIncrement(activeDimension-1);
            }
        }
    }

    int[] newIndex() {
        return new int[shape.length];
    }
    
}
