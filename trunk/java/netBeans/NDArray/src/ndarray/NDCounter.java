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
public class NDCounter extends NDEntity {
    protected int[] nDIndex;
    protected int   index;

    public NDCounter() {
    }

    
    public NDCounter(int[] shape) {
        super(shape);
        initNDIndex();
    }
    
    protected final void initNDIndex(){
        nDIndex = new int[nDimensions()];
        for (int i=0;i<nDIndex.length;i++){
            nDIndex[i]=0;
        }
    }
    
    public final boolean hasNext(){
        return index<(nElements);
    }
    private boolean dimensionHasNext(int dimension){
        return nDIndex[dimension]<(shape[dimension]-1);
    }
    
    public final int[] next() throws NoSuchElementException{
        int dimCounter = nDimensions();
        int[] idx = getCurrentIndex();
        
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
    
    protected int[] getCurrentIndex(){
        return idxCopy(nDIndex);
    }
    
    
}
