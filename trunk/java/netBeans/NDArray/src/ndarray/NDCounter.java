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
    
    protected int[] nDCounterPosition;
    
    protected int   counterPosition;

    public NDCounter() {
        super();
    }

    public NDCounter(int[] shape) {
        super(shape);
        initCounterPosition();
    }
    
    protected final void initCounterPosition(){
        nDCounterPosition = new int[nDimensions()];
        for (int i=0;i<nDCounterPosition.length;i++){
            nDCounterPosition[i]=0;
        }
    }
    
    
    public final boolean hasNext(){
        return counterPosition<(nElements);
    }
    private boolean dimensionHasNext(int dimension){
        return nDCounterPosition[dimension]<(shape[dimension]-1);
    }
    
    /*
     * get the next index from the counter and increment it. 
     * returns - idx, and array of integers representing an n-dimensional counter 
     */
    public final int[] next() throws NoSuchElementException{
        int dimCounter = nDimensions();
        int[] idx = getCurrentIndex();
        
        if (hasNext()){
            recursiveIncrement(dimCounter-1);
            counterPosition++;
        }
        else{
            throw new NoSuchElementException("this iterator is exhausted");
        }
        return idx;
    }
    
    public final int nextFlat() throws NoSuchElementException{
        return this.flattenIndex(this.next());
    }
    
    
    /*
     * gives the internal position of the counter, which may differ from the 
     * index it returns with next
     */
    public final int[] getCounterPosition(){
        return idxCopy(nDCounterPosition);
    }
    
    /*
     * tells the current index into whatever the counter is counting for; 
     * by default it's the same as the counter position
     */
    protected int[] getCurrentIndex(){
        return idxCopy(nDCounterPosition);
    }
    
    
    private void recursiveIncrement(int activeDimension){
        if (dimensionHasNext(activeDimension)){
            nDCounterPosition[activeDimension]+=1;
        }
        else{
            nDCounterPosition[activeDimension] = 0;
            if (activeDimension>0){
                recursiveIncrement(activeDimension-1);
            }
        }
    }
    
    
    
    
    
}
