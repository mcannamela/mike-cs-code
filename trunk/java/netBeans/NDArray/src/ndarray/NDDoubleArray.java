/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.ArrayList;
import java.util.NoSuchElementException;

/**
 *
 * @author mcannamela
 */
public class NDDoubleArray extends NDEntity{
    
    private double[] array;
    private int nElements;
    private int[] strides; 

    public NDDoubleArray(int[] shape) {
        this.shape = idxCopy(shape);
        initNElements();
        array = new double[nElements];
        initStrides();
    }

    public NDDoubleArray(NDDoubleArray original) {
        throw new UnsupportedOperationException("Not yet implemented");
    }
    
    public int[] shape(){
        return idxCopy(shape);
    }
    public int nElements(){
        return nElements;
    }
    
    
    
    public int nDTo1D(int[] nDIndex){
        int idx = nDIndex[0];
        for(int i=1;i<strides.length;i++){
            idx+=strides[i]*nDIndex[i];
        }
        return idx;
    }
    public double getElement(int index){
        return array[index];
    }
    public double getElement(int[] nDIndex){
        return getElement(nDTo1D(nDIndex));
    }
    public void setElement(int index, double value){
        array[index]=value;
    }
    public void setElement(int[] nDIndex, double value){
        setElement(nDTo1D(nDIndex), value);
    }
    public NDDoubleArray getSlice(ArrayList<Slice> rawSlices) {
        SliceList slices = prepSliceList(rawSlices);
        NDIterator iterator = slices.iterator();
        NDDoubleArray slicedArray = new NDDoubleArray(slices.shape());
        int[] idx = new int[shape.length];
        int idx1d;
        int cnt = 0;
        while (iterator.hasNext()){
            idx = iterator.next();
            idx1d = nDTo1D(slices.getSlicePosition(idx));
            slicedArray.setElement(cnt, getElement(idx1d));
        }
        return slicedArray;
    }
    
    private SliceList prepSliceList(ArrayList<Slice> rawSlices) throws NoSuchElementException{
        assert rawSlices.size()==shape.length:"must provide a slice for each dimension";
        ArrayList<Slice> cookedSlices = new ArrayList<>(rawSlices.size());
        Slice cookedSlice;
        int dimension = 0;
        for (Slice slice:rawSlices){
            cookedSlice = slice.getFixedEndpointSlice(shape[dimension]);
            cookedSlices.add(cookedSlice);
        }
        for (Slice slice:cookedSlices){
            validateSlice(slice, dimension);
        }
        return new SliceList(cookedSlices);
    }
    private void validateSlice(Slice slice, int dimension)  throws NoSuchElementException{
        if (slice.start()<0){
            throw new NoSuchElementException("slice start is "+ slice.start()
                    +"  in dimension "+dimension);
        }
        if (slice.stop()>shape[dimension]){
            throw new NoSuchElementException("slice stop is "+ slice.stop()
                    +"  in dimension "+dimension+" where the array has size "+shape[dimension]);
        }
        if (slice.nIndices()==0){
            throw new NoSuchElementException("slice has no elements in dimension "
                        +dimension+": zero size arrays unsupported...for now");
        }
    }
    
    
    
    
    
    
}
