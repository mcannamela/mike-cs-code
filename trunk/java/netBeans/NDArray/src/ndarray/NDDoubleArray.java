/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 *
 * @author mcannamela
 */
public class NDDoubleArray extends NDEntity{
    private double[] array;
    
    public static double[] arrayCopy(double[] original){
        double[] copy = new double[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    
    public NDDoubleArray(int[] shape) {
        super(shape);
        this.shape = idxCopy(shape);
        initNElements();
        array = new double[nElements];
        initStrides();
    }

    public NDDoubleArray(NDDoubleArray original) {
        super(original.shape);
        array = arrayCopy(original.getArray());
    }
    

    private void operate(elementOperator operator, 
                            NDCounter thisCounter, 
                            NDCounter otherCounter, 
                            NDCounter resultCounter){
        assert nElements()==thisCounter.nElements():"mismatched number of elements";
        assert thisCounter.nElements()==resultCounter.nElements():"mismatched number of elements";
        assert thisCounter.nElements()==otherCounter.nElements():"mismatched number of elements";
        
        int[] thisIndex, otherIndex, resultIndex;      
        
        if (otherCounter==null){
            otherIndex = null;
            while (thisCounter.hasNext()){
            thisIndex = thisCounter.next();
            resultIndex = resultCounter.next();
            operator.operate(thisIndex, otherIndex, resultIndex);
            }
        }
        else if (resultCounter==null){
            resultIndex=null;
            while (thisCounter.hasNext()){
            thisIndex = thisCounter.next();
            otherIndex = otherCounter.next();
            
            operator.operate(thisIndex, otherIndex, resultIndex);
            }
        }
        else {
            while (thisCounter.hasNext()){
            thisIndex = thisCounter.next();
            otherIndex = otherCounter.next();
            resultIndex = thisCounter.next();
            operator.operate(thisIndex, otherIndex, resultIndex);
            }
        }
        
        
    }
    private class AssignmentOperator extends elementOperator{
        public AssignmentOperator(NDDoubleArray other) {
            this.other = other;
            this.result= null;
        }
        
        @Override
        public void operate(int[] ownIndex, int[] otherIndex, int[] resultIndex) {
            array[nDTo1D(ownIndex)]=other.getElement(otherIndex);
        }
    }
    private class AccessOperator extends elementOperator{
        public AccessOperator(NDDoubleArray result) {
            this.result = result;
            this.other=null;
        }
        
        @Override
        public void operate(int[] thisIndex, int[] otherIndex, int[] resultIndex) {
            result.setElement(thisIndex, array[nDTo1D(thisIndex)]);
        }
    }
    private void assign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other){
        AssignmentOperator operator = new AssignmentOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }
    
    public NDDoubleArray get(ArrayList<Slice> rawSlices) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        NDDoubleArray slicedArray = new NDDoubleArray(ownCounter.shape());
        AccessOperator operator = new AccessOperator(slicedArray);
        operate(operator, ownCounter, null, slicedArray.getNewCounter());
        return slicedArray;
    }
    public void assign(ArrayList<Slice> rawSlices, NDDoubleArray other){
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        if (!Arrays.equals(ownCounter.shape(),other.shape())){
            throw new DimensionMismatchException("shape of slices does not match shape of array");
        }
        assign(ownCounter, other.getNewCounter(), other);
    }
    
    public void assign(ArrayList<Slice> rawSlices, NDDoubleArray other, ArrayList<Slice> rawOtherSlices){
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        SliceCounter otherCounter = other.getSliceCounter(rawOtherSlices);
        if (!Arrays.equals(ownCounter.shape(),otherCounter.shape())){
            throw new DimensionMismatchException("shape of slices do not match");
        }
        assign(ownCounter, otherCounter, other);
    }
    
    
    public final SliceCounter getSliceCounter(ArrayList<Slice> rawSlices) throws NoSuchElementException{
        assert rawSlices.size()==shape.length:"must provide a slice for each dimension";
        ArrayList<Slice> cookedSlices = new ArrayList<>(rawSlices.size());
        Slice cookedSlice;
        int dimension = 0;
        for (Slice slice:rawSlices){
            if (slice.isRaw()){
                cookedSlice = slice.getCookedSlice(shape[dimension]);
            }
            else{
                cookedSlice = slice;
            }
            cookedSlices.add(cookedSlice);
        }
        for (Slice slice:cookedSlices){
            validateSlice(slice, dimension);
        }
        return new SliceCounter(cookedSlices);
    }
    
    private void validateSlice(Slice slice, int dimension)  throws NoSuchElementException{
        if (slice.start()<0||slice.start()>shape[dimension]){
            throw new NoSuchElementException("slice start is "+ slice.start()
                    +"  in dimension "+dimension+" where the array has size "+shape[dimension]);
        }
        if (slice.stop()<0||slice.stop()>shape[dimension]){
            throw new NoSuchElementException("slice stop is "+ slice.stop()
                    +"  in dimension "+dimension+" where the array has size "+shape[dimension]);
        }
        if (slice.nIndices()==0){
            throw new NoSuchElementException("slice has no elements in dimension "
                        +dimension+": zero size arrays unsupported...for now");
        }
    }
    
    
    
    public NDCounter getNewCounter(){
        return new NDCounter(shape());
    }
    
    public double getElement(int index){
        return array[index];
    }
    public double getElement(int[] nDIndex){
        return getElement(nDTo1D(nDIndex));
    }

    public double[] getArray() {
        return array;
    }
    
    public void setElement(int index, double value){
        array[index]=value;
    }
    public void setElement(int[] nDIndex, double value){
        setElement(nDTo1D(nDIndex), value);
    }
    
}