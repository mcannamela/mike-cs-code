/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 *
 * @author mcannamela
 */
public class NDDoubleArray extends NDArraySkeleton{
    private double[] array;
    
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
    
    private class AssignmentOperator extends elementOperator<NDDoubleArray>{
        public AssignmentOperator(NDDoubleArray other) {
            this.other = other;
            this.result= null;
        }
        
        @Override
        public void operate(int ownIndex, int otherIndex, int resultIndex) {
            array[ownIndex]=other.getElement(otherIndex);
        }
    }
    
    private class AccessOperator extends elementOperator<NDDoubleArray>{
        public AccessOperator(NDDoubleArray result) {
            this.result = result;
            this.other=null;
        }
        
        @Override
        public void operate(int thisIndex, int otherIndex, int resultIndex) {
            result.setElement(thisIndex, array[thisIndex]);
        }
    }
    
    private void assign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other){
        AssignmentOperator operator = new AssignmentOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }
    
    private void setBroadcasting(NDCounter ownCounter, NDCounter otherCounter){
        if (!Arrays.equals(ownCounter.shape(),otherCounter.shape())){
            int[] broadcastMask = makeBroadcastMask(otherCounter.shape);
            otherCounter.setBroadcasting(broadcastMask);
        }
    }
    public void assign(Slice[] rawSlices, NDDoubleArray other){
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        NDCounter otherCounter = other.getNewCounter();
        setBroadcasting(ownCounter, otherCounter);
        
        assign(ownCounter, other.getNewCounter(), other);
    }
    public void assign(Slice[] rawSlices, NDDoubleArray other, Slice[] rawOtherSlices){
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        SliceCounter otherCounter = other.getSliceCounter(rawOtherSlices);
        setBroadcasting(ownCounter, otherCounter);
        assign(ownCounter, otherCounter, other);
    }
    
    public NDDoubleArray get(Slice[] rawSlices) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        NDDoubleArray slicedArray = new NDDoubleArray(ownCounter.shape());
        AccessOperator operator = new AccessOperator(slicedArray);
        operate(operator, ownCounter, null, slicedArray.getNewCounter());
        return slicedArray;
    }
    
    public double getElement(int index){
        return array[index];
    }
    public double getElement(int[] nDIndex){
        return getElement(flattenIndex(nDIndex));
    }

    public double[] getArray() {
        return array;
    }
    
    public void setElement(int index, double value){
        array[index]=value;
    }
    public void setElement(int[] nDIndex, double value){
        setElement(flattenIndex(nDIndex), value);
    }
    public void setElement(int[] nDIndex, double value, IndexFlattener flattener){
        setElement(flattener.flatten(nDIndex), value);
    }
}
