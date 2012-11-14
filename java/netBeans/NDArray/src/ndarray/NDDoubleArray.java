/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.Arrays;

/**
 *
 * @author mcannamela
 */
public class NDDoubleArray extends NDArraySkeleton {

    private double[] array;

    public NDDoubleArray(int[] shape) {
        super(shape);
        this.shape = arrayCopy(shape);
        initNElements();
        array = new double[nElements];
        initStrides();
    }

    public NDDoubleArray(NDDoubleArray original) {
        super(original.shape);
        array = arrayCopy(original.getArray());
    }

    private class AssignmentOperator extends elementOperator<NDDoubleArray> {

        public AssignmentOperator(NDDoubleArray other) {
            this.other = other;
            this.result = null;
        }

        @Override
        public void operate(int ownIndex, int otherIndex, int resultIndex) {
            array[ownIndex] = other.getElement(otherIndex);
        }
    }

    private class TimesEqualsOperator extends elementOperator<NDDoubleArray> {

        public TimesEqualsOperator(NDDoubleArray other) {
            this.other = other;
            this.result = null;
        }

        @Override
        public void operate(int ownIndex, int otherIndex, int resultIndex) {
            array[ownIndex] *= other.getElement(otherIndex);
        }
    }

    private class DivideEqualsOperator extends elementOperator<NDDoubleArray> {

        public DivideEqualsOperator(NDDoubleArray other) {
            this.other = other;
            this.result = null;
        }

        @Override
        public void operate(int ownIndex, int otherIndex, int resultIndex) {
            array[ownIndex] /= other.getElement(otherIndex);
        }
    }

    private class PlusEqualsOperator extends elementOperator<NDDoubleArray> {

        public PlusEqualsOperator(NDDoubleArray other) {
            this.other = other;
            this.result = null;
        }

        @Override
        public void operate(int ownIndex, int otherIndex, int resultIndex) {
            array[ownIndex] += other.getElement(otherIndex);
        }
    }

    private class MinusEqualsOperator extends elementOperator<NDDoubleArray> {

        public MinusEqualsOperator(NDDoubleArray other) {
            this.other = other;
            this.result = null;
        }

        @Override
        public void operate(int ownIndex, int otherIndex, int resultIndex) {
            array[ownIndex] -= other.getElement(otherIndex);
        }
    }

    private class AccessOperator extends elementOperator<NDDoubleArray> {

        public AccessOperator(NDDoubleArray result) {
            this.result = result;
            this.other = null;
        }

        @Override
        public void operate(int thisIndex, int otherIndex, int resultIndex) {
            result.setElement(thisIndex, array[thisIndex]);
        }
    }

    private void assign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other) {
        AssignmentOperator operator = new AssignmentOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }
    private void timesAssign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other) {
        TimesEqualsOperator operator = new TimesEqualsOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }
    private void divideAssign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other) {
        DivideEqualsOperator operator = new DivideEqualsOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }
    private void plusAssign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other) {
        PlusEqualsOperator operator = new PlusEqualsOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }
    private void minusAssign(NDCounter ownCounter, NDCounter otherCounter, NDDoubleArray other) {
        MinusEqualsOperator operator = new MinusEqualsOperator(other);
        operate(operator, ownCounter, otherCounter, null);
    }

    private static void setAssignmentBroadcasting(NDCounter ownCounter, NDCounter otherCounter) {
        if (!Arrays.equals(ownCounter.shape(), otherCounter.shape())) {
            int[] broadcastMask = makeBroadcastMask(ownCounter.shape(), otherCounter.shape);
            otherCounter.setBroadcasting(broadcastMask);
        }
    }

    private static void setBinaryBroadcasting(NDCounter ownCounter, NDCounter otherCounter) {
        int[] targetShape, ownMask, otherMask;
        if (!Arrays.equals(ownCounter.shape(), otherCounter.shape())) {
            targetShape = negotiateTargetShape(ownCounter.shape(), otherCounter.shape);
            ownMask = makeBroadcastMask(targetShape, ownCounter.shape);
            otherMask = makeBroadcastMask(targetShape, otherCounter.shape);
            ownCounter.setBroadcasting(ownMask);
            otherCounter.setBroadcasting(otherMask);
        }
    }

    public void assign(NDDoubleArray other) {
        NDCounter ownCounter = getNewCounter();
        NDCounter otherCounter = other.getNewCounter();
        setAssignmentBroadcasting(ownCounter, otherCounter);
        assign(ownCounter, otherCounter, other);
    }

    public void assign(Slice[] rawSlices, NDDoubleArray other) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        NDCounter otherCounter = other.getNewCounter();
        setAssignmentBroadcasting(ownCounter, otherCounter);
        assign(ownCounter, otherCounter, other);
    }

    public void assign(Slice[] rawSlices, NDDoubleArray other, Slice[] rawOtherSlices) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        SliceCounter otherCounter = other.getSliceCounter(rawOtherSlices);
        setAssignmentBroadcasting(ownCounter, otherCounter);
        assign(ownCounter, otherCounter, other);
    }
    public void plusAssign(NDDoubleArray other) {
        NDCounter ownCounter = getNewCounter();
        NDCounter otherCounter = other.getNewCounter();
        setAssignmentBroadcasting(ownCounter, otherCounter);
        plusAssign(ownCounter, otherCounter, other);
    }

    public void plusAssign(Slice[] rawSlices, NDDoubleArray other) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        NDCounter otherCounter = other.getNewCounter();
        setAssignmentBroadcasting(ownCounter, otherCounter);
        plusAssign(ownCounter, otherCounter, other);
    }

    public void plusAssign(Slice[] rawSlices, NDDoubleArray other, Slice[] rawOtherSlices) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        SliceCounter otherCounter = other.getSliceCounter(rawOtherSlices);
        setAssignmentBroadcasting(ownCounter, otherCounter);
        plusAssign(ownCounter, otherCounter, other);
    }
    

    public NDDoubleArray get(Slice[] rawSlices) {
        SliceCounter ownCounter = getSliceCounter(rawSlices);
        NDDoubleArray slicedArray = new NDDoubleArray(ownCounter.shape());
        AccessOperator operator = new AccessOperator(slicedArray);
        operate(operator, ownCounter, null, slicedArray.getNewCounter());
        return slicedArray;
    }

    public double getElement(int index) {
        return array[index];
    }

    public double getElement(int[] nDIndex) {
        return getElement(flattenIndex(nDIndex));
    }

    public double[] getArray() {
        return array;
    }

    public void setElement(int index, double value) {
        array[index] = value;
    }

    public void setElement(int[] nDIndex, double value) {
        setElement(flattenIndex(nDIndex), value);
    }

    public void setElement(int[] nDIndex, double value, IndexFlattener flattener) {
        setElement(flattener.flatten(nDIndex), value);
    }
}
