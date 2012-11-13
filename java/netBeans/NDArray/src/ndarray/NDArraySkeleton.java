/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.NoSuchElementException;

/**
 *
 * @author Michael
 */
public abstract class NDArraySkeleton extends NDEntity {



    public NDArraySkeleton(int[] shape) {
        super(shape);
    }
    
    public static double[] arrayCopy(double[] original) {
        double[] copy = new double[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    public static int[] arrayCopy(int[] original) {
        int[] copy = new int[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    public static boolean[] arrayCopy(boolean[] original) {
        boolean[] copy = new boolean[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }

    public NDCounter getNewCounter() {
        return new NDCounter(shape());
    }

    public final SliceCounter getSliceCounter(Slice[] rawSlices) throws NoSuchElementException {
        assert rawSlices.length == shape.length : "must provide a slice for each dimension";
        Slice[] cookedSlices = new Slice[rawSlices.length];
        Slice cookedSlice;
        Slice slice;
        int dimension = 0;
        for (int i = 0; i < rawSlices.length; i++) {
            slice = rawSlices[i];
            if (slice.isRaw()) {
                cookedSlice = slice.getCookedSlice(shape[dimension]);
            } else {
                cookedSlice = slice;
            }
            cookedSlices[i] = cookedSlice;
        }
        for (Slice cooked : cookedSlices) {
            validateSlice(cooked, dimension);
        }
        return new SliceCounter(cookedSlices);
    }

    /*
     * set up loops to apply operator element by element to the array
     * //<editor-fold defaultstate="collapsed" desc=" private void operate(...) ">
     * thisCounter - the counter that iterates over this array
     * otherCounter - the counter that iterates over the other array in a
     *                  for a binary operation
     * resultCounter - the counter that iterates over the result array
     *                  into which the result of an operation that makes a new
     *                  array will be placed
     */
    protected void operate(elementOperator operator, NDCounter thisCounter, NDCounter otherCounter, NDCounter resultCounter) {
        assert nElements() == thisCounter.nElements() : "mismatched number of elements";
        assert thisCounter.nElements() == resultCounter.nElements() : "mismatched number of elements";
        assert thisCounter.nElements() == otherCounter.nElements() : "mismatched number of elements";
        int thisIndex;
        int otherIndex;
        int resultIndex;
        if (otherCounter == null) {
            otherIndex = -1;
            while (thisCounter.hasNext()) {
                thisIndex = thisCounter.nextFlat();
                resultIndex = resultCounter.nextFlat();
                operator.operate(thisIndex, otherIndex, resultIndex);
            }
        } else if (resultCounter == null) {
            resultIndex = -1;
            while (thisCounter.hasNext()) {
                thisIndex = thisCounter.nextFlat();
                otherIndex = otherCounter.nextFlat();
                operator.operate(thisIndex, otherIndex, resultIndex);
            }
        } else {
            while (thisCounter.hasNext()) {
                thisIndex = thisCounter.nextFlat();
                otherIndex = otherCounter.nextFlat();
                resultIndex = thisCounter.nextFlat();
                operator.operate(thisIndex, otherIndex, resultIndex);
            }
        }
    }
    //</editor-fold>

    protected void validateSlice(Slice slice, int dimension) throws NoSuchElementException {
        if (slice.start() < 0 || slice.start() > shape[dimension]) {
            throw new NoSuchElementException("slice start is " + slice.start() + "  in dimension " + dimension + " where the array has size " + shape[dimension]);
        }
        if (slice.stop() < 0 || slice.stop() > shape[dimension]) {
            throw new NoSuchElementException("slice stop is " + slice.stop() + "  in dimension " + dimension + " where the array has size " + shape[dimension]);
        }
        if (slice.nIndices() == 0) {
            throw new NoSuchElementException("slice has no elements in dimension " + dimension + ": zero size arrays unsupported...for now");
        }
    }
    
}
