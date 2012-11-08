/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

/**
 *
 * @author mcannamela
 */
public class NDEntity {

    protected int nElements;
    protected int[] shape;
    protected int[] strides;

    public NDEntity() {
        this.shape = new int[]{1};
        initNElements();
        initStrides();
    }

    public NDEntity(int[] shape) {
        this.shape = shape;
        initNElements();
        initStrides();
    }
    
    public static int[] idxCopy(int[] original) {
        int[] copy = new int[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }

    public static String idxToString(int[] idx) {
        if (idx.length == 0) {
            return "()";
        }
        if (idx.length == 1) {
            return "(" + idx[0] + ")";
        }

        String idxString = "(" + idx[0];
        for (int i = 1; i < idx.length; i++) {
            idxString += ", " + idx[i];
        }
        return idxString + ")";
    }

    public final int[] shape() {
        return idxCopy(shape);
    } 
    public final int nDimensions() {
        return shape.length;
    }
    public final int nElements() {
        return nElements;
    }
    public int nDTo1D(int[] nDIndex) {
        int idx = nDIndex[0];
        for (int i = 1; i < strides.length; i++) {
            idx += strides[i] * nDIndex[i];
        }
        return idx;
    }

    protected final void initNElements() {
        nElements = 1;
        for (int n : shape) {
            nElements *= n;
        }
    }
    protected final void initStrides() {
        strides = new int[nDimensions()];
        strides[0] = 1;
        strides[1] = shape[0];
        for (int i = 2; i < nDimensions(); i++) {
            strides[i] = strides[i - 1] * shape[i - 1];
        }
    }

    public final void setStrideAt(int index, int stride) {
        strides[index] = stride;
    }
    public final void setStrides(int[] strides) {
        this.strides = strides;
    }
    public final int getStrideAt(int index) {
        return strides[index];
    }
    public final int[] getStrides() {
        return strides;
    }

    int[] newIndex() {
        return new int[nDimensions()];
    }
}
