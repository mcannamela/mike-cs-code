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

    public int getStrideAt(int index) {
        return strides[index];
    }

    public int[] getStrides() {
        return strides;
    }

    protected static final int[] idxCopy(int[] original) {
        int[] copy = new int[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    

    protected void initNElements() {
        nElements = 1;
        for (int n : shape) {
            nElements *= n;
        }
    }

    protected void initStrides() {
        strides = new int[shape.length];
        strides[0] = 1;
        strides[1] = shape[0];
        for (int i = 2; i < shape.length; i++) {
            strides[i] = strides[i - 1] * shape[i - 1];
        }
    }

    public int nDimensions() {
        return shape.length;
    }
    
}
