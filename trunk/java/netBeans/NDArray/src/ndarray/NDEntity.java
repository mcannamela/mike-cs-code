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
public class NDEntity {

    protected int nElements;
    protected int[] shape;
    protected int[] strides;
    
    protected IndexFlattener flattener;

    public NDEntity() {
        this.shape = new int[]{1};
        initNElements();
        initStrides();
        flattener = getDefaultIndexFlattener();
    }

    public NDEntity(int[] shape) {
        this.shape = shape;
        initNElements();
        initStrides();
        flattener = getDefaultIndexFlattener();
    }
    
    public static int[] arrayCopy(int[] original) {
        int[] copy = new int[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    public static double[] arrayCopy(double[] original) {
        double[] copy = new double[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    public static boolean[] arrayCopy(boolean[] original) {
        boolean[] copy = new boolean[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }
    static int[] negotiateTargetShape(int[] firstShape, int[] secondShape){
        int[] bigShape, smallShape;
        if (firstShape.length>=secondShape.length){
            bigShape = firstShape;
            smallShape = secondShape;
        }
        else{
            bigShape = secondShape;
            smallShape = firstShape;
        }
        int[] targetShape = new int[bigShape.length];
        
        int offset = bigShape.length - smallShape.length;
        
        System.arraycopy(bigShape, 0, targetShape, 0, offset);
        for (int i = offset; i<bigShape.length;i++){
            if (bigShape[i]==1){
                targetShape[i] = smallShape[i-offset];
            }
            else if (smallShape[i-offset]==1){
                targetShape[i] = bigShape[i];
            }
            else{
                assert bigShape[i]==smallShape[i-offset]:"target shape negotiation failed";
                targetShape[i] = bigShape[i];
            }
        }
        return targetShape;
    }
    static int[] makeBroadcastMask(int[] targetShape, int[] broadcastingShape) {
        int[] mask = new int[targetShape.length];
               
        if (targetShape.length<broadcastingShape.length){
            throw new UnsupportedOperationException(
                "broadcasting to arrays of lesser dimension not yet supported");
        }
        else if (targetShape.length==broadcastingShape.length){
            int offset = (targetShape.length-broadcastingShape.length);
            
            //prepend with broadcasting dimensions
            for (int i=0;i<offset;i++){
                mask[i] = 0;
            }
            
            //singleton dimensions must also be broadcast
            for (int i=offset;i<targetShape.length;i++){
                if (broadcastingShape[i-offset]==1){
                    mask[i] = 0;
                }
                else{
                    assertBroadcastDimensionEqual(i, targetShape[i],broadcastingShape[i-offset]);
                    mask[i]=1;
                }
            }
        }
        return mask;
    }
    private static void assertBroadcastDimensionEqual(int idx, int thisShape, int thatShape){
        assert thisShape==thatShape:
            String.format("broadcast mismatch in dimension "+
                          "%d: passed size was %d but this size is %d",
                          idx, thatShape, thisShape);
    }
    
    private class DefaultIndexFlattener extends IndexFlattener{

        public DefaultIndexFlattener() {
            flattenerStrides = strides;
        }
    }
    
    /*
     * //<editor-fold defaultstate="collapsed" desc=" class BroadcastingFlattener extends IndexFlattener ">
     */
    class BroadcastingFlattener extends IndexFlattener{
        
        public BroadcastingFlattener(){
            int[] mask = new int[shape.length];
            for (int i = 0;i<nDimensions();i++){
                if (shape[i]==1){
                    mask[i] = 0;
                }
                else{
                    mask[i] = 1;
                }
            }
            setBroadcastMask(mask);
        }
    
        public BroadcastingFlattener(int[] broadcastMask){
            setBroadcastMask(broadcastMask);
        }
        
        private void setBroadcastMask(int[] mask){
            assert mask.length>shape.length:"broadcasting to array of lower"+
                                            "dimension unsupported";
            validateMask(mask);
            
            int cnt = 0;
            flattenerStrides = new int[mask.length];
            for (int i = 0;i<mask.length;i++){
                if (mask[i]==0){
                    flattenerStrides[i]=0;
                }
                else{
                    while (cnt<shape.length && shape[cnt]==1){
                        cnt++;
                    }
                    flattenerStrides[i] = strides[cnt];
                    cnt++;
                }
            }
        }
    }
    //</editor-fold>
    private void validateMask(int[] mask){
        int nMask = 0;
        int nShape=0;
        for (int m:mask){
            if (m>0){
                nMask++;
            }
        }
        for (int s:shape){
            if (s>1){
                nShape++;
            }
        }
        assert nMask==nShape: 
            String.format("broadcast dimension mismatch: "+
                "broadcast mask expects %d non-singleton dimensions but "
                +"shape has %d", nMask, nShape);            
        
    }
    protected final void initNElements() {
        nElements = 1;
        for (int n : shape) {
            nElements *= n;
        }
    }
    protected final void initStrides() {
        assert nDimensions()>=1:"0-D arrays currently unsupported";
        strides = new int[nDimensions()];
        strides[0] = 1;
        if (nDimensions()==1){
            return;
        }
        strides[1] = shape[0];
        for (int i = 2; i < nDimensions(); i++) {
            strides[i] = strides[i - 1] * shape[i - 1];
        }
    }
    
    int[] newIndex() {
        return new int[nDimensions()];
    }

    public final int[] shape() {
        return arrayCopy(shape);
    } 
    public final int nDimensions() {
        return shape.length;
    }
    public final int nSingletonDimensions() {
        int cnt = 0;
        for(int n:shape){
            if (n==1){
                cnt++;
            }
        }
        return cnt;
    }
    public final int nNonSingletonDimensions(){
        return nDimensions()-nSingletonDimensions();
    }
    public final int nElements() {
        return nElements;
    }
    
    public int flattenIndex(int[] nDIndex) {
        return flattener.flatten(nDIndex);
    }
    public void reshape(int[] shape){
        int p = 1;
        for (int s: shape){
            p*=s;
        }
        assert p==nElements():"new shape must have same product as old shape";
        
        this.shape = arrayCopy(shape);
        initNElements();
        initStrides();
        flattener = getDefaultIndexFlattener();
        
    }
    public void squeeze(){
        int[] squeezedShape = new int[nNonSingletonDimensions()];
        int cnt = 0;
        for(int i = 0; i<nDimensions();i++){
            if (shape[i]!=1){
                squeezedShape[cnt] = shape[i];
                cnt++;
            }
        }
        reshape(squeezedShape);
    }
    
    public final IndexFlattener getDefaultIndexFlattener(){
        return new DefaultIndexFlattener();
    }
    BroadcastingFlattener getBroadcastFlattener(int[] broadcastMask){
        return new BroadcastingFlattener(broadcastMask);
    }
    
    public void setBroadcasting(int[] broadcastMask){
        flattener = getBroadcastFlattener(broadcastMask);
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
    
    public final int[] getFlattenerStrides(){
        return flattener.getFlattenerStrides();
    }
    
}
