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
    
    public static int[] idxCopy(int[] original) {
        int[] copy = new int[original.length];
        System.arraycopy(original, 0, copy, 0, original.length);
        return copy;
    }

    protected int[] makeBroadcastMask(int[] shape) {
        throw new UnsupportedOperationException("Not yet implemented");
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
            int cnt = 0;
            flattenerStrides = new int[mask.length];
            for (int i = 0;i<mask.length;i++){
                if (mask[i]==0){
                    flattenerStrides[i]=0;
                    if (shape[i]==1){
                        cnt++;
                    }
                }
                else{
                    flattenerStrides[i] = strides[cnt];
                    cnt++;
                }
            }
            assert (cnt)==(nDimensions()):
                        String.format("expected %d nonzero entries in broadcastMask, got only %d", 
                        (nDimensions()-nSingletonDimensions()), cnt);
            
        }
        
    }
    //</editor-fold>
    
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
        return idxCopy(shape);
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
    public final int nElements() {
        return nElements;
    }
    
    public int flattenIndex(int[] nDIndex) {
        return flattener.flatten(nDIndex);
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
    
    public final  int[] getFlattenerStrides(){
        return flattener.getFlattenerStrides();
    }
    
}
