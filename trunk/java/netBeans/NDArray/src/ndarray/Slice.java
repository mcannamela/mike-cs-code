/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

/**
 *
 * @author mcannamela
 */
public class Slice{
    private int start, stop, step;
    private int[] indices;
    private static int UNCOOKED_STOP_VALUE = -1;
    private static int DEFAULT_START_VALUE = 0;
    private static int DEFAULT_STEP_VALUE = 1;

    public Slice() {
        this.start = DEFAULT_START_VALUE;
        this.stop = UNCOOKED_STOP_VALUE;
        this.step = DEFAULT_STEP_VALUE;
        validate();
    }
    public Slice(int start) {
        this.start = start;
        this.stop = UNCOOKED_STOP_VALUE;
        this.step = DEFAULT_STEP_VALUE;
        validate();
    }
    public Slice(int start, int stop) {
        this.start = start;
        this.stop = stop;
        this.step = DEFAULT_STEP_VALUE;
        validate();
    }
    public Slice(int start, int stop, int step) {
        this.start = start;
        this.stop = stop;   
        this.step = step;
        validate();
    }
    
    private void validate(){
        validateStart();
        validateStop();
        validateStep();
    }
    
    private void validateStart(){
        assert start>=0:"start must be positive";
    }
    private void validateStop(){
        assert (stop>=0||stop==UNCOOKED_STOP_VALUE):"stop must be positive, or -1 for end element";
    }
    private void validateStep(){
        if (start<stop){
            assert step>=0:"this is an incrementing slice, step must be positive";
        }
        else if (stop<start){
            assert step<0:"this is a decrementing slice, step must be negative";
        }
        else{
            assert step==0:"step is neither 0, positive, or negative. something's rotten";
            assert step!=0:"step must be nonzero";
        }
    }
    
    private void setIndices(){
        indices = new int[nIndices()];
        for (int i=0; i<nIndices(); i++){
            indices[i] = start+step*i;
        }
    }
    
    public final boolean isRaw() {
        return stop==UNCOOKED_STOP_VALUE;
    }
    
    public final int nIndices(){
        int w = Math.abs(stop-start);
        int n = Math.abs(w/step);
        if ((w-n*Math.abs(step))>0){
            n++;
        }
        return n;
    }
    
    public final Slice getCookedSlice(int stop){
        Slice slice = new Slice(start, stop, step);
        return slice;
    }
    public final int start() {
        return start;
    }
    public final int stop() {
        return stop;
    }
    public final int step() {
        return step;
    }
    
    public final int get(int index){
        return indices[index];
    }

    
    
}
