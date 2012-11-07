/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.Iterator;

/**
 *
 * @author mcannamela
 */
public class Slice{
    private int start, stop, step;
    private int[] indices;

    public Slice(int start, int stop, int step) {
        this.start = start;
        this.stop = stop;   
        this.step = step;
    }
    public Slice(int start, int stop) {
        this.start = start;
        this.stop = stop;
        this.step = 1;
    }
    public Slice(int start) {
        this.start = start;
        this.stop = -1;
        this.step = 1;
    }
    public Slice() {
        this.start = 0;
        this.stop = -1;
        this.step = 1;
    }
    
    public int nIndices(){
        if (stop<=start){
            return 0;
        }
        int w = stop-start;
        int n = w/step;
        if ((w-n*step)>0){
            n++;
        }
        return n;
    }
    
    private void setIndices(){
        indices = new int[nIndices()];
        indices[0] = start;
        for (int i=1; i<nIndices(); i++){
            indices[i] = indices[i-1]+1;
        }
    }
    public boolean hasFixedEndpoint(){
        return stop>0;
    }
    public Slice getFixedEndpointSlice(int stop){
        Slice slice = new Slice(start, stop, step);
        return slice;
    }
    public int start() {
        return start;
    }
    public int stop() {
        return stop;
    }
    public int step() {
        return step;
    }
    
    public int get(int index){
        return indices[index];
    }
    
}
