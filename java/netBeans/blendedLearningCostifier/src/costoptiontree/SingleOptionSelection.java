/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import java.util.ArrayList;

/**
 *
 * @author wichtelwesen
 */
public class SingleOptionSelection implements OptionSelectionInterface<Integer>{
    private int nOptions;
    private int selection;

    public SingleOptionSelection() {
        nOptions = 1;
        selection = 0;
    }
    
    public SingleOptionSelection(SingleOptionSelection other) {
        nOptions = other.nOptions;
        selection = other.selection;
    }
    
    public SingleOptionSelection(int nOptions) {
        assert nOptions>0: "must have at least one option";
        this.nOptions = nOptions;
        this.selection = 0;
    }

    public SingleOptionSelection(int nOptions, int selection) {
        assert nOptions>0: "must have at least one option";
        this.nOptions = nOptions;
        this.selection = selection;
    }

    public int getSelection() {
        return selection;
    }

    @Override
    public int nOptions() {
        return nOptions;
    }

    public void setnOptions(int nOptions) {
        this.nOptions = nOptions;
    }

    @Override
    public void select(Integer selection) {
        this.selection = selection;
    }
    
    

    @Override
    public ArrayList<Double> getOptionBlendingFactors() {
        ArrayList<Double> optionWeights = new ArrayList<>(nOptions);
        
        for (int i = 0; i<nOptions;i++){
            optionWeights.add((double) 0);
            if (i==selection){
                optionWeights.set(i, (double) 1);
            }
        }
        return optionWeights;
    }
    
    
}
