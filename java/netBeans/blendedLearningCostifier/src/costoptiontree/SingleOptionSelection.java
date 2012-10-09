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
public class SingleOptionSelection implements OptionSelectionInterface{
    private int nOptions;
    private int selectedOptionIndex;

    public SingleOptionSelection() {
        nOptions = 0;
        selectedOptionIndex = 0;
    }

    public SingleOptionSelection(int nOptions, int selectedOptionIndex) {
        this.nOptions = nOptions;
        this.selectedOptionIndex = selectedOptionIndex;
    }

    public int getSelectedOptionIndex() {
        return selectedOptionIndex;
    }

    public int getnOptions() {
        return nOptions;
    }

    public void setnOptions(int nOptions) {
        this.nOptions = nOptions;
    }

    public void setSelectedOptionIndex(int selectedOptionIndex) {
        this.selectedOptionIndex = selectedOptionIndex;
    }
    
    

    @Override
    public ArrayList<Double> getOptionBlendingFactors() {
        ArrayList<Double> optionWeights = new ArrayList<>(nOptions);
        for (int i = 0; i<nOptions;i++){
            optionWeights.set(i, (double) 0);
            if (i==selectedOptionIndex){
                optionWeights.set(i, (double) 1);
            }
        }
        return optionWeights;
    }
    
    
}
