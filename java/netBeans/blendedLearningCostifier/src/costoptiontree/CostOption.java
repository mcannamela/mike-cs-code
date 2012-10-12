/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.AbstractScalingLaw;

/**
 *
 * @author mcannamela
 */
public class CostOption implements Comparable<CostOption> {
    private String label;
    private String description;
    
    private int minCost;
    private int maxCost;
    private int selectedCost;
    
    private AbstractScalingLaw scalingLaw;

    public CostOption() {
        minCost = 0;
        maxCost = 100;
        selectedCost = 50;
        
    }
    
    public CostOption(CostOption other){
        this.label = other.label;
        this.description = other.description;
        this.minCost = other.minCost;
        this.maxCost = other.maxCost;
        this.selectedCost = other.selectedCost;
        this.scalingLaw = other.scalingLaw;
    }
    

    public CostOption(String label, String description, int minCost, int maxCost,
                        int selectedCost, AbstractScalingLaw scalingLaw) {
        this.label = label;
        this.description = description;
        assert minCost<maxCost:"minCost "+minCost+" must be less than maxCost "+maxCost;
        this.minCost = minCost;
        this.maxCost = maxCost;
        if (selectedCost<minCost || selectedCost>maxCost){
            this.selectedCost = (minCost+maxCost)/2;
        }
        else{
            this.selectedCost = selectedCost;
        }
        this.scalingLaw = scalingLaw;
    }
    
    
    public int getScaledCost(){
        return scalingLaw.scale(selectedCost);
    }

    public String getLabel() {
        return label;
    }
    public String getDescription() {
        return description;
    }
    
    public int getMeanCost(){
        return (getMinCost()+getMaxCost())/2;
    }
    
    public int getMaxCost() {
        return maxCost;
    }
    public int getMinCost() {
        return minCost;
    }
    public int getScaledMeanCost(){
        return scalingLaw.scale(getMeanCost());
    }
    public int getScaledMinCost(){
        return scalingLaw.scale(getMinCost());
    }
    public int getScaledMaxCost(){
        return scalingLaw.scale(getMaxCost());
    }
    public int getSelectedCost() {
        return selectedCost;
    }
    public AbstractScalingLaw getScalingLaw() {
        return scalingLaw;
    }

    public void setLabel(String label) {
        this.label = label;
    }
    public void setDescription(String description) {
        this.description = description;
    }
    public void setMinCost(int minCost) {
        this.minCost = minCost;
    }
    public void setMaxCost(int maxCost) {
        this.maxCost = maxCost;
    }
    public void setSelectedCost(int selectedCost) {
        this.selectedCost = selectedCost;
    }
    public void setScalingLaw(AbstractScalingLaw scalingLaw) {
        this.scalingLaw = scalingLaw;
    }
    
    
    /**
     *implements Comparable interface so that CostIntervals can be compared.
     * they are compared first by mean, cost, then by minimum cost, then by 
     * maximum cost. 
     * @param other another cost interval to compare to this one
     * @return -1,0,1 depending upon whether other is >, ==, < this
     */
    @Override
    public int compareTo(CostOption other){
               
        int retVal=0;
        
        if (getScaledMeanCost()>other.getScaledMeanCost()){retVal= 1;}
        else if (getScaledMeanCost()<other.getScaledMeanCost()){retVal= -1;}
        else if (getScaledMeanCost()==other.getScaledMeanCost()){
            if (getScaledMinCost()>other.getScaledMinCost()){retVal= 1;}
            else if (getScaledMinCost()<other.getScaledMinCost()){ retVal= -1;}
            else if (getScaledMinCost()==other.getScaledMinCost()){
                if (getScaledMaxCost()>other.getScaledMaxCost()){retVal= 1;}
                else if (getScaledMaxCost()<other.getScaledMaxCost()){retVal= -1;}
                else if (getScaledMaxCost()==other.getScaledMaxCost()){retVal= 0;}
            }
        }
        return retVal;
    }

    @Override
    public String toString() {
        return getLabel()+", ("+getMinCost()+", "+getMaxCost()+", "+getSelectedCost()
                +") - "+getDescription();
    }
    
}
