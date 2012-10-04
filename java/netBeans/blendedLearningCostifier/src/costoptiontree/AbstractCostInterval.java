/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

/**
 *
 * @author mcannamela
 */
public abstract class AbstractCostInterval implements Comparable {
    private int minCost;
    private int maxCost;
    private int selectedCost;
    private ScalingLaw scalingLaw;

    public AbstractCostInterval() {
        minCost = 0;
        maxCost = 100;
        selectedCost = 50;
    }
    public AbstractCostInterval(int minCost, int maxCost, int selectedCost, ScalingLaw scalingLaw) {
        this.minCost = minCost;
        this.maxCost = maxCost;
        this.selectedCost = selectedCost;
        this.scalingLaw = scalingLaw;
    }
    
    public int getScaledCost(){
        return scalingLaw.scale(selectedCost);
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
    public int getSelectedCost() {
        return selectedCost;
    }
    public ScalingLaw getScalingLaw() {
        return scalingLaw;
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
    public void setScalingLaw(ScalingLaw scalingLaw) {
        this.scalingLaw = scalingLaw;
    }
    
    
    /**
     *implements Comparable interface so that cost intervals can be compared.
     * they are compared first by mean, cost, then by minimum cost, then by 
     * maximum cost. 
     * @param other another cost interval to compare to this one
     * @return -1,0,1 depending upon whether other is >, ==, < this
     */
    public int compareTo(AbstractCostInterval other){
        int retVal=0;
        if (getMeanCost()>other.getMeanCost()){retVal= 1;}
        else if (getMeanCost()<other.getMeanCost()){retVal= -1;}
        else if (getMeanCost()==other.getMeanCost()){
            if (getMinCost()>other.getMinCost()){retVal= 1;}
            else if (getMinCost()<other.getMinCost()){ retVal= -1;}
            else if (getMinCost()==other.getMinCost()){
                if (getMaxCost()>other.getMaxCost()){retVal= 1;}
                else if (getMaxCost()<other.getMaxCost()){retVal= -1;}
                else if (getMaxCost()==other.getMaxCost()){retVal= 0;}
            }
        }
        return retVal;
    }
}
