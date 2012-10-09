/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import java.util.ArrayList;

/**
 *
 * @author mcannamela
 */
public class CostOptionNode {
    private CostOptionNode parent;
    private ArrayList<CostOptionNode> children;
    private ArrayList<CostOption> costOptions;
    private OptionSelectionInterface optionSelection;
    private String name;
    private String description;

    public CostOptionNode() {
        parent = null;
        children = new ArrayList<>();
        costOptions = new ArrayList<>();
        optionSelection = new SingleOptionSelection();
        name = "aCostOptionNode";
        description = "the default node is pretty boring";
        
    }

    public CostOptionNode(CostOptionNode parent, 
            ArrayList<CostOptionNode> children, 
            ArrayList<CostOption> costOptions, 
            OptionSelectionInterface optionSelection, 
            String name, String description) {
        this.parent = parent;
        this.children = children;
        this.costOptions = costOptions;
        this.optionSelection = optionSelection;
        this.name = name;
        this.description = description;
    }
    public CostOptionNode(CostOptionNode parent, 
            ArrayList<CostOption> costOptions, 
            OptionSelectionInterface optionSelection, 
            String name, String description) {
        this.parent = parent;
        this.children = new ArrayList<>();
        this.costOptions = costOptions;
        this.optionSelection = optionSelection;
        this.name = name;
        this.description = description;
    }

    
    
    
    
    public int nCostOptions(){
        return costOptions.size();
    }
    public int nChildren(){
        return children.size();
    }
    public boolean  isLeaf(){
        return nChildren()==0;
    }
    public void addChild(CostOptionNode childNode){
        children.add(childNode);
    }
    public CostOptionNode removeChild(){
        CostOptionNode child = children.remove(children.size()-1);
        return child;
    }
    public CostOptionNode removeChild(int index){
        CostOptionNode child = children.remove(index);
        return child;
    }
    
    public ArrayList<Double> getOptionBlendingFactors(){
        return optionSelection.getOptionBlendingFactors();
    }
    
    public ArrayList<Integer> selectedCosts(){
        ArrayList<Integer> costs = new ArrayList<>(nCostOptions());
        for (CostOption costOption : costOptions){
            costs.add(costOption.getScaledCost());
        }
        return costs;
    }
    public int ownCost(){
        ArrayList<Double> weights = getOptionBlendingFactors();
        ArrayList<Integer> costs = selectedCosts();
        assert weights.size()==costs.size() : "weights and costs must have same size";
        int cost = 0;
        for(int i=0; i<weights.size();i++){
            cost+= (int) ((double) weights.get(i))*((double)costs.get(i));
        }
        return cost;
    }
    public ArrayList<Integer> childCosts(){
        ArrayList<Integer> costs = new ArrayList<>(nChildren());
        for (CostOptionNode child : children){
            costs.add(child.cost());
        }
        return costs;
    }
    public int cost(){
        ArrayList<Integer> childCosts = childCosts();
        int cost = ownCost();
        for (Integer childCost : childCosts){
            cost+=childCost;
        }
        return cost;
    }

    public ArrayList<CostOptionNode> getChildren() {
        return children;
    }
    public ArrayList<CostOption> getCostOptions() {
        return costOptions;
    }
    public CostOptionNode getParent() {
        return parent;
    }
    public OptionSelectionInterface getOptionSelection() {
        return optionSelection;
    }

    public void setChildren(ArrayList<CostOptionNode> children) {
        this.children = children;
    }
    public void setCostOptions(ArrayList<CostOption> costOptions) {
        this.costOptions = costOptions;
    }
    public void setOptionSelection(OptionSelectionInterface optionSelection) {
        this.optionSelection = optionSelection;
    }
    public void setParent(CostOptionNode parent) {
        this.parent = parent;
    }
    
    
}
