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
    private ArrayList<CostOptionNode> children = new ArrayList<>();
    private ArrayList<CostOption> costOptions= new ArrayList<>();
    private OptionSelectionInterface optionSelection;
    private String name;
    private String description;
    private int treeLevel;
    private boolean hasParent;

    public CostOptionNode() {
        parent = null;
        hasParent = false;
        children = new ArrayList<>();
        costOptions = new ArrayList<>();
        optionSelection = new SingleOptionSelection();
        name = "aCostOptionNode";
        description = "the default node is pretty boring";
        treeLevel = 0;
        
    }

    public CostOptionNode(CostOptionNode parent, 
            ArrayList<CostOptionNode> children, 
            ArrayList<CostOption> costOptions, 
            OptionSelectionInterface optionSelection, 
            String name, String description) {
        this.parent = parent;
        hasParent = true;
        this.children = children;
        this.costOptions = costOptions;
        this.optionSelection = optionSelection;
        this.name = name;
        this.description = description;
        setTreeLevel();
    }
    public CostOptionNode(CostOptionNode parent, 
            ArrayList<CostOption> costOptions, 
            OptionSelectionInterface optionSelection, 
            String name, String description) {
        this.parent = parent;
        hasParent = true;
        this.costOptions = costOptions;
        this.optionSelection = optionSelection;
        this.name = name;
        this.description = description;
        setTreeLevel();
    }
    
    public CostOptionNode( 
            ArrayList<CostOption> costOptions, 
            OptionSelectionInterface optionSelection, 
            String name, String description) {
        hasParent = false;
        this.children = new ArrayList<>();
        this.costOptions = costOptions;
        this.optionSelection = optionSelection;
        this.name = name;
        this.description = description;
        setTreeLevel();
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
    public boolean hasParent(){
        return hasParent;
    }
    public void addChild(CostOptionNode childNode){
        children.add(childNode);
        childNode.setParent(this);
    }
    public CostOptionNode removeChild(){
        CostOptionNode child = children.remove(children.size()-1);
        child.setParent(null);
        return child;
    }
    public CostOptionNode removeChild(int index){
        CostOptionNode child = children.remove(index);
        return child;
    }
    
    public ArrayList<Double> getOptionBlendingFactors(){
        System.out.println(optionSelection.nOptions() +"selection options");
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
    public CostOption getCostOption(int index){
        return costOptions.get(index);
    }
    public CostOptionNode getParent() {
        return parent;
    }
    public OptionSelectionInterface getOptionSelection() {
        return optionSelection;
    }
    public String getName() {
        return name;
    }
    public String getDescription() {
        return description;
    }
    public int getTreeLevel() {
        return treeLevel;
    }
    
    public void setChildren(ArrayList<CostOptionNode> children) {
        this.children = children;
        for(CostOptionNode child: this.children){
            child.setParent(this);
        }
    }
    public void setCostOptions(ArrayList<CostOption> costOptions) {
        this.costOptions = costOptions;
    }
    public void setOptionSelection(OptionSelectionInterface optionSelection) {
        this.optionSelection = optionSelection;
    }
    public void setParent(CostOptionNode parent) {
        this.parent = parent;
        if (parent==null){
            hasParent = false;
        }
        else{
            hasParent = true;
        }
        setTreeLevel();
    }
    public void setDescription(String description) {
        this.description = description;
    }
    public void setName(String name) {
        this.name = name;
    }    
    public final void setTreeLevel(){
        if (hasParent){
            treeLevel = parent.getTreeLevel()+1;
        }
        else{
            treeLevel = 0;
        }
    }

    @Override
    public String toString() {
        int nTabs = getTreeLevel();
        String tabString ="";
        for(int i = 0; i<getTreeLevel();i++){
            tabString+="    ";
        }
        String nodeString = "\n";
        nodeString+=tabString + getName()+" - "+getDescription();
        if (nCostOptions()>0){
            for (CostOption option : costOptions){
                nodeString+="\n"+tabString+option;
            }
        }
        if (nChildren()>0){
            for (CostOptionNode child: children){
                nodeString+= child;
            }
        }
        return nodeString;
        
    }
    
    
    
}
