/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import java.util.ArrayList;
import java.util.Collections;

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
        hasParent = false;
        setParent(null);
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
        hasParent = false;
        setParent(parent);
        setChildren(children);
        
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
        hasParent = false;
        setParent(parent);
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
    
    public void addCostOption(CostOption option){
        costOptions.add(option);
        Collections.sort(costOptions);
    }
    public void removeCostOption(CostOption option){
        costOptions.remove(option);
    }
    public void addChild(CostOptionNode childNode){
        children.add(childNode);
        childNode.setParent(this);
    }
    public boolean isChild(CostOptionNode node){
        boolean isChild = false;
        for(CostOptionNode child: children){
            if (node==child){
                isChild = true;
            }
        }
        return isChild;
    }
    public CostOptionNode removeChild(){
        CostOptionNode child = children.remove(children.size()-1);
        child.setParent(null);
        return child;
    }
    public CostOptionNode removeChild(CostOptionNode node){
        children.remove(node);
        return node;
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
    
    public final void setChildren(ArrayList<CostOptionNode> children) {
        this.children = new ArrayList<>(children);
        for(CostOptionNode child: this.children){
            child.setParent(this);
        }
    }
    public void setCostOptions(ArrayList<CostOption> costOptions) {
        this.costOptions = new ArrayList<>(costOptions);
        Collections.sort(this.costOptions);
    }
    public void setOptionSelection(OptionSelectionInterface optionSelection) {
        this.optionSelection = optionSelection;
    }
    public final void setParent(CostOptionNode parent) {
        if (hasParent() && this.parent!=parent){
            this.parent.removeChild(this);
        }
        this.parent = parent;
        if (parent==null){
            hasParent = false;
        }
        else{
            hasParent = true;
            if (!parent.isChild(this)){
                parent.addChild(this);
            }
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
        String nodeString = "\n------node start from level "+getTreeLevel()+"------\n";
        nodeString+=tabString + getName()+" - "+getDescription();
        if (nCostOptions()>0){
            for (CostOption option : costOptions){
                nodeString+="\n  "+tabString+option;
            }
        }
        if (nChildren()>0){
            for (CostOptionNode child: children){
                nodeString+= child.toString();
            }
        }
        nodeString+="\n------node end at level "+getTreeLevel()+"------\n";
        return nodeString;
        
    }
    
    
    
}
