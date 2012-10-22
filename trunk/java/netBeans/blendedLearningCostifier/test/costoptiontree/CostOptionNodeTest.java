/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.AbstractScalingLaw;

import blendedlearningprogram.BaseBlendedLearningModel;
import blendedlearningprogram.ConstantScalingLaw;
import blendedlearningprogram.LinearInStudentsScalingLaw;
import blendedlearningprogram.LinearInTeachersScalingLaw;
import blendedlearningprogram.ProgramSize;
import blendedlearningprogram.StandardBlendedLearningModel;
import java.util.ArrayList;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Michael
 */
public class CostOptionNodeTest {
    private static ArrayList<CostOption> costOptions = new ArrayList<>(3);
    private static ProgramSize programSize;
    private static BaseBlendedLearningModel blendedLearningModel;
    private static int nrStudents;
    private static int nrPeriods;
    
    private static ArrayList<String> labels = new ArrayList<>(3);
    private static ArrayList<String> descriptions = new ArrayList<>(3);
    private static int[] minCosts = {10, 20, 30};
    private static int[] maxCosts = {100, 200, 300};
    private static int[] selectedCosts = {50, 90, 150};
    private static ArrayList<AbstractScalingLaw> scalingLaws = new ArrayList<>(3);
    
    private CostOptionNode parent;
    
    private OptionSelectionInterface optionSelection;
    private CostOptionNode instance;
    private ArrayList<CostOptionNode> children = new ArrayList<>(3);
    private String name;
    private String optionNodeDescription;
    
    public CostOptionNodeTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        blendedLearningModel = new StandardBlendedLearningModel();
        nrStudents = 300;
        nrPeriods = 6;
        programSize = new ProgramSize(blendedLearningModel, nrStudents, nrPeriods);
        
        for (int i = 0; i<3;i++){
            labels.add("option"+i);
            descriptions.add("description"+i);
        }
        scalingLaws.add(new ConstantScalingLaw(programSize));
        scalingLaws.add(new LinearInStudentsScalingLaw(programSize));
        scalingLaws.add(new LinearInTeachersScalingLaw(programSize));
        
        for (int i=0;i<3;i++){
            costOptions.add( new CostOption( labels.get(i),  descriptions.get(i), 
                    minCosts[i],  maxCosts[i], selectedCosts[i],  
                    scalingLaws.get(i)));
        }
        
        
        
        
    }
    
    @Before
    public void setUp() {
        
        parent = new CostOptionNode();
        
        optionSelection = new SingleOptionSelection(costOptions.size(),0);
//        System.out.println("there are "+costOptions.size() +" options");
        
        SingleOptionSelection copyOptionSelection = new SingleOptionSelection((SingleOptionSelection) optionSelection);
//        System.out.println("there are "+copyOptionSelection.nOptions() +" selection options");
        name = "nodeName";
        optionNodeDescription = "a cost option node";
        instance = new CostOptionNode( parent, 
            costOptions, 
            optionSelection, 
            name, optionNodeDescription);
        
        children = new ArrayList<>(3);
        for (int i=0;i<3;i++){
            children.add(new CostOptionNode(copyCostOptions(costOptions), 
                                new SingleOptionSelection((SingleOptionSelection) optionSelection), 
                                "childNode"+i, 
                                "description of childNode"+i)
                        );
        }
            
    }
    
    private ArrayList<CostOption> copyCostOptions(ArrayList<CostOption> options){
        ArrayList<CostOption> copyOfOptions = new ArrayList<>(options.size()); 
        for (CostOption x:options){
            copyOfOptions.add(new CostOption(x));
        }
        return copyOfOptions;
    }
    
    @Test
    public void testArrayListEquals(){
        ArrayList<Double> a = new ArrayList<>(2);
        ArrayList<Double> b = new ArrayList<>(2);
        ArrayList<Double> c = new ArrayList<>();
        
        for(int i = 0;i<2;i++){
            a.add((double)i);
            b.add((double)i);
            c.add((double)i);
        }
        assertEquals(a,b);
        assertEquals(a,c);
        assertEquals(b,c);
    }

    @Test
    public void testNCostOptions() {
        System.out.println("nCostOptions");
        
        int expResult = 3;
        int result = instance.nCostOptions();
        assertEquals(expResult, result);
        
        
    }

    @Test
    public void testNChildren() {
        System.out.println("nChildren");
        
        assertEquals(0,instance.nChildren());
        
        int childCnt = 0;
        for(CostOptionNode child:children){
            instance.addChild(child);
            childCnt++;
            assertEquals(childCnt,instance.nChildren());
        }
        for(CostOptionNode child:children){
            instance.removeChild();
            childCnt--;
            assertEquals(childCnt,instance.nChildren());
        }
        
    }

    @Test
    public void testIsLeaf() {
        System.out.println("isLeaf");
        
        assertEquals(true, instance.isLeaf());
        int childCnt = 0;
        for(CostOptionNode child:children){
            instance.addChild(child);
            childCnt++;
            assertEquals(false, instance.isLeaf());
            assertEquals(childCnt,instance.nChildren());
        }
        
    }

    @Test
    public void testAddChild() {
        System.out.println("addChild");

        instance.addChild(children.get(0));
        assertEquals(1,instance.nChildren());
        assertEquals(children.get(0),instance.getChildren().get(0));
        assertTrue(instance.isChild(children.get(0)));

    }
    
    @Test
    
    public void testIsChild(){
        System.out.println("isChild");
        
        for(CostOptionNode node: children){
            assertFalse(instance.isChild(node));
        }
        
        instance.setChildren(children);
        for(CostOptionNode node: children){
            assertTrue(instance.isChild(node));
        }
        for(CostOptionNode node: children){
            instance.removeChild(node);
        }
        
        for(CostOptionNode node: children){
            assertFalse(instance.isChild(node));
        }
        
        for(CostOptionNode node: children){
            instance.addChild(node);
        }
        
        for(CostOptionNode node: children){
            assertTrue(instance.isChild(node));
        }
        
    }

    @Test
    public void testRemoveChild_0args() {
        System.out.println("removeChild 0args");
        instance.addChild(children.get(0));
        CostOptionNode child = instance.removeChild();
        assertEquals(children.get(0), child);
    }
    
    @Test
    public void testRemoveChild_CostOptionNode() {
        System.out.println("removeChild node");
        for (CostOptionNode node: children){
            instance.addChild(node);
        }
        for (CostOptionNode node: children){
            instance.removeChild(node);
        }
        
        instance.addChild(children.get(0));
        instance.addChild(children.get(1));
        instance.removeChild(children.get(1));
        assertEquals(1, instance.nChildren());
        assertEquals(children.get(0), instance.getChildren().get(0));
        
    }

    @Test
    public void testRemoveChild_int() {
        System.out.println("removeChild int");
        
        int childCnt = 0;
        for(CostOptionNode child:children){
            instance.addChild(child);
            childCnt++;
            assertEquals(childCnt,instance.nChildren());
        }
        
        CostOptionNode child = instance.removeChild(1);
        assertEquals(children.get(1), child);
        
    }

    @Test
    public void testGetOptionBlendingFactors() {
        System.out.println("getOptionBlendingFactors");
        
        ArrayList<Double> expResult = new ArrayList<>();
        double[] expected = {1.0, 0.0, 0.0};
        for (int i=0;i<3;i++){
            expResult.add(expected[i]);
        }
        ArrayList result = instance.getOptionBlendingFactors();
        assertEquals(expResult, result);
        
        
    }

    @Test
    public void testSelectedCosts() {
        System.out.println("selectedCosts");
        
        int[] expResult = {50, 27000, 300};
        ArrayList result = instance.selectedCosts();
        for (int i=0;i<3;i++){
            assertEquals(expResult[i], result.get(i));
        }
        
        
        
        
    }

    @Test
    
    public void testOwnCost() {
        System.out.println("ownCost");
        
        assertEquals(50, instance.ownCost());
        
        ((SingleOptionSelection)optionSelection).select(1);
        assertEquals(27000, instance.ownCost());
        
        ((SingleOptionSelection)optionSelection).select(2);
        assertEquals(300, instance.ownCost());
    }

    @Test
    
    public void testChildCosts() {
        System.out.println("childCosts");
        
        assertEquals(0, instance.childCosts().size());
        
        int childCnt = 0;
        for(CostOptionNode child:children){
            instance.addChild(child);
            childCnt++;
            assertEquals(childCnt,instance.nChildren());
        }
        
        int[] expResult = {50, 50,50};
        ArrayList<Integer> result = instance.childCosts();
        for (int i =0;i<3;i++){
            assertEquals((int)expResult[i], (int)result.get(i));
        }
        ArrayList<CostOptionNode> theChildren = instance.getChildren();
        theChildren.get(1).setOptionSelection(new SingleOptionSelection(3,1));
        theChildren.get(2).setOptionSelection(new SingleOptionSelection(3,2));
        
        expResult[1] = 27000;
        expResult[2] = 300;
        
        result = instance.childCosts();
        for (int i =0;i<3;i++){
            assertEquals((int)expResult[i], (int)result.get(i));
        }
        

    }

    @Test
    public void testCost() {
        System.out.println("cost");
        
        assertEquals(50, instance.cost());
        
        instance.setChildren(children);
        assertEquals(200, instance.cost());
    }

    @Test
    @Ignore
    public void testGetChildren() {
        System.out.println("getChildren");
        
        ArrayList expResult = null;
        ArrayList result = instance.getChildren();
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetCostOptions() {
        System.out.println("getCostOptions");
        
        ArrayList expResult = null;
        ArrayList result = instance.getCostOptions();
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetCostOption() {
        System.out.println("getCostOption");
        int index = 0;
        
        CostOption expResult = null;
        CostOption result = instance.getCostOption(index);
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetParent() {
        System.out.println("getParent");
        
        CostOptionNode expResult = null;
        CostOptionNode result = instance.getParent();
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetOptionSelection() {
        System.out.println("getOptionSelection");
        
        OptionSelectionInterface expResult = null;
        OptionSelectionInterface result = instance.getOptionSelection();
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetName() {
        System.out.println("getName");
        
        String expResult = "";
        String result = instance.getName();
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetDescription() {
        System.out.println("getDescription");
        
        String expResult = "";
        String result = instance.getDescription();
        assertEquals(expResult, result);
        
        
    }

    @Test
    @Ignore
    public void testGetTreeLevel() {
        System.out.println("getTreeLevel");
        
        int expResult = 0;
        int result = instance.getTreeLevel();
        assertEquals(expResult, result);
        
        
    }

    @Test
    public void testSetChildren() {
        System.out.println("setChildren");
        
        instance.setChildren(children);
        assertEquals(3,instance.nChildren());
        assertTrue(instance.getChildren().get(0).hasParent());
        assertTrue(instance.getChildren().get(1).hasParent());
        assertTrue(instance.getChildren().get(2).hasParent());
        assertEquals(2,instance.getChildren().get(0).getTreeLevel());
        assertEquals(2,instance.getChildren().get(1).getTreeLevel());
        assertEquals(2,instance.getChildren().get(2).getTreeLevel());
        
        
        
    }

    @Test
    @Ignore
    public void testSetCostOptions() {
        System.out.println("setCostOptions");
        ArrayList<CostOption> costOptions = null;
        
        instance.setCostOptions(costOptions);
        
        
    }

    @Test
    @Ignore
    public void testSetOptionSelection() {
        System.out.println("setOptionSelection");
        OptionSelectionInterface optionSelection = null;
        
        instance.setOptionSelection(optionSelection);
        
        
    }

    @Test
    public void testSetParent() {
        System.out.println("setParent");
        
        assertTrue(parent.isChild(instance));
        instance.setParent(null);
        assertFalse(parent.isChild(instance));
        assertFalse(instance.hasParent());
        assertEquals(0,instance.getTreeLevel());
        
        instance.setParent(parent);
        assertTrue(instance.hasParent());
        assertTrue(parent.isChild(instance));
        assertEquals(1,instance.getTreeLevel());
        
        CostOptionNode otherParent = new CostOptionNode();
        assertFalse(otherParent==parent);
        instance.setParent(otherParent);
        assertTrue(instance.hasParent());
        assertFalse(parent.isChild(instance));
        assertTrue(otherParent.isChild(instance));
        
    }

    @Test
    @Ignore
    public void testSetDescription() {
        System.out.println("setDescription");
        String description = "";
        
        instance.setDescription(description);

    }

    @Test
    @Ignore
    public void testSetName() {
        System.out.println("setName");
        String name = "";
        
        instance.setName(name);
        
        
    }

    @Test
    @Ignore
    public void testSetTreeLevel() {
        System.out.println("setTreeLevel");
        
        instance.setTreeLevel();
        
        
    }

    @Test
    public void testToString() {
        System.out.println("toString");
        System.out.println("from parent:");
        System.out.println(parent);
        
        assertEquals(0,instance.nChildren());
        instance.setChildren(children);
        
        System.out.println("with children:");
        System.out.println(parent);
        
        
        
        
        
    }
}
