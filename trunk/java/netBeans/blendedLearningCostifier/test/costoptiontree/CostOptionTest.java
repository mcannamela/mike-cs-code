/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.AbstractScalingLaw;
import blendedlearningprogram.BaseBlendedLearningModel;
import blendedlearningprogram.LinearInStudentsScalingLaw;
import blendedlearningprogram.LinearInTeachersScalingLaw;
import blendedlearningprogram.ProgramSize;
import blendedlearningprogram.StandardBlendedLearningModel;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Michael
 */
public class CostOptionTest {
    private CostOption instance;
    private static ProgramSize programSize;
    private static BaseBlendedLearningModel blendedLearningModel;
    private static int nrStudents;
    private static int nrPeriods;
    
    private String label;
    private String description;
    private int minCost;
    private int maxCost;
    private int selectedCost;
    private AbstractScalingLaw scalingLaw;
    
    
    public CostOptionTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        blendedLearningModel = new StandardBlendedLearningModel();
        nrStudents = 300;
        nrPeriods = 6;
        programSize = new ProgramSize(blendedLearningModel, nrStudents, nrPeriods);
    }
    
    @Before
    public void setUp() {
        label = "the label";
        description = "this is the description";
        minCost = 10;
        maxCost = 100;
        selectedCost = 50;
        scalingLaw = new LinearInStudentsScalingLaw(programSize);
        
        
       
        instance = new CostOption( label,  description,  minCost,  maxCost,
                         selectedCost,  scalingLaw);
    }
    
    

    /**
     * Test of getScaledCost method, of class CostOption.
     */
    @Test
    public void testGetScaledCost() {
        System.out.println("getScaledCost");
        
        int expResult = 300*50;
        int result = instance.getScaledCost();
        assertEquals(expResult, result);
        
    }

    /**
     * Test of getLabel method, of class CostOption.
     */
    @Test
    public void testGetLabel() {
        System.out.println("getLabel");
        
        String expResult = label;
        String result = instance.getLabel();
        assertEquals(expResult, result);
        
    }

    /**
     * Test of getDescription method, of class CostOption.
     */
    @Test
    public void testGetDescription() {
        System.out.println("getDescription");
        
        String expResult = description;
        String result = instance.getDescription();
        assertEquals(expResult, result);
        
    }

    /**
     * Test of getMeanCost method, of class CostOption.
     */
    @Test
    public void testGetMeanCost() {
        System.out.println("getMeanCost");
        
        int expResult = 55;
        assertEquals(expResult, instance.getMeanCost());
        
        //check rounding behavior
        instance.setMinCost(11);
        assertEquals(11, instance.getMinCost());
        assertEquals(expResult, instance.getMeanCost());
    }

    /**
     * Test of getMaxCost method, of class CostOption.
     */
    @Test
    public void testGetMaxCost() {
        System.out.println("getMaxCost");
        
        int expResult = maxCost;
        int result = instance.getMaxCost();
        assertEquals(expResult, result);
    }

    /**
     * Test of getMinCost method, of class CostOption.
     */
    @Test
    public void testGetMinCost() {
        System.out.println("getMinCost");
        
        int expResult = minCost;
        int result = instance.getMinCost();
        assertEquals(expResult, result);
    }

    /**
     * Test of getSelectedCost method, of class CostOption.
     */
    @Test
    public void testGetSelectedCost() {
        System.out.println("getSelectedCost");
        
        int expResult = selectedCost;
        int result = instance.getSelectedCost();
        assertEquals(expResult, result);
    }

    /**
     * Test of getScalingLaw method, of class CostOption.
     */
    @Test
    public void testGetScalingLaw() {
        System.out.println("getScalingLaw");
        
        AbstractScalingLaw expResult = scalingLaw;
        AbstractScalingLaw result = instance.getScalingLaw();
        assertEquals(expResult, result);
    }

    /**
     * Test of setLabel method, of class CostOption.
     */
    @Test
    public void testSetLabel() {
        System.out.println("setLabel");
        label = "newLabel";
        
        instance.setLabel(label);
        assertEquals(label, instance.getLabel());
    }

    /**
     * Test of setDescription method, of class CostOption.
     */
    @Test
    public void testSetDescription() {
        System.out.println("setDescription");
        description = "newDescription";
        
        instance.setDescription(description);
        assertEquals(description, instance.getDescription());
    }

    /**
     * Test of setMinCost method, of class CostOption.
     */
    @Test
    public void testSetMinCost() {
        System.out.println("setMinCost");
        minCost = 0;
        
        instance.setMinCost(minCost);
        assertEquals(minCost, instance.getMinCost());
    }

    /**
     * Test of setMaxCost method, of class CostOption.
     */
    @Test
    public void testSetMaxCost() {
        System.out.println("setMaxCost");
        maxCost = 0;
        
        instance.setMaxCost(maxCost);
        assertEquals(maxCost, instance.getMaxCost());
    }

    /**
     * Test of setSelectedCost method, of class CostOption.
     */
    @Test
    public void testSetSelectedCost() {
        System.out.println("setSelectedCost");
        selectedCost = 20;
        
        instance.setSelectedCost(selectedCost);
        assertEquals(selectedCost, instance.getSelectedCost());
    }

    /**
     * Test of setScalingLaw method, of class CostOption.
     */
    @Test
    public void testSetScalingLaw() {
        System.out.println("setScalingLaw");
        LinearInTeachersScalingLaw linearInTeachersScalingLaw = new LinearInTeachersScalingLaw(programSize);
        
        instance.setScalingLaw(linearInTeachersScalingLaw);
        assertEquals(linearInTeachersScalingLaw, instance.getScalingLaw());
    }

    /**
     * Test of compareTo method, of class CostOption.
     */
    @Test
    public void testCompareTo() {
        System.out.println("compareTo");
        CostOption other = new CostOption( label,  description,  minCost,  maxCost,
                         selectedCost,  scalingLaw);
        
        //min, max are equal, should be equal
        assertEquals(0, instance.compareTo(other));
        
        //instance<other
        other.setMinCost(20);
        assertEquals(-1, instance.compareTo(other));
        
        //instance>other
        other.setMinCost(5);
        assertEquals(1, instance.compareTo(other));
        
        //same mean values, but tie broken by max value
        other.setMinCost(10);
        other.setMaxCost(101);
        assertEquals(-1, instance.compareTo(other));
        assertEquals(instance.getMeanCost(), other.getMeanCost());
        
        //same mean values, tie broken by min value
        other.setMinCost(11);
        other.setMaxCost(100);
        assertEquals(-1, instance.compareTo(other));
        assertEquals(instance.getMeanCost(), other.getMeanCost());
        
    }

    /**
     * Test of toString method, of class CostOption.
     */
    @Test
    public void testToString() {
        System.out.println("toString");
        System.out.println(instance);
        
    }
}
