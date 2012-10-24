/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.ProgramSize;
import blendedlearningprogram.StandardBlendedLearningModel;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.apache.commons.configuration.ConfigurationException;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Michael
 */
public class CostOptionNodeFactoryTest {
    private static final StandardBlendedLearningModel blendedLearningModel = new StandardBlendedLearningModel();
    private static final int nrStudents = 300;
    private static final int nrPeriods = 6;
    private static final Path testConfigurationPath = Paths.get(".","testConfigurationPath");
    private static ProgramSize programSize;
    
    private CostOptionNodeFactory instance;
    
    public CostOptionNodeFactoryTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        programSize = new ProgramSize(blendedLearningModel, nrStudents, nrPeriods);
    }
    
    @Before
    public void setUp() {
        instance = new CostOptionNodeFactory(programSize);
    }

    @Test
    public void testNOptions() {
        System.out.println("nOptions");
        
        
        int expResult = 3;
        int result = CostOptionNodeFactory.nOptions(testConfigurationPath);
        assertEquals(expResult, result);
        
    }

    @Test
    public void testNChildren() {
        System.out.println("nChildren");
                
        int expResult = 2;
        int result = CostOptionNodeFactory.nChildren(testConfigurationPath);
        assertEquals(expResult, result);
    }

    @Test
    public void testMakeCostOptionNode() throws ConfigurationException {
        System.out.println("makeCostOptionNode");
        CostOptionNode parent = null;
        CostOptionNode root = instance.makeCostOptionNode(testConfigurationPath, parent);
        
        System.out.println(root);
//        assertEquals(2, root.nChildren());
        assertEquals(3, root.nCostOptions());
        assertEquals(2, root.getChildren().get(0).nCostOptions());
        assertEquals(2, root.getChildren().get(1).nCostOptions());
        
        
    }
}
