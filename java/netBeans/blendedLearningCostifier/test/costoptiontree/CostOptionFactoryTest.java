/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.ProgramSize;
import blendedlearningprogram.StandardBlendedLearningModel;
import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import org.apache.commons.configuration.Configuration;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Michael
 */
public class CostOptionFactoryTest {
    private static final StandardBlendedLearningModel blendedLearningModel = new StandardBlendedLearningModel();
    private static final int nrStudents = 300;
    private static final int nrPeriods = 6;

    private static final Path testConfigurationPath = Paths.get(".","testConfigurationPath");
    
    private static ProgramSize programSize;
    private static ArrayList<Path> optionsConfigPaths = new ArrayList<>(3);

    private CostOptionFactory optionFactory;    
    private ArrayList<CostOption> costOptions = new ArrayList<>(3);
    private PropertiesConfiguration config; 
    
    public CostOptionFactoryTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        programSize = new ProgramSize(blendedLearningModel, nrStudents, nrPeriods);
    }
    
    @Before
    public void setUp() {
        optionFactory = new CostOptionFactory(programSize);    
        
        optionsConfigPaths.add(Paths.get( "optionOne.option"));
        optionsConfigPaths.add(Paths.get( "optionTwo.option"));
        optionsConfigPaths.add(Paths.get( "optionThree.option"));
    }
    
    @After
    public void tearDown() {
        optionsConfigPaths.clear();
        costOptions.clear();
    }

    @Test 
    public void testMakeCostOption_0args(){
        System.out.println("makeCostOption_0args");
        CostOptionFactory defaultFactory = new CostOptionFactory();
        CostOption option = defaultFactory.makeCostOption();
        System.out.println(option);
    }
    @Test
    public void testMakeCostOption() {
        System.out.println("makeCostOption_config");
        for (Path path : optionsConfigPaths){
            try{
                config = new PropertiesConfiguration(testConfigurationPath.resolve(path).toFile());
                costOptions.add(optionFactory.makeCostOption(config));
            }
            catch(ConfigurationException e){
                System.out.println(e);
            }
        }
            int[] minCosts = {100,200,300};
            int[] maxCosts = {1000,2000,3000};
            int[] selectedCosts = {500, 1100, 1000};
            int[] scaledCosts = {500, 330000, 2000};
            String[] labels = {"optionOne", "optionTwo", "optionThree"};
            
            for (int i=0;i<3;i++){
                assertEquals(minCosts[i], costOptions.get(i).getMinCost());
                assertEquals(maxCosts[i], costOptions.get(i).getMaxCost());
                assertEquals(selectedCosts[i], costOptions.get(i).getSelectedCost());
                assertEquals(scaledCosts[i], costOptions.get(i).getScaledCost());
                assertEquals(labels[i], costOptions.get(i).getLabel());
            }
            
        
    }
}
