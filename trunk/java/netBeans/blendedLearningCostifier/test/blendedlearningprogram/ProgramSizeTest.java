/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Ignore;

/**
 *
 * @author wichtelwesen
 */
public class ProgramSizeTest {
    private static ProgramSize defaultInstance;
    private ProgramSize instance;
    private BlendedLearningModelInterface blendedLearningModel;
    private int nrStudents;
    private int nrPeriods;
    
    public ProgramSizeTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        defaultInstance = new ProgramSize();
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp(){
        blendedLearningModel = new StandardBlendedLearningModel();
        nrStudents = 300;
        nrPeriods = 6;
        instance = new ProgramSize(blendedLearningModel, nrStudents, nrPeriods);
    }

    /**
     * Test of getNrPeriods method, of class ProgramSize.
     */
    @Test
    public void testGetNrPeriods() {
        System.out.println("getNrPeriods");
        assertEquals(8, defaultInstance.getNrPeriods());
        assertEquals(nrPeriods, instance.getNrPeriods());
        
    }

    /**
     * Test of getNrStudents method, of class ProgramSize.
     */
    @Test
    public void testGetNrStudents() {
        System.out.println("getNrStudents");
        assertEquals(100, defaultInstance.getNrStudents());
        assertEquals(nrStudents, instance.getNrStudents());
    }

    /**
     * Test of getBlendedLearningModel method, of class ProgramSize.
     */
    @Test
    @Ignore
    public void testGetBlendedLearningModel() {
        System.out.println("getBlendedLearningModel");
        ProgramSize instance = new ProgramSize();
        BlendedLearningModelInterface expResult = null;
        BlendedLearningModelInterface result = instance.getBlendedLearningModel();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setNrStudents method, of class ProgramSize.
     */
    @Test
    public void testSetNrStudents() {
        System.out.println("setNrStudents");
        int newNrStudents = 80;
        instance.setNrStudents(newNrStudents);
        assertEquals(newNrStudents, instance.getNrStudents());
    }

    /**
     * Test of setNrPeriods method, of class ProgramSize.
     */
    @Test
    @Ignore
    public void testSetNrPeriods() {
        System.out.println("setNrPeriods");
        int nrPeriods = 0;
        ProgramSize instance = new ProgramSize();
        instance.setNrPeriods(nrPeriods);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setBlendedLearningModel method, of class ProgramSize.
     */
    @Test
    @Ignore
    public void testSetBlendedLearningModel() {
        System.out.println("setBlendedLearningModel");
        BlendedLearningModelInterface blendedLearningModel = null;
        ProgramSize instance = new ProgramSize();
        instance.setBlendedLearningModel(blendedLearningModel);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getNrTeachers method, of class ProgramSize.
     */
    @Test
    public void testGetNrTeachers() {
        System.out.println("getNrTeachers");
        assertEquals(1, defaultInstance.getNrTeachers());
        assertEquals(2, instance.getNrTeachers());
        instance.setNrStudents(301);
        assertEquals(3, instance.getNrTeachers());
    }

    /**
     * Test of getNrSimultaneousSections method, of class ProgramSize.
     */
    @Test
    public void testGetNrSimultaneousSections() {
        System.out.println("getNrSimultaneousSections");
        assertEquals(1, defaultInstance.getNrSimultaneousSections());
        assertEquals(2, instance.getNrSimultaneousSections());
        instance.setNrStudents(301);
        assertEquals(3, instance.getNrTeachers());
    }
}
