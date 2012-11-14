/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.Arrays;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author mcannamela
 */
public class NDDoubleArrayTest {
    NDDoubleArray instance;
    int[] shape = {2, 3, 4};
    public NDDoubleArrayTest() {
    }
    
    @Before
    public void setUp() {
        instance = NDDoubleArray.range(shape);
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of assign method, of class NDDoubleArray.
     */
    @Test
    public void testAssign_NDDoubleArray() {
        System.out.println("assign - NDDoubleArray");
        NDDoubleArray other = NDDoubleArray.ones(shape);
        instance.assign(other);
        assertTrue(instance.equals(other));
        other = NDDoubleArray.zeros(shape);
        assertFalse(instance.equals(other));
        instance.assign(other);
        assertTrue(instance.equals(other));
    }

    /**
     * Test of assign method, of class NDDoubleArray.
     */
    @Ignore
    @Test
    public void testAssign_SliceArr_NDDoubleArray() {
        System.out.println("assign");
        Slice[] rawSlices = null;
        NDDoubleArray other = null;
        NDDoubleArray instance = null;
        instance.assign(rawSlices, other);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of assign method, of class NDDoubleArray.
     */
    
    @Ignore
    @Test
    public void testAssign_3args() {
        System.out.println("assign");
        Slice[] rawSlices = null;
        NDDoubleArray other = null;
        Slice[] rawOtherSlices = null;
        NDDoubleArray instance = null;
        instance.assign(rawSlices, other, rawOtherSlices);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of plusAssign method, of class NDDoubleArray.
     */
    @Ignore
    @Test
    public void testPlusAssign_NDDoubleArray() {
        System.out.println("plusAssign");
        NDDoubleArray other = null;
        NDDoubleArray instance = null;
        instance.plusAssign(other);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of plusAssign method, of class NDDoubleArray.
     */
    @Ignore
    @Test
    public void testPlusAssign_SliceArr_NDDoubleArray() {
        System.out.println("plusAssign");
        Slice[] rawSlices = null;
        NDDoubleArray other = null;
        NDDoubleArray instance = null;
        instance.plusAssign(rawSlices, other);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of plusAssign method, of class NDDoubleArray.
     */
    @Ignore
    @Test
    public void testPlusAssign_3args() {
        System.out.println("plusAssign");
        Slice[] rawSlices = null;
        NDDoubleArray other = null;
        Slice[] rawOtherSlices = null;
        NDDoubleArray instance = null;
        instance.plusAssign(rawSlices, other, rawOtherSlices);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of get method, of class NDDoubleArray.
     */
    @Ignore
    @Test
    public void testGet() {
        System.out.println("get");
        Slice[] rawSlices = null;
        NDDoubleArray instance = null;
        NDDoubleArray expResult = null;
        NDDoubleArray result = instance.get(rawSlices);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getElement method, of class NDDoubleArray.
     */
    
    @Test
    public void testGetElement_int() {
        System.out.println("getElement");
        int index;
        NDCounter counter = instance.getNewCounter();
        while (counter.hasNext()){
            index = counter.nextFlat();
            assertEquals((double) index, instance.getElement(index), 1e-15);
        }
        
    }

    /**
     * Test of getElement method, of class NDDoubleArray.
     */
    @Test
    public void testGetElement_intArr() {
        System.out.println("getElement");
        int[] nDIndex;
        int index;
        NDCounter counter = instance.getNewCounter();
        while (counter.hasNext()){
            nDIndex = counter.next();
            index = counter.flattenIndex(nDIndex);
            assertEquals(index, instance.getElement(nDIndex), 1e-15);
            assertEquals(instance.getElement(index),instance.getElement(nDIndex), 1e-15);
            System.out.println("nidx, idx, value = "+Arrays.toString(nDIndex)+", "+index+", "+instance.getElement(nDIndex));
        }
        
        for (int i = 0;i<instance.nElements();i++){
            assertEquals(i,instance.getElement(i), 1e-15);
        }
        
    }

    /**
     * Test of getArray method, of class NDDoubleArray.
     */
    @Test
    public void testGetArray() {
        System.out.println("getArray");
        int p = 1;
        for(int i:shape){
            p*=i;
        }
        assertEquals(p, instance.nElements());
        for (int i = 0;i<p;i++){
            assertEquals((double) i, instance.getElement(i), 1e-15);
        }
    }

    /**
     * Test of setElement method, of class NDDoubleArray.
     */
    @Test
    public void testSetElement_int_double() {
        System.out.println("setElement");
        
        double value = 30.1;
        for (int i = 0; i<instance.nElements();i++){
            instance.setElement(i, value);
            assertEquals(value, instance.getElement(i), 1e-15);
        } 
        
    }

    /**
     * Test of setElement method, of class NDDoubleArray.
     */
    @Test
    public void testSetElement_intArr_double() {
        System.out.println("setElement");
        double value = 20.1;
        NDCounter counter = instance.getNewCounter();
        while (counter.hasNext()){
            instance.setElement(counter.next(), value);
        }
        for (double x: instance.getArray()){
            assertEquals(x, value, 1e-15);
        }
    }

    /**
     * Test of setElement method, of class NDDoubleArray.
     */
    @Ignore
    @Test
    public void testSetElement_3args() {
        System.out.println("setElement");
        int[] nDIndex = null;
        double value = 0.0;
        IndexFlattener flattener = null;
        NDDoubleArray instance = null;
        instance.setElement(nDIndex, value, flattener);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}
