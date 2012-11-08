/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import org.junit.After;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author mcannamela
 */
public class SliceTest {
    private Slice instance;
    
    
    public SliceTest() {
    }
    
    @Before
    public void setUp() {
        instance = new Slice(0,10,1);
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of isRaw method, of class Slice.
     */
    @Test
    public void testIsRaw() {
        System.out.println("isRaw");
        
        assertFalse(instance.isRaw());
        
        instance = new Slice();
        assertTrue(instance.isRaw());
        
        instance = new Slice(0,-1,1);
        assertTrue(instance.isRaw());
        
    }

    /**
     * Test of nIndices method, of class Slice.
     */
    @Test
    public void testNIndices() {
        System.out.println("nIndices");
        
        assertEquals(10, instance.nIndices());
        
        instance = new Slice(1,10,1);
        assertEquals(9, instance.nIndices());
        
        instance = new Slice(0,1,1);
        assertEquals(1, instance.nIndices());
        
        instance = new Slice(10,0,-1);
        assertEquals(10, instance.nIndices());
        
        instance = new Slice(2,12,1);
        assertEquals(10, instance.nIndices());
        
        instance = new Slice(2,12,2);
        assertEquals(5, instance.nIndices());
        instance = new Slice(2,13,2);
        assertEquals(6, instance.nIndices());
        instance = new Slice(2,14,2);
        assertEquals(6, instance.nIndices());
        
        instance = new Slice(2,11,3);
        assertEquals(3, instance.nIndices());
        
        instance = new Slice(2,12,3);
        assertEquals(4, instance.nIndices());
        
        instance = new Slice(2,13,3);
        assertEquals(4, instance.nIndices());
        instance = new Slice(2,14,3);
        assertEquals(4, instance.nIndices());
        instance = new Slice(2,15,3);
        assertEquals(5, instance.nIndices());
    }

    /**
     * Test of getCookedSlice method, of class Slice.
     */
    @Test
    public void testGetCookedSlice() {
        System.out.println("getCookedSlice");
        int stop = 10;
        instance = new Slice();
        
        Slice result = instance.getCookedSlice(stop);
        assertEquals(stop, result.stop());
        assertEquals(0, result.start());
        assertEquals(1, result.step());
        
        instance = new Slice(1,-1,3);
        
        result = instance.getCookedSlice(stop);
        assertEquals(stop, result.stop());
        assertEquals(1, result.start());
        assertEquals(3, result.step());
        
    }

    /**
     * Test of start method, of class Slice.
     */
    @Test
    public void testStart() {
        System.out.println("start");
        int expResult = 0;
        int result = instance.start();
        assertEquals(expResult, result);
    }

    /**
     * Test of stop method, of class Slice.
     */
    @Test
    public void testStop() {
        System.out.println("stop");
        int expResult = 10;
        int result = instance.stop();
        assertEquals(expResult, result);
    }

    /**
     * Test of step method, of class Slice.
     */
    @Test
    public void testStep() {
        System.out.println("step");
        int expResult = 1;
        int result = instance.step();
        assertEquals(expResult, result);
    }

    /**
     * Test of get method, of class Slice.
     */
    @Test
    public void testGet() {
        System.out.println("get");
        int[] indices, results;
        
        instance = new Slice(0,10,1);
        indices = new int[]{0,1,2,3,4,5,6,7,8,9};
        results = new int[indices.length];
        assertEquals(indices.length, instance.nIndices());
        for (int i=0;i<indices.length;i++){
//            System.out.println(""+i);
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
        
        instance = new Slice(0,10,2);
        indices = new int[]{0,2,4,6,8};
        results = new int[indices.length];
        for (int i=0;i<indices.length;i++){
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
        
        instance = new Slice(1,10,2);
        indices = new int[]{1,3,5,7,9};
        results = new int[indices.length];
        for (int i=0;i<indices.length;i++){
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
        
        instance = new Slice(2,13,3);
        indices = new int[]{2,5,8,11};
        results = new int[indices.length];
        for (int i=0;i<indices.length;i++){
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
        
        instance = new Slice(2,14,3);
        indices = new int[]{2,5,8,11};
        results = new int[indices.length];
        for (int i=0;i<indices.length;i++){
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
        
        instance = new Slice(2,15,3);
        indices = new int[]{2,5,8,11,14};
        results = new int[indices.length];
        for (int i=0;i<indices.length;i++){
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
        
        instance = new Slice(15,2,-3);
        indices = new int[]{15,12,9,6,3};
        results = new int[indices.length];
        for (int i=0;i<indices.length;i++){
            results[i] = instance.get(i);
        }
        assertArrayEquals(indices, results);
    }
}
