/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author Michael
 */
public class NDCounterTest {
    NDCounter instance;
    int[] defaultShape = {2,3,4};
    public NDCounterTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @Before
    public void setUp() {
        instance = new NDCounter(defaultShape);
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSetStrides() {
        System.out.println("setStrides");
        int[] strides = {2, 1, 0};
        instance.setStrides(strides);
        for (int i=0;i<3;i++){
            assertEquals(strides[i], instance.getStrides()[i]);
        }
    }

    @Test
    public void testSetStrideAt() {
        System.out.println("setStrideAt");
        int index = 2;
        int stride = 0;
        instance.setStrideAt(index, stride);
        assertEquals(stride, instance.getStrideAt(index));
    }

    @Test
    public void testHasNext() {
        System.out.println("hasNext");
        int nElements = 1;
        for (int i:defaultShape){
            nElements*=i;
        }
        
        int elCount = 0;
        boolean expResult;
        for (int i=0;i<nElements;i++){
            if (i<(nElements)){
                assertTrue(instance.hasNext());
                instance.next();
                elCount++;
            }
            else{
                assertFalse(instance.hasNext());
            }
            
        }
        assertEquals(elCount, instance.nElements());
    }

    @Test
    public void testNext() {
        System.out.println("next");
        int[] shape = {2,2,2};
        instance = new NDCounter(shape);
        int[][] expResults = {
                              {0,0,0},
                              {0,0,1},
                              {0,1,0},
                              {0,1,1},
                              {1,0,0},
                              {1,0,1},
                              {1,1,0},
                              {1,1,1}
        };
        int[] idx;
        int cnt=0;
        while (instance.hasNext()){
            idx = instance.next();
            System.out.println(NDEntity.idxToString(idx));
            for(int i=0;i<shape.length;i++){
                assertEquals(expResults[cnt][i], idx[i]);
            }
            cnt++;
        }
    }

    @Test
    @Ignore
    public void testNewIndex() {
        System.out.println("newIndex");
        NDCounter instance = null;
        int[] expResult = null;
        int[] result = instance.newIndex();
        assertArrayEquals(expResult, result);
        fail("The test case is a prototype.");
    }
}
