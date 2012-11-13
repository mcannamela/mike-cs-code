/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.Arrays;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Michael
 */
public class NDEntityTest {
    private NDEntity instance;
    private int[] defaultShape = {2,3,4};
    public NDEntityTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @Before
    public void setUp() {
        instance = new NDEntity(defaultShape);
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGetStrideAt() {
        System.out.println("getStrideAt");
        int[] shape = {3, 4, 5};
        int[] expResult = {1, 3, 12};
        int[] results = new int[shape.length];
        instance = new NDEntity(shape);
        for (int i=0;i<shape.length;i++){
            results[i] = instance.getStrideAt(i);
        }
        assertArrayEquals(expResult, results);
        
    }

    @Test
    public void testGetStrides() {
        System.out.println("getStrides");
        int[] expResult = {1,2,6};
        int[] result = instance.getStrides();
        assertArrayEquals(expResult, result);
        
    }

    @Test
    public void testIdxCopy() {
        System.out.println("idxCopy");
        int[] original = {0,1,2};
        int[] expResult = {0,1,2};
        int[] result = NDEntity.idxCopy(original);
        assertArrayEquals(expResult, result);
        
    }

    @Test
    public void testIdxToString() {
        System.out.println("idxToString");
        int[] idx = {0,1,2};
        String expResult = "[0, 1, 2]";
        String result = Arrays.toString(idx);
        assertEquals(expResult, result);
        
        int[] idx1d = {1};
        assertEquals("[1]", Arrays.toString(idx1d));
        
        int[] idx0d = {};
        assertEquals("[]", Arrays.toString(idx0d));
    }

    @Test
    public void testInitNElements() {
        System.out.println("initNElements");
        int nElements = 1;
        for(int i: defaultShape){
            nElements*=i;
        }
        assertEquals(nElements, instance.nElements());
    }

    @Test
    public void testInitStrides() {
        System.out.println("initStrides");
        int[] shape = {2, 3, 4};
        int[] expStrides = {1, 2, 6};
        instance = new NDEntity(shape);
        assertArrayEquals(expStrides, instance.getStrides());
        
    }

    @Test
    public void testNDimensions() {
        System.out.println("nDimensions");
        
        int expResult = 3;
        int result = instance.nDimensions();
        assertEquals(expResult, result);
        
    }

    @Test
    public void testNElements() {
        System.out.println("nElements");
        int[] shape = {3,4,5};
        instance = new NDEntity(shape);
        int expResult = 3*4*5;
        int result = instance.nElements();
        assertEquals(expResult, result);
        
    }
    
    @Test
    public void testFlattenIndex(){
        System.out.println("flattenIndex");
        int[][] testIndices = {
                                {0,0,0},
                                {1,0,0},
                                {0,1,0},
                                {1,1,0},
                                {0,2,0},
                                {0,0,1},
                                {0,1,1} 
        };
        int[] expecteds = {0,1,2,3,4,6,8};
        int[] results = new int[testIndices.length];
        for (int i=0;i<testIndices.length;i++){
            results[i] = instance.flattenIndex(testIndices[i]);
        }
        assertArrayEquals(expecteds, results);
    }
    
    @Test
    public void testBroadcastFlattenIndex(){
        System.out.println("BroadcastFlattenIndex");
        int[] shape = {1,2,3};
        instance = new NDEntity(shape);
        int[][] testIndices = {
                                {0,0,0},
                                {0,1,0},
                                {0,0,1},
                                {0,1,1},
                                {0,0,2},
                                {0,1,2} ,
                                {2,0,0},
                                {2,1,0},
                                {2,0,1},
                                {2,1,1},
                                {2,0,2},
                                {2,1,2}
        };
        System.out.println("flattener strides are:");
        System.out.println(Arrays.toString(instance.getFlattenerStrides()));
        
        instance.setBroadcasting(new int[] {0,1,1});
        System.out.println("flattener strides are now:");
        System.out.println(Arrays.toString(instance.getFlattenerStrides()));
        int[] expecteds = {0,1,2,3,4,5,0,1,2,3,4,5};
        int[] results = new int[testIndices.length];
        for (int i=0;i<testIndices.length;i++){
            results[i] = instance.flattenIndex(testIndices[i]);
        }
        System.out.println(Arrays.toString(results));
        assertArrayEquals(expecteds, results);
    }
    
}
