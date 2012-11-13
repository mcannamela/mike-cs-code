/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

import java.util.Arrays;
import org.junit.After;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 * @author mcannamela
 */
public class SliceCounterTest {
    
    SliceCounter instance;
    int[] defaultShape = {2,3,4};
    public SliceCounterTest() {
    }
    
    
    @Before
    public void setUp() {
        instance = new SliceCounter(defaultShape);
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
        instance = new SliceCounter(shape);
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
        int[] idx,counterPosition;
        int cnt=0;
        while (instance.hasNext()){
            counterPosition = instance.getCounterPosition();
            idx = instance.next();
            
            assertArrayEquals(expResults[cnt], idx);
            cnt++;
            
            System.out.println("1-D index: "+instance.nDTo1D(counterPosition)+
                                ", position: "+Arrays.toString(counterPosition)+
                                ", value:"+Arrays.toString(idx));
        }
    }
    
    @Test
    public void testNextWithSkips() {
        System.out.println("nextWithSkips");
        Slice[] slices = new Slice[]{
                            new Slice(0,4,2),
                            new Slice(1,7,3),
                            new Slice(2,0,-1)
                            };
                
        instance = new SliceCounter(slices);
        System.out.println(String.format("slice counter has %d dimensions", instance.nDimensions()));
        
        int[][] expResults = {
                              {0,1,2},
                              {0,1,1},
                              {0,4,2},
                              {0,4,1},
                              {2,1,2},
                              {2,1,1},
                              {2,4,2},
                              {2,4,1}
                             };
        
        int[] idx, counterPosition;
        int cnt=0;
        while (instance.hasNext()){
            counterPosition = instance.getCounterPosition();
            idx = instance.next();
            
            assertArrayEquals(expResults[cnt], idx);
            cnt++;
            
//            System.out.println("1-D index: "+instance.nDTo1D(counterPosition)+
//                                ", position: "+Arrays.toString(counterPosition)+
//                                ", value:"+Arrays.toString(idx));
        }
    }
    
    @Ignore
    @Test
    public void testNextWithBroadcast() {
        System.out.println("nextWithBroadcast");
        Slice[] slices = new Slice[]{
                            new Slice(0,4,2),
                            new Slice(1,7,3),
                            new Slice(2,0,-1)
                            };
                
        instance = new SliceCounter(slices);
        instance.setStrideAt(1, 0);
        System.out.println(String.format("slice counter has %d dimensions", instance.nDimensions()));
        System.out.println("Strides are "+Arrays.toString(instance.getStrides()));
        
        int[][] expResults = {
                              {0,1,2},
                              {0,1,1},
                              {0,1,2},
                              {0,1,1},
                              {2,1,2},
                              {2,1,1},
                              {2,1,2},
                              {2,1,1}
                             };
        
        int[] idx, counterPosition;
        int cnt=0;
        while (instance.hasNext()){
            counterPosition = instance.getCounterPosition();
            idx = instance.next();
            
            assertArrayEquals(expResults[cnt], idx);
            cnt++;
            
            System.out.println("1-D index: "+instance.nDTo1D(counterPosition)+
                                ", position: "+Arrays.toString(counterPosition)+
                                ", value:"+Arrays.toString(idx));
        }
    }
}
