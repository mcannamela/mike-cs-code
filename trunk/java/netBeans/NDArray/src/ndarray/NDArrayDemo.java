/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ndarray;

/**
 *
 * @author mcannamela
 */
public class NDArrayDemo {
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        int[] shape = {2, 3, 4};
        int[] idx = new int[shape.length];
        String idxString;
        NDCounter iterator = new NDCounter(shape);
        
        while (iterator.hasNext()){
            idxString = "";
            idx = iterator.next();
            idxString+=idx[0];
            for(int i=1;i<idx.length;i++){
                idxString+=", "+idx[i];
            }
            System.out.println(idxString);
        }
        
        for (int i=1;i<1;i++){
            System.out.println("loop executed "+i+" times");
        }
    }
}
