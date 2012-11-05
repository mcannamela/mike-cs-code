/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package earthmoverdistancesample;

import java.util.ArrayList;


/**
 *
 * @author wichtelwesen
 */
public class Simple2DDoubleArray implements Simple2DArrayInterface<Double>{
    private double[][] valueArray;
    private int nRows;
    private int nColumns;
    
    public Simple2DDoubleArray(int nRows, int nColumns, Double val){
        this.nRows = nRows;
        this.nColumns = nColumns;
        setConstantValue(val);
    }
    
    public Simple2DDoubleArray(int nRows, int nColumns){
        this.nRows = nRows;
        this.nColumns = nColumns;
        setConstantValue(0.0);   
    }
    
    private void setConstantValue(double val){
        valueArray = new double[nRows][nColumns];
        for(int i=0; i<nRows;i++){
            for(int j=0; j<nColumns;j++){
                valueArray[i][j]=val;
            }
        }
    }
    
    public static Simple2DDoubleArray zeros(int nRows, int nColumns){
        return new Simple2DDoubleArray(nRows,nColumns, 0.0);
    }
    public static Simple2DDoubleArray ones(int nRows, int nColumns){
        return new Simple2DDoubleArray(nRows,nColumns, 1.0);
    }
    public static double[] ArrayListToDoubleArray(ArrayList<Double> arr){
        double[] A = new double[arr.size()];
        for (int i=0;i<arr.size();i++){
            A[i] = arr.get(i);
        }
        return A;
    }
    public Simple2DDoubleArray zeros(){
        return zeros(nRows, nColumns);
    }
    @Override
    public Simple2DDoubleArray onesInRow(int row){
        Simple2DDoubleArray arr = zeros();
        for(int j=0;j<nColumns;j++){
            arr.setValueAt(row,j,1.0);
        }
        return arr;
    }
    @Override
    public Simple2DDoubleArray onesInColumn(int column){
        Simple2DDoubleArray arr = zeros();
        for(int i=0;i<nRows;i++){
            arr.setValueAt(i,column,1.0);
        }
        return arr;
    }

    @Override
    public Simple2DArrayInterface ravel(ArrayList<Double> flatArray, int nRows, int nColumns) {
        assert flatArray.size()==nRows*nColumns:"bad ravel size, flatArray is size:"+flatArray.size()+
                   ", while nRows, nColumns is "+nRows +", "+nColumns;
        
        Simple2DDoubleArray raveledArray = new Simple2DDoubleArray(nRows, nColumns);
        int cnt=0;
        for(int i=0; i<nRows;i++){
            for(int j=0; j<nColumns;j++){
                raveledArray.setValueAt(i,j, flatArray.get(cnt));
                cnt++;
            }
        }
        return raveledArray;
        
    }

    @Override
    public Simple2DArrayInterface ravel(ArrayList<Double> flatArray) {
        return ravel(flatArray, shape()[0], shape()[1]);
    }
    
        

    @Override
    public int[] shape(){
        final int[] s = {nRows, nColumns};
        return s;
    }
    @Override
    public Double valueAt(int row, int column){
        return valueArray[row][column];
    }
    @Override
    public void setValueAt(int row, int column, Double value ){
        valueArray[row][column]=value;
    }
    
    @Override
    public Simple2DArrayInterface elementWiseOperate(Simple2DArrayInterface<Double> other,
                                    BinaryOperatorInterface<Double> operator){
        Simple2DDoubleArray result = new Simple2DDoubleArray(nRows,nColumns);
        for(int i=0; i<nRows;i++){
            for(int j=0; j<nColumns;j++){
                result.setValueAt(i,j,operator.operate(valueAt(i,j),other.valueAt(i,j)));
            }
        }
        return result;  
    }
   
    @Override
    public Simple2DArrayInterface elementWiseAdd(Simple2DArrayInterface other){
        doubleAdder adder = new doubleAdder();
        return elementWiseOperate(other, adder);
    }
    @Override
    public Simple2DArrayInterface elementWiseMultiply(Simple2DArrayInterface other){
        doubleMultiplier multiplier = new doubleMultiplier();
        return elementWiseOperate(other, multiplier);
    }
    
    @Override
    public ArrayList<Double> flatten(){
        ArrayList<Double> flatValueArray = new ArrayList<>(nRows*nColumns);
        int cnt=0;
        for(int i=0; i<nRows;i++){
            for(int j=0; j<nColumns;j++){
                flatValueArray.set(cnt, valueArray[i][j]);
                cnt++;
            }
        }
        return flatValueArray;
    }
    
    @Override
    public Double sum(){
        ArrayList<Double> flatValueArray = flatten();
        double s = 0;
        for(double element:flatValueArray){
            s+=element;
            }
        return s;
    }
    
    @Override
    public ArrayList<Double> sumOfColumns(){
        ArrayList<Double> s = new ArrayList<>();
        for(int i=0; i<nRows;i++){
            s.set(i, 0.0);
            for(int j=0; j<nColumns;j++){
                s.set(i, s.get(i)+valueArray[i][j]);
            }
        }
        return s;
    }
    
    @Override
    public ArrayList<Double> sumOfRows(){
        ArrayList<Double> s = new ArrayList<>();
        for(int j=0; j<nColumns;j++){
            s.set(j,0.0);
            for(int i=0; i<nRows;i++){
                s.set(j, s.get(j)+valueArray[i][j]);
            }
        }
        return s;
    }
    
}
class doubleAdder implements BinaryOperatorInterface<Double>{
    @Override
    public Double operate(Double x, Double y){
        return x+y;
    }
}

class doubleMultiplier implements BinaryOperatorInterface<Double>{
    @Override
    public Double operate(Double x, Double y){
        return x*y;
    }
}