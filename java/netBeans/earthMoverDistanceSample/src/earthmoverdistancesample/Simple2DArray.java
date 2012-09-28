/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package earthmoverdistancesample;

/**
 *
 * @author wichtelwesen
 */
public class Simple2DArray implements Simple2DArrayInterface {
    private double[][] valueArray;
    private int nRows;
    private int nColumns;
    
    public Simple2DArray(int nRows, int nColumns){
        valueArray = new double[nRows][nColumns];
        for(int i=0; i<nRows;i++){
            for(int j=0; j<nColumns;j++){
                valueArray[i][j]=0;
            }
        }
    }
        

    @Override
    public int[] shape(){
        final int[] s = {nRows, nColumns};
        return s;
    }
    @Override
    public double valueAt(int row, int column){
        return valueArray[row][column];
    }
    @Override
    public void setValueAt(int row, int column, double value ){
        valueArray[row][column]=value;
    }
    
    @Override
    public Simple2DArrayInterface elementWiseOperate(Simple2DArrayInterface other,
                                    doubleOperatorInterface operator){
        Simple2DArray result = new Simple2DArray(nRows,nColumns);
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
    public double[] flatten(){
        double[] flatValueArray = new double[nRows*nColumns];
        int cnt=0;
        for(int i=0; i<nRows;i++){
            for(int j=0; j<nColumns;j++){
                flatValueArray[cnt]=valueArray[i][j];
            }
        }
        return flatValueArray;
    }
    
    @Override
    public double sum(){
        double[] flatValueArray = flatten();
        double s = 0;
        for(double element:flatValueArray){
            s+=element;
            }
        return s;
    }
    
    @Override
    public double[] sumOfColumns(){
        double[] s = new double[nRows];
        for(int i=0; i<nRows;i++){
            s[i]=0;
            for(int j=0; j<nColumns;j++){
                s[i]+=valueArray[i][j];
            }
        }
        return s;
    }
    
    @Override
    public double[] sumOfRows(){
        double[] s = new double[nColumns];
        for(int j=0; j<nColumns;j++){
            s[j]=0;
            for(int i=0; i<nRows;i++){
                s[j]+=valueArray[i][j];
            }
        }
        return s;
    }
    
}

class doubleAdder implements doubleOperatorInterface{
    @Override
    public double operate(double x, double y){
        return x+y;
    }
}

class doubleMultiplier implements doubleOperatorInterface{
    @Override
    public double operate(double x, double y){
        return x*y;
    }
}