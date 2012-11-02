/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package earthmoverdistancesample;

/**
 *
 * @author mcannamela
 */
public class TransportProblem implements Runnable{
    private AbstractSimple2DArray costMatrix;
    private Signature supplySignature;
    private Signature demandSignature;

    
    public TransportProblem(Signature supplySignature, Signature demandSignature) {    
        this.supplySignature = supplySignature;
        this.demandSignature = demandSignature;
        
        setCostMatrix();
    }
    
    @Override
    public void run() {
        
    }
    
    private void setCostMatrix(){
        costMatrix = supplySignature.distanceMatrix(demandSignature);
    }
    
    

    public void setDemandSignature(Signature demandSignature) {
        this.demandSignature = demandSignature;
        setCostMatrix();
    }

    public void setSupplySignature(Signature supplySignature) {
        this.supplySignature = supplySignature;
        setCostMatrix();
    }
    
    

    
    
    
    
    
    
    
    
    
    
}
