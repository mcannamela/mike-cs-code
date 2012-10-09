/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;

import blendedlearningprogram.ProgramSize;

/**
 *
 * @author wichtelwesen
 */
public class CostOptionNodeFactory {
    private ProgramSize programSize;

    public CostOptionNodeFactory(ProgramSize programSize) {
        this.programSize = programSize;
    }
    
    public CostOptionNode makeCostOptionNode(){
        return new CostOptionNode();
    }
}
