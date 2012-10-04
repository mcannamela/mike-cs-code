/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;
import blendedlearningprogram.ProgramSize;
/**
 *
 * @author mcannamela
 */
public abstract class ScalingLaw {
    private ProgramSize programSize;
    
    public abstract int scale(int cost);
    
}
