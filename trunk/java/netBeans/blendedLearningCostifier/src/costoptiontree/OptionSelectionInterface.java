/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;
import java.util.ArrayList;
/**
 *
 * @author wichtelwesen
 */
public interface OptionSelectionInterface<T> {
    public ArrayList<Double> getOptionBlendingFactors();
    
    public void select(T selection);
    public int nOptions();
    
}
