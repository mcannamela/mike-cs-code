/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;
/**
 *
 * @author mcannamela
 */
public abstract class AbstractScalingLaw {
    protected ProgramSize programSize;

    public AbstractScalingLaw() {
        programSize = new ProgramSize();
    }

    public AbstractScalingLaw(ProgramSize programSize) {
        this.programSize = programSize;
    }
    
    
    public abstract int scale(int cost);

    public ProgramSize getProgramSize() {
        return programSize;
    }

    public void setProgramSize(ProgramSize programSize) {
        this.programSize = programSize;
    }
    
    
}
