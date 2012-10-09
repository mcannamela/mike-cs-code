/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package blendedlearningprogram;

/**
 *
 * @author wichtelwesen
 */
public class LinearInTeachersScalingLaw extends AbstractScalingLaw{

    public LinearInTeachersScalingLaw(ProgramSize programSize) {
        super(programSize);
    }
    
    @Override
    public int scale(int cost) {
        return cost*programSize.getNrTeachers();
    }
}
