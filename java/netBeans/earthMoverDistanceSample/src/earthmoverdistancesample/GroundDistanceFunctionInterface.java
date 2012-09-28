/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package earthmoverdistancesample;

/**
 *
 * @author wichtelwesen
 */
public interface GroundDistanceFunctionInterface<T> {
    double distance(T firstFeature, T secondFeature);
}
