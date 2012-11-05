package earthmoverdistancesample;
import java.util.ArrayList;
/**
 *
 * @author wichtelwesen
 */
public class Signature {
    private ArrayList<BaseFeature> features;
    private ArrayList<Double> weights;
    
    Signature(ArrayList<BaseFeature> features, ArrayList<Double> weights){ 
        assert features.size()==weights.size(): assertString(features.size(), 
                                                            weights.size());
        this.features = new ArrayList<>(features);
        this.weights = new ArrayList<>(weights);
    }
    
    public int nFeatures(){
        return features.size();
    }
    public int nWeights(){
        return weights.size();
    }
    
    public ArrayList<BaseFeature> getFeatures(){
        return new ArrayList<>(features);
    }
    public ArrayList<Double> getWeights(){
        return new ArrayList<>(weights);
    }
    
    public double getWeightsNorm(){
        double s = 0;
        for (Double x: weights){
            s+=x;
        }
        return s;
    }
    
    public BaseFeature featureAt(int i){
        return features.get(i);
    }
    
    public double weightAt(int i){
        return weights.get(i);
    }
    
    public void setFeatures(ArrayList<BaseFeature> features){
        assert features.size()==nWeights(): assertString(features.size(),
                                                        nWeights());
        this.features = new ArrayList<>(features);
    }
    public void setWeights(ArrayList<Double> weights){
        assert nFeatures()==weights.size(): assertString(nFeatures(),
                                                        weights.size());
        this.weights = new ArrayList<>(weights);
    }
    
    public Simple2DArrayInterface distanceMatrix( Signature  other){
        assert nFeatures()==nWeights(): assertString(nFeatures(), nWeights());
        assert other.nFeatures()==other.nWeights(): "other "+assertString(
                                           other.nFeatures(), other.nWeights());
        Simple2DDoubleArray distanceArray = new Simple2DDoubleArray(nFeatures(), other.nFeatures());
        for(int i=0;i<nFeatures();i++){
            for(int j=0;j<other.nFeatures();j++){
                BaseFeature thisFeature = featureAt(i);
                BaseFeature thatFeature = other.featureAt(j);
                distanceArray.setValueAt(i, j, thisFeature.distance(thatFeature));
            }
        }
        return distanceArray;
    }
    
    
    public void normalizeWeights(){
        normalizeWeights(1.0);
    }
    public void normalizeWeights(double normValue){
        double weightSum = getWeightsNorm();
        double newVal;
        
        for (int i = 0; i<nWeights();i++){
            newVal = normValue*weights.get(i)/weightSum;
            weights.set(i, newVal);
        }
    }
    
    private String assertString(int nFeatures, int nWeights){
        return "features has "+ nFeatures+ "elements, weights has"
                                   +nWeights +"elements";
    }
}
