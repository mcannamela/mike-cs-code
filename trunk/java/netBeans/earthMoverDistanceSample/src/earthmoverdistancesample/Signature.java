package earthmoverdistancesample;
import java.util.ArrayList;
/**
 *
 * @author wichtelwesen
 */
class Signature {
    private ArrayList<Feature> features;
    private ArrayList<Double> weights;
    
    Signature(ArrayList<Feature> features, ArrayList<Double> weights){ 
        assert features.size()==weights.size(): assertString(features.size(), 
                                                            weights.size());
        this.features = features;
        this.weights = weights;
    }
    
    public int nFeatures(){
        return features.size();
    }
    public int nWeights(){
        return weights.size();
    }
    
            
    public ArrayList<Feature> getFeatures(){
        return features;
    }
    public ArrayList<Double> getWeights(){
        return weights;
    }
    
    public Feature featureAt(int i){
        return features.get(i);
    }
    
    public double weightAt(int i){
        return weights.get(i);
    }
    
    public void setFeatures(ArrayList<Feature> features){
        assert features.size()==nWeights(): assertString(features.size(),
                                                        nWeights());
        this.features = features;
    }
    public void setWeights(ArrayList<Double> weights){
        assert nFeatures()==weights.size(): assertString(nFeatures(),
                                                        weights.size());
        this.weights = weights;
    }
    
    public Simple2DArrayInterface distanceMatrix( Signature  other){
        assert nFeatures()==nWeights(): assertString(nFeatures(), nWeights());
        assert other.nFeatures()==other.nWeights(): "other "+assertString(
                                           other.nFeatures(), other.nWeights());
        Simple2DArray distanceArray = new Simple2DArray(nFeatures(), other.nFeatures());
        for(int i=0;i<nFeatures();i++){
            for(int j=0;j<nFeatures();j++){
                Feature thisFeature = featureAt(i);
                Feature thatFeature = other.featureAt(j);
                distanceArray.setValueAt(i, j, thisFeature.distance(thatFeature));
            }
        }
        
        return distanceArray;
                
    }
    
    public void normalizeWeights(double normValue){
        
    }
    
    private String assertString(int nFeatures, int nWeights){
        return "features has "+ nFeatures+ "elements, weights has"
                                   +nWeights +"elements";
    }
    
    
    
    
}
