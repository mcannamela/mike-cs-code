package earthmoverdistancesample;

/**
 *
 * @author wichtelwesen
 */
public abstract class Feature<T> {
    private T value;
    private GroundDistanceFunctionInterface<Feature<T>> distanceFunction;
    
    Feature(T value,GroundDistanceFunctionInterface<Feature<T>> distanceFunction ){
        this.value = value;
        this.distanceFunction = distanceFunction;
    }
    
    double distance(Feature<T> other){
        return distanceFunction.distance(this, other);
    }
    
    T getValue(){
        return value;
    }
}


