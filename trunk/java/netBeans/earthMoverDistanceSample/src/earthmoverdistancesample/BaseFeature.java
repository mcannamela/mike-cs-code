package earthmoverdistancesample;

/**
 *
 * @author wichtelwesen
 */
public class BaseFeature<T> {
    private T value;
    private GroundDistanceFunctionInterface<BaseFeature<T>> distanceFunction;
    
    BaseFeature(T value,GroundDistanceFunctionInterface<BaseFeature<T>> distanceFunction ){
        this.value = value;
        this.distanceFunction = distanceFunction;
    }
    
    public final double distance(BaseFeature<T> other){
        return distanceFunction.distance(this, other);
    }
    
    T getValue(){
        return value;
    }
}


