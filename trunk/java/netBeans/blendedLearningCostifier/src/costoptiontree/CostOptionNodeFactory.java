/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;
import blendedlearningprogram.ProgramSize;
import java.io.File;
import java.io.IOException;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.Collections;
import org.apache.commons.configuration.*;


/**
 *
 * @author wichtelwesen
 */
public class CostOptionNodeFactory {
    private ProgramSize programSize;
    
    public static final String NODE_DESCRIPTION_KEY = "description";
    public static final String NODE_CONFIG_FILENAME = "node.config";
    public static final String OPTION_SUFFIX = "option";

    public CostOptionNodeFactory(ProgramSize programSize) {
        this.programSize = programSize;
    }
    
    public static String reconstituteCommaConfigString(String[] arrayOfStrings){
        String description = "";
        if (arrayOfStrings.length>0){
            for (String d: arrayOfStrings){
                description+=d+", ";
            }
            description = description.substring(0, description.length()-2);
        }
        
        return description;
    }
    
    public static void assertExists(Path path){
        assert Files.exists(path) && !Files.notExists(path): 
                "file DNE or can't be accessed:"+path;
    }
    static ArrayList<Path> getCostOptionPaths(Path nodePath){
        assertExists(nodePath);
        
        ArrayList<Path> costOptionPaths = new ArrayList<>();
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(nodePath, "*.{"+OPTION_SUFFIX+"}")) {
            for (Path path: stream) {
                costOptionPaths.add(path);
            }
        } catch (IOException | DirectoryIteratorException x) {
            System.err.println(x);
        }
        return costOptionPaths;
    }
    static int nOptions(Path costOptionsPath){
        return getCostOptionPaths(costOptionsPath).size();
    }
    
    static ArrayList<Path> getChildNodePaths(Path parentPath){
        assertExists(parentPath);
         ArrayList<Path> childPaths = new ArrayList<>();
         
         try (DirectoryStream<Path> stream = Files.newDirectoryStream(parentPath)) {
            for (Path path: stream) {
                if (Files.isDirectory(path)){
                    childPaths.add(path);
                }
            }
        } catch (IOException | DirectoryIteratorException x) {
            System.err.println(x);
        }
         return childPaths;
    }
    
    public static int nChildren(Path parentPath){
        return getChildNodePaths(parentPath).size();
    }
    
    private void setChildren(Path nodePath, CostOptionNode parent) throws ConfigurationException{
        CostOptionNodeFactory childFactory = new CostOptionNodeFactory(programSize);
        
        ArrayList<Path> childPaths = getChildNodePaths(nodePath);
        ArrayList<CostOptionNode> children = new ArrayList<>();
        for (Path path: childPaths) {
//                    System.out.println(path);
                    children.add(childFactory.makeCostOptionNode(path, parent));
        }
    }
    
    ArrayList<CostOption> getCostOptions(Path nodePath) throws ConfigurationException {
        CostOptionFactory optionFactory = new CostOptionFactory(programSize);
        ArrayList<CostOption> costOptions = new ArrayList<>();
        ArrayList<Path> costOptionPaths = getCostOptionPaths(nodePath);
        PropertiesConfiguration optionConfig;
        
             
        for (Path path: costOptionPaths) {
            optionConfig = new PropertiesConfiguration(path.toFile());
            costOptions.add(optionFactory.makeCostOption(optionConfig));
        }
        
        Collections.sort(costOptions);
        return costOptions;
    }
    
    private OptionSelectionInterface getOptionSelection(AbstractConfiguration config, Path nodePath){
        //later, can support choice and configuration of optionSelection type via the config file
        SingleOptionSelection optionSelection = new SingleOptionSelection(nOptions(nodePath), 0);
        return optionSelection;
    }
    public CostOptionNode makeCostOptionNode(Path nodePath, CostOptionNode parent ) throws ConfigurationException{
        CostOptionNode node;
        ArrayList<CostOption> costOptions = getCostOptions( nodePath);
        String name = nodePath.getFileName().toString();       
        
        assert !"".equals(name): "name string is empty or null!";
        
        File configFile;
        PropertiesConfiguration config; 
        OptionSelectionInterface optionSelection; 
        
        configFile= nodePath.resolve(Paths.get(NODE_CONFIG_FILENAME)).toFile();
        config= new PropertiesConfiguration(configFile);
        optionSelection = getOptionSelection(config, nodePath);
        
        String[] descriptions = config.getStringArray(NODE_DESCRIPTION_KEY);
        String description = reconstituteCommaConfigString(descriptions);
            
        node = new CostOptionNode(parent, costOptions, 
                                optionSelection, name, description);
             
        setChildren(nodePath,node);
             
        return node;
    }
    
}
