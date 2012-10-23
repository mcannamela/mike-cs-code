/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package costoptiontree;
import blendedlearningprogram.ProgramSize;
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
    public int nOptions(Path costOptionsPath){
        assert Files.exists(costOptionsPath) && !Files.notExists(costOptionsPath): 
                "file DNE or can't be accessed:"+costOptionsPath;
        int nOptions = 0;
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(costOptionsPath, "*.{config}")) {
            //iterate over all config files in the nodeCostOptions directory
            for (Path file: stream) {
                nOptions++;
            }
        } catch (IOException | DirectoryIteratorException x) {
            System.err.println(x);
        }
        return nOptions;
    }
    
    public int nChildren(Path childDirectory){
        int nChildren = 0;
        
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(childDirectory)) {
            for (Path file: stream) {
                assert Files.isDirectory(file): "child directory must contain "+
                        "only costOptionNode directories:"+file;
                nChildren++;
            }
        } catch (IOException | DirectoryIteratorException x) {
            System.err.println(x);
        }
        return nChildren;
    }
    private void setChildren(Path nodePath, CostOptionNode parent) throws ConfigurationException{
        Path childDirectory = Paths.get(nodePath.toString(), "children");
        assert Files.isDirectory(childDirectory): "nodePath must have subdirectory called children:"+nodePath;
        CostOptionNodeFactory childFactory = new CostOptionNodeFactory(programSize);
        ArrayList<CostOptionNode> children = new ArrayList<>();
        if (nChildren(childDirectory)>0){
            try (DirectoryStream<Path> stream = Files.newDirectoryStream(childDirectory)) {
                for (Path file: stream) {
                    System.out.println(file);
                    children.add(childFactory.makeCostOptionNode(file, parent));
                }
            } catch (IOException | DirectoryIteratorException x) {
                System.err.println(x);
            }
        }
        
//        for (CostOptionNode child : children){
//            parent.addChild(child);
//        }
        
    }
    private ArrayList<CostOption> getCostOptions(Path costOptionsPath) throws ConfigurationException {
        assert Files.exists(costOptionsPath) && !Files.notExists(costOptionsPath): 
                "file DNE or can't be accessed:"+costOptionsPath;
        CostOptionFactory optionFactory = new CostOptionFactory(programSize);
        ArrayList<CostOption> costOptions = new ArrayList<>();
        PropertiesConfiguration optionConfig;
        //open the sub-directory that contains the costOptions config files
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(costOptionsPath)) {
            //iterate over all config files in the nodeCostOptions directory
            for (Path file: stream) {
                //open the option config file and add it to the list of costOptions for this node
                optionConfig = new PropertiesConfiguration(file.toFile());
                costOptions.add(optionFactory.makeCostOption(optionConfig));
            }
        } catch (IOException | DirectoryIteratorException x) {
            System.err.println(x);
        }
//        System.out.println("Before sorting:");
//        System.out.println(costOptions);
        Collections.sort(costOptions);
//        System.out.println("After sorting:");
//        System.out.println(costOptions);
        return costOptions;
    }
    private OptionSelectionInterface getOptionSelection(AbstractConfiguration config, Path costOptionsPath){
        //later, can support choice and configuration of optionSelection type via the config file
        SingleOptionSelection optionSelection = new SingleOptionSelection(nOptions(costOptionsPath), 0);
        return optionSelection;
    }
    public CostOptionNode makeCostOptionNode(Path nodePath, CostOptionNode parent ) throws ConfigurationException{
        
        String name = nodePath.getFileName().toString();       
        
        assert !"".equals(name): "name string is empty or null!";
        Path costOptionsPath = Paths.get(nodePath.toString(), "nodeCostOptions");
        ArrayList<CostOption> costOptions = getCostOptions( costOptionsPath);
        
        
        //open the config file for this node
        CostOptionNode node;
        PropertiesConfiguration config = 
                    new PropertiesConfiguration(nodePath.resolve(Paths.get("node.config")).toFile());

        OptionSelectionInterface optionSelection = 
                   getOptionSelection(config, costOptionsPath);

        
        String[] descriptions = config.getStringArray("description");
        
        String description = reconstituteCommaConfigString(descriptions);
            
        node = new CostOptionNode(parent, costOptions, 
                                optionSelection, name, description);
             
        if (Files.exists(nodePath.resolve("children"))){
                 setChildren(nodePath,node);
             }

               
        return node;
    }
    
}
