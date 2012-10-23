/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import blendedlearningprogram.ProgramSize;
import blendedlearningprogram.ProgramSizeFactory;
import costoptiontree.CostOptionNode;
import costoptiontree.CostOptionNodeFactory;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JDesktopPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.PropertiesConfiguration;

/**
 *
 * @author mcannamela
 */
public class BlendedLearningCostifier extends JFrame
                               implements ActionListener {
    private JDesktopPane desktop;
    private JSplitPane hSplitPane, vSplitPane;
    private ProgramSizePanel programSizePanel = new ProgramSizePanel();
    private CostSummaryPanel costSummaryPanel = new CostSummaryPanel();
    
    private ProgramSize programSize;
    private JPanel jPanel_rootNodes = new JPanel();
    private JPanel jPanel_rightSide = new JPanel();
    private JMenu jMenu = new JMenu("File");
    
    private ArrayList<RootNodePanel> rootNodePanels= new ArrayList<>();
    private Path rootPath = Paths.get("C:\\Users\\Michael\\Dropbox\\timewise_blendedLearningEvaluator\\aBlendedLearningProgram");   
    private CostOptionNode rootNode;
    
    private static final String ACTION_OPEN = "open";
    private static final String ACTION_CLOSE = "quit";
    private static final String ACTION_DEFAULT_DEMO = "demo";
    
    private JButton dButton_close = new JButton();
    
    private boolean hasOpenConfiguration = false;
    
    public BlendedLearningCostifier() throws ConfigurationException  {
        super("BlendedLearningCostifier");
        
        
        dButton_close.setVisible(false);
        dButton_close.setActionCommand(ACTION_CLOSE);
        dButton_close.addActionListener(this);
        
        
        int inset = 50;
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        setBounds(2*inset, 0,
                  screenSize.width  - inset*3,
                  screenSize.height - inset*0);

        desktop = new JDesktopPane();
       
        vSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, jPanel_rootNodes, desktop);
        
        jPanel_rightSide.setLayout(new javax.swing.BoxLayout(jPanel_rightSide, javax.swing.BoxLayout.Y_AXIS));
        jPanel_rightSide.add(programSizePanel);
        jPanel_rightSide.add(costSummaryPanel);
        
        hSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, vSplitPane, jPanel_rightSide);
        hSplitPane.setResizeWeight(1.0);
        
        setContentPane(hSplitPane);
        setJMenuBar(createMenuBar());
        
        desktop.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);
    }
    
    protected JMenuBar createMenuBar() {
        JMenuBar menuBar = new JMenuBar();

        //Set up the lone menu.
        jMenu.setMnemonic(KeyEvent.VK_F);
        menuBar.add(jMenu);

        //Set up the first menu item.
        JMenuItem menuItem = new JMenuItem("Open Configuration");
        menuItem.setMnemonic(KeyEvent.VK_N);
        menuItem.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_N, ActionEvent.ALT_MASK));
        menuItem.setActionCommand(ACTION_OPEN);
        menuItem.addActionListener(this);
        jMenu.add(menuItem);
        
        
        menuItem = new JMenuItem("Close Configuration");
        menuItem.setMnemonic(KeyEvent.VK_L);
        menuItem.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_L, ActionEvent.ALT_MASK));
        menuItem.setActionCommand(ACTION_CLOSE);
        menuItem.addActionListener(this);
        jMenu.add(menuItem);
        
        
        menuItem = new JMenuItem("Show Demo Configuration");
        menuItem.setMnemonic(KeyEvent.VK_D);
        menuItem.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_D, ActionEvent.ALT_MASK));
        menuItem.setActionCommand(ACTION_DEFAULT_DEMO);
        menuItem.addActionListener(this);
        jMenu.add(menuItem);

        return menuBar;
    }

    public void readProgramSize() throws ConfigurationException {
        PropertiesConfiguration sizeConfig = new PropertiesConfiguration(rootPath.resolve("programSize.config").toFile());
        programSize = new ProgramSizeFactory().makeProgramSize(sizeConfig);
        
        programSizePanel.setProgramSize(programSize);
        programSizePanel.addProgramSizeChangedListener(this);
    }
    
    public void readRootNode() throws ConfigurationException{
        CostOptionNodeFactory factory = new CostOptionNodeFactory(programSize);
        rootNode = factory.makeCostOptionNode(rootPath, null);
    }
    
    public void openNodeDirectory() throws ConfigurationException{
        if (hasOpenConfiguration){
            dButton_close.doClick();            
        }
        readProgramSize();
        readRootNode();
        createRootNodePanels();
        costSummaryPanel.displayCost(rootNode.cost());
        validate();
        vSplitPane.resetToPreferredSizes();
        hasOpenConfiguration = true;
    }
    
    
    //React to menu selections.
    public void actionPerformed(ActionEvent e) {
        String command = e.getActionCommand();
        if (ACTION_OPEN.equals(command)) { //new
            JFileChooser chooser = new JFileChooser();
            File file;
            chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            int returnVal = chooser.showOpenDialog(this);
            if(returnVal == JFileChooser.APPROVE_OPTION) {
               file = chooser.getSelectedFile();
               System.out.println("You chose to open this directory: " +
                    file.getName());
               rootPath = file.toPath();
               try {
                openNodeDirectory();
                } catch (ConfigurationException ex) {
                    Logger.getLogger(BlendedLearningCostifier.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
        }
        else if (ACTION_DEFAULT_DEMO.equals(command)){
            try {
                openNodeDirectory();
            } catch (ConfigurationException ex) {
                Logger.getLogger(BlendedLearningCostifier.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        else if (ProgramSizePanel.ACTION_PROGRAM_SIZE_CHANGED.equals(command)){
            System.out.println("rootNodePanels heard that program size changed");
            for (RootNodePanel panel : rootNodePanels){
                panel.displayCost();
            }
            for (JInternalFrame frame : desktop.getAllFrames()){
                ((CostOptionNodeInternalFrame)frame).programSizeChanged();
            }
        }
        else if (RootNodePanel.ACTION_COST_CHANGED.equals(command)){
            costSummaryPanel.displayCost(rootNode.cost());
        }
        else if (ACTION_CLOSE.equals(command)) { //quit
            System.out.println("closing...");
            for(JInternalFrame frame:desktop.getAllFrames()){
                frame.setVisible(false);
                frame.dispose();
            }
            jPanel_rootNodes.removeAll();
            rootNodePanels.clear();
            
            validate();
            jPanel_rootNodes.validate();
            vSplitPane.resetToPreferredSizes();
            hasOpenConfiguration = false;
        }
    }

    //Create a new internal frame.
    protected void createRootNodePanels() {
        for (CostOptionNode node : rootNode.getChildren()){
            RootNodePanel panel = new RootNodePanel();
            panel.setNode(node);
            panel.setDesktop(desktop);
            jPanel_rootNodes.add(panel); 
            rootNodePanels.add(panel);
            panel.addCostChangedListener(this);
//            System.out.println(node);
        }
    }
    

    //Quit the application.
    protected void quit() {
        System.exit(0);
    }
    private int nFrames(){
        return desktop.getAllFrames().length;
    }

    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    private static void createAndShowGUI() throws ConfigurationException {
        //Make sure we have nice window decorations.
//        JFrame.setDefaultLookAndFeelDecorated(true);

        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(ChildNodeView.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //Create and set up the window.
        BlendedLearningCostifier frame = new BlendedLearningCostifier();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Display the window.
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                try {
                    createAndShowGUI();
                } catch (ConfigurationException ex) {
                    Logger.getLogger(BlendedLearningCostifier.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        });
    }
}
