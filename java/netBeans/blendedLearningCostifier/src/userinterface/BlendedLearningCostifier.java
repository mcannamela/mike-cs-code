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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JDesktopPane;
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
    
    private ProgramSize programSize;
    private JPanel jPanel_rootNodes = new JPanel();
    private JPanel jPanel_rightSide = new JPanel();
    
    private ArrayList<RootNodePanel> rootNodePanels= new ArrayList<>();
    private static final Path rootPath = Paths.get("C:\\Users\\Michael\\Dropbox\\timewise_blendedLearningEvaluator\\aBlendedLearningProgram");   
    private CostOptionNode rootNode;
    
    public BlendedLearningCostifier() throws ConfigurationException  {
        super("BlendedLearningCostifier");
        
        PropertiesConfiguration sizeConfig = new PropertiesConfiguration(rootPath.resolve("programSize.config").toFile());
        programSize = new ProgramSizeFactory().makeProgramSize(sizeConfig);
        
        
        CostOptionNodeFactory factory = new CostOptionNodeFactory(programSize);
        rootNode = factory.makeCostOptionNode(rootPath, null);
        
        programSizePanel.setProgramSize(programSize);
        programSizePanel.addProgramSizeChangedListener(this);
        
        
        //Make the big window be indented 50 pixels from each edge
        //of the screen.
        int inset = 50;
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        setBounds(2*inset, 0,
                  screenSize.width  - inset*3,
                  screenSize.height - inset*0);

        desktop = new JDesktopPane(); //a specialized layered pane
       
        vSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, jPanel_rootNodes, desktop);
        
        jPanel_rightSide.setLayout(new javax.swing.BoxLayout(jPanel_rightSide, javax.swing.BoxLayout.Y_AXIS));
        jPanel_rightSide.add(programSizePanel);
        hSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, vSplitPane, jPanel_rightSide);
        hSplitPane.setResizeWeight(1.0);
        
        createRootNodePanels(); //create first "window"
        setContentPane(hSplitPane);
        
        
        setJMenuBar(createMenuBar());

        //Make dragging a little faster but perhaps uglier.
        desktop.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);
    }

    protected JMenuBar createMenuBar() {
        JMenuBar menuBar = new JMenuBar();

        //Set up the lone menu.
        JMenu menu = new JMenu("File");
        menu.setMnemonic(KeyEvent.VK_D);
        menuBar.add(menu);

        //Set up the first menu item.
        JMenuItem menuItem = new JMenuItem("New");
        menuItem.setMnemonic(KeyEvent.VK_N);
        menuItem.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_N, ActionEvent.ALT_MASK));
        menuItem.setActionCommand("new");
        menuItem.addActionListener(this);
        menu.add(menuItem);

        //Set up the second menu item.
        menuItem = new JMenuItem("Quit");
        menuItem.setMnemonic(KeyEvent.VK_Q);
        menuItem.setAccelerator(KeyStroke.getKeyStroke(
                KeyEvent.VK_Q, ActionEvent.ALT_MASK));
        menuItem.setActionCommand("quit");
        menuItem.addActionListener(this);
        menu.add(menuItem);

        return menuBar;
    }

    //React to menu selections.
    public void actionPerformed(ActionEvent e) {
        String command = e.getActionCommand();
        if ("new".equals(command)) { //new
            createRootNodePanels();
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
        else if ("quit".equals(command)) { //quit
            quit();
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
