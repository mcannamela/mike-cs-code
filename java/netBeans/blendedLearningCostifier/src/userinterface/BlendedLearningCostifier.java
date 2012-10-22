/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import blendedlearningprogram.ProgramSize;
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
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import org.apache.commons.configuration.ConfigurationException;

/**
 *
 * @author mcannamela
 */
public class BlendedLearningCostifier extends JFrame
                               implements ActionListener {
    private JDesktopPane desktop;
    private JSplitPane hSplitPane, vSplitPane;
    private ProgramSizePanel programSizePanel = new ProgramSizePanel();;
    private JPanel rootPanelArea = new JPanel();
    
    private static final Path rootPath = Paths.get("C:\\Users\\Michael\\Dropbox\\timewise_blendedLearningEvaluator\\aBlendedLearningProgram");   
    private CostOptionNode rootNode;

    public BlendedLearningCostifier() throws ConfigurationException {
        super("BlendedLearningCostifier");
        
        CostOptionNodeFactory factory = new CostOptionNodeFactory(new ProgramSize());
        rootNode = factory.makeCostOptionNode(rootPath, null);
        
        //Make the big window be indented 50 pixels from each edge
        //of the screen.
        int inset = 50;
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        setBounds(2*inset, 0,
                  screenSize.width  - inset*3,
                  screenSize.height - inset*0);

        //Set up the GUI.
        desktop = new JDesktopPane(); //a specialized layered pane
//        desktop.setMinimumSize(new Dimension(900, 800));
//        desktop.setLayout(new javax.swing.BoxLayout(desktop, javax.swing.BoxLayout.LINE_AXIS));
//        desktop.setLayout(new FlowLayout());
        
        vSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, rootPanelArea, desktop);
//        vSplitPane.setResizeWeight(1.0);
        
//        programSizePanel.setMaximumSize(programSizePanel.getPreferredSize());
        hSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, vSplitPane, programSizePanel);
        hSplitPane.setResizeWeight(1.0);
        
        createFrame(); //create first "window"
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
        if ("new".equals(e.getActionCommand())) { //new
            createFrame();
        } else { //quit
            quit();
        }
    }

    //Create a new internal frame.
    protected void createFrame() {
//        CostOptionNodeInternalFrame frame = new CostOptionNodeInternalFrame();
        for (CostOptionNode node : rootNode.getChildren()){
            RootNodePanel panel = new RootNodePanel();
            panel.setNode(node);
            panel.setDesktop(desktop);
    //        panel.setVisible(true);
    //        panel.setLocation(0,0);
            rootPanelArea.add(panel); 
        }
        
        
        
//        frame.setNode(rootNode);
//        frame.setVisible(true); //necessary as of 1.3
//        frame.setLocation(0,0);
//        addFrame(frame);
    }
    
    public void addFrame(CostOptionNodeInternalFrame frame){
        desktop.add(frame);
        try {
            frame.setSelected(true);
        } catch (java.beans.PropertyVetoException e) {}
    }

    //Quit the application.
    protected void quit() {
        System.exit(0);
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
