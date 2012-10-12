/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import blendedlearningprogram.ProgramSize;
import costoptiontree.CostOptionNode;
import costoptiontree.CostOptionNodeFactory;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import javax.swing.Box;
import org.apache.commons.configuration.*;

/**
 *
 * @author Michael
 */
public class ChildNodeSummaryTest extends javax.swing.JFrame{
    private ArrayList<ChildNodeSummaryView> childNodeSummaryViews = new ArrayList<>();
    private int nOptionViews;
    
    private CostOptionNode rootNode;
    private ProgramSize programSize = new ProgramSize();
    private static final Path rootPath = Paths.get("C:\\Users\\Michael\\Dropbox\\timewise_blendedLearningEvaluator\\testConfigurationPath");   
    
    public ChildNodeSummaryTest() throws ConfigurationException{
        
        CostOptionNodeFactory factory = new CostOptionNodeFactory(programSize);
        
        
        rootNode = factory.makeCostOptionNode(rootPath, null);
        
        
        
        ChildNodeSummaryView view;
        for (CostOptionNode node: rootNode.getChildren()){
            view = new ChildNodeSummaryView();
            childNodeSummaryViews.add(view);
            view.setChildNode(node);
        }
        
        initComponents();
    }

    
    @SuppressWarnings("unchecked")
    private void initComponents() {
        
        
        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        getContentPane().setLayout(new javax.swing.BoxLayout(getContentPane(), javax.swing.BoxLayout.Y_AXIS));
        
        for(ChildNodeSummaryView view : childNodeSummaryViews){
            getContentPane().add(view);
            getContentPane().add(Box.createVerticalGlue());
        }

        pack();
    }

    public static void main(String args[]) {
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
            java.util.logging.Logger.getLogger(ChildNodeSummaryView.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                try{
                    new ChildNodeSummaryTest().setVisible(true);
                }
                catch (ConfigurationException e){
                    System.out.println(e);
                }
            }
        });
    }    
}
