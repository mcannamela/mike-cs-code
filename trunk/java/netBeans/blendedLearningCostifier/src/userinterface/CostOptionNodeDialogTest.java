/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import blendedlearningprogram.ProgramSize;
import costoptiontree.CostOptionNode;
import costoptiontree.CostOptionNodeFactory;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.nio.file.Path;
import java.nio.file.Paths;
import javax.swing.JButton;
import org.apache.commons.configuration.ConfigurationException;

/**
 *
 * @author mcannamela
 */
public class CostOptionNodeDialogTest extends javax.swing.JFrame implements ActionListener{
    private static final Path rootPath = Paths.get("C:\\Users\\Michael\\Dropbox\\timewise_blendedLearningEvaluator\\testConfigurationPath");   
    private CostOptionNode rootNode;
    private JButton showButton = new JButton();
    private static final String showDialogAction = "show";
    
    public CostOptionNodeDialogTest()  throws ConfigurationException {
        CostOptionNodeFactory factory = new CostOptionNodeFactory(new ProgramSize());
        rootNode = factory.makeCostOptionNode(rootPath, null);
          
        initComponents();
    }
    
    @SuppressWarnings("unchecked")
    private void initComponents() {
        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        getContentPane().setLayout(new javax.swing.BoxLayout(getContentPane(), javax.swing.BoxLayout.Y_AXIS));
        
        showButton.setText("show me");
        showButton.setActionCommand(showDialogAction);
        
        getContentPane().add(showButton);
        
        showButton.addActionListener(this);
        
        setLocation(500,500);

        pack();
    }
    
    
    @Override
    public void actionPerformed(ActionEvent evt) {
        String command = evt.getActionCommand();
        if (showDialogAction.equals(command)) {
            CostOptionNodeDialog dialog = new CostOptionNodeDialog(this, false);
            dialog.setNode(rootNode);
            dialog.setVisible(true);
        }
        
    }
    public static void main(String args[]) throws ConfigurationException {
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
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                try{
                    new CostOptionNodeDialogTest().setVisible(true);
                }
                catch (ConfigurationException e){
                    System.out.println(e);
                }
            }
        });
        
    }
}
