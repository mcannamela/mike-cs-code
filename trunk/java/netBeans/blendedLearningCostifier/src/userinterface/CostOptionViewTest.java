/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOption;
import costoptiontree.CostOptionFactory;
import java.util.ArrayList;
import javax.swing.Box;

/**
 *
 * @author Michael
 */
public class CostOptionViewTest extends javax.swing.JFrame{
    
    private ArrayList<CostOptionView> costOptionViews = new ArrayList<>();
    private int nOptionViews;
    
    private CostOption option= (new CostOptionFactory()).makeCostOption();
    
    public CostOptionViewTest() {
        
        nOptionViews = 3;
        CostOptionView view;
        for (int i=0;i<nOptionViews;i++){
            System.out.println("creating option view");
            view = new CostOptionView();
            costOptionViews.add(view);
        }
        costOptionViews.get(0).setOption(option);
        initComponents();
    }

    
    @SuppressWarnings("unchecked")
    private void initComponents() {
        
        
        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        getContentPane().setLayout(new javax.swing.BoxLayout(getContentPane(), javax.swing.BoxLayout.Y_AXIS));
        
        for(CostOptionView view: costOptionViews){
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
            java.util.logging.Logger.getLogger(CostOptionViewTest.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        java.awt.EventQueue.invokeLater(new Runnable() {
            @Override
            public void run() {
                new CostOptionViewTest().setVisible(true);
            }
        });
    }

}
