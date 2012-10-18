/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOption;
import costoptiontree.CostOptionNode;
import costoptiontree.OptionSelectionInterface;
import java.awt.Container;
import java.awt.Dialog;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JRootPane;
import javax.swing.JSlider;

/**
 *
 * @author Michael
 */
public class ChildNodeView extends javax.swing.JPanel implements ActionListener {
    
    private CostOptionNode node;
    
    private JButton button_dummyCostChanged= new JButton();
    private JButton button_dummyNodeDoubleClicked= new JButton();
    
    public static final String ACTION_COST_CHANGED = "costChanged";
    public static final String ACTION_NODE_DOUBLE_CLICKED = "nodeDoubleClicked";
    /**
     * Creates new form ChildNodeView
     */
    public ChildNodeView() {
        initComponents();
        
        button_dummyCostChanged.setVisible(false);
        button_dummyCostChanged.setActionCommand(ACTION_COST_CHANGED);
        button_dummyNodeDoubleClicked.setVisible(false);
        button_dummyNodeDoubleClicked.setActionCommand(ACTION_NODE_DOUBLE_CLICKED);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel_name = new javax.swing.JLabel();
        jLabel_selectedOptionLabel = new javax.swing.JLabel();
        jSlider_selectedOptionCost = new javax.swing.JSlider();
        jLabel_scaledMin = new javax.swing.JLabel();
        jLabel_scaledMax = new javax.swing.JLabel();
        jLabel_scaledSelected = new javax.swing.JLabel();
        jLabel_scaled = new javax.swing.JLabel();

        jLabel_name.setFont(new java.awt.Font("Tahoma", 0, 18)); // NOI18N
        jLabel_name.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_name.setText("childName");

        jLabel_selectedOptionLabel.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel_selectedOptionLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_selectedOptionLabel.setText("selectedOptionLabel");

        jSlider_selectedOptionCost.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSlider_selectedOptionCostStateChanged(evt);
            }
        });

        jLabel_scaledMin.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_scaledMin.setText("scaledMin");

        jLabel_scaledMax.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_scaledMax.setText("scaledMax");

        jLabel_scaledSelected.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_scaledSelected.setText("scaledSelected");

        jLabel_scaled.setFont(new java.awt.Font("Tahoma", 0, 12)); // NOI18N
        jLabel_scaled.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_scaled.setText("$,scaled");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel_name)
                .addGap(18, 18, 18)
                .addComponent(jLabel_selectedOptionLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel_scaledMin)
                        .addGap(18, 18, Short.MAX_VALUE)
                        .addComponent(jLabel_scaledSelected)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(jLabel_scaledMax))
                    .addComponent(jSlider_selectedOptionCost, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 10, Short.MAX_VALUE)
                .addComponent(jLabel_scaled)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel_scaledMax)
                    .addComponent(jLabel_scaledSelected)
                    .addComponent(jLabel_scaled)
                    .addComponent(jLabel_scaledMin))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jSlider_selectedOptionCost, javax.swing.GroupLayout.PREFERRED_SIZE, 26, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(jLabel_selectedOptionLabel)
                        .addComponent(jLabel_name)))
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    public void addCostChangedListener(ActionListener listener){
        button_dummyCostChanged.addActionListener(listener);
    }
    void addNodeDoubleClickedListener(ActionListener listener) {
        button_dummyNodeDoubleClicked.addActionListener(listener);
    }
    
    private void jSlider_selectedOptionCostStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSlider_selectedOptionCostStateChanged
        CostOption option = getSelectedOption();
        int value = jSlider_selectedOptionCost.getValue();

        JSlider source = (JSlider)evt.getSource();
        if (!source.getValueIsAdjusting()) {
            int selectedCost = (int)source.getValue();
            option.setSelectedCost(selectedCost);
            jLabel_scaledSelected.setText(new Integer(option.getScaledCost()).toString());
            button_dummyCostChanged.doClick();
        }
        
    }//GEN-LAST:event_jSlider_selectedOptionCostStateChanged
    
    @Override
    public void actionPerformed(ActionEvent e) {
        if ("delete".equals(e.getActionCommand())) {
            delete();
        } 
    }
    
    public void delete(){
        node.getParent().removeChild(node);
        this.getParent().remove(this);
    }
    
    public CostOption getSelectedOption(){
        OptionSelectionInterface optionSelection = node.getOptionSelection();
        int selection = optionSelection.getSummaryOptionIndex();
        return node.getCostOption(selection);
    }
    public void setChildNode(CostOptionNode node) {
        this.node = node;
        jLabel_name.setText(node.getName());
        selectedOptionChanged();
        
    }

    public CostOptionNode getNode() {
        return node;
    }
    
    public void selectedOptionChanged(){
        CostOption option = getSelectedOption();
        jLabel_selectedOptionLabel.setText(option.getDescription());
        
        jLabel_scaledMax.setText(new Integer(option.getScaledMaxCost()).toString());
        jLabel_scaledMin.setText(new Integer(option.getScaledMinCost()).toString());
        
        jSlider_selectedOptionCost.setMaximum(option.getMaxCost());
        jSlider_selectedOptionCost.setMinimum(option.getMinCost());
        
        jSlider_selectedOptionCost.setValue(option.getSelectedCost());
        
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel_name;
    private javax.swing.JLabel jLabel_scaled;
    private javax.swing.JLabel jLabel_scaledMax;
    private javax.swing.JLabel jLabel_scaledMin;
    private javax.swing.JLabel jLabel_scaledSelected;
    private javax.swing.JLabel jLabel_selectedOptionLabel;
    private javax.swing.JSlider jSlider_selectedOptionCost;
    // End of variables declaration//GEN-END:variables

    void refresh() {
        System.out.println("    refreshing view");
        selectedOptionChanged();
    }

    
}
