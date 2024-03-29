/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOption;
import costoptiontree.CostOptionNode;
import costoptiontree.SingleOptionSelection;
import java.awt.Dialog;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyVetoException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;

/**
 *
 * @author mcannamela
 */
public class CostOptionNodeInternalFrame extends javax.swing.JInternalFrame implements ActionListener {
    private CostOptionNode node;
    private CostOptionViewList costOptionViewList = new CostOptionViewList();
    private ChildNodeViewList childNodeViewList= new ChildNodeViewList();
    
    public static final String ACTION_LIST_SELECTION_CHANGED = SelectableComponentList.ACTION_SELECTION_CHANGED;
    
    private JButton button_dummySelectionChanged;
    public static final String ACTION_SELECTION_CHANGED = "nodeSelectionChanged";
    
    private JButton button_dummyCostChanged;
    public static final String ACTION_COST_CHANGED = "nodeCostChanged";
    /**
     * Creates new form CostOptionNodeInternalFrame
     */
    public CostOptionNodeInternalFrame() {
        initComponents();
        initComponentLists();
        setResizable(true);
    }
    public void addCostChangedListener(ActionListener listener){
        button_dummyCostChanged.addActionListener(listener);
    }
    public void addSelectionChangedListener(ActionListener listener){
        button_dummySelectionChanged.addActionListener(listener);
    }
    
    private void initComponentLists(){
        
        button_dummySelectionChanged = DummyButtonFactory.makeDummyButton(ACTION_SELECTION_CHANGED);
        button_dummyCostChanged = DummyButtonFactory.makeDummyButton(ACTION_COST_CHANGED);
        
        jScrollPane_optionList.setViewportView(costOptionViewList);
        jScrollPane_childChoices.setViewportView(childNodeViewList);
        
        costOptionViewList.addCostChangedListener(this);
        costOptionViewList.addSelectionChangedListener(this);
        childNodeViewList.addCostChangedListener(this);
        childNodeViewList.addNodeExpandListener(this);        
    }
    
    public void setNode(CostOptionNode node){
        this.node = node;
        jLabel_nodeName.setText(node.getName());
        jTextArea_nodeDescription.setText(node.getDescription());
        
        for (CostOption option: node.getCostOptions()){
            costOptionViewList.addOption(option);
        }
        
        for (CostOptionNode child : node.getChildren()){
            childNodeViewList.addChildNodeSummary(child);
        }
        
        int selection = ((SingleOptionSelection)node.getOptionSelection()).getSelection();
        costOptionViewList.select(selection);
        
        if (node.getTreeLevel()==0){
            setLocation(100,20);
        }
       
        displayCost();   
    }
    
    public void displayCost(){
        jLabel_scaledCost.setText(((Integer)node.cost()).toString());
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel_nodeName = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        jTextArea_nodeDescription = new javax.swing.JTextArea();
        jLabel_scaledCost = new javax.swing.JLabel();
        jLabel_scaledCostLabel = new javax.swing.JLabel();
        jScrollPane_optionList = new javax.swing.JScrollPane();
        jLabel_options = new javax.swing.JLabel();
        jLabel_childChoices = new javax.swing.JLabel();
        jScrollPane_childChoices = new javax.swing.JScrollPane();
        jButton_newChildChoice = new javax.swing.JButton();
        jButton_newOption = new javax.swing.JButton();
        jButton_OK = new javax.swing.JButton();

        jLabel_nodeName.setFont(new java.awt.Font("Tahoma", 0, 24)); // NOI18N
        jLabel_nodeName.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_nodeName.setText("nodeName");

        jTextArea_nodeDescription.setEditable(false);
        jTextArea_nodeDescription.setColumns(20);
        jTextArea_nodeDescription.setLineWrap(true);
        jTextArea_nodeDescription.setRows(5);
        jTextArea_nodeDescription.setWrapStyleWord(true);
        jScrollPane1.setViewportView(jTextArea_nodeDescription);

        jLabel_scaledCost.setFont(new java.awt.Font("Tahoma", 0, 24)); // NOI18N
        jLabel_scaledCost.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        jLabel_scaledCost.setText("cost");

        jLabel_scaledCostLabel.setFont(new java.awt.Font("Tahoma", 0, 24)); // NOI18N
        jLabel_scaledCostLabel.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_scaledCostLabel.setText("scaled cost, $");

        jLabel_options.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel_options.setText("OPTIONS");

        jLabel_childChoices.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jLabel_childChoices.setText("CHILD CHOICES");

        jButton_newChildChoice.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jButton_newChildChoice.setText("+");
        jButton_newChildChoice.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton_newChildChoiceActionPerformed(evt);
            }
        });

        jButton_newOption.setFont(new java.awt.Font("Tahoma", 0, 14)); // NOI18N
        jButton_newOption.setText("+");
        jButton_newOption.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton_newOptionActionPerformed(evt);
            }
        });

        jButton_OK.setText("OK");
        jButton_OK.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton_OKActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(20, 20, 20)
                .addComponent(jLabel_childChoices)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jButton_newChildChoice)
                .addContainerGap(460, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(20, 20, 20)
                                .addComponent(jLabel_options)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(jButton_newOption)
                                .addGap(0, 0, Short.MAX_VALUE))
                            .addGroup(layout.createSequentialGroup()
                                .addContainerGap()
                                .addComponent(jLabel_nodeName)
                                .addGap(48, 48, 48)
                                .addComponent(jButton_OK)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(jLabel_scaledCostLabel)))
                        .addGap(10, 10, 10)
                        .addComponent(jLabel_scaledCost))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jScrollPane_optionList, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jScrollPane1, javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jScrollPane_childChoices))))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(jLabel_nodeName)
                            .addComponent(jLabel_scaledCost)
                            .addComponent(jLabel_scaledCostLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 20, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(jButton_OK)
                        .addGap(18, 18, 18)))
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 61, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel_options)
                    .addComponent(jButton_newOption))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jScrollPane_optionList, javax.swing.GroupLayout.DEFAULT_SIZE, 303, Short.MAX_VALUE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel_childChoices)
                    .addComponent(jButton_newChildChoice))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane_childChoices, javax.swing.GroupLayout.PREFERRED_SIZE, 130, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jButton_newChildChoiceActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton_newChildChoiceActionPerformed
        System.out.println("would add a new child choice");
    }//GEN-LAST:event_jButton_newChildChoiceActionPerformed

    private void jButton_newOptionActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton_newOptionActionPerformed
        System.out.println("would add a new option");
    }//GEN-LAST:event_jButton_newOptionActionPerformed

    private void jButton_OKActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton_OKActionPerformed
        setVisible(false);
        dispose();
    }//GEN-LAST:event_jButton_OKActionPerformed

    @Override
    public void actionPerformed(ActionEvent evt) {
        String command = evt.getActionCommand();
        System.out.println("\nAction in CostOptionInternalFrame "+node.getName()+": "+command);
        
        if (CostOptionViewList.ACTION_COST_CHANGED.equals(command) || 
                ChildNodeViewList.ACTION_COST_CHANGED.equals(command)||
                CostOptionNodeInternalFrame.ACTION_COST_CHANGED.equals(command)){
//            System.out.println("Action in CostOptionNodeInternalFrame: "+command);
            displayCost();
            button_dummyCostChanged.doClick();
            
        }
        else if (ACTION_LIST_SELECTION_CHANGED.equals(command)){
            node.setOptionSelection(costOptionViewList.getOptionSelection());
            displayCost();
            button_dummySelectionChanged.doClick();
        }
        else if (ACTION_SELECTION_CHANGED.equals(command)){
            System.out.println("will now refresh child views");
            childNodeViewList.refresh();
            displayCost();
            button_dummyCostChanged.doClick();
        }
        else if (ChildNodeViewList.ACTION_EXPAND_NODE.equals(command)){
            CostOptionNodeInternalFrame frame = new CostOptionNodeInternalFrame();
            frame.setNode(childNodeViewList.getClickedNode());
            
            frame.addCostChangedListener(this);
            frame.addSelectionChangedListener(this);
            frame.setVisible(true);
            
            
            getParent().add(frame);
            frame.setLocation(getLocation().x+100, getLocation().y + 50);
            frame.moveToFront();
            try {
                frame.setSelected(true);
            } catch (PropertyVetoException ex) {
                Logger.getLogger(CostOptionNodeInternalFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    public void programSizeChanged(){
        childNodeViewList.refresh();
        costOptionViewList.refresh();
        displayCost();
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jButton_OK;
    private javax.swing.JButton jButton_newChildChoice;
    private javax.swing.JButton jButton_newOption;
    private javax.swing.JLabel jLabel_childChoices;
    private javax.swing.JLabel jLabel_nodeName;
    private javax.swing.JLabel jLabel_options;
    private javax.swing.JLabel jLabel_scaledCost;
    private javax.swing.JLabel jLabel_scaledCostLabel;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane_childChoices;
    private javax.swing.JScrollPane jScrollPane_optionList;
    private javax.swing.JTextArea jTextArea_nodeDescription;
    // End of variables declaration//GEN-END:variables
}
