/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOptionNode;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;

/**
 *
 * @author mcannamela
 */
public class ChildNodeSummaryViewList extends BaseComponentListWithDeleteButtons{
    
    private JButton button_dummyCostChanged;
    public static final String ACTION_COST_CHANGED = ChildNodeSummaryView.ACTION_COST_CHANGED;
    
    @Override
    protected void initComponents() {
        super.initComponents();
        button_dummyCostChanged = new JButton();
        button_dummyCostChanged.setVisible(false);
        button_dummyCostChanged.setActionCommand(ACTION_COST_CHANGED);
    }
    
    public void addCostChangedListener(ActionListener listener){
        button_dummyCostChanged.addActionListener(listener);
    }
    public void addChildNodeSummary(CostOptionNode child){
        ChildNodeSummaryView view = new ChildNodeSummaryView();
        view.setChildNode(child);
        addComponent(view);
    }
    
    @Override
    public void actionPerformed(ActionEvent evt) {
        super.actionPerformed(evt);
        if (evt.getActionCommand().equals(ACTION_COST_CHANGED)){
            System.out.println("Action in ChildNodeSummaryViewList: "+ACTION_COST_CHANGED);
            button_dummyCostChanged.doClick();
        }
    }
}
