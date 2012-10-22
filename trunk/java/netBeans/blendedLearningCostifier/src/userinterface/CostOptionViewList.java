/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOption;
import costoptiontree.SingleOptionSelection;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;

/**
 *
 * @author Michael
 */
public class CostOptionViewList extends SelectableComponentList implements ActionListener{
    private JButton button_dummyCostChanged;
    public static final String ACTION_COST_CHANGED = CostOptionView.ACTION_COST_CHANGED;
    
    @Override
    protected void initComponents() {
        super.initComponents();
        button_dummyCostChanged = new JButton();
        button_dummyCostChanged.setVisible(false);
        button_dummyCostChanged.setActionCommand(CostOptionViewList.ACTION_COST_CHANGED);
    }
    
    public void addCostChangedListener(ActionListener listener){
        button_dummyCostChanged.addActionListener(listener);
    }
    public void addOption(CostOption option){
        CostOptionView view = new CostOptionView();
        view.setOption(option);
        view.addCostChangedListener(this);
        addComponent(view);
    }
    
    public SingleOptionSelection getOptionSelection(){
            int nOptions = nComponents();
            SingleOptionSelection optionSelection;
            optionSelection= new SingleOptionSelection(nOptions, selection);
            return optionSelection;
    }

    @Override
    public void actionPerformed(ActionEvent evt) {
        super.actionPerformed(evt);
        if (evt.getActionCommand().contains(SELECTION_ACTION_PREFIX)){
         
        }
        if (evt.getActionCommand().equals(CostOptionViewList.ACTION_COST_CHANGED)){
            System.out.println("Action in CostOptionViewList: "+CostOptionViewList.ACTION_COST_CHANGED);
            button_dummyCostChanged.doClick();
        }
    }
    
    public void refresh(){
        for (Component view : componentList){
            ((CostOptionView) view).refresh();
        }
    }


    
    
    
}
