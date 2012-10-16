/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOption;

/**
 *
 * @author Michael
 */
public class CostOptionViewList extends SelectableComponentList{
    public void addOption(CostOption option){
        CostOptionView view = new CostOptionView();
        view.setOption(option);
        addComponent(view);
    }
}
