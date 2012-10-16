/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOptionNode;

/**
 *
 * @author mcannamela
 */
public class ChildNodeSummaryViewList extends BaseComponentListWithDeleteButtons{
    public void addChildNodeSummary(CostOptionNode child){
        ChildNodeSummaryView view = new ChildNodeSummaryView();
        view.setChildNode(child);
        addComponent(view);
    }
}
