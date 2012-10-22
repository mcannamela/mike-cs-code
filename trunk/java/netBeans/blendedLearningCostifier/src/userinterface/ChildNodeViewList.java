/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import costoptiontree.CostOptionNode;
import java.awt.Component;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;

/**
 *
 * @author mcannamela
 */
public class ChildNodeViewList extends BaseComponentListWithDeleteButtons implements ActionListener{
    
    private CostOptionNode clickedNode;
    
    private JButton button_dummyCostChanged;
    public static final String ACTION_COST_CHANGED = ChildNodeView.ACTION_COST_CHANGED;
    
    private JButton button_dummyNodeClicked;
    public static final String ACTION_EXPAND_NODE = "expandNode";
    
    @Override
    protected void initComponents() {
        super.initComponents();
        
        addMouseListener(new java.awt.event.MouseAdapter() {
            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                clickEvent(evt);
            }
        });
        button_dummyCostChanged = DummyButtonFactory.makeDummyButton(ACTION_COST_CHANGED);
        button_dummyNodeClicked = DummyButtonFactory.makeDummyButton(ACTION_EXPAND_NODE);
    }
    
    public void addCostChangedListener(ActionListener listener){
        button_dummyCostChanged.addActionListener(listener);
    }
    public void addNodeExpandListener(ActionListener listener){
        button_dummyNodeClicked.addActionListener(listener);
    }
    public void addChildNodeSummary(CostOptionNode child){
        ChildNodeView view = new ChildNodeView();
        view.setChildNode(child);
        view.addCostChangedListener(this);
        addComponent(view);
    }

    public CostOptionNode getClickedNode() {
        return clickedNode;
    }
    
    @Override
    public void actionPerformed(ActionEvent evt) {
        super.actionPerformed(evt);
        if (evt.getActionCommand().equals(ACTION_COST_CHANGED)){
            System.out.println("Action in ChildNodeSummaryViewList: "+ACTION_COST_CHANGED);
            button_dummyCostChanged.doClick();
        }
    }
    
    private void clickEvent(java.awt.event.MouseEvent evt) {                                         
        
        Point p = evt.getPoint();
        
//        System.out.println("the point clicked is "+p);
        int cnt = 0;
        for (Component c: componentList){
            System.out.println(c.getParent().getBounds());
            if (c.getParent().getBounds().contains(p)){
//                System.out.println("this is the "+cnt+"th node");
                clickedNode = ((ChildNodeView) c).getNode();
            }
            cnt++;
        }
        if (clickedNode!=null){
        
            if (evt.getClickCount()==1){
                System.out.println("Mouse clicked on node "+ clickedNode.getName());
            }
            else if (evt.getClickCount()==2){
                System.out.println("Mouse double-clicked on node "+ clickedNode.getName());
                button_dummyNodeClicked.doClick();

            }
        }
        else{
            System.out.println("clicked node is null");
        }
    }

    void refresh() {
//        System.out.println("refreshing views");
        for (Component c : componentList){
            ((ChildNodeView) c).refresh();
        }
    }
}
