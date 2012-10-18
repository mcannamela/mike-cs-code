/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 *
 * @author Michael
 */
public class SelectableComponentList extends BaseComponentListWithDeleteButtons{
    
    
    private ButtonGroup selectionButtonGroup = new ButtonGroup();
    private ArrayList<JRadioButton> selectionButtonList = new ArrayList<>();
    protected int selection = 0;
    
    private JButton button_dummySelectionChanged;
    public static final String ACTION_SELECTION_CHANGED = "SelectableComponentList:selectionChanged";
    public static final String SELECTION_ACTION_PREFIX = "selection";

    @Override
    protected void initComponents() {
        super.initComponents();
        button_dummySelectionChanged = new JButton();
        button_dummySelectionChanged.setVisible(false);
        button_dummySelectionChanged.setActionCommand(ACTION_SELECTION_CHANGED);
    }

    
    @Override
    protected JPanel assembleComponentPanel(Component component, 
        JButton deleteButton, JPanel componentPanel) {
        JRadioButton selectionButton = new JRadioButton();
        Integer nButtons = selectionButtonGroup.getButtonCount();
        
        selectionButton.setActionCommand(makeActionCommand(nButtons));
        selectionButton.addActionListener(this);
        selectionButtonGroup.add(selectionButton);
        selectionButtonList.add(selectionButton);
        
        componentPanel = super.assembleComponentPanel(component, 
                                    deleteButton, componentPanel);
        componentPanel.add(selectionButton,2);
//        componentPanel.add(Box.createHorizontalGlue(),3);        
        
//        if (nComponents()==1){
//            System.out.println("first component added, it will be selected");
//            selectionButton.setSelected(true);
//        }
        System.out.println("Selectable: assembled component number "+nComponents());
        return componentPanel;
    }
    
    @Override
    public void removeComponent(int index){
        JRadioButton selectionButton = selectionButtonList.get(index);
        selectionButtonList.remove(selectionButton);
        selectionButtonGroup.remove(selectionButton);
        super.removeComponent(index);
        
    }
    
    protected String makeActionCommand(int i){
        return SELECTION_ACTION_PREFIX+"_"+i;
    }
    
    @Override
    public void actionPerformed(ActionEvent evt) {
        super.actionPerformed(evt);
        assert selectionButtonGroup.getButtonCount()==nComponents(): "should be as many buttons as components";
        for (int i=0;i<nComponents();i++){
            if ( makeActionCommand(i).equals(evt.getActionCommand()) ){
                selection = i;
                button_dummySelectionChanged.doClick();
//                System.out.println("option "+i+" is selected");
            }
        }
    }
    
    public void addSelectionChangedListener(ActionListener listener){
        button_dummySelectionChanged.addActionListener(listener);
    }

    public int getSelection() {
        return selection;
    }
    
    public void select(int i){
        JRadioButton button = selectionButtonList.get(i);
        button.doClick();
        
    }
    
    
}
