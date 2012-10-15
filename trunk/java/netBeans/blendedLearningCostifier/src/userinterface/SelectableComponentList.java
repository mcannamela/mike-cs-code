/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 *
 * @author Michael
 */
public class SelectableComponentList extends BaseComponentListWithDeleteButtons{
    
    private JButton button_dummySelectionChanged;
    private ButtonGroup selectionButtonGroup = new ButtonGroup();
    private int selection = 0;
    
    public static final String actionSelectionChanged = "selectionChanged";

    @Override
    protected void initComponents() {
        super.initComponents();
        button_dummySelectionChanged = new JButton();
        button_dummySelectionChanged.setVisible(false);
        button_dummySelectionChanged.setActionCommand(actionSelectionChanged);
    }

    
    @Override
    protected JPanel assembleComponentPanel(Component component, 
        JButton deleteButton, JPanel componentPanel) {
        JRadioButton selectionButton = new JRadioButton();
        Integer nButtons = selectionButtonGroup.getButtonCount();
        
        selectionButton.setActionCommand(nButtons.toString());
        selectionButton.addActionListener(this);
        selectionButtonGroup.add(selectionButton);
        
        componentPanel = super.assembleComponentPanel(component, 
                                    deleteButton, componentPanel);
        componentPanel.add(selectionButton,2);
//        componentPanel.add(Box.createHorizontalGlue(),3);        
        return componentPanel;
    }
    
    @Override
    public void actionPerformed(ActionEvent evt) {
        super.actionPerformed(evt);
        for (int i=0;i<selectionButtonGroup.getButtonCount();i++){
            if ( ((Integer)i).toString().equals(evt.getActionCommand()) ){
                selection = i;
                button_dummySelectionChanged.doClick();
                System.out.println("option "+i+" is selected");
                
            }
        }
    }
    
    public void addSelectionChangedListener(ActionListener listener){
        button_dummySelectionChanged.addActionListener(listener);
    }

    public int getSelection() {
        return selection;
    }
    
    
}
