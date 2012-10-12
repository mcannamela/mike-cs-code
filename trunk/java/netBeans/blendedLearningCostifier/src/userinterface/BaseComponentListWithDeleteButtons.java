/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;

/**
 *
 * @author Michael
 */
abstract public class BaseComponentListWithDeleteButtons extends JPanel{
    
    private ArrayList<Component> componentList = new ArrayList<>();
    private ArrayList<JButton> deleteButtonList = new ArrayList<>();
    
    private static final String actionDeleteComponent = "deleteComponent";
    private static final String deleteButtonText = "x";
    
    public BaseComponentListWithDeleteButtons() {
        initComponents();
    }
    
    
    
    private JPanel makeNewComponentPanel(){
        JPanel panel = new JPanel();
        BoxLayout layout = new BoxLayout(this, BoxLayout.X_AXIS);
               
        panel.setLayout(layout);
        
        
        return panel;
    }
    
    private javax.swing.JButton makeNewButton(){
        JButton button = new JButton();
        
        return button;
    }
    private javax.swing.JButton makeNewDeleteButton(){
        JButton button = makeNewButton();
        button.setText(deleteButtonText);
        button.setActionCommand(actionDeleteComponent);
        return button;
    }

    @SuppressWarnings("unchecked")
    private void initComponents() {
        BoxLayout layout = new BoxLayout(this, BoxLayout.Y_AXIS);
        this.setLayout(layout);
    }
    
    public void addComponent(Component component){
       JPanel componentPanel = makeNewComponentPanel();
       JButton button = makeNewButton();
       
       
       
       componentPanel.add(button);
       componentPanel.add(Box.createHorizontalGlue());
       componentPanel.add(component);
       
       deleteButtonList.add(button);
       componentList.add(component);
       
    }
    
    
    
    public void removeComponent(int index){
        componentList.remove(index);
        deleteButtonList.remove(index);
        remove(index);
    }
    public void removeComponent(Component component){
        int index = componentList.indexOf(component);
        assert index>=0:"component not found in componentList";
        removeComponent(index);
    }
    public void removeComponent(JButton pushButton){
        int index = deleteButtonList.indexOf(pushButton);
        assert index>=0:"button not found in pushButtonList";
        removeComponent(index);
    }
    
    public void actionPerformed(ActionEvent evt) {
        if (actionDeleteComponent.equals(evt.getActionCommand())) {
            JButton button = (JButton)evt.getSource();
            removeComponent(button);
        }
} 

        
}
