/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;

/**
 *
 * @author Michael
 */
public class BaseComponentListWithDeleteButtons extends JPanel implements ActionListener{
    
    protected ArrayList<Component> componentList = new ArrayList<>();
    protected ArrayList<JButton> deleteButtonList = new ArrayList<>();
    
    public static final String actionDeleteComponent = "deleteComponent";
    private static final String deleteButtonText = "x";
    
    public BaseComponentListWithDeleteButtons() {
        initComponents();
    }
    
    public int nComponents(){
        return componentList.size();
    }
    
    protected JPanel makeNewComponentPanel(){
        JPanel panel = new JPanel();
        BoxLayout layout = new BoxLayout(panel, BoxLayout.X_AXIS);
               
        panel.setLayout(layout);
        
        
        return panel;
    }
    
    protected javax.swing.JButton makeNewButton(){
        JButton button = new JButton();
        
        return button;
    }
    protected javax.swing.JButton makeNewDeleteButton(){
        JButton button = new JButton();
        button.setText(deleteButtonText);
        button.setActionCommand(actionDeleteComponent);
        button.addActionListener(this);
        return button;
    }

    @SuppressWarnings("unchecked")
    protected void initComponents() {
        BoxLayout layout = new BoxLayout(this, BoxLayout.Y_AXIS);
        this.setLayout(layout);
    }
    
    protected JPanel assembleComponentPanel(Component component, 
            JButton deleteButton, JPanel componentPanel){
        
       
        componentPanel.add(deleteButton);
        componentPanel.add(Box.createHorizontalGlue());
        componentPanel.add(component);
       
        deleteButtonList.add(deleteButton);
        componentList.add(component);
       
        deleteButton.addActionListener((ActionListener) component);
        System.out.println("Base: assembled component number "+nComponents());
        return componentPanel;
    }
    
    public JButton addComponent(Component component){
        JPanel componentPanel = makeNewComponentPanel();
        JButton deleteButton = makeNewDeleteButton();
        
        componentPanel.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
        add(assembleComponentPanel(component, deleteButton, componentPanel));
        validate();
        return deleteButton;
       
    }
    
    public void removeComponent(int index){
        
        componentList.remove(index);
        deleteButtonList.remove(index);
        remove(index);
        validate();
//        repaint(getBounds());
        System.out.println("component should be removed!");
    }
    public void removeComponent(Component component){
        int index = componentList.indexOf(component);
        assert index>=0:"component not found in componentList";
        removeComponent(index);
    }
    public void removeComponent(JButton pushButton){
        int index = deleteButtonList.indexOf(pushButton);
        assert index>=0:"button not found in pushButtonList";
        System.out.println("index of button is: "+index);
        removeComponent(index);
    }
    
    public Component getListComponent(int i){
        return componentList.get(i);
    }
    
    @Override
    public void actionPerformed(ActionEvent evt) {
//        System.out.println("actionPerformed in BaseComponentListWithDeleteButtons: "+evt.getActionCommand());
        if (actionDeleteComponent.equals(evt.getActionCommand())) {
            JButton button = (JButton)evt.getSource();
            removeComponent(button);
        }
} 

        
}
