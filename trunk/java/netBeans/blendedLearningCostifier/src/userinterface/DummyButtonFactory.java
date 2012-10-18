/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import javax.swing.JButton;

/**
 *
 * @author Michael
 */
public class DummyButtonFactory {
    
    public static JButton makeDummyButton(String actionCommand){
        JButton button = new JButton();
        button.setActionCommand(actionCommand);
        button.setVisible(false);
        
        return button;
    }
    
}
