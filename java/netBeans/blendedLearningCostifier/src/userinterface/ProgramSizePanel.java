/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package userinterface;

import blendedlearningprogram.BlendedLearningModelEnum;
import blendedlearningprogram.BlendedLearningModelFactory;

import blendedlearningprogram.BaseBlendedLearningModel;
import blendedlearningprogram.ProgramSize;
import java.awt.event.ActionListener;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JSpinner;

/**
 *
 * @author mcannamela
 */
public class ProgramSizePanel extends javax.swing.JPanel {
    ProgramSize programSize = new ProgramSize();
    BlendedLearningModelFactory blendedLearningModelFactory = new BlendedLearningModelFactory() ;
    
    private JButton dButton_programSizeChanged = new JButton();
    
    public static final String ACTION_PROGRAM_SIZE_CHANGED = "programSizeChanged";
    
    /**
     * Creates new form ProgramSizePanel
     */
    public ProgramSizePanel() {
        initComponents();
        jComboBox_blendedLearningModel.setModel(new ComboBoxModel_blendedLearningModels());
        dButton_programSizeChanged.setActionCommand(ACTION_PROGRAM_SIZE_CHANGED);
    }
    
    public void addProgramSizeChangedListener(ActionListener listener){
        dButton_programSizeChanged.addActionListener(listener);
    }
    
    public void setProgramSize(ProgramSize programSize){
        this.programSize = programSize;
        jSpinner_nrStudents.setValue(programSize.getNrStudents());
        jComboBox_blendedLearningModel.setSelectedItem(programSize.getBlendedLearningModel().getType());
        programSizeChanged();
    }
    
    

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jSpinner_nrStudents = new javax.swing.JSpinner();
        jLabel_programSize = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        jPanel1 = new javax.swing.JPanel();
        jLabel5 = new javax.swing.JLabel();
        jLabel_displayStudentTeacherRatio = new javax.swing.JLabel();
        jComboBox_blendedLearningModel = new javax.swing.JComboBox();
        jLabel4 = new javax.swing.JLabel();
        jButton1 = new javax.swing.JButton();
        jScrollPane1 = new javax.swing.JScrollPane();
        jTextArea_description = new javax.swing.JTextArea();
        jLabel7 = new javax.swing.JLabel();
        jLabel_displayNrTeachers = new javax.swing.JLabel();
        jLabel6 = new javax.swing.JLabel();
        jSpinner_nrPeriods = new javax.swing.JSpinner();

        jSpinner_nrStudents.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(100), Integer.valueOf(1), null, Integer.valueOf(10)));
        jSpinner_nrStudents.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSpinner_nrStudentsStateChanged(evt);
            }
        });

        jLabel_programSize.setFont(new java.awt.Font("Tahoma", 0, 18)); // NOI18N
        jLabel_programSize.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_programSize.setText("Program Size");

        jLabel3.setText("nrStudents:");

        jPanel1.setBorder(javax.swing.BorderFactory.createEtchedBorder(0));

        jLabel5.setText("Student/Teacher Ratio:");

        jLabel_displayStudentTeacherRatio.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_displayStudentTeacherRatio.setText("##");
        jLabel_displayStudentTeacherRatio.setBorder(javax.swing.BorderFactory.createBevelBorder(1));

        jComboBox_blendedLearningModel.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jComboBox_blendedLearningModelActionPerformed(evt);
            }
        });

        jLabel4.setText("Blended Learning Model:");

        jButton1.setText("+");

        jTextArea_description.setColumns(20);
        jTextArea_description.setLineWrap(true);
        jTextArea_description.setRows(5);
        jTextArea_description.setWrapStyleWord(true);
        jScrollPane1.setViewportView(jTextArea_description);

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel4)
                    .addComponent(jLabel5)
                    .addGroup(jPanel1Layout.createSequentialGroup()
                        .addGap(10, 10, 10)
                        .addComponent(jLabel_displayStudentTeacherRatio, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jComboBox_blendedLearningModel, javax.swing.GroupLayout.PREFERRED_SIZE, 107, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jButton1))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel4)
                .addGap(1, 1, 1)
                .addComponent(jComboBox_blendedLearningModel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton1)
                .addGap(20, 20, 20)
                .addComponent(jLabel5)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel_displayStudentTeacherRatio)
                .addContainerGap())
        );

        jLabel7.setText("nrTeachers:");

        jLabel_displayNrTeachers.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel_displayNrTeachers.setText("##");
        jLabel_displayNrTeachers.setBorder(javax.swing.BorderFactory.createBevelBorder(1));

        jLabel6.setText("nrPeriods:");

        jSpinner_nrPeriods.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(100), Integer.valueOf(1), null, Integer.valueOf(1)));
        jSpinner_nrPeriods.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSpinner_nrPeriodsStateChanged(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel_programSize, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(jLabel7)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                    .addComponent(jLabel_displayNrTeachers, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addGroup(layout.createSequentialGroup()
                                    .addComponent(jLabel3)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                    .addComponent(jSpinner_nrStudents, javax.swing.GroupLayout.PREFERRED_SIZE, 63, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jLabel6)
                                .addGap(18, 18, 18)
                                .addComponent(jSpinner_nrPeriods, javax.swing.GroupLayout.PREFERRED_SIZE, 63, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel_programSize)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jSpinner_nrStudents, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jSpinner_nrPeriods, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel6))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel7)
                    .addComponent(jLabel_displayNrTeachers))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(16, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void jComboBox_blendedLearningModelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jComboBox_blendedLearningModelActionPerformed
        BaseBlendedLearningModel model;
        model = blendedLearningModelFactory.makeBlendedLearningModel( (String) ((JComboBox)evt.getSource()).getSelectedItem());
        programSize.setBlendedLearningModel(model);
        jTextArea_description.setText(model.getDescription());
        programSizeChanged();
        
    }//GEN-LAST:event_jComboBox_blendedLearningModelActionPerformed

    private void jSpinner_nrStudentsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSpinner_nrStudentsStateChanged
        programSize.setNrStudents((int)((JSpinner) evt.getSource()).getValue());
        programSizeChanged();
    }//GEN-LAST:event_jSpinner_nrStudentsStateChanged

    private void jSpinner_nrPeriodsStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSpinner_nrPeriodsStateChanged
        programSize.setNrPeriods((int)((JSpinner) evt.getSource()).getValue());
        programSizeChanged();
    }//GEN-LAST:event_jSpinner_nrPeriodsStateChanged
    
    private void programSizeChanged(){
        displayNrTeachers();
        displayRatio();
        dButton_programSizeChanged.doClick();
        
    }
    private void displayNrTeachers(){
        jLabel_displayNrTeachers.setText(""+programSize.getNrTeachers());
    }
    private void displayRatio(){
        jLabel_displayStudentTeacherRatio.setText(""
                +programSize.getBlendedLearningModel().getStudentToTeacherRatio());
    }
    
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jButton1;
    private javax.swing.JComboBox jComboBox_blendedLearningModel;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel_displayNrTeachers;
    private javax.swing.JLabel jLabel_displayStudentTeacherRatio;
    private javax.swing.JLabel jLabel_programSize;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JSpinner jSpinner_nrPeriods;
    private javax.swing.JSpinner jSpinner_nrStudents;
    private javax.swing.JTextArea jTextArea_description;
    // End of variables declaration//GEN-END:variables

    private class ComboBoxModel_blendedLearningModels extends DefaultComboBoxModel<String>{
        public ComboBoxModel_blendedLearningModels(){
            super();
            for (BlendedLearningModelEnum model : BlendedLearningModelEnum.values()){
                this.addElement(model.toString());
            }
            
        }
        
    }
}
