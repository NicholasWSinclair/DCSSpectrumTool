import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,
    QCheckBox, QLineEdit, QLabel
)
from PyQt5.QtCore import QTimer
from xoppylib.xoppy_xraylib_util import descriptor_kind_index

class MirrorRow(QWidget):
    def __init__(self, MirrorStripe="", Angle=2.1, density='?',parent=None):
        super(MirrorRow, self).__init__(parent)
         
        self.layout = QHBoxLayout()
 
        self.checkbox = QCheckBox()
        self.layout.addWidget(self.checkbox)

        #    fields for MirrorStripe, Angle, and Density
        if MirrorStripe:
            self.MirrorStripe_edit = QLineEdit(MirrorStripe)
        else:
            self.MirrorStripe_edit = QLineEdit()
        self.MirrorStripe_edit.setPlaceholderText("MirrorStripe")
        self.layout.addWidget(self.MirrorStripe_edit)
        
        
        if self.isvalidvalue(Angle):
            self.Angle_edit = QLineEdit(str(Angle))
        else:
            self.Angle_edit = QLineEdit()
        self.Angle_edit.setPlaceholderText("Angle (mrad)")
        self.layout.addWidget(self.Angle_edit)

        if self.isvalidvalue(density):
            self.density_edit = QLineEdit(str(density))
        else:
            self.density_edit = QLineEdit()
        self.density_edit.setPlaceholderText("Density (g/cc)")
        self.layout.addWidget(self.density_edit)
 
        self.delete_button = QPushButton("Delete")
        self.delete_button.clicked.connect(self.delete_row)
        self.layout.addWidget(self.delete_button)
        
         
        self.MirrorStripe_edit.editingFinished.connect(lambda: parent.update_lists())
        self.Angle_edit.editingFinished.connect(lambda: parent.update_lists())
        self.density_edit.editingFinished.connect(lambda: parent.update_lists())
        self.checkbox.toggled.connect(lambda: parent.update_lists())
        
        self.setLayout(self.layout)
        
    def isvalidvalue(self,value):
        if value=='?':
            return 1
        try:
            if float(value)>0:
                return 1
            else:
                return 0
        except:
            return 0
                
        
    def get_values(self):
        if self.checkbox.isChecked():
            return (
                self.MirrorStripe_edit.text(),
                self.Angle_edit.text(),
                self.density_edit.text()
            )
        return None
    def delete_row(self):
        self.parent().remove_row(self)
        
    def blink_red(self,widgettoblink): 
        widgettoblink.setStyleSheet("QLineEdit { background-color: red; }")
 
        QTimer.singleShot(500, lambda: self.revert_color(widgettoblink))  #  set color back after 500 ms

    def revert_color(self,widgettoblink):
        widgettoblink.setStyleSheet("")

class MirrorManager(QWidget):
    def __init__(self,parent_obj):
        super().__init__()
        self.parent_obj = parent_obj
        self.setWindowTitle("Mirror Manager")
        self.layout = QVBoxLayout()
        
        
        self.add_button = QPushButton("Add Mirror")
        self.add_button.clicked.connect(self.add_Mirror)
        self.layout.addWidget(self.add_button)


        self.Mirrors_container = QVBoxLayout()
        
        
        self.layout.addLayout(self.Mirrors_container)

        self.setLayout(self.layout)
        self.activeMirrorStripes = []
        self.activeAngles = []
        self.activeDensities = []
        

    def add_Mirror(self, MirrorStripe="", Angle=0.0, density=0.0):
        new_row = MirrorRow(MirrorStripe, Angle, density, self)
        self.Mirrors_container.addWidget(new_row)
        
    def isvalidcompound(self,compound):
        
        try:
            kind = descriptor_kind_index(compound)
        except:
            return 0         
        if not (kind==-1):
            return 1
            

    def remove_row(self, row):
        row.deleteLater()  
        self.Mirrors_container.removeWidget(row)  
        row.setParent(None)  
        
    def update_lists(self):
        
        self.activeMirrorStripes.clear()
        self.activeAngles.clear()
        self.activeDensities.clear()

        
        for i in range(self.Mirrors_container.count()):
            row = self.Mirrors_container.itemAt(i).widget()
            values = row.get_values()
            badval = 0
            if values:
                MirrorStripe, Angle, density = values
                if MirrorStripe and Angle and density:
                    if self.isvalidvalue(Angle) and self.isvalidvalue(density) and self.isvalidcompound(MirrorStripe):
                        
                        
                        self.activeMirrorStripes.append(MirrorStripe)
                        self.activeAngles.append(float(Angle))
                        if density=='?':
                            self.activeDensities.append(density)
                        else:
                            self.activeDensities.append(float(density))
                    else:
                        badval=1
                else:
                    badval=1
                if badval:
                    row.checkbox.setCheckState(0)
                    if not self.isvalidvalue(Angle):
                        row.blink_red(row.Angle_edit)
                    if not self.isvalidvalue(density):
                        row.blink_red(row.density_edit)  
                    if not self.isvalidcompound(MirrorStripe):
                        row.blink_red(row.MirrorStripe_edit)
            else:
                badval=1
            
                          
        i = 1        
        for stripe,density,angle in zip(self.activeMirrorStripes,self.activeDensities,self.activeAngles):
            if density =='?':
                d = ''
            else:
                d = f'({float(density):.2f} g/cc)'
            print(f'Active Mirror {i}: {stripe}{d} @ {angle:.2f} mrad')
            i+=1
    
        
    def isvalidvalue(self,value):
        if value=='?':
            return 1
        try:
            if float(value)>0:
                return 1
            else:
                return 0
        except:
            return 0
        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MirrorManager(app)
    window.resize(600, 200)
    
    # Example: Adding some initial rows with values
    window.add_Mirror("Pt", 2.1, '?')

    window.show()
    sys.exit(app.exec_())
