import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,
    QCheckBox, QLineEdit, QLabel
)
from PyQt5.QtCore import QTimer
from xoppylib.xoppy_xraylib_util import descriptor_kind_index

class filterRow(QWidget):
    def __init__(self, compound="", thickness=0.0, density=0.0,parent=None):
        super(filterRow, self).__init__(parent)
        
        self.layout = QHBoxLayout()

        self.checkbox = QCheckBox()
        self.layout.addWidget(self.checkbox)

        # fields for Compound, Thickness, and Density
        if compound:
            self.compound_edit = QLineEdit(compound)
        else:
            self.compound_edit = QLineEdit()
        self.compound_edit.setPlaceholderText("Compound")
        self.layout.addWidget(self.compound_edit)
        
        
        if self.isvalidvalue(thickness):
            self.thickness_edit = QLineEdit(str(thickness))
        else:
            self.thickness_edit = QLineEdit()
        self.thickness_edit.setPlaceholderText("Thickness (mm)")
        self.layout.addWidget(self.thickness_edit)

        if self.isvalidvalue(density):
            self.density_edit = QLineEdit(str(density))
        else:
            self.density_edit = QLineEdit()
        self.density_edit.setPlaceholderText("Density (g/cc)")
        self.layout.addWidget(self.density_edit)

        self.delete_button = QPushButton("Delete")
        self.delete_button.clicked.connect(self.delete_row)
        self.layout.addWidget(self.delete_button)
        
        self.compound_edit.editingFinished.connect(lambda: parent.update_lists())
        self.thickness_edit.editingFinished.connect(lambda: parent.update_lists())
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
                self.compound_edit.text(),
                self.thickness_edit.text(),
                self.density_edit.text()
            )
        return None
    def delete_row(self):
        self.parent().remove_row(self)
        
    def blink_red(self,widgettoblink):
        widgettoblink.setStyleSheet("QLineEdit { background-color: red; }")
    
        # timer to revert color back
        QTimer.singleShot(500, lambda: self.revert_color(widgettoblink))  # Change color back after 500 ms

    def revert_color(self,widgettoblink):
        # Revert to original style
        widgettoblink.setStyleSheet("")

class FilterManager(QWidget):
    def __init__(self,parent_obj):
        super().__init__()
        self.parent_obj = parent_obj
        self.setWindowTitle("Filter Manager")
        self.layout = QVBoxLayout()
        
        self.add_button = QPushButton("Add filter")
        self.add_button.clicked.connect(self.add_filter)
        self.layout.addWidget(self.add_button)


        self.filters_container = QVBoxLayout()
        
        self.layout.addLayout(self.filters_container)

        self.setLayout(self.layout)
        self.activeCompounds = []
        self.activeThicknesses = []
        self.activeDensities = []
        

    def add_filter(self, compound="", thickness=0.0, density=0.0,enable=False):
        new_row = filterRow(compound, thickness, density, self)
        self.filters_container.addWidget(new_row)
        if enable:
            new_row.checkbox.setCheckState(2)
            self.update_lists()
            

        
    def isvalidcompound(self,compound):
        
        try:
            kind = descriptor_kind_index(compound)
        except:
            return 0         
        if not (kind==-1):
            return 1
            

    def remove_row(self, row):
        row.deleteLater()  
        self.filters_container.removeWidget(row)  
        row.setParent(None)  
        
    def update_lists(self):
        self.activeCompounds.clear()
        self.activeThicknesses.clear()
        self.activeDensities.clear()

        # go through rows and update lists based on checked boxes
        for i in range(self.filters_container.count()):
            row = self.filters_container.itemAt(i).widget()
            values = row.get_values()
            badval = 0
            if values:
                compound, thickness, density = values
                if compound and thickness and density:
                    if self.isvalidvalue(thickness) and self.isvalidvalue(density) and self.isvalidcompound(compound):
                        
                        
                        self.activeCompounds.append(compound)
                        self.activeThicknesses.append(float(thickness))
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
                    if not self.isvalidvalue(thickness):
                        row.blink_red(row.thickness_edit)
                    if not self.isvalidvalue(density):
                        row.blink_red(row.density_edit)  
                    if not self.isvalidcompound(compound):
                        row.blink_red(row.compound_edit)
            else:
                badval=1
            
        # print the lists 
        print("Active Compounds:", self.activeCompounds)
        print("Active Thicknesses:", self.activeThicknesses)
        print("Active Densities:", self.activeDensities)        
        
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
    window = FilterManager(app)
    window.resize(600, 200)
    
    # Example: Adding some initial rows with values
    window.add_filter("Al", 1, 2.7)

    window.show()
    sys.exit(app.exec_())
