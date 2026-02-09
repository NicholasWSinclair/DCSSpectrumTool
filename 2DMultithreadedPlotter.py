import sys
import numpy as np
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel, QLineEdit, QPushButton, 
    QProgressBar, QApplication
)
from PyQt5.QtCore import QThreadPool, QRunnable, pyqtSignal, QObject
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class WorkerSignals(QObject):
    progress = pyqtSignal(int, int, float)  # x, y coordinates and pixel value
    finished = pyqtSignal()

class PixelCalculationTask(QRunnable):
    def __init__(self, x_range, y_range, signals):
        super().__init__()
        self.x_range = x_range  # Range of x-coordinates to calculate
        self.y_range = y_range  # Range of y-coordinates to calculate
        self.signals = signals  # Signals to communicate progress
    
    def run(self):
        # Perform the calculation for each pixel in the assigned range
        for x in self.x_range:
            for y in self.y_range:
                # Example pixel calculation (e.g., pattern or gradient)
                value = (np.sin(x * 0.1) + np.cos(y * 0.1)) * 127.5 + 127.5
                self.signals.progress.emit(x, y, value)  # Emit progress for each pixel
                
            
        self.signals.finished.emit()  # Emit finished signal when done

class PlotterWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        
        # Main layout
        main_layout = QVBoxLayout(self)
        
        # Setup matplotlib figure and canvas for plotting
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(self.canvas)
        
        # GroupBox for parameter inputs
        params_groupbox = QGroupBox("Image Parameters")
        params_layout = QVBoxLayout()
        
        # Image width and height inputs
        width_layout = QHBoxLayout()
        width_layout.addWidget(QLabel("Width:"))
        self.width_edit = QLineEdit("100")
        width_layout.addWidget(self.width_edit)
        
        height_layout = QHBoxLayout()
        height_layout.addWidget(QLabel("Height:"))
        self.height_edit = QLineEdit("100")
        height_layout.addWidget(self.height_edit)
        
        # Add layouts to the groupbox
        params_layout.addLayout(width_layout)
        params_layout.addLayout(height_layout)
        
        # Button for calculation and plotting
        self.plot_button = QPushButton("Generate Image")
        self.plot_button.clicked.connect(self.start_calculation)
        params_layout.addWidget(self.plot_button)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        params_layout.addWidget(self.progress_bar)
        
        # Set layout for the groupbox and add to main layout
        params_groupbox.setLayout(params_layout)
        main_layout.addWidget(params_groupbox)

        # Initialize thread pool
        self.thread_pool = QThreadPool()
        self.thread_pool.setMaxThreadCount(4)  # Set number of threads in the pool
        self.image_data = None
        self.total_pixels = 0
        self.calculated_pixels = 0
    
    def start_calculation(self):
        # Retrieve values from input fields
        try:
            width = int(self.width_edit.text())
            height = int(self.height_edit.text())
        except ValueError:
            print("Invalid input; please enter integer values.")
            return

        # Initialize image data array (grayscale image)
        self.image_data = np.zeros((height, width), dtype=np.uint8)
        
        # Set up the plot
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.image_display = self.ax.imshow(self.image_data, cmap="gray", vmin=0, vmax=255)
        self.ax.axis('off')  # Hide axis for a cleaner view
        
        # Setup progress bar
        self.total_pixels = width * height
        self.calculated_pixels = 0
        self.progress_bar.setMaximum(self.total_pixels)
        self.progress_bar.setValue(0)

        # Divide image into chunks and start worker tasks
        chunk_size = width // self.thread_pool.maxThreadCount()
        for i in range(self.thread_pool.maxThreadCount()):
            x_range = range(i * chunk_size, (i + 1) * chunk_size)
            y_range = range(height)
            
            # Setup signals for each worker
            signals = WorkerSignals()
            signals.progress.connect(self.update_image_pixel)
            signals.finished.connect(self.finish_calculation)
            
            # Create and start the task
            task = PixelCalculationTask(x_range, y_range, signals)
            self.thread_pool.start(task)
    
    def update_image_pixel(self, x, y, value):
        # Update the pixel in the image data
        self.image_data[y, x] = value  # Note: numpy uses (y, x) indexing
        self.image_display.set_data(self.image_data)  # Update image data
        
        # Update progress and display
        self.calculated_pixels += 1
        self.progress_bar.setValue(self.calculated_pixels)
        
        if np.mod(self.calculated_pixels,100)==0:
            self.canvas.draw()
            
        

    def finish_calculation(self):
        # Check if all pixels are calculated
        if self.calculated_pixels >= self.total_pixels:
            print("Image generation finished.")
            
            
    # Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PlotterWidget()
    window.setWindowTitle("Image Generator with Multithreaded Pixel Calculation")
    window.resize(400, 400)
    window.show()
    sys.exit(app.exec_())
