import sys
import time
from helpers import *
import numpy as np
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.uic import loadUi
from vispy import app

app.use_app('pyqt5')
class userInterface(QMainWindow):
    def __init__(self):
        '''
        Loading the external user interface file
        Please note that the .ui file was created via a program called
        QtDesigner which is a part of the pyt5-tools module, 
        while not necessary to run the code, 
        it is recommended for editing the ui.
        
        Another option would be to convert the .ui file into a py file using pyqt5,
        and then manually making changes to the resulting .py file
        '''
        super(userInterface,self).__init__()
        #QMainWindow.__init__(self)
        loadUi("userInterface.ui",self)
        self.setWindowTitle("Electrodynamics Simulator")
        '''
        Setup a Threadpool to simultaneously run gui and computation
        This threadpool is needed because otherwise while running the computation,
        the GUI will become unresponsive until the computation is finished,
        while this might be acceptable for simpler computation work,
        it is not the case for any computation that takes greater than roughly hundred milliseconds
        The 100ms is to my knowledge the average response time of the human eye.
        '''
        self.threadpool=QThreadPool()
        self.PlotpushButton.clicked.connect(self.plotGraph)

    def plotGraph(self):
        '''
        Grabbing the parameters from the sliders and text of GUI
        Some parameters are manually definied such as epsilon, mu, iterations
        This would be troubling but my skills at making a decent GUI are limited,
        it might be improved at a future update but for manually changing the 
        paramters over here is the only option. Apologies.
        '''
        chargePosExpr=self.PositionLine.text()
        chargeExpr=self.ChargeLine.text()
        currentPosExpr=self.CurrentPositionLine.text()
        currentExpr=self.CurrentLine.text()
        booleanNormalize=self.NormalizecheckBox.isChecked()
        booleanColour=self.ColourcheckBox.isChecked()
        stepSize=int(self.StepSizeSlider.value())/10
        vectorLength=self.LengthSlider.value()
        if booleanNormalize:
                vectorLength/=10
        epsilon=1
        mu=1
        iterations=8
        meshSide=3
        '''
        Starting the computation based on the above paramters
        '''
        mesh=generateMesh(-meshSide,meshSide,-meshSide,meshSide,-meshSide,meshSide,stepSize)
        self.graphWidget.graph.updateParameters(mesh,stepSize,0,0.1,
        vectorLength,booleanNormalize,booleanColour)
        self.graphWidget.graph.updateFields(chargePosExpr,chargeExpr,currentPosExpr,currentExpr)
        self.Exprs=(chargePosExpr,chargeExpr,currentPosExpr,currentExpr)
        self.graphWidget.graph.updateData(epsilon,mu,iterations)
        self.variables=(epsilon,mu,iterations)
        self.graphWidget.graph.processData()
        self.graphWidget.graph.show()
        self.update()

    def update(self):
        '''
        Function that updates the vectors by re-running the computation and then plotting them on the enviornment
        '''
        updater=TimeStepper(self)
        self.threadpool.start(updater)

class TimeStepper(QRunnable):
    def __init__(self,plotterSelf):
        super(TimeStepper,self).__init__()
        self.plotter=plotterSelf
    @pyqtSlot()
    def run(self):
        for i in range(10):
            Exprs=self.plotter.Exprs
            epsilon,mu,iterations=self.plotter.variables
            self.plotter.graphWidget.graph.updateTime()
            self.plotter.graphWidget.graph.updateFields(Exprs[0],Exprs[1],Exprs[2],
                    Exprs[3])
            self.plotter.graphWidget.graph.updateData(epsilon,mu,iterations)
            self.plotter.graphWidget.graph.processData()
            self.plotter.graphWidget.graph.show()

app=QApplication(['Electrodynamics'])
window=userInterface()
window.show()
sys.exit(app.exec_())
