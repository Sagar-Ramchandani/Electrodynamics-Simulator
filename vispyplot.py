from PyQt5.QtWidgets import *
from helpers import *
import numpy as np
from vispy import app, visuals, scene
from pprint import pprint
app.use_app('pyqt5')

class Plot():
    def __init__(self):
        '''
        Initializes the Plotting enviornment
        '''
        self.setupCanvas()
        self.initialPlot=True
    def setupCanvas(self):
        '''
        Does basic setup inculding the plot enviornment and the XYZ Axis lines,
        and the camera and its parameters such as POV and distance
        '''
        self.PlotAxis=scene.visuals.create_visual_node(visuals.XYZAxisVisual)
        self.canvas=scene.SceneCanvas(keys='interactive', title='plot3d',always_on_top=True)
        self.view=self.canvas.central_widget.add_view()
        self.view.camera='turntable'
        self.view.camera.fov=45
        self.view.camera.distance=10
        axisLength=100
        self.PlotAxis(pos=np.array([
                [-axisLength,0,0],[axisLength,0,0],
                [0,-axisLength,0],[0,axisLength,0],
                [0,0,-axisLength],[0,0,axisLength]]),parent=self.view.scene)

    def updateParameters(self,mesh,stepSize,time,timeStep,scale,booleanNormalize,booleanColour):
        '''
        Does initialization of the parameters of a certain plot, 
        It has to be called everytime the parameters are changed.

        Parameters:
        mesh:
            Pass on the mesh generated by np.meshgrid()
            The region in the mesh is the only region where computation will be done
        stepSize:
            The step size of the mesh.
            One assumes that the stepSize in all directions is the same
        time:
            Pass on the current time of the system
        timeStep:
            Pass on the time step of the system to be simulated. The units are arbritary
        scale:
            Adjusts the size of the vectors by scaling them with this constant
        booleanNormalize:
            Flag to normalize the vectors if set to True
        booleanColour:
            Flag to assign colour to the vectors based on their magnitude if set to True

        Returns:
        None 
        But updates the self object with the required parameters
        '''
        self.time=time
        self.timeStep=timeStep
        self.scale=scale
        self.booleanNormalize=booleanNormalize
        self.booleanColour=booleanColour
        self.mesh=mesh
        self.stepSize=stepSize
        self.Potential=np.zeros(mesh[0].shape)
        self.MagneticPotential=np.asarray([np.zeros(mesh[0].shape),
            np.zeros(mesh[0].shape),np.zeros(mesh[0].shape)])

    def updateFields(self,chargePosExpr,chargeExpr,currentPosExpr,currentExpr):
        '''
        Generates the matrices that hold info on charges and currents,
        from the expressions given as input and stores them in self object

        Paramters:
        chargePosExpr:
            The expression that defines 3D position of charges
        chargeExpr:
            The expression that defines the strength of charges in the entire mesh
        currentPosExpr:
            The expression defining the 3D Position where a current is present
        currentExpr:
            The expression defining the strength of current in the entire mesh

        Returns:
        None
        But stores the generated matrices in self object
        '''
        self.chargeMatrix=generateChargeMatrix(self.mesh,self.time,chargePosExpr,chargeExpr)
        self.currentMatrix=generateCurrentMatrix(self.mesh,self.time,currentPosExpr,currentExpr)

    def updateData(self,epsilon,mu,iterations):
        '''
        Generates the final co-ordinate data and vector data based on the charge and current Matrices,
        and then stores them in the self object

        Parameters:
        epsilon:
            Physical constant for charges, set to 1 for most demonstration purposes
        mu:
            Physical constant for currents, set to 1 for most demonstration purposes
        iterations:
            Number of iterations to complete in the computation process for the fields

        Returns:
        None
        But stores co-ordinate data in xMatrix,yMatrix,zMatrix
        and vector data in uMatrix,vMatrix,wMatrix of the self object
        '''
        #Generating the Potentials
        Potential=generatePotentialMatrix(self.mesh,self.stepSize,epsilon,self.chargeMatrix,iterations)
        MagneticPotential=generateMagneticPotentialMatrix(
                self.mesh,self.stepSize,mu,epsilon,self.currentMatrix,iterations)

        #Generating the fields
        xe,ye,ze,ue,ve,we=generateElectricField(self.mesh,Potential)
        xm,ym,zm,um,vm,wm=generateMagneticField(self.mesh,self.stepSize,MagneticPotential)
        
        #An attempt at finding the induced Electric Field
        xei,yei,zei,uei,vei,wei=inducedElectricField(self.mesh,self.timeStep,self.MagneticPotential,MagneticPotential)

        #Storing the Magnetic potential for current time to be reused for induced Electric Field calcs later
        self.MagneticPotential=MagneticPotential

        #One can and should simply use the co-ordinate data from a single field
        #This is because all the fields are computed on the same mesh at present
        self.xMatrix=xe# or xm
        self.yMatrix=ye# or ym
        self.zMatrix=ze# or zm
        self.uMatrix=(ue+um+uei)*self.scale
        self.vMatrix=(ve+vm+vei)*self.scale
        self.wMatrix=(we+wm+wei)*self.scale

    def processData(self):
        '''
        Processes the generated co-ordinate and vector data to be able to plot them on the plotting enviornment
        '''
        #Computing the colour of the vector based on its magnitude
        vectorMags=np.sqrt(np.square(self.uMatrix)+np.square(self.vMatrix)+np.square(self.wMatrix))
        vectorMagsNormalized=(vectorMags.flatten()-vectorMags.min())/ vectorMags.ptp()
        '''
        The scaling constant of vectorMagsNormalized selects the range of Hues in the HSV format to use

        By default it should be set from 0 to 255 but can be selected from 0 to 360
        '''
        scaledVectors=vectorMagsNormalized*255
        colourPairs=[]
        for vector in scaledVectors:
                colourPairs.append(vector)
                colourPairs.append(vector)
        colourPairs=np.asarray(colourPairs)
        
        #We convert the colours from HSV system to RGB because the plotting library supports RGB
        self.colourVectors=list(map(lambda x: HSVToRGB(x,1,1),scaledVectors))
        #Legacy code commented out
        #self.colourPairs=np.repeat(self.colourVectors,2)
        self.colourPairs=list(map(lambda x: HSVToRGB(x,1,1),colourPairs))
        
        #Normalizes the vectors and scales them up appropriately
        #Suggestion: Could possibly use a higher dimensional matrix to store all co-ordinate data
        if self.booleanNormalize:
                self.uMatrix=self.uMatrix*self.scale/vectorMags
                self.vMatrix=self.vMatrix*self.scale/vectorMags
                self.wMatrix=self.wMatrix*self.scale/vectorMags

        #Computes the shifted co-ordinates as a result of the vectors
        xPrime=self.xMatrix+self.uMatrix
        yPrime=self.yMatrix+self.vMatrix
        zPrime=self.zMatrix+self.wMatrix
        #Computes the start and end co-ordinates of the vectors
        vectorStart=np.c_[self.xMatrix.flatten(),self.yMatrix.flatten(),self.zMatrix.flatten()]
        vectorEnd=np.c_[xPrime.flatten(),yPrime.flatten(),zPrime.flatten()]
        #Changes the format of the vector data to match that of the plotting library
        self.line=[]
        self.arrows=[]
        for i in range(vectorStart.shape[0]):
                self.line.append(vectorStart[i])
                self.line.append(vectorEnd[i])
                self.arrows.append(np.array([vectorStart[i],vectorEnd[i]]).reshape((1,6)))
        self.line=np.asarray(self.line)
        self.arrows=np.asarray(self.arrows)
        self.arrows=self.arrows.reshape((self.arrows.shape[0],6))

        #Filter out arrows below certain magnitude threshold to clean up the plot
        for i in range(len(vectorMagsNormalized)):
            #A Threshold of x cuts off any vectors with magnitude of x% of longest vector
            if np.sqrt(abs(np.dot(vectorMagsNormalized[i],vectorMagsNormalized[i])))<0.05:
                        self.arrows[i]=None
                        if self.booleanNormalize:
                                self.line[2*i]=None
                                self.line[2*i-1]=None
    def update(self):
        '''
        Since the vectors are already created and number of vectors won't change for given mesh.
        We simply change the values of the vectors in memory instead of recreating them for performance.

        This function is typically meant to be used to do time-dependent simulations and update the vectors with time
        '''
        self.PlotArrows3D._arrow_color = self.colourVectors
        self.PlotArrows3D._arrows_changed = True
        self.PlotArrows3D.set_data(pos=self.line, color=self.colourPairs, width=5, connect='segments',arrows=self.arrows)
    def updateTime(self):
        '''
        Simply updates the internal time of the system by the timeStep entered during initialization
        '''
        self.time+=self.timeStep
    def show(self):
        '''
        This function plots the vectors on the created enviornment 
        '''
        if self.booleanColour!=True:
                self.colourVectors=[1,1,1]
                self.colourPairs=[1,1,1]
        if self.initialPlot:
                self.PlotArrows3D=scene.visuals.Arrow(pos=self.line,color=self.colourPairs, width=5, method='gl',
                                        arrow_type="curved",arrow_size=5.0,
                                        arrows=self.arrows,connect='segments',
                                        arrow_color=self.colourVectors,
                                        antialias=True,parent=self.view.scene)
                self.initialPlot=False
        else:
                self.update()

class VispyPlot(QWidget):
    def  __init__(self,parent=None):
        '''
        Initializes the plotting system so it can be used in a pyqt5 GUI
        '''
        QWidget.__init__(self,parent)
        self.graph=Plot()
        self.canvas=self.graph.canvas
        self.verticalLayout=QVBoxLayout()
        self.verticalLayout.addWidget(self.canvas.native)
        self.setLayout(self.verticalLayout)
