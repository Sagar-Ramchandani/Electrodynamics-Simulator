import numpy as np
import parser
import re

def generateMesh(x1,x2,y1,y2,z1,z2,step):
    '''
    Creates the required mesh for all future computations.
    The commented out line is legacy code 

    arange vs linspace

    arange: It has a flaw that if the distance between x1 and x2 is not divisible by the step size,
    then the mesh will leave out some region close to x2 which is the remainder,
    and that would not be included in the computation

    linspace: with the way linspace is setup over here, it will always create equally spaced points between 
    the endpoints and there won't be a remainder left which is why arange was replaced by linspace
    '''
    #return np.meshgrid(np.arange(x1,x2,step),np.arange(y1,y2,step),np.arange(z1,z2,step))
    return np.meshgrid(np.linspace(x1,x2,int((x2-x1)/step)),
            np.linspace(y1,y2,int((y2-y1)/step)),
            np.linspace(z1,z2,int((z2-z1)/step)))

def HSVToRGB(h,s,v):
    '''
    I link the following answer from computerscience stack exchange which explains in great detail 
    how the HSV to RGB conversion works, the code below is a translation of the pseudocode written in the answer

    https://cs.stackexchange.com/questions/64549/convert-hsv-to-rgb-colors
    
    Should the above link not work at a future date, I am also linking the respective wikipedia article

    https://en.wikipedia.org/wiki/HSL_and_HSV#From_HSV

    Assuming that no matter when you reference the code, a copy of wikipedia will survive.

    The need for this conversion is because
    We want to assign a particular colour or Hue to the magnitude of a vector
    This naturally means we are working in a HSV format where we don't care about the Saturation and the Value
    I set these to 1 during the processData function for simplicity.

    The Hue is what sets the colour and we simply map the magnitudes to the range of possible values of Hue

    The conversion is necessary because the plotting library to the best of my knowledge only supports RGB.
    '''
    maxValue=v
    chroma=s*v
    minValue=maxValue-chroma
    if h>=300:
        hPrime=(h-360)/60
    else:
        hPrime=(h/60)
    if hPrime>=-1 and hPrime<1:
        if hPrime<0:
            r=maxValue
            g=minValue
            b=g-hPrime*chroma
        else:
            r=maxValue
            b=minValue
            g=b+hPrime*chroma
    elif hPrime<3:
        if hPrime<2:
            g=maxValue
            b=minValue
            r=b-(hPrime-2)*chroma
        else:
            r=minValue
            g=maxValue
            b=r+(hPrime-2)*chroma
    elif hPrime<5:
        if hPrime<4:
            r=minValue
            b=maxValue
            g=r-(hPrime-4)*chroma
        else:
            g=minValue
            b=maxValue
            r=g+(hPrime-4)*chroma
    else:
        r,g,b=1,1,1
    
    return (r,g,b)
    
def inverseLaplacian(laplacianValue,field,step):
    '''
    The following function performs a single inverse laplacian routine,
    for the most accuracy one would ideally have to perform the iteration infinitely many times,
    but for most cases it converges pretty quickly.

    It is simply an algebric inverse of the defintion of the discrete numerical laplacian operation.

    The use of numpy roll function is to ease the computation, it simply shifts the matrix in a given axis by a given step
    or if you imagine a 2D matrix on a paper, then it rolls the paper into a cylinder and rotates it by a certain step
    
    The use of the roll function does create an issue at the boundaries because it creates periodic boundaries, 
    this would not be ideal for the vast majority of use cases and thus it is recommended to not place charges 
    very close to the boundary of the mesh.

    Parameters:
    laplacianValue:
        It is the known value of the laplacian which has to be inverted, in our case it is the potentials.
    field:
        The initial field assumption, which will be iterated multiple times to get closer to actual field value,
        usually in our case we simply use a zero field, this is because of the fact that for any physical case
        where there is a finite amount of charge or current the field will be zero at infinity.
    step:
        The step size used when creating the mesh needs to be passed on to the computing function.

    Returns:
        The result of a single iteration of the inverseLaplacian operation on the given field.
    '''
    iteratedField=(
    np.roll(field, -1, axis=0)+np.roll(field, 1, axis=0)+
    np.roll(field, -1, axis=1)+np.roll(field, 1, axis=1)+
    np.roll(field, -1, axis=2)+np.roll(field, 1, axis=2)-
    step**2*laplacianValue)*(1/8)
    return iteratedField

def generatePotentialMatrix(mesh,step,epsilon,chargeDensityMatrix,iterations):
    '''
    Generates the potential from the charges
    It simply generates an initial assumption of zero for the field,
    uses the charge data as the laplacianValue 
    and iterates using the inverseLaplacian function a set number of times

    Parameters:
    mesh:
        The mesh on which the computation will occur
    step:
        The step size using which the mesh was created
    epsilon:
        Physical constant for electric charges
    chargeDensityMatrix:
        Matrix containing information on the strength of charges over the mesh
        is generated via the use of generateChargeMatrix function.
    iterations:
        Number of iterations to be used for the inverseLaplacian
        
    Returns:
    The Potential field for a given set of charges

    See also:
    generateMagneticPotentialMatrix()
    '''
    potentialField=np.zeros(mesh[0].shape)
    laplacianValue=(-1/epsilon)*chargeDensityMatrix
    for i in range(iterations):
        potentialField=inverseLaplacian(
            laplacianValue,potentialField,step)
    return potentialField

def generateMagneticPotentialMatrix(mesh,step,mu,epsilon,
        currentDensity,iterations,inducedField=0):
    '''
    Generates a potential from the currents
    It simply generates an initial assumption of zero for the field,
    the inital field is a vector field which is different as compared to the 
    scalar field used in generatePotentialMatrix()
    It then uses the current data as the laplacianValue
    and iteraties using the inverseLaplacian seperately on all 3 dimension
    of the vector field a set number of times.

    Paramaters:
    mesh:
        The mesh on which the computation will occur
    step:
        The step size using which the mesh was created
    mu:
        Physical constant for currents
    epsilon:
        Physical constant for electric charges
    currentDensity:
        Matrix containing information of the strength of currents over the mesh
        is generated using the generateCurrentMatrix function
    iterations:
       Number of iterations to be used for the inverseLaplacian function
    inducedField:
        An attempt at account for the induced fields, set by default to zero
    
    Returns:
    The potential field for a given set of currents

    See also:
    generatePotentialMatrix()
    '''
    meshShape=mesh[0].shape
    magneticPotentialField=np.asarray([
    np.zeros(meshShape),np.zeros(meshShape),np.zeros(meshShape)])
    laplacianValue=(-1*mu)*currentDensity+(mu*epsilon)*inducedField
    for i in range(iterations):
        result=[]
        for dimension in range(3):
            result.append(inverseLaplacian(
            laplacianValue[dimension],magneticPotentialField[dimension],step))
        magneticPotentialField=np.asarray(result)
    return magneticPotentialField

def generateElectricField(mesh,potentialField):
    '''
    A function to generate the ElectricField from a given Potential Field,
    It works by simply taking the gradient of the potential field over the given mesh
    to get the Electric Field

    See also:
    generateMagneticField()
    '''
    x,y,z=mesh
    u,v,w=np.gradient(potentialField)
    #Changing around the vectors and their directions to make it compatible
    u,v=v,u
    u=-u
    v=-v
    w=-w
    return (x,y,z,u,v,w)

def curl(field,step):
    '''
    A function that simply takes the curl of the given field with a given step size
    '''
    field_u,field_v,field_w=field

    du_dy=(np.roll(field_u, -1, axis=1)-np.roll(field_u, 1, axis=1))/(2*step)
    du_dz=(np.roll(field_u, -1, axis=2)-np.roll(field_u, 1, axis=2))/(2*step)
    dv_dx=(np.roll(field_v, -1, axis=0)-np.roll(field_v, 1, axis=0))/(2*step)
    dv_dz=(np.roll(field_v, -1, axis=2)-np.roll(field_v, 1, axis=2))/(2*step)
    dw_dx=(np.roll(field_w, -1, axis=0)-np.roll(field_w, 1, axis=0))/(2*step)
    dw_dy=(np.roll(field_w, -1, axis=1)-np.roll(field_w, 1, axis=1))/(2*step)

    curl=np.asarray([(dw_dy-dv_dz),(du_dz-dw_dx),(dv_dx-du_dy)])
    return curl

def generateMagneticField(mesh,step,potentialField):
    '''
    A function to generate the MagneticField from a given Potential Field,
    It works by simply taking the curl of the potential field over the given mesh
    to get the Magnetic Field

    See also:
    generateElectricField()
    '''

    x,y,z=mesh
    u,v,w=curl(potentialField,step)
    u,v=v,u
    u=-u
    v=-v
    w=-w
    return (x,y,z,u,v,w)

def inducedElectricField(mesh,timeStep,magneticPotentialInitial,
        magneticPotentialFinal):
    '''
    An attempt at account for the induced Electric field due to changing magnetic fields over time
    It works by computing the the numerical derivative of the field over each co-ordinate seperately 
    '''
    x,y,z=mesh
    u=(magneticPotentialInitial[0]-magneticPotentialFinal[0])/timeStep
    v=(magneticPotentialInitial[1]-magneticPotentialFinal[1])/timeStep
    w=(magneticPotentialInitial[2]-magneticPotentialFinal[2])/timeStep
    u,v=v,u
    return (x,y,z,u,v,w) 

def inducedMagneticField(mesh,timeStep,electricFieldInitial,electricFieldFinal):
    '''
    An incomplete attempt at accounting for the induced Magnetic field due to changing electric fields
    '''
    #result=(electricFieldFinal-electricFieldInitial)/timeStep
    #return result
    pass

def chargeObjectsFromExpression(posExpr,chargeExpr):
    '''
    This function takes the description of the position and strength of charges 
    in the form of a mathematical expression and 
    filters out the expression to a suitable format to generate a chargeMatrix
    '''
    expression=re.sub(r'\^','**',posExpr)
    individualObjects=re.split(r'\s*or\s*',expression)
    individualCharges=re.split(r'\s*or\s*',chargeExpr)
    andObjects=[]
    orObjects=[]
    andObjCharge=[]
    orObjCharge=[]
    for obj in range(len(individualObjects)):
        if re.search(r'\s*and\s*',individualObjects[obj]):
            andObjects.append(re.split(r'\s*and\s*',individualObjects[obj]))
            andObjCharge.append(individualCharges[obj])
        else:
            orObjects.append(individualObjects[obj])
            orObjCharge.append(individualCharges[obj])
    return (orObjects,andObjects,orObjCharge,andObjCharge)

def currentObjectsFromExpression(posExpr,currentExpr):
    '''
    This function takes the description of the position and strength of currents
    in the form of a mathematical expression and
    filters out the expression to a suitable format to generate a currentMatrix
    '''
    expression=re.sub(r'\^','**',posExpr)
    individualObjects=re.split(r'\s*or\s*',expression)
    individualCharges=re.split(r'\s*or\s*',currentExpr)
    andObjects=[]
    orObjects=[]
    orObjCurrent=[]
    andObjCurrent=[]
    for obj in individualObjects:
        if re.search(r'\s*and\s*',obj):
            andObjects.append(re.split(r'\s*and\s*',obj))
        else:
            orObjects.append(obj)
    return (orObjects,andObjects)

def generateChargeMatrix(mesh,time,posExpr,chargeExpr):
    '''
    This function generates a chargeMatrix that contains information on the 
    position and strength of charges based on the information from the
    mathematical expression. 
    It internally calls the function chargeObjectsFromExpression() for 
    formatting of the expression.

    orObjects are objects defined by a single expression for their position.

    andObjects are definied by more than one expression for the position,
    in essence by doing an and or an intersection in a 3D Venn Diagram
    
    in process it does use eval() and unfortunately to the best of my current 
    knowledge at the time of the writing I did not find a better way,

    eval() should be avoided because it allows for the execution of code entered by the user.
    Although in this particular case where the user has complete control over the source code,
    it would be acceptable, but it cannot be readily deloyed on something like a server without risks.
    '''
    x,y,z=mesh

    r=np.sqrt(np.square(x)+np.square(y)+np.square(z))
    s=np.sqrt(np.square(x)+np.square(y))

    phi=np.arctan2(y,x)*180/np.pi+180*np.ones(x.shape)
    theta=np.arctan2(s,z)*180/np.pi

    t=time*np.ones(x.shape)
    result=np.zeros(x.shape)

    orObjects,andObjects,orObjCharge,andObjCharge=chargeObjectsFromExpression(
            posExpr,chargeExpr)
    if orObjects!=['']:
        for i in range(len(orObjects)):
            objectPosition=np.vectorize(int)(eval(parser.expr(orObjects[i]).compile()))
            result+=objectPosition*eval(parser.expr(orObjCharge[i]).compile())

    if andObjects!=[]:
        for i in range(len(andObjects)):
            objectPosition=np.ones(x.shape)
            for j in range(len(andObjects[i])):
                objectPosition=np.logical_and(objectPosition,
                eval(parser.expr(andObjects[i][j]).compile()))
            result+=objectPosition*eval(parser.expr(andObjCharge[i]).compile())
    return result
 
def generateCurrentMatrix(mesh,time,posExpr,currentExpr):
    '''
    This function generates a currentMatrix that contains information on the 
    position and strength of currents based on the information from the
    mathematical expression. 
    It internally calls the function currentObjectsFromExpression() for 
    formatting of the expression.

    orObjects are objects defined by a single expression for their position.

    andObjects are definied by more than one expression for the position,
    in essence by doing an and or an intersection in a 3D Venn Diagram
    
    in process it does use eval() and unfortunately to the best of my current 
    knowledge at the time of the writing I did not find a better way,

    eval() should be avoided because it allows for the execution of code entered by the user.
    Although in this particular case where the user has complete control over the source code,
    it would be acceptable, but it cannot be readily deloyed on something like a server without risks.
    '''

    x,y,z=mesh

    r=np.sqrt(np.square(x)+np.square(y)+np.square(z))
    s=np.sqrt(np.square(x)+np.square(y))

    phi=np.arctan2(y,x)*180/np.pi+180*np.ones(x.shape)
    theta=np.arctan2(s,z)*180/np.pi

    t=time*np.ones(x.shape)
    result=np.asarray([np.zeros(x.shape),np.zeros(y.shape),np.zeros(z.shape)])

    orObjects,andObjects,orObjCurrent,andObjCurrent=chargeObjectsFromExpression(
            posExpr,currentExpr)
    orObjCurrent=list(map(lambda string: re.split(r'\s*,\s*',string),orObjCurrent))
    andObjCurrent=list(map(lambda string: re.split(r'\s*,\s*',string),andObjCurrent))
    if orObjects!=['']:
        for i in range(len(orObjects)):
            objectPosition=np.vectorize(int)(eval(parser.expr(orObjects[i]).compile()))
            result[0]+=objectPosition*eval(parser.expr(orObjCurrent[i][0]).compile())
            result[1]+=objectPosition*eval(parser.expr(orObjCurrent[i][1]).compile())			
            result[2]+=objectPosition*eval(parser.expr(orObjCurrent[i][2]).compile())
    if andObjects!=[]:
        for i in range(len(andObjects)):
            objectPosition=np.ones(x.shape)
            for j in range(len(andObjects[i])):
                objectPosition=np.logical_and(objectPosition,
                eval(parser.expr(andObjects[i][j]).compile()))
            result[0]+=objectPosition*eval(parser.expr(andObjCurrent[i][0]).compile())
            result[1]+=objectPosition*eval(parser.expr(andObjCurrent[i][1]).compile())
            result[2]+=objectPosition*eval(parser.expr(andObjCurrent[i][2]).compile())
    return result

