import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

pi = np.pi
app = QtGui.QApplication([])

I = 1.0j

class Component(object):
    def __init__(self,attribute_parameters, attribute_names = None):
        self.attribute_parameters = attribute_parameters
        if attribute_names == None:
            attribute_names = attribute_parameters.keys()
        self.attribute_names  = attribute_names
        
        self.attributes = {}
        for name, values in self.attribute_parameters.items():
            self.setAttribute(name, values[2])
        
        self.outputs = {}
        self.outputs_names = []
        
    def setAttribute(self,name,x):
        self.attributes[name] = x
    
    def getAttribute(self,name):
        return self.attributes[name]
        
    
class DoubleBusRing(Component):
    def __init__(self): 
        
        attribute_parameters = {"r1"  : (0. ,1., 0.98), 
                                     "r2"  : (0. ,1., 0.98),
                                     "L"   : (50.,500.,70.),
                                     "tau" : (0.5,1.0,0.985)} # name : min, max, default
        
        attribute_names = ["r1", "r2", "tau", "L"]  # To have a list (ordered) of parameters
        super(self.__class__,self).__init__(attribute_parameters,attribute_names)
        
        # Define some constants
        self.ng = 4.2 # group index
        
        # Define the output functions:
        self.outputs = {"Through I" : self.getThroughIntensity,
                        "Drop I" : self.getDropIntensity,
                        "Through Phase" : self.getThroughPhase,
                        "Drop Phase" : self.getDropPhase
                        }
        self.outputs_names = ["Through I","Drop I","Through Phase","Drop Phase" ]
      
    def theta(self,lbda):
        L = self.getAttribute("L")
        return 2*pi*L*self.ng/lbda
    
    def getThrough(self,lbda):
        r1 = self.getAttribute("r1")
        r2 = self.getAttribute("r2")
        tau = self.getAttribute("tau")
        eth = np.exp(I*self.theta(lbda))
        amp =(r1-r2*tau*eth)/(1-r1*r2*tau*eth)
        return amp
    
    def getThroughPhase(self, lbda):
        amp = self.getThrough(lbda)
        return np.angle(amp)

    def getDropPhase(self, lbda):
        amp = self.getDrop(lbda)
        return np.angle(amp)    
    
    def getThroughIntensity(self, lbda):
        amp = self.getThrough(lbda)
        return np.abs(amp)**2

    def getDropIntensity(self, lbda):
        amp = self.getDrop(lbda)
        return np.abs(amp)**2    
    
    def getDrop(self,lbda):
        r1 = self.getAttribute("r1")
        r2 = self.getAttribute("r2")
        t1 = np.sqrt(1-r1**2)
        t2 = np.sqrt(1-r2**2)
        
        tau = self.getAttribute("tau")
        th = self.theta(lbda)
        eth = np.exp(I*th)
        amp =  -t1*t2*np.sqrt(tau)*np.exp(I*th/2)/(1-r1*r2*tau*eth)    
        return amp        
        
class SingleBusCROW(Component):
    def __init__(self): 
        attribute_parameters = {"r1"  : (0. ,1., 0.98), 
                                     "r2"  : (0. ,1., 0.98),
                                     "L1"   : (50.,500.,70.),
                                     "L2"   : (50.,500.,70.),
                                     "tau1" : (0.5,1.0,0.985), # name : min, max, default
                                     "tau2" : (0.5,1.0,0.985)} # name : min, max, default
        
        attribute_names = ["r1", "tau1" , "L1",  "r2", "tau2" , "L2"]  # To have a list (ordered) of parameters
 
        super(self.__class__,self).__init__(attribute_parameters,attribute_names)
        
        self.ng = 4.2
        
        # Define the output functions:
        self.outputs = {"Through I" : self.getThroughIntensity,
                        "Through Phase" : self.getThroughPhase,
                        "FOM" : self.getFOM,
                        }
        self.outputs_names = ["Through I","Through Phase","FOM"]
     
    
    def getThrough(self,lbda):
        r1 = self.getAttribute("r1")
        r2 = self.getAttribute("r2")

        tau1 = self.getAttribute("tau1")
        tau2 = self.getAttribute("tau2")

        t2 = np.sqrt(1-r2**2)
        t1 = np.sqrt(1-r1**2)

        L1 = self.getAttribute("L1")
        L2 = self.getAttribute("L2")

        g1 = tau1*np.exp(I*2*pi*L1*self.ng/lbda)
        g2 = tau2*np.exp(I*2*pi*L2*self.ng/lbda)
        
        A = g1 * (r2-g2)/(1-g2*r2)
        res = (r1 - A ) / (1-r1*A)
        return res 
    
    def getThroughPhase(self, lbda):
        amp = self.getThrough(lbda)
        return np.angle(amp)

    def getThroughIntensity(self, lbda):
        amp = self.getThrough(lbda)
        return 10*np.log10(np.abs(amp)**2)
    
    def getFOM(self,lbda):
        amp = self.getThrough(lbda)
        return  np.abs(amp)**2 * self.getThroughPhase(lbda) 

    
class Manipulate():
    def __init__(self, component , xs):
        self.xs = xs
        self.component = component
        self.plots = {}
        self.curves = {}
        # For now the default strategy is one plot per output. TODO: May want to have multiple curves on same plot for some components. 
        for name in self.component.outputs_names:
            out_func = self.component.outputs[name]
            plot = pg.PlotWidget(title = name)
            
            ys = out_func(xs)
            plot.setXRange(min(xs), max(xs), padding=0, update=True)
            plot.setYRange(min(ys), max(ys), padding=0, update=True) # 
            self.plots[name] =  plot
            self.curves[name] = plot.plot(xs,ys, pen=(255,0,0),fillLevel=0. )    
        
        self.setupWidget()
    
    def setupWidget(self):
        layout = QtGui.QGridLayout()
        i = 0
        j = 0
        N = len(self.plots)
        i0 = 4
        if N % 3 == 0:
            i0 = 3
        for name in self.component.outputs_names:
            plot = self.plots[name]
            layout.addWidget(plot,j,i)
            i += 1
            if i == i0:
                i = 0
                j += 1
                
        # Widgets 
        # Spin widgets
        
        def valueChanged_default(sb):
            self.update()
            pass
        
        spinlayout=QtGui.QGridLayout() 
        spins_widget=QtGui.QWidget()  
        spins_widget.setLayout(spinlayout)
        
        i = 0
        j = 0        
        funcs = {}
        for var in self.component.attribute_names:
            mini,maxi,start_value = self.component.attribute_parameters[var]                       
            curr_spin = pg.SpinBox(value=start_value, suffix='',step= (maxi - mini) / 100, minStep=(maxi - mini) / 10000)        
            curr_spin.setMinimum(mini)    
            curr_spin.setMaximum(maxi)

            funcs[var] = lambda sb,value, var = var : self.component.setAttribute(var,value)
            curr_spin.sigValueChanging.connect(funcs[var])
            curr_spin.sigValueChanged.connect(valueChanged_default)  
            
            label = QtGui.QLabel(var)
            spinlayout.addWidget(label,i,0+2*j)
            spinlayout.addWidget(curr_spin,i,1+2*j)
            i += 1
            if i == 4:
                i = 0
                j += 1
            layout.addWidget(spins_widget,2,0)   
        
        self.layout = layout
    
    def update(self):       
        for name,out_func in self.component.outputs.items():
            ys = out_func(self.xs)
            self.plots[name].setXRange(min(self.xs), max(self.xs), padding=0, update=True)
            self.plots[name].setYRange(min(ys), max(ys), padding=0, update=True)
            self.curves[name].setData(self.xs,ys )
        

def main():
    print "Starting main"
    
    win = QtGui.QMainWindow()    
    
    cw = QtGui.QWidget()
    
    #component = DoubleBusRing()
    component = SingleBusCROW()
    lbda_step = 0.00001
    lbdas = np.arange(1.5455,1.5485+lbda_step,lbda_step)
    manip = Manipulate(component, lbdas)
    
    cw.setLayout(manip.layout)
    win.setCentralWidget(cw)
    win.show()    
    
    win.setCentralWidget(cw)
   
    pg.setConfigOptions(antialias=True)
    win.setGeometry(100,100,1200,600)
    win.setWindowTitle('Analysor')   
         
    import sys
    pg.setConfigOptions(useWeave=False)
    print "Hello"
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        print "Running"
        QtGui.QApplication.instance().exec_()
        print "Done"
if __name__ == '__main__':
    print "Launching..."
    main()


    