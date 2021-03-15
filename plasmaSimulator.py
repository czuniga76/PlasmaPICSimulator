import numpy as np
import matplotlib.pyplot as plt

# Christian Zuniga, March 2021. Based on Master's thesis "Particle In Cell Simulation
# of Cusp-Confined Plasma With Monte Carlo Collisions".
# Substituted some functions with Python functions from Numpy and Scipy

# define constants
pi = np.pi
u0 = 4*pi*10**(-7)  # magnetic permeability free space (henrys/ m)
etome = 1.7588196*10**11  # electron |q|/m  TODO use modified charge distribution
md = 1591.05350		#e*e/2 me e0  
miontoelec = 256

cr = 461341     # check

### 03/25/2019 TODO check movement 1 particle without E field, only magnetic field
# obtained oscillations



class PlasmaChamber():
    def __init__(self):
        return None


         
    
class Cusp(PlasmaChamber):
    def __init__(self,I=600,Lx=0.01,Ly=0.01):
        self.LX = Lx		# Lx, Ly are the size of rectangle for simulation (m)
        self.LY = Ly
        self.I = I      # Amperes current in wires
        self.bconst1 = u0*self.I/(2*pi)
        self.Nx = 100	# number of grid points
        self.Ny = 100
        self.dx = Lx/self.Nx
        self.dy = Ly/self.Ny
    
    def getMagneticFieldAtPosition(self,x,y):   # TODO check if extension to case x,y arrays
        r = x**2 + y**2
        Bx = self.bconst1*(-y/r + y/((x-self.LX)**2 + y**2))      # Magnetic field of 2 wires in z, with current +z for wire at 0, -z wire at Lx
        By = self.bconst1*(x/r + (self.LX-x)/((x-self.LX)**2 + y**2))
        return (Bx,By)   # TODO check if it's better to return matrix. Yes, return 2D array or save as internal
    
    def getMagneticRotMatrixAtPosition(self,P,qtomdt): # receive position matrices and charge/mass ratio x dt/2
        
        # TODO Check if B field values are sensible
        numCols = 3*P.shape[1]
        Htot = np.zeros((3,numCols)) # check
        #print " Htot shape ", Htot.shape
        # ### for p in Position:
        #     ##x,y = p[0],p[1]
        Bx,By = self.getMagneticFieldAtPosition(P[0,:],P[1,:])
        print (" magnetic field Bx, By ", Bx, By)
        #     Form H and add to Htot
        a = qtomdt*By
        b = qtomdt*Bx
        a2 = a**2
        b2 = b**2
        
        #print " a2 , b2 shapes ", a2.shape, b2.shape
        onesM = np.ones(a.shape)
        S = a2 + b2  + onesM 
        #print " S shape ", S.shape
        H11 = (b2 - a2 + onesM)/S
        H12 = 2*a*b/S
        H13 = -2*a/S
        H21 = H12
        H22 =  (onesM + a2 - b2)/S
        H23 = 2 *b/S
        H31 = -H13
        H32 = -H23
        H33 = (onesM -a2 - b2)/S
        Htot[0,0:numCols:3] = H11
        Htot[0,1:numCols:3] = H12
        Htot[0,2:numCols:3] = H13
        Htot[1,0:numCols:3] = H21
        Htot[1,1:numCols:3] = H22
        Htot[1,2:numCols:3] = H23
        Htot[2,0:numCols:3] = H31
        Htot[2,1:numCols:3] = H32
        Htot[2,2:numCols:3] = H33
        
        #print Htot
        #print H.shape
        return Htot # return matrix of  H matrices, size 3, 3Np
        
    def calculateElectricField(self,chargeDistribution): # TODO solve Poisson's Equation on grid
        return None
    
    def checkBoundaries(self):
        #TODO go through each particle
        # check crossing on bottom, sides, and top
        # if bottom, replace partical randomly, others, reflect particle?
        return None
        
class PICSimulation():
    def __init__(self,T=30*(10**-9),Np=1000):
        self.Np = Np    # number of particles Np electrons, Np ions
        self.simTime = T   # simulation time (sec)
        self.dt = 10**(-7) #(-11)   # timestep
        self.Nt = int(self.simTime/self.dt)
        print(" Number of time steps ", self.Nt, " time step (s) ", self.dt, " time ", self.simTime)
        self.I = 600   # current for B field (seems excessive)
        self.Lx = 0.01  # width (m)
        self.Ly = 0.01  # height simulation area (m)
        self.qtomdte = -etome*self.dt/2
        self.qtomdti = etome*self.dt/(2*miontoelec)
        self.cusp = Cusp(self.I,self.Lx,self.Ly)   # init chamber
        self.plasma = Plasma(self.Np,self.cusp.Nx,self.cusp.Ny) # init plasma
        self.trajectory = np.zeros((self.Nt,2)) # electron trajectories
        self.trajectoryi = np.zeros((self.Nt,2)) # ion trajectories
        radiusL = self.plasma.vthermale/(1.76*2.4*10**8)
        print (" sample radius ", radiusL)
        periodL = 2*pi*radiusL/self.plasma.vthermale
        print (" sample period (s)" , periodL)
        
        dx = self.cusp.dx
        dy = self.cusp.dy
        self.trajectory[0,:] = self.plasma.loceMatrix[0,:]
        self.trajectoryi[0,:] = self.plasma.lociMatrix[0,:]
        # rotate velocities by half time step
        #HMe = self.cusp.getMagneticRotMatrixAtPosition(dx*self.plasma.loceMatrix,-self.qtomdte)
        #HMi = self.cusp.getMagneticRotMatrixAtPosition(dx*self.plasma.lociMatrix,-self.qtomdti)
        #HMe = HMe.reshape(self.Np,3,3)  # TODO check if this op can be done earlier
        #HMi = HMi.reshape(self.Np,3,3)
        # once for elec, once for ions
        Ve_neg = self.plasma.veleMatrix    
        Vi_neg = self.plasma.veliMatrix
        #self.plasma.veleMatrix = np.sum(np.transpose(HMe,(0,2,1)).reshape(self.Np,3,3,1)*Ve_neg.reshape(self.Np,3,1,1),-3)
        #self.plasma.veliMatrix = np.sum(np.transpose(HMi,(0,2,1)).reshape(self.Np,3,3,1)*Vi_neg.reshape(self.Np,3,1,1),-3)
        #self.plasma.veleMatrix = self.plasma.veleMatrix.reshape(3,self.Np)
        #self.plasma.veliMatrix = self.plasma.veliMatrix.reshape(3,self.Np) 
        # Add diagnostic variables for potential, sample trajectories
        teMatrix = np.zeros((Np,3))
        seMatrix = np.zeros((Np,3))
        tiMatrix = np.zeros((Np,3))
        siMatrix = np.zeros((Np,3))
        
        Bxe,Bye =  self.cusp.getMagneticFieldAtPosition(dx*self.plasma.loceMatrix[:,0],dy*self.plasma.loceMatrix[:,1])
        Bxi,Byi =  self.cusp.getMagneticFieldAtPosition(dx*self.plasma.lociMatrix[:,0],dy*self.plasma.lociMatrix[:,1])
        teMatrix[:,0] = -self.qtomdte*Bxe
        teMatrix[:,1] = -self.qtomdte*Bye
        tsq = np.sum(teMatrix*teMatrix,axis=1)
        dene = np.ones((Np,)) + tsq
        dene = dene.reshape(Np,1)
        seMatrix = 2*teMatrix/dene
        
        tiMatrix[:,0] = -self.qtomdti*Bxi
        tiMatrix[:,1] = -self.qtomdti*Byi
        tsqi = np.sum(tiMatrix*tiMatrix,axis=1)
        deni = np.ones((Np,)) + tsqi
        deni = deni.reshape(Np,1)
        siMatrix = 2*tiMatrix/deni
        
        velepMatrix = self.plasma.veleMatrix + np.cross(self.plasma.veleMatrix,teMatrix)
        velePosMatrix = self.plasma.veleMatrix + np.cross(velepMatrix,seMatrix)
        self.plasma.veleMatrix = velePosMatrix
        velipMatrix = self.plasma.veliMatrix + np.cross(self.plasma.veliMatrix,tiMatrix)
        veliPosMatrix = self.plasma.veliMatrix + np.cross(velipMatrix,siMatrix)
        self.plasma.veliMatrix = veliPosMatrix
        self.vesq = []  # proxy for kinetic energy
        
        
    def runSimulation(self):
        # TODO
        # Calculate electric field
        # Advance particles to new positions with half acceleration
        # Rotate velocities according to B field
        # Advance particles another half step
        # Check if Boundaries have been crossed
        # Add MC collisions
        # Store diagnostics
        k1e = -etome*self.dt/2
        k1i = k1e/miontoelec
        dx = self.cusp.dx
        dy = self.cusp.dy
        Np = self.Np
        teMatrix = np.zeros((Np,3))
        seMatrix = np.zeros((Np,3))
        tiMatrix = np.zeros((Np,3))
        siMatrix = np.zeros((Np,3))
        for k in range(1,self.Nt):
            # advance particles
            #print " vele matrix shape ", self.plasma.veleMatrix.shape
            self.plasma.loceMatrix += self.dt*self.plasma.veleMatrix[:,:2]
            self.plasma.lociMatrix += self.dt*self.plasma.veliMatrix[:,:2]
            #print self.plasma.loceMatrix
            # store sample trajectory
            self.trajectory[k,:]  = self.plasma.loceMatrix[0,:]  # .reshape(1,2)
            self.trajectoryi[k,:] = self.plasma.lociMatrix[0,:]
            # calculate V- with half acceleration. Only Vx,Vy impacted by Efield
            #HM = np.random.rand(self.Np,3,3) # TODO generate HM from particle positions
            #HMe = self.cusp.getMagneticRotMatrixAtPosition(dx*self.plasma.loceMatrix,self.qtomdte)
            #HMi = self.cusp.getMagneticRotMatrixAtPosition(dx*self.plasma.lociMatrix,self.qtomdti)
            #HMe = HMe.reshape(self.Np,3,3)  # TODO check if this op can be done earlier
            #HMi = HMi.reshape(self.Np,3,3)
            #print " H at " ,k , HMe
            # once for elec, once for ions
            Ve_neg = self.plasma.veleMatrix
            self.vesq.append(np.sum(Ve_neg[0]**2))
            #print "vel before H ", np.sum(Ve_neg[0]**2)
            #Ve_neg[:2,:] += k1e*self.plasma.electricField
            Vi_neg = self.plasma.veliMatrix
            #Vi_neg[:2,:] += k1i*self.plasma.electricField
            # Vectorization thanks to https://jameshensman.wordpress.com/2010/06/14/multiple-matrix-multiplication-in-numpy/
            #Ve_pos = np.sum(np.transpose(HMe,(0,2,1)).reshape(self.Np,3,3,1)*Ve_neg.reshape(self.Np,3,1,1),-3)
            #Vi_pos = np.sum(np.transpose(HMi,(0,2,1)).reshape(self.Np,3,3,1)*Vi_neg.reshape(self.Np,3,1,1),-3)
            # TODO first skip electric field effect, just magnetifc
            # calculate new velocities by half acceleration
            #### test cross product instead
            Bxe,Bye =  self.cusp.getMagneticFieldAtPosition(dx*self.plasma.loceMatrix[:,0],dy*self.plasma.loceMatrix[:,1])
            Bxi,Byi =  self.cusp.getMagneticFieldAtPosition(dx*self.plasma.lociMatrix[:,0],dy*self.plasma.lociMatrix[:,1])
            teMatrix[:,0] = -self.qtomdte*Bxe
            teMatrix[:,1] = -self.qtomdte*Bye
            tsq = np.sum(teMatrix*teMatrix,axis=1)
            dene = np.ones((Np,)) + tsq
            dene = dene.reshape(Np,1)
            seMatrix = 2*teMatrix/dene
        
            tiMatrix[:,0] = -self.qtomdti*Bxi
            tiMatrix[:,1] = -self.qtomdti*Byi
            tsqi = np.sum(tiMatrix*tiMatrix,axis=1)
            deni = np.ones((Np,)) + tsqi
            deni = deni.reshape(Np,1)
            siMatrix = 2*tiMatrix/deni
            
            velepMatrix = self.plasma.veleMatrix + np.cross(self.plasma.veleMatrix,teMatrix)
            velePosMatrix = self.plasma.veleMatrix + np.cross(velepMatrix,seMatrix)
            self.plasma.veleMatrix = velePosMatrix
            velipMatrix = self.plasma.veliMatrix + np.cross(self.plasma.veliMatrix,tiMatrix)
            veliPosMatrix = self.plasma.veliMatrix + np.cross(velipMatrix,siMatrix)
            self.plasma.veliMatrix = veliPosMatrix
            
            #self.plasma.veleMatrix = Ve_pos #+ k1e*self.plasma.electricField # TODO calc electric field for elec and for ions
            #self.plasma.veliMatrix = Vi_pos #+ k1i*self.plasma.electricField
            
            #self.plasma.veleMatrix = self.plasma.veleMatrix.reshape(3,self.Np)
            #print ("vel after H "), np.sum(self.plasma.veleMatrix[0]**2)
            #self.vesq.append(np.sum(self.plasma.veleMatrix[0]**2))
            #self.plasma.veliMatrix = self.plasma.veliMatrix.reshape(3,self.Np) 
            
            self.cusp.checkBoundaries
            
            # TODO  calculate electric fields from new charge positions
    
   
    def plotTrajectory(self):
        print (" Number of time steps ", self.Nt)
        #plt.plot(self.trajectory[:,0]) # ,self.trajectory[:,1])
        #plt.figure()
        #plt.plot(self.trajectory[:,1])
        #plt.figure()
        print ("initial e pos ", self.cusp.dx*self.trajectory[0,:])
        print (" final e pos ", self.cusp.dx*self.trajectory[-1,:])
        print ("initial i pos ", self.cusp.dx*self.trajectoryi[0,:])
        print (" final i pos ", self.cusp.dx*self.trajectoryi[-1,:])
        plt.plot(self.cusp.dx*self.trajectory[:,0],self.cusp.dy*self.trajectory[:,1])
        plt.figure()
        plt.plot(self.cusp.dx*self.trajectoryi[:,0],self.cusp.dy*self.trajectoryi[:,1])
        #plt.plot(self.vesq)
        plt.show()
   
    def runSimulationCollsions(self):
        # same as above but with collisions
        return None         
########

class Plasma():
    def __init__(self,Np=1000,Nx=256,Ny=256,vthermale=343400,vthermali=42775): 
        # thermal vel in m/s
        self.Np = Np
        print( "Nx ", Nx)
        self.vthermale = vthermale
        self.vthermali = vthermali
        self.loadParticles(Np,Nx,Ny)
        self.electricField = np.zeros((2,Np))
        
  
    
    def loadParticles(self,Np,Nx,Ny):
        # initialize velocities from a Maxwell distribution
        # initialize positions to be in upper quarter of array
        vthermale = self.vthermale
        vthermali = self.vthermali
        self.veliMatrix = np.zeros((Np,3))
        self.veleMatrix = np.zeros((Np,3))
        self.loceMatrix = np.zeros((Np,2))
        self.lociMatrix = np.zeros((Np,2))
        print ("elec vel and pos shape ", self.veleMatrix.shape, self.loceMatrix.shape)
        sqt2 = np.sqrt(2)
        
        for i in range(3):
            ve = sqt2*vthermale*np.random.normal(size=Np)
            vi = sqt2*vthermali*np.random.normal(size=Np)
            self.veleMatrix[:,i] = ve
            self.veliMatrix[:,i] = vi
        if self.veleMatrix[0,1] > 0:
            self.veleMatrix[0,1] = -self.veleMatrix[0,1]
        if self.veleMatrix[0,0] > 0:
            self.veleMatrix[0,0] = -self.veleMatrix[0,0]
        self.veleMatrix[0,2] = 0.0
        self.loceMatrix[:,0] = np.random.randint(0,Nx,size=Np)
        self.loceMatrix[:,1] = np.random.randint(0.75*Ny,Ny,size=Np)
        
        self.lociMatrix[:,0] = np.random.randint(0,Nx,size=Np)
        self.lociMatrix[:,1] = np.random.randint(0.75*Ny,Ny,size=Np)
        
   
    
   
 
def main():
    print("test sim")
    
    sim = PICSimulation(T=(10**-3),Np=1)
    sim.runSimulation()
    sim.plotTrajectory()
    
if __name__ == "__main__":
    main()
    