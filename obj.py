class wellGrid:
    def __init__(self,nRows,nCols, dt, wellR, pipeR, shoeDepth, wellDepth, qBit,
                 kMud=1.2*0.5778,denMud=10, cpHatMud = 0.3824, tSurf=85,gtGradient=1.1/100):
        import numpy as np
        self.nRows = nRows
        self.nCols = nCols
        self.wellR = wellR
        self.shoeDepth = shoeDepth
        self.dr = wellR/(nCols-1)
        self.dz = wellDepth/(nRows-1)
        self.dt = dt
        self.tSurf = tSurf
        self.qBit = qBit
        self.kMud = kMud
        self.gtGradient = gtGradient
        self.matrix = np.zeros((nRows,nCols))
        self.velocities = self._calcVelocities(np.zeros((nRows,nCols)), self.dr, wellR, pipeR)
        idx = np.arange(0,nRows,1)*self.dz*gtGradient + tSurf
        self.matrix = np.repeat(idx.reshape(-1,1),nCols, axis=1)
        self.alpha = 3.7E-5 * kMud / (denMud * cpHatMud)
        self.lambda0 = ((self.dr)**2.)/(self.alpha*dt)
        self.lambda1, self.lambda2 = self._calcLambdas(self.velocities, self.dr, self.dz, self.alpha)
        
        

    def _calcVelocities(self, matrix, dr, wellR, pipeR):
        import numpy as np
        radii=np.arange(self.nCols)*dr
        dummyRow = np.zeros((self.nCols,))
        for i, r in enumerate(radii):
            if r<pipeR:
                dummyRow[i] = 0.32 * (1 - r/pipeR) ** 2.
            else:
                dummyRow[i] = 175.7 * (1 - (r / wellR) ** 2. - ((1 - (pipeR / wellR) ** 2.) / np.log(1 / (pipeR / wellR))) * np.log(wellR / r))
        return np.repeat(dummyRow.reshape(1,-1), self.nRows, axis=0)
        
        
    def _calcLambdas(self, velsMatrix, dr, dz, alpha):
        import numpy as np
        radii=np.arange(self.nCols)*dr
        lambda1 = np.zeros(self.nCols)
        lambda2 = np.zeros(self.nCols)
        for i, r in enumerate(radii):
            lambda1[i] = velsMatrix[1,i]*(dr/12)**2/(2*alpha*dz)
            if r==0:
                lambda2[i]==0
            else:
                lambda2[i] = dr/(2*r)
        
        lambda1 = np.repeat(lambda1.reshape(1,-1), self.nRows, axis=0)  
        lambda2 = np.repeat(lambda2.reshape(1,-1), self.nRows, axis=0)
        return lambda1, lambda2
    
    def calcTemperatures(self, timeFrame=300):
        import numpy as np
        for k in range(100):
            for i in range(self.nRows-1, -1, -1):
                for j in range(0, self.nCols, 1):
                    if j==0:
                        if i==0:
                            self.matrix[i,j]=self.tSurf
                        elif i==self.nRows-1:
                            self.matrix[i,j]=(1/self.lambda0+2)*(self.lambda0*self.matrix[i,j] + 2*self.matrix[i,j+1] + self.lambda1[i,j]*(2*self.dz*self.qBit/self.kMud))
                        else:
                            self.matrix[i,j]=(1/self.lambda0+2)*(self.lambda0*self.matrix[i,j] + 2*self.matrix[i,j+1] - self.lambda1[i,j]*(self.matrix[i-1,j]-self.matrix[i+1,j]))
                    elif j==self.nCols-1:
                        if i==0:
                            self.matrix[i,j]=self.tSurf
                        elif i==self.nRows-1:
                            dtInf = ((i*self.dz*self.gtGradient+self.tSurf)-self.matrix[i,j])/(2*self.wellR)
                            self.matrix[i,j]=(1/self.lambda0+2)*(self.lambda0*self.matrix[i,j] + 2*self.matrix[i,j-1] + 2*self.dr*dtInf*(1+self.lambda2[i,j]) + self.lambda1[i,j]*(2*self.dz*self.qBit/self.kMud))
                        elif i*self.dz>self.shoeDepth:
                            dtInf = ((i*self.dz*self.gtGradient+self.tSurf)-self.matrix[i,j])/(2*self.wellR)
                            self.matrix[i,j]=(1/self.lambda0+2)*(self.lambda0*self.matrix[i,j] + 2*self.matrix[i,j-1] + 2*self.dr*dtInf*(1+self.lambda2[i,j]) + self.lambda1[i,j]*(self.matrix[i-1,j]-self.matrix[i+1,j]))
                        else:
                            self.matrix[i,j]=(1/self.lambda0+2)*(self.lambda0*self.matrix[i,j] + 2*self.matrix[i,j-1] + self.lambda1[i,j]*(self.matrix[i-1,j]-self.matrix[i+1,j]))
                            

    
    
myWellGrid = wellGrid(nRows=20,nCols=34, dt=1, wellR=3.25, pipeR=2.5, shoeDepth=1000, wellDepth=3000, qBit=926)
myWellGrid.calcTemperatures()