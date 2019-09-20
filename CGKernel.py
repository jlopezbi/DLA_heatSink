import math
import rhinoscriptsyntax as rs

class GKernel:
  gaussKernel = []
  def __init__(self,gStepSize,maxGrowLen,minGrowLen,cutoffDist):
    self.gStepSize = gStepSize
    self.maxGrowLen = maxGrowLen
    self.minGrowLen = minGrowLen
    self.cutoffDist = cutoffDist
    self.gaussKernel = self.createGaussKernel()


  def createGaussKernel(self):
    gStepSize = self.gStepSize
    maxGrowLen = self.maxGrowLen
    minGrowLen = self.minGrowLen
    cutoffDist = self.cutoffDist

    def gaussFunc(x,a,b,c):
      return a*math.exp(-((x-b)**2.0)/(2.0*(c**2.0)))

    gaussKernel = []

    a = maxGrowLen
    b = 0
    c = math.sqrt(-(cutoffDist**2.0)/math.log(minGrowLen/maxGrowLen))

    x = 0
    y = maxGrowLen
    assert(y==gaussFunc(0,a,b,c)), "problems with gaussFunc"
    while(x<cutoffDist):
      y = gaussFunc(x,a,b,c)
      gaussKernel.append(y)

      x +=gStepSize
    return gaussKernel

  def plot(self):
    points = []
    for i in range(len(self.gaussKernel)):
      x = i*self.gStepSize
      y = self.gaussKernel[i]
      points.append([x,y,0])
      #rs.AddPoint(x,y,0)
    rs.AddPolyline(points)
