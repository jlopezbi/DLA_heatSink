# create gaussian kernel
import rhinoscriptsyntax as rs
import random
import math
import time
import Rhino
import scriptcontext
import System.Guid

def plotGaussKernel():
	stepSize = .01
	nIter = 120
	maxGrowLen = 2
	minGrowLen = .001
	cutoffDist = 3
	s = 1.0/(maxGrowLen*math.sqrt(2.0*math.pi)) # calculate s so that gauss(x=0)=maxGrowLen
	u = 0

	gaussKernel = createGaussKernel(stepSize,maxGrowLen,minGrowLen,cutoffDist)
	for i in range(len(gaussKernel)):
		x = i*stepSize
		y = gaussKernel[i]
		rs.AddPoint(x,y,0)



def createGaussKernel(stepSize,maxGrowLen,minGrowLen,cutoffDist):

	def gaussFunc(x,a,b,c):
		return a*math.exp(-((x-b)**2.0)/(2.0*(c**2.0)))

	gaussKernel = []

	a = maxGrowLen
	b = 0
	c = math.sqrt(-(cutoffDist**2.0)/math.log(minGrowLen/maxGrowLen))

	x = 0
	y = maxGrowLen
	assert(y==gaussFunc(0,a,b,c))
	while(x<cutoffDist):
		y = gaussFunc(x,a,b,c)
		gaussKernel.append(y)

		x +=stepSize
	return gaussKernel

plotGaussKernel()