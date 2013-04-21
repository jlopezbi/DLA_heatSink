# testing grow region

import rhinoscriptsyntax as rs
import heapq
import random
import math
import time
import Rhino
import scriptcontext
#import System.Guid

def growRegion():
	filter = Rhino.DocObjects.ObjectType.Mesh
	rc, objRef = Rhino.Input.RhinoGet.GetOneObject("select testMesh",False,filter)
	if not objRef or rc!=Rhino.Commands.Result.Success: return rc
	mesh  = objRef.Mesh()
	if not mesh: return
	mesh.Compact()

	randIdx = int(random.uniform(0,mesh.Vertices.Count-1))
	tVertIdxRoot = mesh.TopologyVertices.TopologyVertexIndex(randIdx)
	vertPnt = mesh.Vertices[randIdx]
	rad = .1
	rs.AddSphere(vertPnt,rad)

	growVerts = []

	stepSize = .01
	maxGrowLen = .5
	minGrowLen = .02
	cutoffDist = .7
	gKernel = GKernel(stepSize,maxGrowLen,minGrowLen,cutoffDist)
	gKernel.plot()


	conVertsIdx = mesh.Vertices.GetConnectedVertices(randIdx)
	print type(conVertsIdx)
	#conVertsIdx is an Array[int] in .NET framework.
	print str(conVertsIdx.Length)
	for i in range(conVertsIdx.Length):
		idx = conVertsIdx[i]
		if(idx != randIdx):
			tVertIdx = mesh.TopologyVertices.TopologyVertexIndex(idx)
			dist = lenBetweenTVerts(tVertIdxRoot,tVertIdx,mesh)
			lookUpIdx = int(round(dist/stepSize))
			distStr = "d:%1.2f,i:%d"%(dist,lookUpIdx)
			rs.AddTextDot(distStr, mesh.Vertices[idx])
			if(dist<cutoffDist):
				growVerts.append([idx,lookUpIdx])
		else:
			growVerts.append([idx,0])

	"""GROW REGION"""
	for i in range(len(growVerts)):
		vertIdx = growVerts[i][0]
		kernelIdx = growVerts[i][1]
		growLength = gKernel.gaussKernel[kernelIdx]

		vert = mesh.Vertices[vertIdx]
		vertNormal = mesh.Normals[vertIdx]
		growVec = vertNormal.Multiply(vertNormal,growLength)
		newLoc = rs.VectorAdd(vert,growVec)
		
		#normalArrow = rs.AddLine(vert,newLoc)
		#rs.CurveArrows(normalArrow,2)
		
		mesh.Vertices.SetVertex(vertIdx,newLoc.X,newLoc.Y,newLoc.Z)
		scriptcontext.doc.Objects.Replace(objRef, mesh)			


def assignGrowLens(vertIdx,gaussKernel):
	pass

def lenBetweenTVerts(tVertIdx1,tVertIdx2,mesh):
		p1 = mesh.TopologyVertices[tVertIdx1]
		p2 = mesh.TopologyVertices[tVertIdx2]
		return p1.DistanceTo(p2)

class GKernel:
	gaussKernel = []
	def __init__(self,stepSize,maxGrowLen,minGrowLen,cutoffDist):
		self.stepSize = stepSize
		self.maxGrowLen = maxGrowLen
		self.minGrowLen = minGrowLen
		self.cutoffDist = cutoffDist
		self.gaussKernel = self.createGaussKernel()


	def createGaussKernel(self):
		stepSize = self.stepSize
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

			x +=stepSize
		return gaussKernel

	def plot(self):
		for i in range(len(self.gaussKernel)):
			x = i*self.stepSize
			y = self.gaussKernel[i]
			rs.AddPoint(x,y,0)



if __name__=="__main__":
	growRegion()

