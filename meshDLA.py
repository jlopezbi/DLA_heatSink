#DLA on a mesh

import rhinoscriptsyntax as rs
import random
import math
import time
import Rhino
import scriptcontext
import System.Guid

def meshDLA_MAIN():
	#mesh = rs.GetObject("select seed mesh", filter=32, preselect=True,)
	filter = Rhino.DocObjects.ObjectType.Mesh
	rc, objRef = Rhino.Input.RhinoGet.GetOneObject("select seedMesh",False,filter)
	if not objRef or rc!=Rhino.Commands.Result.Success: return rc
	print "BEGIN_________________________"
	print "type objRef: " + str(type(objRef))
	mesh  = objRef.Mesh()
	mesh.Compact()
	mesh.Normals.ComputeNormals()
	#mesh.Normals.Flip()
	mesh.Normals.UnitizeNormals()
	nVerts = mesh.Vertices.Count
	if not mesh: return

	# seedMeshID = scriptcontext.doc.Objects.Find(objRef.ObjectId)
	# seedMesh = seedMeshID.Geometry
	# seedVerts = seedMesh.Vertices

	filter = Rhino.DocObjects.ObjectType.PolysrfFilter
	rc, boundBoxRef = Rhino.Input.RhinoGet.GetOneObject("select boundingBox",False,filter)
	if not boundBoxRef or rc!=Rhino.Commands.Result.Success: return rc
	boxID = scriptcontext.doc.Objects.Find(boundBoxRef.ObjectId)
	boundBox = rs.BoundingBox(boxID.Geometry)

	# pnt1 = rs.GetObject("select test vert")
	# if not pnt1: return
	# isInside = rs.IsPointInSurface(boundBoxID,pnt1)
	# print "isInside: " + str(isInside)

	# if boundBox:
	# 	for i, point in enumerate(boundBox):
	# 		rs.AddTextDot(i,point)

	particles = []

	xMin = boundBox[0].X
	xMax = boundBox[1].X
	xRange = [boundBox[0].X,boundBox[1].X]
	yMin = boundBox[0].Y
	yMax = boundBox[2].Y
	yRange = [boundBox[0].Y,boundBox[2].Y]
	z = boundBox[4].Z

	radius = .05
	nParticles =20
	timeSteps = 1000
	growLength = .001
	feedLength = .01
	for i in range(nParticles):
		#x = random.uniform(xMin,xMax)
		#y = random.uniform(yMin,yMax)
		#posVec = rs.VectorCreate([x,y,z],[0,0,0])
		p = Particle(radius)
		p.setToSpawnLoc(xRange,yRange,z)
		p.drawParticle(i)
		particles.append(p)


	
	#print "nParticles = " + str(len(particles))

	
	for t in range(timeSteps):
		#time.sleep(0.01*10**-7)
		#Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()
		for i in range(len(particles)):

			p = particles[i]
			p.moveParticle(speed = .01)

			if not p.inBounds(boxID):
				p.setToSpawnLoc(xRange,yRange,z)

			stickIdx = p.didStick(mesh)
			if stickIdx >=0:
				growVertice(objRef,mesh,stickIdx,p.posVec,.006)
				checkVertNeighborEdges(objRef,mesh,stickIdx,[.05,.11])
				mesh.Normals.ComputeNormals()
				mesh.Normals.UnitizeNormals()
				#displayFeedNormals(mesh,feedLength)
				#displayVertices(mesh)
				p.setToSpawnLoc(xRange,yRange,z)
				
				
						
			p.clearParticle()
			p.drawParticle(i)

def growVertice(objRef, mesh,idx,foodVec,growLength):
	vert = mesh.Vertices[idx]
	vertNormal = mesh.Normals[idx]
	growVec = vertNormal.Multiply(vertNormal,growLength)
	newLoc = rs.VectorAdd(vert,growVec)
	normalArrow = rs.AddLine(vert,newLoc)
	rs.CurveArrows(normalArrow,2)
	#newLoc = Rhino.Geometry.Vector3f.Add()
	#print type(newLoc)
	mesh.Vertices.SetVertex(idx,newLoc.X,newLoc.Y,newLoc.Z)
	scriptcontext.doc.Objects.Replace(objRef, mesh)



def checkVertNeighborEdges(objRef, mesh,idx,lengthRange):
	minLength = lengthRange[0]
	maxLength = lengthRange[1]

	vert = mesh.Vertices[idx]
	#rs.AddTextDot("v",vert)
	connectedVertsIdx = mesh.Vertices.GetConnectedVertices(idx)	
	

	tVertIdx = mesh.TopologyVertices.TopologyVertexIndex(idx)
	tVert = mesh.TopologyVertices[tVertIdx]
	assert (tVert==vert), "topolgy vert and vert not the same!"


	for neighVertIdx in connectedVertsIdx:
		if(neighVertIdx !=idx):

			tCenterVertIdx = mesh.TopologyVertices.TopologyVertexIndex(idx)
			tCenterVert = mesh.TopologyVertices[tCenterVertIdx]

			tNeighVertIdx = mesh.TopologyVertices.TopologyVertexIndex(neighVertIdx)
			tNeighVert = mesh.TopologyVertices[tNeighVertIdx]
			#rs.AddSphere(tNeighVert,r)
			dist = rs.Distance(tCenterVert,tNeighVert)
			strDist = "%.2f" % dist
			if(dist >= maxLength):

				foundEdgeIdx = mesh.TopologyEdges.GetEdgeIndex(tCenterVertIdx,tNeighVertIdx)
				mesh.TopologyEdges.SplitEdge(foundEdgeIdx,.5)

			elif (dist<=minLength):

				foundEdgeIdx = mesh.TopologyEdges.GetEdgeIndex(tCenterVertIdx,tNeighVertIdx)
				mesh.TopologyEdges.CollapseEdge(foundEdgeIdx)
				#print "FuseEdge dist: " + strDist
				#rs.AddTextDot(strDist,tNeighVert)

			# foundEdge = mesh.TopologyEdges.GetEdgeIndex(tCenterVertIdx,tNeighVertIdx)
			# mesh.TopologyEdges.SplitEdge(foundEdge,.5)
	scriptcontext.doc.Objects.Replace(objRef,mesh)


def displayFeedNormals(mesh,feedLength):
	for i in range(mesh.Vertices.Count):
		vertNormal = mesh.Normals[i]
		feedVec = vertNormal.Multiply(vertNormal,feedLength)
		vert = mesh.Vertices[i]
		newLoc = rs.VectorAdd(vert,feedVec)
		feedLine = rs.AddLine(vert,newLoc)

def displayVertices(mesh):
	for i in range(mesh.Vertices.Count):
		vert = mesh.Vertices[i]
		rs.AddPoint(vert)

	



class Particle:
	geom = None
	textDot = None
	def __init__(self, radius):
		self.posVec = [0,0,0]
		self.radius = radius
		self.geom = [0,0] #idx 0 => point, idx 1 => sphere

	def moveParticle(self,speed):
		inclination = random.triangular(0,math.pi+.001,math.pi)
		azimuth = random.uniform(0,math.pi*2.0)

		velX = speed*math.sin(azimuth)*math.cos(inclination)
		velY = speed*math.sin(azimuth)*math.sin(inclination)
		velZ = speed*math.cos(inclination)
		vel = rs.VectorCreate([velX,velY,velZ],[0,0,0])

		self.posVec  = rs.VectorAdd(self.posVec,vel)
		#rs.MoveObjects(self.geom,vel)
		#rs.MoveObject(obj,vel)
	

	def inBounds(self,boundBoxID):
		pnt = self.geom[0]
		return rs.IsPointInSurface(boundBoxID,pnt)

	def setToSpawnLoc(self,xRange,yRange,z):
		x = random.uniform(xRange[0],xRange[1])
		y = random.uniform(yRange[0],yRange[1])
		self.posVec =  rs.VectorCreate([x,y,z],[0,0,0])

	def didStick(self,mesh):
		nVerts = mesh.Vertices.Count
		stickIdx = -1
		for i in range(nVerts):
			vert = mesh.Vertices[i]
			dist = rs.Distance(self.posVec,vert)
			#take the first stick distance found
			if(dist<=self.radius):
				stickIdx = i
				return stickIdx
		return stickIdx
				


	def drawParticle(self,i):
		#print "geom[0]:" + str(self.geom[0])
		self.geom[0] = rs.AddPoint(self.posVec)
		self.geom[1] = rs.AddSphere(self.posVec,self.radius)
		"""
		if self.textDot:
			rs.TextDotPoint(self.textDot,self.posVec)
		else:
			self.textDot = rs.AddTextDot(i,self.posVec)
		"""

	def clearParticle(self):
		rs.DeleteObjects(self.geom)


if __name__=="__main__":
	meshDLA_MAIN()












