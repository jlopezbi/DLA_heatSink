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

	radius = .2
	nParticles = 20
	for i in range(nParticles):
		#x = random.uniform(xMin,xMax)
		#y = random.uniform(yMin,yMax)
		#posVec = rs.VectorCreate([x,y,z],[0,0,0])
		p = Particle(radius)
		p.setToSpawnLoc(xRange,yRange,z)
		p.drawParticle(i)
		particles.append(p)

	nVerts = mesh.Vertices.Count
	#print "nParticles = " + str(len(particles))

	timeSteps = 100
	for t in range(timeSteps):
		#time.sleep(0.01*10**-7)
		#Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()
		for i in range(len(particles)):
			#print "particle#: " + str(i)
			p = particles[i]
			p.moveParticle(speed = .08)
			if not p.inBounds(boxID):
				p.setToSpawnLoc(xRange,yRange,z)
			stickIdx = p.didStick(mesh,nVerts)
			if stickIdx >=0:
				#rs.AddSphere(p.posVec,p.radius)
				growVertice(objRef,mesh,stickIdx,p.posVec)
				checkVertNeighborEdges(objRef,mesh,stickIdx,[.1,.7])
				p.setToSpawnLoc(xRange,yRange,z)				
			p.clearParticle()
			p.drawParticle(i)

def growVertice(objRef, mesh,idx,foodVec):
	#scriptcontext.doc.Objects.Delete(objRef,True)
	vert = mesh.Vertices[idx]
	newLoc = rs.VectorSubtract(vert,foodVec)
	mesh.Vertices.SetVertex(idx,foodVec.X,foodVec.Y,foodVec.Z)
	scriptcontext.doc.Objects.Replace(objRef, mesh)
	#print type(scriptcontext.doc.Objects)
	#scriptcontext.doc.Views.Redraw()
	# if scriptcontext.doc.Objects.AddMesh(mesh)!=System.Guid.Empty:
	# 	scriptcontext.doc.Views.Redraw()
	# 	return Rhino.Commands.Result.Success
	# return Rhino.Commands.Result.Failure

def checkVertNeighborEdges(objRef, mesh,idx,lengthRange):

	minLength = lengthRange[0]
	maxLength = lengthRange[1]

	vert = mesh.Vertices[idx]
	#rs.AddTextDot("v",vert)
	connectedVertsIdx = mesh.Vertices.GetConnectedVertices(idx)
	# for i in connectedVertsIdx:
	# 	if(i!=idx):
	# 		v = mesh.Vertices[i]
	# 		rs.AddSphere(v,.01)
	# print type(connectedVertsIdx)
	# print str(connectedVertsIdx)
	

	tVertIdx = mesh.TopologyVertices.TopologyVertexIndex(idx)
	tVert = mesh.TopologyVertices[tVertIdx]
	assert (tVert==vert), "topolgy vert and vert not the same!"
	vertsRepesent = mesh.TopologyVertices.MeshVertexIndices(tVertIdx)

	tConnectedVerts = mesh.TopologyVertices.ConnectedTopologyVertices(tVertIdx)

	#rs.AddPoint(tVert)
	# print "meshType in neighbor check: " + str(type(mesh))
	for neighIdx in connectedVertsIdx:
		if(neighIdx !=idx):
			neighVert = mesh.Vertices[neighIdx]
			dist = rs.Distance(vert,neighVert)

			distStr = "%.2f" %dist
			#textDot = rs.AddTextDot(distStr,neighVert)

			tNeighIdx = mesh.TopologyVertices.TopologyVertexIndex(neighIdx)
			edgeIdx = mesh.TopologyEdges.GetEdgeIndex(tVertIdx,tNeighIdx)

			if dist >= maxLength:
				mesh.TopologyEdges.SplitEdge(edgeIdx,0.5)
				scriptcontext.doc.Objects.Replace(objRef, mesh)
				#scriptcontext.doc.Objects.AddMesh(mesh)
				#rs.TextObjectStyle(textDot,Italic)
				print("long edge: %.2f" %dist)
				#find that edge subdivide it
				
			elif dist <= minLength:
				#find that edge collapse it

				print("short edge: %.2f" %dist)
				#rs.TextObjectStyle(textDot,Italic)






	



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

	def didStick(self,mesh,nVerts):
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












