#DLA on a mesh
#http://wiki.mcneel.com/developer/rhinocommonsamples/closestpoint?s[]=rtree
#use R-tree

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
	rc, boxRef = Rhino.Input.RhinoGet.GetOneObject("select boundingBox",False,filter)
	if not boxRef or rc!=Rhino.Commands.Result.Success: return rc
	boxID = scriptcontext.doc.Objects.Find(boxRef.ObjectId)
	boundBoxID = rs.BoundingBox(boxID.Geometry)

	world = World(boxRef)
	world.setSpawnPlane()

	particles = []

	radius = .05
	nParticles =25
	timeSteps = 2000
	growLength = .005
	feedLength = .01
	maxEdgeLength = .11
	minEdgeLen = .001

	for i in range(nParticles):
		p = Particle(radius)
		p.setToSpawnLoc(world)
		p.drawParticle(i)
		particles.append(p)
	
	for t in range(timeSteps):
		#time.sleep(0.01*10**-8)
		#Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()
		for i in range(len(particles)):

			p = particles[i]
			p.moveParticle(speed = .01)

			if not p.inBounds(world):
				p.setToSpawnLoc(world)
			p.clearParticle()
			p.drawParticle(i)

		growVerts = searchMesh(world,mesh,particles)

		for gVertIdx in growVerts:
			growVertice(objRef,mesh,gVertIdx,growLength)
			subdivideLongNeighbors(objRef,mesh,gVertIdx,[.05,.11])

		collapseShortEdges(mesh,minEdgeLen)



def searchMesh(world,mesh,particles):
	tree = Rhino.Geometry.RTree()

	for i,vertex in enumerate(mesh.Vertices):
		tree.Insert(vertex, i)

	def SearchCallback(sender, data):
		sData = data.Tag
		sData.vertices.add(data.Id)
		sData.addedVert = True

	class SearchData:

		def __init__(self):
			self.vertices = set()
			self.addedVert = False

	sData = SearchData()
			
	for particle in particles:
		x = particle.posVec.X
		y = particle.posVec.Y
		z = particle.posVec.Z
		pnt = Rhino.Geometry.Point3d(x,y,z)
		sphere = Rhino.Geometry.Sphere(pnt,particle.radius)
		tree.Search(sphere, SearchCallback, sData)
		if(sData.addedVert):
			particle.setToSpawnLoc(world)
		sData.addedVert = False
		


	return sData.vertices



def growVertice(objRef, mesh,idx,growLength):
	vert = mesh.Vertices[idx]
	vertNormal = mesh.Normals[idx]
	growVec = vertNormal.Multiply(vertNormal,growLength)
	newLoc = rs.VectorAdd(vert,growVec)
	"""
	normalArrow = rs.AddLine(vert,newLoc)
	rs.CurveArrows(normalArrow,2)
	"""
	#newLoc = Rhino.Geometry.Vector3f.Add()
	#print type(newLoc)
	mesh.Vertices.SetVertex(idx,newLoc.X,newLoc.Y,newLoc.Z)
	scriptcontext.doc.Objects.Replace(objRef, mesh)

def collapseShortEdges(mesh,minEdgeLen):
	collapsedAnEdge = False
	for i in range(mesh.TopologyEdges.Count):
		edgeLine = mesh.TopologyEdges.EdgeLine(i)
		length = edgeLine.Length
		if(length<minEdgeLen):
			mesh.TopologyEdges.CollapseEdge(i)
			collapsedAnEdge = True
	if collapsedAnEdge:
		scriptcontext.doc.Objects.Replace(objRef,mesh)




def subdivideLongNeighbors(objRef, mesh,idx,lengthRange):
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

	scriptcontext.doc.Objects.Replace(objRef,mesh)

#def collapseSmallEdges():

			#elif (dist<=minLength):

				#foundEdgeIdx = mesh.TopologyEdges.GetEdgeIndex(tCenterVertIdx,tNeighVertIdx)
				#mesh.TopologyEdges.CollapseEdge(foundEdgeIdx)
				#print "FuseEdge dist: " + strDist
				#rs.AddTextDot(strDist,tNeighVert)

			# foundEdge = mesh.TopologyEdges.GetEdgeIndex(tCenterVertIdx,tNeighVertIdx)
			# mesh.TopologyEdges.SplitEdge(foundEdge,.5)


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


class World:
	spawnXRange = None
	spawnYRange = None
	spawnZ = None

	def __init__(self,boxRef):
		self.boxID = boxRef.ObjectId
		self.boundBox = rs.BoundingBox(self.boxID)
		print "world created"
		print "world.boundBox type:" + str(type(self.boundBox))
		print "world.boxID type:" + str(type(self.boxID))

	
	def setSpawnPlane(self):
		xMin = self.boundBox[0].X
		xMax = self.boundBox[1].X
		self.spawnXRange = [xMin,xMax]

		yMin = self.boundBox[0].Y
		yMax = self.boundBox[2].Y
		self.spawnYRange = [yMin,yMax]

		self.spawnZ = self.boundBox[4].Z-.01


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
	

	def inBounds(self,world):
		pnt = self.geom[0]
		return rs.IsPointInSurface(world.boxID,pnt)

	def setToSpawnLoc(self,world):
		xRange = world.spawnXRange
		yRange = world.spawnYRange
		z = world.spawnZ

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












