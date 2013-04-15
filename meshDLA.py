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
	if not mesh: return
	mesh.Compact()
	mesh.Normals.ComputeNormals()
	#mesh.Normals.Flip()
	mesh.Normals.UnitizeNormals()
	

	
	filter = Rhino.DocObjects.ObjectType.PolysrfFilter
	rc, boxRef = Rhino.Input.RhinoGet.GetOneObject("select boundingBox",False,filter)
	if not boxRef or rc!=Rhino.Commands.Result.Success: return rc
		
	

	pRadius = .05
	nParticles =25
	timeSteps = 100
	growLength = .001
	feedLength = .01
	maxEdgeLength = .07
	minEdgeLen = .023
	sideMult = 2
	topMult = 4

	world = World(boxRef,pRadius,sideMult,topMult)
	world.getSpawnPlane()
	#world.draw()

	particles = []
	growVerts = set()
	highestPoint = 0

	for i in range(nParticles):
		p = Particle(pRadius)
		p.setToSpawnLoc(world)
		p.drawParticle(i)
		particles.append(p)
	
	for t in range(timeSteps):
		#time.sleep(0.01*10**-8)
		#Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()

		#MOVE PARTICLES
		for i in range(len(particles)):

			p = particles[i]
			p.moveParticle(speed = .01)

			if not p.inBounds(world):
				p.setToSpawnLoc(world)
			p.clearParticle()
			p.drawParticle(i)

		#SEARCH FOR INTERSECTIONS
		growVerts = searchMesh(world,mesh,particles)

		#GROW THE VERTS THAT INTERSECTED
		for gVertIdx in growVerts:
			newPnt = growVertice(objRef,mesh,gVertIdx,growLength)
			if(newPnt >highestPoint):
				highestPoint = newPnt
			subdivideLongNeighbors(objRef,mesh,gVertIdx,maxEdgeLength)
			mesh.Normals.ComputeNormals()
			mesh.Normals.UnitizeNormals()

		#COLLAPSE SHORT EDGES
		didCollapse = collapseShortEdges(mesh,minEdgeLen)
		if(didCollapse):
			mesh.Normals.ComputeNormals()
			mesh.Normals.UnitizeNormals()

		#MOVE TOP OF BOUND BOX


		scriptcontext.doc.Objects.Replace(objRef, mesh)
		

	displayFeedNormals(mesh,growLength)



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
		sphere = particle.sphere
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
	
	normalArrow = rs.AddLine(vert,newLoc)
	rs.CurveArrows(normalArrow,2)
	
	#newLoc = Rhino.Geometry.Vector3f.Add()
	#print type(newLoc)
	mesh.Vertices.SetVertex(idx,newLoc.X,newLoc.Y,newLoc.Z)
	return newLoc.Z
	#scriptcontext.doc.Objects.Replace(objRef, mesh)

def collapseShortEdges(mesh,minEdgeLen):
	collapsedAnEdge = False
	for i in range(mesh.TopologyEdges.Count):
		edgeLine = mesh.TopologyEdges.EdgeLine(i)
		length = edgeLine.Length
		if(length<minEdgeLen):
			mesh.TopologyEdges.CollapseEdge(i)
			collapsedAnEdge = True
	if collapsedAnEdge:
		print "collapsed and edge!"
	return collapsedAnEdge
		#scriptcontext.doc.Objects.Replace(objRef,mesh)




def subdivideLongNeighbors(objRef, mesh,idx,maxEdgeLength):
	#minLength = lengthRange[0]
	#maxLength = lengthRange[1]

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

			if(dist >= maxEdgeLength):

				foundEdgeIdx = mesh.TopologyEdges.GetEdgeIndex(tCenterVertIdx,tNeighVertIdx)
				mesh.TopologyEdges.SplitEdge(foundEdgeIdx,.5)

	#scriptcontext.doc.Objects.Replace(objRef,mesh)


def displayFeedNormals(mesh,displayLength):
	for i in range(mesh.Vertices.Count):
		vertNormal = mesh.Normals[i]
		feedVec = vertNormal.Multiply(vertNormal,displayLength)
		vert = mesh.Vertices[i]
		newLoc = rs.VectorAdd(vert,feedVec)
		feedLine = rs.AddLine(vert,newLoc)

def displayVertices(mesh):
	for i in range(mesh.Vertices.Count):
		vert = mesh.Vertices[i]
		rs.AddPoint(vert)






#---------------------------WORLD----------------------------------
class World:
	spawnXRange = None
	spawnYRange = None
	spawnZ = None
	corners = None
	lineIDs = None

	def __init__(self,inputBoxRef,pRadius,sideMult,topMult):
		self.inputBoxID = inputBoxRef.ObjectId
		#self.inputBoxRef = inputBoxRef.Surface().GetBoundingBox(True)
		# self.boundBox = rs.BoundingBox(self.boxID)

		# minX = self.boundBox[0].X
		# maxX = self.boundBox[1].X
		# minY = self.boundBox[0].Y
		# maxY = self.boundBox[2].Y
		# minZ = self.boundBox[0].Z
		# maxZ = self.boundBox[4].Z

		rs.HideObject(self.inputBoxID)
		#self.boundBoxBetter = Rhino.Geometry.BoundingBox(minX,minY,minZ,maxX,maxY,maxZ)
		self.boundBoxBetter = inputBoxRef.Surface().GetBoundingBox(True)
		self.box = Rhino.Geometry.Box(self.boundBoxBetter)
	
		self.boxBrepID = scriptcontext.doc.Objects.AddBrep(self.box.ToBrep())
		#self.boundBoxBetter.Inflate(pRadius*sideMult,pRadius*sideMult,pRadius*topMult)
		self.corners = self.boundBoxBetter.GetCorners()
		self.lineIDs = []

		print "world created"
		#print "world.boundBox type:" + str(type(self.boundBox))
		#print "world.boxID type:" + str(type(self.boxID))
		print "world.boundBoxBetter type:" + str(type(self.boundBoxBetter))

	
	def getSpawnPlane(self):
		xMin = self.boundBoxBetter.Min.X
		xMax = self.boundBoxBetter.Max.X
		yMin = self.boundBoxBetter.Min.Y
		yMax = self.boundBoxBetter.Max.Y

		self.spawnXRange = [xMin,xMax]
		self.spawnYRange = [yMin,yMax]
		self.spawnZ = self.boundBoxBetter.Max.Z

	def reDraw(self):
		boxBrep = self.box.ToBrep()
		self.boxBrepID = scriptcontext.doc.Objects.Replace(self.boxBrepID, boxBrep)


		
			




#---------------------------Particle--------------------------------
class Particle:
	geom = None
	textDot = None
	sphereID = None

	def __init__(self, radius):
		self.posVec = [0,0,0]
		point = Rhino.Geometry.Point3d(0,0,0)
		self.sphere = Rhino.Geometry.Sphere(point,radius)
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
		self.sphere.Translate(vel)
		#rs.MoveObjects(self.geom,vel)
		#rs.MoveObject(obj,vel)
	

	def inBounds(self,world):
		#pnt = self.geom[0]
		#return rs.IsPointInSurface(world.inputBoxID,pnt)
		pnt = self.sphere.Center
		bbox = world.boundBoxBetter
		return bbox.Contains(pnt)

	def setToSpawnLoc(self,world):
		xRange = world.spawnXRange
		yRange = world.spawnYRange
		z = world.spawnZ

		x = random.uniform(xRange[0],xRange[1])
		y = random.uniform(yRange[0],yRange[1])
		self.posVec =  rs.VectorCreate([x,y,z],[0,0,0])
		spawnPnt = Rhino.Geometry.Point3d(self.posVec)
		self.sphere.Center = spawnPnt

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

		#self.geom[0] = rs.AddPoint(self.posVec)
		#elf.geom[1] = rs.AddSphere(self.posVec,self.radius)
		if(self.sphereID):
			scriptcontext.doc.Objects.Delete(self.sphereID,False)

		self.sphereID = scriptcontext.doc.Objects.AddSphere(self.sphere)


	def clearParticle(self):
		#rs.DeleteObjects(self.geom)
		return



if __name__=="__main__":
	meshDLA_MAIN()












