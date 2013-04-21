#DLA on a mesh
#http://wiki.mcneel.com/developer/rhinocommonsamples/closestpoint?s[]=rtree
#use R-tree to search for points close enough to vertices

import rhinoscriptsyntax as rs
import random
import math
import time
import Rhino
import scriptcontext
import System.Guid

def meshDLA_MAIN():

	filter = Rhino.DocObjects.ObjectType.Mesh
	rc, objRef = Rhino.Input.RhinoGet.GetOneObject("select seedMesh",False,filter)
	if not objRef or rc!=Rhino.Commands.Result.Success: return rc
	mesh  = objRef.Mesh()
	if not mesh: return
	mesh.Compact()
	
	filter = Rhino.DocObjects.ObjectType.PolysrfFilter
	rc, boxRef = Rhino.Input.RhinoGet.GetOneObject("select boundingBox (polySrf)",False,filter)
	if not boxRef or rc!=Rhino.Commands.Result.Success: return rc
		
	

	pRadius = .4
	speed = pRadius*.333
	nParticles =25
	timeSteps = 3000
	growLength = .07
	feedLength = .01
	maxEdgeLength = .6355
	minEdgeLen = .1
	thresMult = 1.3
	peakInclination = (4.0/5.0)*math.pi

	print "INPUT PARAMS________________________"
	print "pRadius = %1.2f in." % pRadius 
	print "nParticles = " +str(nParticles)
	print "timeSteps(ts) = " +str(timeSteps)
	print "speed = %0.2f in./ts" % speed
	print "growLength = %0.2f in." % growLength
	print "maxEdgeLength = %0.2f in." % maxEdgeLength
	print "minEdgeLen = %0.2f in." % minEdgeLen
	print "thresMult = " + str(thresMult)
	peakIncDeg = peakInclination*180.0/math.pi
	print "peakInclination = %1.1fdeg" % peakIncDeg
	print "____________________________________"

	world = World(boxRef,pRadius)
	world.getSpawnPlane()
	#world.draw()
	coral = Coral(objRef,mesh)

	particles = []
	growVerts = set()
	prevHighestPoint = 0;
	highestPoint = 0

	#INITIALIZE PARTICLES
	for i in range(nParticles):
		p = Particle(pRadius)
		p.setToSpawnLoc(world)
		p.drawParticle(i)
		particles.append(p)

	#RUN SIMIULATION
	for t in range(timeSteps):
		#time.sleep(0.01*10**-8)
		#Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()

		#MOVE PARTICLES
		for i in range(len(particles)):

			p = particles[i]
			p.moveParticle(speed,peakInclination)

			if not p.inBounds(world):
				p.setToSpawnLoc(world)
			p.clearParticle()
			p.drawParticle(i)

		#SEARCH FOR INTERSECTIONS
		growVerts = coral.verticesThatAte(world,particles)
		#growVerts = verticesThatAte(world,mesh,particles)

		#GROW THE VERTS THAT INTERSECTED
		for gVertIdx in growVerts:
			newPnt = coral.growVertice(gVertIdx,growLength)
			if(newPnt >highestPoint):
				highestPoint = newPnt

			coral.subdivideLongEdges(maxEdgeLength)
			coral.updateNormals()
			

		#COLLAPSE SHORT EDGES
		didCollapse = coral.collapseShortEdges(minEdgeLen)
		if(didCollapse):
			coral.updateNormals()

		#MOVE TOP OF BOUND BOX`
		if(highestPoint>prevHighestPoint):
			world.moveTop(highestPoint,pRadius,thresMult)
			world.reDraw()
			world.getSpawnPlane()
		prevHighestPoint = highestPoint

		coral.reDraw()

		

	#displayGrowNormals(mesh,growLength)
	for particle in particles:
		scriptcontext.doc.Objects.Hide(particle.sphereID,True)




#---------------------------CORAL----------------------------------
class Coral:
	def __init__(self,objRef,mesh):
		mesh.Normals.ComputeNormals()
		mesh.Normals.UnitizeNormals()
		mesh.Compact()

		self.objRef = objRef
		self.mesh = mesh

	def verticesThatAte(self, world, particles):
		mesh = self.mesh
		tree = Rhino.Geometry.RTree()

		#populate Rtree with mesh vertices
		for i,vertex in enumerate(mesh.Vertices):
			tree.Insert(vertex, i)

		#function that runs when a sphere intersects with a vert
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
			sphere = particle.sphere
			tree.Search(sphere, SearchCallback, sData)
			if(sData.addedVert):
				particle.setToSpawnLoc(world)
			sData.addedVert = False

		return sData.vertices


	def growVertice(self,idx,growLength):
		mesh = self.mesh
		vert = mesh.Vertices[idx]
		vertNormal = mesh.Normals[idx]
		growVec = vertNormal.Multiply(vertNormal,growLength)
		newLoc = rs.VectorAdd(vert,growVec)
		
		#normalArrow = rs.AddLine(vert,newLoc)
		#rs.CurveArrows(normalArrow,2)
		
		mesh.Vertices.SetVertex(idx,newLoc.X,newLoc.Y,newLoc.Z)
		return newLoc.Z
		#scriptcontext.doc.Objects.Replace(objRef, mesh)


	def growRegion(self,idxCenter,growLen,maxDist):

		mesh = self.mesh
		searchedVerts = []

		#???def searchNeighbors(idxCenter)

		def lenBetweenTVerts(tVertIdx1,tVertIdx2):
			p1 = mesh.TopologyVertices(tVertIdx1)
			p2 = mesh.TopologyVertices(tVertIdx2)
			return p1.DistanceTo(p2)


	def collapseShortEdges(self,minEdgeLen):
		mesh = self.mesh
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

	def subdivideLongEdges(self,maxEdgeLength):
		#iterate through all vertices of mesh and subdivide if too long. slower than
		#subdividLongNeighbors, but easier to write
		mesh = self.mesh
		edges = mesh.TopologyEdges
		nEdges = mesh.TopologyEdges.Count
		for i in range(nEdges):
			tVerts = mesh.TopologyEdges.GetTopologyVertices(i)
			p1 = mesh.TopologyVertices[tVerts.I]
			p2 = mesh.TopologyVertices[tVerts.J]
			lenEdge =  p1.DistanceTo(p2)
			if(lenEdge >= maxEdgeLength):
				mesh.TopologyEdges.SplitEdge(i,.5) 


	def subdivideLongNeighbors(self,idx,maxEdgeLength):
		mesh = self.mesh
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


	def displayGrowNormals(self,displayLength):
		mesh = self.mesh
		for i in range(mesh.Vertices.Count):
			vertNormal = mesh.Normals[i]
			feedVec = vertNormal.Multiply(vertNormal,displayLength)
			vert = mesh.Vertices[i]
			newLoc = rs.VectorAdd(vert,feedVec)
			feedLine = rs.AddLine(vert,newLoc)

	def displayVertices(self):
		mesh = self.mesh
		for i in range(mesh.Vertices.Count):
			vert = mesh.Vertices[i]
			rs.AddPoint(vert)

	def reDraw(self):
		scriptcontext.doc.Objects.Replace(self.objRef,self.mesh)

	def updateNormals(self):
		self.mesh.Normals.ComputeNormals()
		self.mesh.Normals.UnitizeNormals()
		self.mesh.Compact()

#-----------------------------------------------------------------------


def verticesThatAte(world,mesh,particles):
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
	
	#normalArrow = rs.AddLine(vert,newLoc)
	#rs.CurveArrows(normalArrow,2)
	
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

def displayGrowNormals(mesh,displayLength):
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

	def __init__(self,inputBoxRef,pRadius):
		self.inputBoxID = inputBoxRef.ObjectId
		rs.HideObject(self.inputBoxID)
	
		self.boundBoxBetter = inputBoxRef.Surface().GetBoundingBox(True)
		box = Rhino.Geometry.Box(self.boundBoxBetter)
	
		self.boxBrepID = scriptcontext.doc.Objects.AddBrep(box.ToBrep())
		self.corners = self.boundBoxBetter.GetCorners()
		self.lineIDs = []

		print "world created"
		print "world.boundBoxBetter type:" + str(type(self.boundBoxBetter))

	
	def getSpawnPlane(self):
		xMin = self.boundBoxBetter.Min.X
		xMax = self.boundBoxBetter.Max.X
		yMin = self.boundBoxBetter.Min.Y
		yMax = self.boundBoxBetter.Max.Y

		self.spawnXRange = [xMin,xMax]
		self.spawnYRange = [yMin,yMax]
		self.spawnZ = self.boundBoxBetter.Max.Z

	def moveTop(self,highestPoint,pRadius,thresMult):
		maxX = self.boundBoxBetter.Max.X
		maxY = self.boundBoxBetter.Max.Y
		maxZ = self.spawnZ

		dist = maxZ - highestPoint
		travelZone = pRadius*thresMult
		if(dist < travelZone):
			#print "highest point inside threshold area"
			newMaxZ = highestPoint+travelZone
			newPnt = Rhino.Geometry.Point3d(maxX,maxY,newMaxZ)
			self.boundBoxBetter.Max = newPnt


	def reDraw(self):
		box = Rhino.Geometry.Box(self.boundBoxBetter)
		scriptcontext.doc.Objects.Replace(self.boxBrepID, box.ToBrep())

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

	def moveParticle(self,speed,peakInclination):
		inclination = random.triangular(0,peakInclination,math.pi)
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












