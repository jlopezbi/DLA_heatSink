#DLA on a mesh
#http://wiki.mcneel.com/developer/rhinocommonsamples/closestpoint?s[]=rtree
#use R-tree to search for points close enough to vertices

import rhinoscriptsyntax as rs
import random, math, time
import Rhino
import scriptcontext
import System.Guid
import System.Drawing

def meshDLA_MAIN():
	
	"""GET MESH"""
	filter = Rhino.DocObjects.ObjectType.Mesh
	rc, objRef = Rhino.Input.RhinoGet.GetOneObject("select seedMesh",False,filter)
	if not objRef or rc!=Rhino.Commands.Result.Success: return rc
	mesh  = objRef.Mesh()
	if not mesh: return
	mesh.Compact()

	"""GET POLYSRF BOX
	filter = Rhino.DocObjects.ObjectType.PolysrfFilter
	rc, boxRef = Rhino.Input.RhinoGet.GetOneObject("select boundingBox (polySrf)",False,filter)
	if not boxRef or rc!=Rhino.Commands.Result.Success: return rc
	#scriptcontext.doc.Objects.Hide(boxRef,True)
	"""


	pRadius = .3
	speed = .05
	nParticles = 25
	timeSteps = 300
	ratioMaxMin = .33
	thresMult = 1.3
	peakInclination = (.6)*math.pi
	stepSize = .01
	maxGrowLen = .06
	minGrowLen = .001
	cutoffDist = 1
	tsSave = 50
	showParticles = True
	showWorld = False


	threshDist = pRadius*1.4

	"""INITIALIZE WORLD, CORAL"""
	world = World(mesh)
	if not showWorld: world.hide()
	coral = Coral(objRef,mesh,ratioMaxMin)
	world.resize(coral,threshDist)
	world.getSpawnPlane()
	world.reDraw()

	print "INPUT PARAMS________________________"
	print "pRadius = %1.2fin." % pRadius 
	print "nParticles = " +str(nParticles)
	peakIncDeg = peakInclination*180.0/math.pi
	print "peakInclination = %1.1fdeg" % peakIncDeg
	print "timeSteps(ts) = " +str(timeSteps)
	print "speed = %0.2f in./ts" % speed
	print "	maxGrowLen = %1.2fin." % maxGrowLen
	print "	minGrowLen = %1.2fin." % minGrowLen
	print "	cutoffDist = %1.2fin." % cutoffDist

	print "maxEdgeLength(Avg) = %0.2f in." % coral.maxEdgeLength
	print "minEdgeLen = %0.2f in." % coral.minEdgeLen
	print "thresMult = " + str(thresMult)
	print "____________________________________"

	
	gKernel = GKernel(stepSize,maxGrowLen,minGrowLen,cutoffDist)
	#gKernel.plot()

	particles = []
	growVerts = set()
	prevHighestPoint = 0;
	highestPoint = 0

	"""INITIALIZE PARTICLES"""
	for i in range(nParticles):
		p = Particle(pRadius)
		p.setToSpawnLoc(world)
		if showParticles: p.drawParticle(i)
		particles.append(p)

	"""RUN SIMIULATION"""
	ts = 0 
	for t in range(timeSteps):
		# ts += 1
		# if(ts>=tsSave):
		# 	coral.saveCopy()
		# 	ts = 0

		#time.sleep(0.3)
		Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()

		"""MOVE PARTICLES"""
		for i in range(len(particles)):
			p = particles[i]
			#boundChecks occur within moveParticle()
			p.moveParticle(speed,peakInclination,world)
			if showParticles: p.clearParticle(), p.drawParticle(i)

		"""SEARCH FOR INTERSECTIONS"""		
		centerVerts = coral.verticesThatAte(world,particles)

		"""GROW REGIONS AROUND CENTERVERTS"""
		growVerts = coral.getGrowData(gKernel,centerVerts)
		coral.grow(growVerts) 

		"""SUBDIVIDE LONG EDGES"""
		coral.subdivideLongEdges()

		"""COLLAPSE SHORT EDGES"""
		didCollapse = coral.collapseShortEdges()

		"""UPDATE MESH"""
		coral.updateNormals()
		coral.mesh.Weld(math.pi)
		

		"""UPDATE CORAL AND WORLD """
		coral.reDraw()

		
		#coral.colorVerts(growVerts,gKernel)

		world.resize(coral,threshDist)
		world.getSpawnPlane()
		if showWorld: world.reDraw()
		
	rs.AddTextDot("done", (0,0,0))
	#displayGrowNormals(mesh,growLength)
	#..coral.displayLineage()

	# for particle in particles:
	# 	scriptcontext.doc.Objects.Hide(particle.sphereID,True)
	# scriptcontext.doc.Views.Redraw()

#---------------------------CORAL----------------------------------
class Coral:
	def __init__(self,objRef,mesh,ratioMaxMin):
		mesh.Normals.ComputeNormals()
		mesh.Normals.UnitizeNormals()
		mesh.Compact()
		#mesh.Weld(math.pi)

		self.objRef = objRef
		self.mesh = mesh

		self.avgEdgeLen = self.getAvgEdgeLen()
		self.maxEdgeLength = self.avgEdgeLen
		self.minEdgeLen = self.avgEdgeLen*ratioMaxMin

		self.subdivideLongEdges()
		self.collapseShortEdges()

		self.lineage = []
		self.lineage.append(self.mesh.Duplicate())

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



	def getGrowData(self,gKernel,centerVerts):
		mesh = self.mesh
		stepSize = gKernel.stepSize
		kernelLen = len(gKernel.gaussKernel)

		#tVertIdxRoot = mesh.TopologyVertices.TopologyVertexIndex(idxCenter)
		#conVertsIdx = mesh.Vertices.GetConnectedVertices(idxCenter)
		growVerts = []

		def lenBetweenTVerts(tVertIdx1,tVertIdx2,mesh):
			p1 = mesh.TopologyVertices[tVertIdx1]
			p2 = mesh.TopologyVertices[tVertIdx2]
			return p1.DistanceTo(p2)

		centerColor = System.Drawing.Color.FromArgb(164,223,45)

		for idxCenter in centerVerts:
			tVertIdxRoot = mesh.TopologyVertices.TopologyVertexIndex(idxCenter)
			conVertsIdx = mesh.Vertices.GetConnectedVertices(idxCenter)

			for i in range(conVertsIdx.Length):
				idxN = conVertsIdx[i]
				if(idxN != idxCenter):
					tVertIdx = mesh.TopologyVertices.TopologyVertexIndex(idxN)
					dist = lenBetweenTVerts(tVertIdxRoot,tVertIdx,mesh)
					lookUpIdx = int(round(dist/stepSize))
					
					#distStr = "d:%1.2f,i:%d"%(dist,lookUpIdx)
					#rs.AddTextDot(distStr, mesh.Vertices[idx])
					if(lookUpIdx<kernelLen):
						growLen = gKernel.gaussKernel[lookUpIdx]
						growVerts.append([idxN,growLen])
				elif(idxN == idxCenter):
					#mesh.VertexColors.SetColor(idx,centerColor)
					growVerts.append([idxCenter,gKernel.maxGrowLen])

		return growVerts

	def grow(self,growVerts):
		mesh = self.mesh
		for i in range(len(growVerts)):
			vertIdx = growVerts[i][0]
			growLen = growVerts[i][1]

			vert = mesh.Vertices[vertIdx]
			vertNormal = mesh.Normals[vertIdx]
			growVec = vertNormal.Multiply(vertNormal,growLen)
			newLoc = Rhino.Geometry.Point3d.Add(vert,growVec)
						

			mesh.Vertices.SetVertex(vertIdx,newLoc)

	def colorVerts(self,growVerts,gKernel):
		meshID = self.objRef.ObjectId

		gaussMax = gKernel.maxGrowLen
		gaussMin = gKernel.minGrowLen

		rMax = 44
		gMax = 255
		bMax = 50

		colors = [None]*self.mesh.Vertices.Count
		for i in range(len(colors)): 
			colors[i] = [100,100,100]

		if self.mesh.Vertices.Count != rs.MeshVertexCount(meshID):
			print "ERROR: mesh.vertices.Count != rs.meshVertesCount!"
		if growVerts:
			for i in range(len(growVerts)):
				vertIdx = growVerts[i][0]
				kernelIdx = growVerts[i][1]
				growLen = gKernel.gaussKernel[kernelIdx]
				ratio = (growLen-gaussMin)/(gaussMax-gaussMin)
				if(ratio >1):
					print "ratio: %1.2f" % ratio

				r = ratio*rMax
				g = ratio*gMax
				b = ratio*bMax

				colors[i] = [r,g,b]



		rs.MeshVertexColors( meshID, colors )
	
	def collapseShortEdges(self):
		mesh = self.mesh
		minEdgeLen = self.minEdgeLen
		collapsedAnEdge = False
		for i in range(mesh.TopologyEdges.Count):
			# edgeLine = mesh.TopologyEdges.EdgeLine(i)
			# length = edgeLine.Length
			length = self.getLenEdge(i)
			if(length<minEdgeLen):
				mesh.TopologyEdges.CollapseEdge(i)
				collapsedAnEdge = True
		# if collapsedAnEdge:
		# 	print "collapsed and edge!"
		return collapsedAnEdge

	def subdivideLongEdges(self):
		#iterate through all vertices of mesh and subdivide if too long. slower than
		#subdividLongNeighbors, but easier to write
		mesh = self.mesh
		maxEdgeLength = self.maxEdgeLength
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
		"""OLD CODE"""
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

	def getAvgEdgeLen(self):
		mesh = self.mesh
		totLen = 0
		for i in range(mesh.TopologyEdges.Count):
			totLen += self.getLenEdge(i)
		self.avgEdgeLen = totLen/mesh.TopologyEdges.Count
		return self.avgEdgeLen

	def getLenEdge(self, edgeIdx):
		mesh = self.mesh
		edgeLine = mesh.TopologyEdges.EdgeLine(edgeIdx)
		return edgeLine.Length 



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

	def saveCopy(self):
		#scriptcontext.doc.Objects.AddMesh(self.mesh)
		print "saved copy"
		self.lineage.append(self.mesh.Duplicate())


	def displayLineage(self):
		bboxMesh = self.mesh.GetBoundingBox(False)
		spacing = max(bboxMesh.Max.X-bboxMesh.Min.X,bboxMesh.Max.Y-bboxMesh.Min.Y)
		print "nSaved: %d" % len(self.lineage)

		for i in range(len(self.lineage)):
			tVec = Rhino.Geometry.Vector3d(0,i*spacing,0)
			bloop = self.lineage[i]
			print "type of lineage[%d]:" %i + str(type(bloop))
			#if(i!=0):
			bloop.Translate(tVec)
			if scriptcontext.doc.Objects.AddMesh(bloop) == System.Guid.Empty:
				print "addMesh fail: %d" %i
			#scriptcontext.doc.Objects.AddMesh(bloop)
#---------------------------GAUSS KERNEL----------------------------
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
#---------------------------WORLD----------------------------------
class World:
	spawnXRange = None
	spawnYRange = None
	spawnZ = None
	lineIDs = None

	def __init__(self,mesh):
			
		self.boundBoxBetter = mesh.GetBoundingBox(True)
		minZ = self.boundBoxBetter.Min.Z
		minX = self.boundBoxBetter.Min.X
		minY = self.boundBoxBetter.Min.Y
		newMinPnt = Rhino.Geometry.Point3d(minX,minY,minZ-.01)
		self.boundBoxBetter.Min = newMinPnt
		#scriptcontext.doc.Objects.AddBrep(Rhino.Geometry.Box(self.boundBoxBetter).ToBrep())

		box = Rhino.Geometry.Box(self.boundBoxBetter)
		
		self.boxBrepID = scriptcontext.doc.Objects.AddBrep(box.ToBrep())
		#scriptcontext.doc.Objects.Hide(self.boxBrepID,True)
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


	def resize(self,coral,threshDist):
		mesh = coral.objRef.Mesh() # was not updating bounding box when using coral.mesh
		bboxCoral = mesh.GetBoundingBox(True)
		#scriptcontext.doc.Objects.AddBrep(bboxCoral.ToBrep())
		if not bboxCoral.IsValid:
			print "bbox mesh not Valid!"

		coralMin = bboxCoral.Min
		coralMax = bboxCoral.Max
		worldMin = self.boundBoxBetter.Min
		worldMax = self.boundBoxBetter.Max

		minOffset = Rhino.Geometry.Vector3d(threshDist,threshDist,0)
		maxOffset = Rhino.Geometry.Vector3d(threshDist,threshDist,threshDist)

		threshMin = worldMin + minOffset
		threshMax = worldMax - maxOffset
		
		newWorldMin = worldMin + (coralMin-threshMin)
		newWorldMax = worldMax + (coralMax-threshMax)
		#print "nwMax:" + str(type(newWorldMax))
		self.boundBoxBetter.Min = newWorldMin
		self.boundBoxBetter.Max = newWorldMax



	def checkIntersect(self,mesh):
		boxBrepID = self.boxBrepID
		pass


	def reDraw(self):
		box = Rhino.Geometry.Box(self.boundBoxBetter)
		scriptcontext.doc.Objects.Replace(self.boxBrepID, box.ToBrep())

	def hide(self):
		scriptcontext.doc.Objects.Hide(self.boxBrepID,True)

#---------------------------PARTICLE--------------------------------
class Particle:
	geom = None
	textDot = None
	sphereID = None

	def __init__(self, radius):
		#self.posVec = [0,0,0]
		point = Rhino.Geometry.Point3d(0,0,0)
		self.sphere = Rhino.Geometry.Sphere(point,radius)
		self.radius = radius
		self.geom = [0,0] #idx 0 => point, idx 1 => sphere

	def moveParticle(self,speed,peakInclination,world):
		inclination = random.triangular(0,peakInclination,math.pi)
		azimuth = random.uniform(0,math.pi*2.0)

		velX = speed*math.sin(azimuth)*math.cos(inclination)
		velY = speed*math.sin(azimuth)*math.sin(inclination)
		velZ = speed*math.cos(inclination)
		vel = rs.VectorCreate([velX,velY,velZ],[0,0,0])

		#self.posVec  = rs.VectorAdd(self.posVec,vel)
		#if(self.posVec.X < world.boundBoxBetter.Min.X
		pntPos = Rhino.Geometry.Point3d.Add(self.sphere.Center,vel)
		if(pntPos.Z < world.boundBoxBetter.Min.Z or pntPos.Z > world.boundBoxBetter.Max.Z):
			self.setToSpawnLoc(world)
			return
		elif(pntPos.X < world.boundBoxBetter.Min.X):
			pntPos.X = world.boundBoxBetter.Max.X
		elif(pntPos.X>world.boundBoxBetter.Max.X):
			pntPos.X = world.boundBoxBetter.Min.X
		elif(pntPos.Y < world.boundBoxBetter.Min.Y):
			pntPos.Y = world.boundBoxBetter.Max.Y
		elif(pntPos.Y>world.boundBoxBetter.Max.Y):
			pntPos.Y = world.boundBoxBetter.Min.Y

		self.sphere.Center = pntPos
		#self.sphere.Translate(vel)


	

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












