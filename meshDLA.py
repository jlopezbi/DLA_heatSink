#DLA on a mesh.
#http://wiki.mcneel.com/developer/rhinocommonsamples/closestpoint?s[]=rtree
#use R-tree to search for points close enough to vertices
# Josh Lopez-Binder
# Make a start mesh in rhino, select it, it will grow

import rhinoscriptsyntax as rs
import random, math, time
import Rhino
import scriptcontext
import System.Guid
import System.Drawing
import CParticle
import CWorld
import CCoral
import CGKernel

GKernel = CGKernel.GKernel
Coral = CCoral.Coral
#/bin/bash: s: command not found
Particle = CParticle.Particle

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

	areaMassProps = Rhino.Geometry.AreaMassProperties.Compute(mesh)
	centroid = areaMassProps.Centroid
	centroid.Z += 2.3
	rs.ViewTarget(view=None,target=centroid)

	showParticles = True
	showWorld = False
	debugTime = False
	saveLineage = False
	spin = False

	pRadius = .5
	stepRatio = 0.5
	speed = pRadius*stepRatio #relate to pRadius and some other len??
	nParticles = 50
	if debugTime:
		timeSteps = 50
	else:
		timeSteps = 1000
	ratioMaxMinEdgeLen = .33
	thresMult = 5
	alpha = 2
	beta = 6
	gravFactor = 1
	peakInclination = gravFactor*math.pi
	gStepSize = .01
	maxGrowLen = .25
	minGrowLen = .001
	cutoffDist = 1
	nSave = 10
	tsSave = timeSteps/nSave

	threshDist = pRadius*thresMult




	"""INITIALIZE WORLD, CORAL"""
	world = CWorld.World(mesh)
	if not showWorld: world.hide()
	coral = Coral(objRef,mesh,ratioMaxMinEdgeLen)
	world.resize(coral,threshDist)
	world.getSpawnPlane()
	world.reDraw()
	"""save param string"""
	paramStr = "\nPARAMS_______________"+\
			   "\npRadius: %1.2fin." % pRadius +\
			   "\nstepRatio: %1.2f"%stepRatio +\
			   "\nnParticles: %d" % nParticles +\
			   "\ntimeSteps: %d" % timeSteps +\
			   "\nmaxEdgeLen: %1.2f" % coral.maxEdgeLength +\
			   "\nminEdgeLen: %1.2f" % coral.minEdgeLen +\
			   " (max/minEdgeLen: %.2f)" % ratioMaxMinEdgeLen +\
			   "\nthreshMult: %1.1f" % thresMult +\
			   "\nmaxGrowLen: %.2fin." % maxGrowLen +\
			   "\nminGrowLen: %.2fin." % minGrowLen +\
			   "\ncutoffDist: %1.1fin." % cutoffDist +\
			   "\nused TrianglularDist" +\
			   "\ngravFactor: %.2f" % gravFactor +\
			   "\nnSave: %d" % nSave
	#"\nbetaDist: a=%d, b=%d" %(alpha,beta) +\

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



	gKernel = GKernel(gStepSize,maxGrowLen,minGrowLen,cutoffDist)
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
		if spin:
			rs.RotateView(angle = 1.0)
		ts += 1
		if(ts>=tsSave):
			#print "time:%d"%timeSteps
			if saveLineage: coral.saveToLineage(t)
			ts = 0

		#time.sleep(0.3)
		Rhino.RhinoApp.Wait()
		scriptcontext.escape_test()

		"""MOVE PARTICLES"""
		for i in range(len(particles)):
			p = particles[i]
			#boundChecks occur within moveParticle()
			p.moveParticle(speed,alpha,beta,peakInclination,world)
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
		coral.setNakedVerts()


		"""UPDATE CORAL AND WORLD """
		coral.reDraw()


		#coral.colorVerts(growVerts,gKernel)

		world.resize(coral,threshDist)
		world.getSpawnPlane()
		if showWorld: world.reDraw()

	#displayGrowNormals(mesh,growLength)
	if saveLineage: coral.packLineage(paramStr)
	if not rs.SetUserText(coral.objRef.ObjectId,"params",paramStr,True):
					print "SetUserText failed"


	# if not showParticles:
	# 	for particle in particles:
	# 		#rs.AddPolyline(particle.pnts)
	# 		scriptcontext.doc.Objects.Hide(particle.sphereID,True)
		#scriptcontext.doc.Views.Redraw()



if __name__=="__main__":
	reload(CParticle)
	reload(CWorld)
	reload(CCoral)
	reload(CGKernel)
	meshDLA_MAIN()
