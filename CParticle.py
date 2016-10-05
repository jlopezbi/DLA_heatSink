import rhinoscriptsyntax as rs
import Rhino
import scriptcontext
import random, math, time

class Particle:
  geom = None
  textDot = None
  sphereID = None

  def __init__(self, radius):
    #self.posVec = [0,0,0]
    point = Rhino.Geometry.Point3d(0,0,0)
    self.sphere = Rhino.Geometry.Sphere(point,radius)
    self.radius = radius
    self.pnts = []
    self.geom = [0,0] #idx 0 => point, idx 1 => sphere

  def moveParticle(self,speed,alpha,beta,peakInclination,world):
    """TRIANGULAR DISTRIBUTION"""
    inclination = random.triangular(0,math.pi+.00001,peakInclination)

    """BETA DISTRIBUTION"""
    """
    rand = 1-random.betavariate(alpha,beta)
    inclination = rand*math.pi
    """

    azimuth = random.uniform(0,math.pi*2.0)

    velX = speed*math.sin(azimuth)*math.cos(inclination)
    velY = speed*math.sin(azimuth)*math.sin(inclination)
    velZ = speed*math.cos(inclination) # <----- whattttt??? this error!?
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

    self.pnts.append(self.sphere.Center)
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
    '''
    appears to be unused
    '''
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
