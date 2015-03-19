import Rhino
import scriptcontext

class World():
  spawnXRange = None
  spawnYRange = None
  spawnZ = None
  lineIDs = None

  def __init__(self,mesh):
      
    self.boundBoxBetter = mesh.GetBoundingBox(True)
  
    #scriptcontext.doc.Objects.AddBrep(Rhino.Geometry.Box(self.boundBoxBetter).ToBrep())

    box = Rhino.Geometry.Box(self.boundBoxBetter)
    
    self.boxBrepID = scriptcontext.doc.Objects.AddBrep(box.ToBrep())

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


