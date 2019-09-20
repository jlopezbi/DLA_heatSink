import Rhino
import scriptcontext
import System.Guid

def sharing_function(mesh, polyp, neighbor_level):
    '''computes how much the given polyp should grow by
    Args:
        mesh (rhinocommon mesh)
        polyp (int): vert index in mesh that is to grow
        neighbor_level (int): level of polyp from polyp-that-ate; acts as proxy for distance
    '''
    #TODO: write a linear func of neighbor level to get grow_length for the polyp
    #get function (or perhaps pass in)
    #evaluate function to get grow_length
    pass

class Coral:
  def __init__(self,objRef,mesh,ratioMaxMin):
    mesh.Normals.ComputeNormals()
    mesh.Normals.UnitizeNormals()
    mesh.Compact()
    #mesh.Weld(math.pi)

    self.objRef = objRef
    self.mesh = mesh
    self.groupIdx = scriptcontext.doc.Groups.Add()
    if(self.groupIdx == -1): print "lineage group failed to intialize"

    self.avgEdgeLen = self.getAvgEdgeLen()
    self.maxEdgeLength = self.avgEdgeLen
    self.minEdgeLen = self.avgEdgeLen*ratioMaxMin

    arrNakedBool = self.mesh.GetNakedEdgePointStatus()
    self.nakedVerts = set()
    for i in range(arrNakedBool.Length):
      if (arrNakedBool[i] == True):
        self.nakedVerts.add(i)


    self.subdivideLongEdges()
    self.collapseShortEdges()

    self.lineage = []
    self.lineage.append(self.mesh.Duplicate())

  def setNakedVerts(self):
    arrNakedBool = self.mesh.GetNakedEdgePointStatus()
    for i in range(arrNakedBool.Length):
      if (arrNakedBool[i] == True):
        self.nakedVerts.add(i)

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

  def get_neighbors(self,vertex,nLevels):
      '''
      Returns:
          list[list[int]] where int is vertex index
      '''
      pass

  def getGrowData(self,gKernel,centerVerts):
    '''
    Args:
        gKernel [?]
        centerVerts (list[int]): verts that 'ate'
    '''
    mesh = self.mesh
    gStepSize = gKernel.gStepSize
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

      #first neihbors
      conVertsIdx = mesh.Vertices.GetConnectedVertices(idxCenter)

      #for each neighbor:
          # compute its distance from the tVertIdxRoot (center vert)
          # so need to come up with concept of distance for an arbitrary vert.
            # shortest path across mesh seems valid
            # even simpler: discrete, based on level of neighbor. breaks down a bit if/when edge
            # lengths start varying more. but I like doing the simple thing first.
          # add to growVerts the (vert_idx, growLen)


      neighbors = get_neighbor_verts(start_vert=idxCenter, mesh=mesh, n_levels=1)
      for i, neighborhood in enumerate(neighbors):
          for vert in neighborhood:
              #TODO: plot the gaussKernel to make sure this use is correct
              grow_length = sharing_function(mesh, vert, i)
              growVerts.append([vert, grow_length])



      for i in range(conVertsIdx.Length):
        idxN = conVertsIdx[i]
        if(idxN != idxCenter and idxN not in centerVerts):
          tVertIdx = mesh.TopologyVertices.TopologyVertexIndex(idxN)
          dist = lenBetweenTVerts(tVertIdxRoot,tVertIdx,mesh)
          lookUpIdx = int(round(dist/gStepSize))

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
    nakedVerts = self.nakedVerts

    for i in range(len(growVerts)):

      vertIdx = growVerts[i][0]
      growLen = growVerts[i][1]
      if(vertIdx not in nakedVerts):
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
    #   print "collapsed and edge!"
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

  def saveToLineage(self,ts):
    print "saved copy at timeStep: %d" %ts
    #scriptcontext.doc.Groups.AddToGroup(self.groupIdx,)
    self.lineage.append(self.objRef.Mesh().Duplicate())


  def displayLineage(self):
    """OLD CODE"""
    groupIdx = scriptcontext.doc.Groups.Add()

    bboxMesh = self.mesh.GetBoundingBox(False)
    spacing = max(bboxMesh.Max.X-bboxMesh.Min.X,bboxMesh.Max.Y-bboxMesh.Min.Y)
    print "nSaved: %d" % len(self.lineage)

    for i in range(len(self.lineage)):
      tVec = Rhino.Geometry.Vector3d(0,i*spacing,0)
      bloop = self.lineage[i]
      print "type of lineage[%d]:" %i + str(type(bloop))
      #if(i!=0):
      bloop.Translate(tVec)
      meshGuid = scriptcontext.doc.Objects.AddMesh(bloop)
      if meshGuid == System.Guid.Empty:
        print "addMesh fail: %d" %i
      else:
        scriptcontext.doc.Groups.AddToGroup(groupIdx,meshGuid)
      #scriptcontext.doc.Objects.AddMesh(bloop)

  def packLineage(self,paramStr):
    groupIdx = scriptcontext.doc.Groups.Add()

    bboxMesh = self.mesh.GetBoundingBox(False)
    spacingVec = Rhino.Geometry.Vector3d(0,0,0)
    for i, mesh in enumerate(self.lineage):
      mesh.Translate(spacingVec)
      meshGuid = scriptcontext.doc.Objects.AddMesh(mesh)
      if meshGuid == System.Guid.Empty:
        print "addMesh failed: %d" %i
      else:
        scriptcontext.doc.Groups.AddToGroup(groupIdx,meshGuid)
      if(i!= len(self.lineage)-1):
        bboxCurr = mesh.GetBoundingBox(True)
        bboxNext = self.lineage[i+1].GetBoundingBox(True)
        bboxCurrLen = bboxCurr.Max.Y-bboxCurr.Min.Y
        bboxNextLen = bboxNext.Max.Y-bboxNext.Min.Y

        spacingVec.Y += (bboxCurrLen+bboxNextLen)/2.0+.000
      else:
        if not rs.SetUserText(meshGuid,"params",paramStr,True):
          print "SetUserText failed"
