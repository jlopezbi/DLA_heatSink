#DLA with branch thickening

import rhinoscriptsyntax as rs
import random as rand
import math
import time
import Rhino
import scriptcontext

class node:
    def __init__(self, posVec, parentID, radius):
        self.posVec = posVec
        self.parentID = parentID
        self.radius = radius
        self.geom = []

    def dispNodeInfo(self):
        self.geom.append(rs.AddPoint(self.posVec));
        #nodeInfoStr = str(self.parentID) + 'r:' + str(self.radius)
        #self.geom.append(rs.AddTextDot(nodeInfoStr,self.posVec))

    def drawEdge(self,nodes):
        p1 = self.posVec
        p2 = nodes[self.parentID].posVec
        self.geom.append(rs.AddLine(p1,p2))

        pNormal = rs.VectorSubtract(p2,p1)
        height = rs.VectorLength(pNormal)
        plane = rs.PlaneFromNormal(p1,pNormal)
        radius = self.radius
        self.geom.append(rs.AddCylinder(plane,height,radius))

    def drawCircle(self):
        p1 = self.posVec
        radius = self.radius
        self.geom.append(rs.AddCircle(p1,radius))

    def clearGeom(self):
        for o in self.geom: rs.DeleteObject(o)
        
    def increaseRadius(self,sFactor, nodes):
        self.clearGeom()
        self.radius += sFactor
            
        
        if(self.parentID >= 0):
            self.drawCircle()
            self.drawEdge(nodes)          
            parentNode = nodes[self.parentID]
            parentNode.increaseRadius(sFactor,nodes)
        elif(self.parentID == -1):
            rootNode = nodes[0]
            rootNode.radius += sFactor
            rootNode.clearGeom()
            rootNode.drawCircle()   
            #rootNode.dispNodeInfo()

        #recursively strengthen every node from current node to root


def branchingDLA_MAIN():
    #create root node
    nodes = [node([0,0,0],-1,.1)]
    rNodes = [rs.AddPoint(nodes[0].posVec)]
    speed = .8 #rhino units per iteration
    sRange = [.5,1.5]
    nThrows = 300
    boundRadius = (sRange[0]+sRange[1])*2
    bGeom = rs.AddCircle([0,0,0],boundRadius)
    sFactor = .5

    passParams = [boundRadius,bGeom,sFactor]

    

    for i in range(nThrows):
        resizeThreshold = .7
        stickRange = rand.uniform(sRange[0],sRange[1])
        passParams = throwNode(nodes,speed,stickRange,passParams,resizeThreshold)
        

        
def throwNode(nodes,speed,stickRange,passParams,resizeThreshold):
    boundRadius = passParams[0]
    bGeom = passParams[1]
    #sFactor = passParams[2]
    #print passParams[2]
    #assumes time steps of 1
    startPos = getBoundaryPos(boundRadius);
    endPos = getBoundaryPos(boundRadius);
    
    direction = rs.VectorSubtract(endPos,startPos)
    direction = rs.VectorUnitize(direction)
    vel = rs.VectorScale(direction,speed);

    currentPos = rs.VectorAdd(startPos,vel)

    previewGeom = []

    isTravelling = True
    while(isTravelling):
        scriptcontext.escape_test() #hit escape to quit NOT WORKING
        time.sleep(0.01*10**-7)
        for o in previewGeom: rs.DeleteObjects(o)
        dist = rs.VectorLength(currentPos) #distance to origin
        #check if particle went out of bounds
        if(dist>boundRadius):
            isTravelling = False
       
        else:
            previewGeom.append(drawPos(currentPos,stickRange))
            for i in range(len(nodes)):
                n = nodes[i]
                if(inStickRange(currentPos,n,stickRange)):
                    #GOT STUCK! add a new node at that position
                    newNode = node(currentPos,i,.08) #parent is idx of node stuck too
                    newNode.increaseRadius(.01, nodes)
                    nodes.append(newNode)
                    #rNodes.append(rs.AddPoint(currentPos))
                    if(math.fabs(boundRadius-dist) <= resizeThreshold):
                        #print "boundRadius should resize"
                        rs.DeleteObjects(bGeom)
                        boundRadius += resizeThreshold/2 #arbitrary
                        bGeom = rs.AddCircle([0,0,0],boundRadius)
                        passParams[0] = boundRadius
                        passParams[1] = bGeom

                    isTravelling = False
                    for o in previewGeom: rs.DeleteObjects(o)
                    break
            
            currentPos = rs.VectorAdd(currentPos,vel)
        Rhino.RhinoApp.Wait()
        
    return passParams

def drawPos(currentPos,stickRange):
    tempGeom = []
    tempGeom.append(rs.AddPoint(currentPos))
    tempGeom.append(rs.AddCircle(currentPos,stickRange))
    return tempGeom

def inStickRange(currentPos,node,stickRange):
    assert(stickRange>0), "stickRange must be greater than zero"
    nodePos = node.posVec
    distVec = rs.VectorSubtract(nodePos,currentPos)
    distance = rs.VectorLength(distVec)
    return (distance <= stickRange)

currentPos = rs.VectorCreate([1,0,0],[0,0,0])
node1 = node([0,0,0],0,10)
assert (inStickRange(currentPos,node1,1) == True), "inStickRange error"
assert (inStickRange(currentPos,node1,1.0000001) == True), "inStickRange error"
assert (inStickRange(currentPos,node1,.999999) == False), "inStickRange error"
    
def getBoundaryPos(radius):
    angle = rand.uniform(0,2.0*math.pi)
    x = radius*math.cos(angle)
    y = radius*math.sin(angle)
    z = 0;
    
    pos = rs.VectorCreate([0,0,0],[x,y,z])
    return pos
                                    
if __name__=="__main__":
    branchingDLA_MAIN()
