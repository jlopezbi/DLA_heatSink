#DLA with branch thickening

import rhinoscriptsyntax as rs
import random as rand
import math
import time
import Rhino

class node:
    def __init__(self, posVec, parentID, radius):
        self.posVec = posVec
        self.parentID = parentID
        self.radius = radius

    def drawNode(self):
        rs.AddPoint(self.posVec);
        rs.AddTextDot(self.parentID,self.posVec);
        
    def increaseRadius(self, nodes):
        
        self.radius = self.radius + 1
        if(self.parentID >= 0):
            parentNode = nodes[self.parentID]
            parentNode.increaseRadius(nodes)

def branchingDLA_MAIN():
    #nodes = [ ]
    #create root node
    nodes = [node([0,0,0],0,1)]
    rNodes = [rs.AddPoint(nodes[0].posVec)]
    speed = .5 #rhino units per iteration
    stickRange = 1.5
    nThrows = 300
    boundRadius = 10.0

    rs.AddCircle([0,0,0],boundRadius)

    for i in range(nThrows):
        throwNode(nodes,speed,stickRange,rNodes,boundRadius)
        

def randPntsTest(nodes,nThrows):
    for i in range(nThrows):
        x = rand.uniform(-10,10)
        y = rand.uniform(-10,10)
        posVec = rs.VectorCreate([0,0,0],[x,y,0])
        n = node(posVec,i,i+1)
        n.drawNode()
        
def throwNode(nodes,speed,stickRange,rNodes,boundRadius):
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
        time.sleep(0.0001)
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
                    newNode = node(currentPos,i,1) #parent is idx of node stuck to, r = 1
                    rs.AddLine(n.posVec,currentPos)
                    nodes.append(newNode)
                    #rNodes.append(rs.AddPoint(currentPos))
                    isTravelling = False
                    for o in previewGeom: rs.DeleteObjects(o)
                    break
            
            currentPos = rs.VectorAdd(currentPos,vel)
        Rhino.RhinoApp.Wait()

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
