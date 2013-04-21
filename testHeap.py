# testing heaps

import rhinoscriptsyntax as rs
import heapq
import random
import math
import time
import Rhino
import scriptcontext

def testHeap():
	filter = Rhino.DocObjects.ObjectType.Mesh
	rc, objRef = Rhino.Input.RhinoGet.GetOneObject("select testMesh",False,filter)
	if not objRef or rc!=Rhino.Commands.Result.Success: return rc
	mesh  = objRef.Mesh()
	if not mesh: return
	mesh.Compact()

	randIdx = int(random.uniform(0,mesh.Vertices.Count-1))
	vertPnt = mesh.Vertices[randIdx]
	rad = .1
	rs.AddSphere(vertPnt,rad)

if __name__=="__main__":
	testHeap()

