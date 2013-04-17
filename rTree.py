import Rhino
import rhinoscriptsyntax as rs

 
def RunSearch():
    id = rs.GetObject("select mesh", rs.filter.mesh)
    mesh = rs.coercemesh(id)
    if mesh:
        rs.UnselectObject(id)
        tree = Rhino.Geometry.RTree()
        # I can add a RhinoCommon function that just builds an rtree from the
        # vertices in one quick shot, but for now...
        for i,vertex in enumerate(mesh.Vertices): 
            tree.Insert(vertex, i)
 
        while(True):
            point = rs.GetPoint("test point")
            if not point: break
 
            data = SearchData(mesh, point)
            sphere = Rhino.Geometry.Sphere(point, 2)

            vertices = set()

            def SearchCallback(sender, data):
                vertices = data.Tag
                vertices.add(data.Id)

            if tree.Search(sphere, SearchCallback, vertices):
                for index in vertices:
                    rs.AddPoint(mesh.Vertices[index])
 
if __name__=="__main__":
    RunSearch()