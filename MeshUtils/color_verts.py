import mesh_neighbors
import rhinoscriptsyntax as rs
import System.Drawing.Color as Color

mesh = rs.GetObject("sel mesh", rs.filter.mesh)

col = [255, 0, 255]
one_color_list = [col for i in range(rs.MeshVertexCount(mesh))]

rs.MeshVertexColors(mesh, colors=one_color_list)
