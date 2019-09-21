import mesh_neighbors as mn
import System.Drawing.Color as Color
import scriptcontext

def test_color_a_vert(mesh_guid, mesh):
    start_vert =0
    mesh.VertexColors.CreateMonotoneMesh(Color.FromArgb(0,0,255))
    mesh.VertexColors.SetColor(start_vert,255,0,255)
    scriptcontext.doc.Objects.Replace(mesh_guid, mesh)

if __name__=="__main__":
    mesh_guid, mesh = mn.user_select_mesh()
    test_color_a_vert(mesh_guid, mesh)
