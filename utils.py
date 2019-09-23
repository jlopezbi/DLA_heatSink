import rhinoscriptsyntax as rs
import scriptcontext

def user_select_mesh(message="select a mesh"):
    mesh_guid = rs.GetObject(message,filter=rs.filter.mesh, preselect=True)
    return mesh_guid, get_geom_from_guid(mesh_guid)

def get_geom_from_guid(guid):
    obj = scriptcontext.doc.Objects.Find(guid)
    mesh =  obj.Geometry
    return mesh
