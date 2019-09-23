import utils
import importlib
import CCoral

def test_grow_one_neighborhood():
    ratio = .33 #copied from MeshDLA.py
    mesh_guid, mesh = utils.user_select_mesh()
    coral = CCoral.Coral(mesh_guid, mesh, ratio)
    grow_verts = coral.getGrowData(gKernel=None, centerVerts=[0])
    coral.grow(grow_verts)

if __name__=="__main__":
    test_grow_one_neighborhood()
