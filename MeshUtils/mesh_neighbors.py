import rhinoscriptsyntax as rs
import scriptcontext
import Rhino
# probs should be inmport from meshutils in rhino_unfolder

def user_select_mesh(message="select a mesh"):
    mesh_guid = rs.GetObject(message,filter=32,preselect=True)
    return get_geom_from_guid(mesh_guid)

def get_geom_from_guid(guid):
    obj = scriptcontext.doc.Objects.Find(guid)
    mesh =  obj.Geometry
    print mesh
    return mesh

def get_vert_neighbors(vert,mesh):
    return set(list(mesh.Vertices.GetConnectedVertices(vert)))

def visualize_all(neighbors,mesh):
    shade = 0
    shade_increment = 255/len(neighbors)
    for neighborhood in neighbors:
        visualize_a_neighborhood(neighborhood,mesh,(shade,shade,shade))
        shade += shade_increment
        print shade

def visualize_a_neighborhood(neighborhood,mesh,color):
    for vert in neighborhood:
        geom = show_mesh_vert(vert,mesh)
        rs.ObjectColor(geom,color)

def get_topo_vert(vert,mesh):
    return mesh.TopologyVertices.TopologyVertexIndex(vert)

def show_mesh_vert(vert,mesh):
    '''note this is a re-write of a mesh method probably found 
    in meshutils in rhinounfolder
    '''
    coordinates = mesh.Vertices[vert]
    #point = rs.AddPoint(coordinates)
    dot = rs.AddTextDot(str(vert),coordinates)
    return dot


def get_neighbor_verts(start_vert,mesh,n_levels):
    '''
    for a given vertex on a graph, find the n_separation neighbors 
    (made-up term) for n_levels out. Note this problem has likeley 
    already been solved! however this solutions is rather simple
    '''
    neighbors = []
    seen_verts = set([start_vert])
    first_neighbors = get_vert_neighbors(start_vert,mesh)
    prev_neighbors = first_neighbors
    neighbors.append(first_neighbors)
    seen_verts.update(first_neighbors)
    for i in range(n_levels):
        new_neighbors = set([])
        for vert in prev_neighbors:
            potential_neighbors = get_vert_neighbors(vert,mesh)
            newely_seen = potential_neighbors.difference(seen_verts)
            new_neighbors.update(newely_seen)
        neighbors.append(new_neighbors)
        prev_neighbors = new_neighbors
        seen_verts.update(new_neighbors)
    return neighbors

if __name__=="__main__":
    mesh = user_select_mesh() 
    start_vert = 458
    n_levels = 4
    neighbors = get_neighbor_verts(start_vert,mesh,n_levels)
    print neighbors
    visualize_all(neighbors,mesh)
