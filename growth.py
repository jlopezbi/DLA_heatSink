import rhinoscriptsyntax as rs

MAX_GROW_LEN = .25
MIN_GROW_LEN = .01
MAX_NEIGHBOR_LEVEL = 5

def sharing_function(mesh, polyp, neighbor_level):
    '''computes how much the given polyp should grow by
    Args:
        mesh (rhinocommon mesh)
        polyp (int): vert index in mesh that is to grow
        neighbor_level (int): level of polyp from polyp-that-ate; acts as proxy for distance
    '''

    #get function (or perhaps pass in)
    #evaluate function to get grow_length
    pass

def grow_length(neighbor_level):
    ''' computes how much a polyp at neighbor_level from polyp_that_ate should grow by

    Args:
        neighbor_level (int): number corresponding to 'ring' of neighbor relative to the
        vert-that-ate
    Returns:
        float (grow_length)
    '''
    #TODO: write a linear func of neighbor level to get grow_length for the polyp
    slope = (MAX_GROW_LEN - MIN_GROW_LEN) / (0 - MAX_NEIGHBOR_LEVEL)
    y_intercept = MAX_GROW_LEN
    return slope * neighbor_level + y_intercept

def plot(X, Y):
    points = []
    for x, y in zip(X,Y):
        points.append((x,y,0))
    rs.AddPolyline(points)

if __name__=="__main__":
    # test grow_length
    x = list(range(MAX_NEIGHBOR_LEVEL+1))
    y = [grow_length(d) for d in x]
    print(x)
    print(y)
    assert(abs(y[-1] - MIN_GROW_LEN) < .001)
    plot(x, y)
