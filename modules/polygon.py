import numpy as np
from numba import jit, njit
import numba

@jit(nopython=True)
def is_inside_sm(polygon, point):
    """ Returns whether a point falls within a polygon

    ****** PARAMETERS ******
    1. polygon: defined by at least three points [[x1, y1], [x2, y2], [x3,y3]]
    2. point: a single point [xp, yp]

    ******* RETURNS ********
    1. intersections: point is inside polygon if this variable is greater than 0
    """
    length = len(polygon)-1
    dy2 = point[1] - polygon[0][1]
    intersections = 0
    ii = 0
    jj = 1

    while ii<length:
        dy  = dy2
        dy2 = point[1] - polygon[jj][1]

        # consider only lines which are not completely above/bellow/right from the point
        if dy*dy2 <= 0.0 and (point[0] >= polygon[ii][0] or point[0] >= polygon[jj][0]):

            # non-horizontal line
            if dy<0 or dy2<0:
                F = dy*(polygon[jj][0] - polygon[ii][0])/(dy-dy2) + polygon[ii][0]

                if point[0] > F: # if line is left from the point the ray moving towards left will intersect it
                    intersections += 1
                elif point[0] == F: # point on line
                    return 2

            # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and dx*dx2<=0)
            elif dy2==0 and (point[0]==polygon[jj][0] or (dy==0 and \
                 (point[0]-polygon[ii][0])*(point[0]-polygon[jj][0])<=0)):
                return 2
        ii = jj
        jj += 1
    return intersections & 1


@njit(parallel=True)
def is_inside_sm_parallel(points, polygon):
    """ Iterate whether a list of points are within a polygon

    ****** PARAMETERS ******
    1. polygon: defined by at least three points [[x1, y1], [x2, y2], [x3,y3]]
    2. points: a list of points [[xp1, yp1], [xp2, yp2], ..., [xpn, ypn]]

    ******* RETURNS ********
    1. inside: List containing which points are within a polygon (> 0)
    """
    ln = len(points)
    inside = np.empty(ln, dtype=numba.boolean)
    for i in numba.prange(ln):
        inside[i] = is_inside_sm(polygon, points[i])
    return inside

