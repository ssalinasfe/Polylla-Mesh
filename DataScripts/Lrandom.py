# Generate a PLSG with L-shape with random points inside of
# Input: nPoints: number of points inside the square and in the borders
#        Tolerance: Distance to consirate a point to close to the border, so it must move to the border
# Output: .poly file


import numpy
import random
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import nearest_points
import sys

def square2x2(num_points, tolerance):
    border_original = [(-1.0, 1.0), (-1.0, -1.0), (0.0, -1.0), (0.0, 0.0), (1.0, 0.0), (1.0,1.0)]
    new_border = border_original.copy()

    poly = Polygon(border_original) 
    
    s = set()
    while len(s) + len(new_border) != num_points:
        x = random.uniform(-1.0, 1.0 )
        y = random.uniform(-1.0, 1.0 )
        pt = Point(x,y)
        if(poly.contains(pt)):
            #distance = numpy.round(poly.boundary.distance(pt), len(str(tolerance))-1 )
            if( poly.boundary.distance(pt) < tolerance):
                border_pt = poly.exterior.interpolate(poly.boundary.project(pt))
                #border_pt = nearest_points(poly, pt)
                #s.add((border_pt[0].x, border_pt[0].y))
                insert_border_ccw(new_border, border_pt)
            else:
                s.add((x, y))
    #print(new_border)
    constrai = [[k, k + 1] for k in range(len(new_border)-1)] + [[len(new_border)-1, 0]]
    write_poly(new_border + list(s),constrai)

def insert_border_ccw(border, pt):
    for i in range(len(border)):
        e1 = border[i % len(border)]
        e2 = border[(i+1) % len(border)]
        line = LineString([Point(e1),Point(e2)])
        #print(e1,e2,line.wkt)
        if (pt.within(line)):
            border.insert((i+1)%len(border),(pt.x, pt.y))
            break

def write_poly(points, segments):
    n_points =len(points)
    n_segments = len(segments)
    #print(n_points, 2, 0, 0)
    #for i in range(0, n_points):
    #    print(i + 1, points[i][0], points[i][1])
    #print(n_segments, 0)
    #for i in range(0, n_segments):
    #    print(i + 1, segments[i][0] + 1, segments[i][1] + 1)
        
    f = open('input/RP2x2L_' + str(n_points) + '.poly', 'w')
    f.write("{} 2 0 0\n".format(n_points))
    for i in range(0, n_points):
        f.write('{0} {1} {2}\n'.format(i + 1, points[i][0], points[i][1]))
    f.write("{} 0\n".format(n_segments))
    for i in range(0, n_segments):
        f.write( '{0} {1} {2}\n'.format(i + 1, segments[i][0] + 1, segments[i][1] + 1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    random.seed(1312314)
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    num_points = int(argument_list[0])
    tolerance = float(argument_list[1])
    square2x2(num_points, tolerance)
