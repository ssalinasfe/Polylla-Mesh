# Generate a PLSG in a square of 2x2 with a quartercircle of radio 0.4 in the left botton corner and random points inside.
# Input: nPoints: number of points inside the square and in the borders
#        Tolerance: Distance to consirate a point to close to the border, so it must move to the border
#        h_input: Separation beetwen the points of the quartercircle.
# Output: A .poly file


import numpy
import getopt, sys
from shapely.geometry import Point, Polygon, LineString
import random

def quatercircle2x2(num_P, tolerance, h_input):
    
    r = 0.4 #radio
    h = h_input #separation
    n = int( (numpy.pi/2)/ h)
    points = numpy.array(
        [
            [r*numpy.cos(alpha), r*numpy.sin(alpha)]
            for alpha in numpy.linspace(numpy.pi/2 - h, 0.0 , n + 1, endpoint=False)
        ]
    )
    points = numpy.concatenate((points, [[0.4,0.0],  [2.0,0.0], [2.0,2.0], [0, 2.0], [0.0, 0.4]]), axis=0)

    border_original = points.tolist()
    poly = Polygon(border_original) 
   # num_P = 20
    
        #if(numpy.sqrt(x*x + y*y) > 0.4):
        #    s.add(move_point(2.0, x, y))
    s = set()
    random.seed(138)
    while len(s) + len(border_original) != num_P:
        x = random.uniform(0.0, 2.0 )
        y = random.uniform(0.0, 2.0 )
        pt = Point(x,y)
        if(poly.contains(pt) and numpy.sqrt(x*x + y*y) > 0.4):
            #distance = numpy.round(poly.boundary.distance(pt), len(str(tolerance))-1 )
            if( poly.boundary.distance(pt) < tolerance):
                border_pt = poly.exterior.interpolate(poly.boundary.project(pt))
                #border_pt = nearest_points(poly, pt)
                #s.add((border_pt[0].x, border_pt[0].y))
                insert_border_ccw(border_original, border_pt)
            else:
                s.add((x, y))
    #print(new_border)
    constrai = [[k, k + 1] for k in range(len(border_original)-1)] + [[len(border_original)-1, 0]]
    write_poly(border_original + list(s),constrai, tolerance, h_input)

def insert_border_ccw(border, pt):
    for i in range(len(border)):
        e1 = border[i % len(border)]
        e2 = border[(i+1) % len(border)]
        line = LineString([Point(e1),Point(e2)])
        #print(e1,e2,line.wkt)
        if (pt.within(line)):
            border.insert((i+1)%len(border),(pt.x, pt.y))
            break



def write_poly(points, segments, tolerance, h_input):
    n_points =len(points)
    n_segments = len(segments)
    #print(n_points, 2, 0, 0)
    #for i in range(0, n_points):
    #    print(i + 1, points[i][0], points[i][1])
    #print(n_segments, 0)
    #for i in range(0, n_segments):
    #    print(i + 1, segments[i][0] + 1, segments[i][1] + 1)
        
    f = open('RP2x2quartercircleV2_' + str(n_points) + '.poly', 'w')
    f.write("#{} {} {} {}\n".format(n_points, n_segments, tolerance, h_input))
    f.write("{} 2 0 0\n".format(n_points))
    for i in range(0, n_points):
        f.write('{0} {1} {2}\n'.format(i + 1, points[i][0], points[i][1]))
    f.write("{} 0\n".format(n_segments))
    for i in range(0, n_segments):
        f.write( '{0} {1} {2}\n'.format(i + 1, segments[i][0] + 1, segments[i][1] + 1))
    f.write('\n')
    f.close()
    print("write " + "RP2x2quartercircleV2_" + str(n_points) + '.poly')

if __name__ == "__main__":
    #test_quater_disk()
   # test_rectangle()
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    nPoints = int(argument_list[0])
    tolerance = float(argument_list[1])
    h_input = float(argument_list[2])
    
    quatercircle2x2(nPoints, tolerance, h_input)