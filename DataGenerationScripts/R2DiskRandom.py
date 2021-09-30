# Generate a PLSG of a disk of radio 2 with random points inside
# Input: nPoints: number of points inside the square and in the borders
#        Tolerance: Distance to consirate a point to close to the border, so it must move to the border
#        h_input: Separation beetwen the points of the disk.
# Output: A .poly file



import numpy
import getopt, sys
import pygalmesh
import random

def quatercircle2x2(num_P, h_O, tolerance):
    
    h = h_O
    n = int(2 * numpy.pi / h)
    
    points = numpy.array(
        [
            [numpy.cos(alpha), numpy.sin(alpha)]
            for alpha in numpy.linspace(0.0, 2 * numpy.pi, n + 1, endpoint=False)
        ]
    )
    constraints = [[k, k + 1] for k in range(n)] + [[n, 0]]
    random.seed(138)
    s = set([(random.uniform(0.00000001, 2.0000000 ), random.uniform(0.00000001, 2.0000000 ))])

    #print(num_P, len(s) + len(points))
    while len(s) + len(points) != num_P:
        x = random.uniform(-1.0, 1.0 )
        y = random.uniform(-1.0, 1.0 )
        if(numpy.sqrt(x*x + y*y) < (1 - tolerance)*(1 - tolerance) ):
            s.add( ( x, y) )
    
    putos = numpy.concatenate((points, list(s)))

    #print(len(putos), len(constrai))
    write_poly(putos,constraints)


def write_poly(points, segments):
    n_points =len(points)
    n_segments = len(segments)
    #print(n_points, 2, 0, 0)
    #for i in range(0, n_points):
    #    print(i + 1, points[i][0], points[i][1])
    #print(n_segments, 0)
    #for i in range(0, n_segments):
    #    print(i + 1, segments[i][0] + 1, segments[i][1] + 1)
        
    f = open('input/2x2RPDisk_' + str(n_points) + '.poly', 'w')
    f.write("{} 2 0 0\n".format(n_points))
    for i in range(0, n_points):
        f.write('{0} {1} {2}\n'.format(i + 1, points[i][0], points[i][1]))
    f.write("{} 0\n".format(n_segments))
    for i in range(0, n_segments):
        f.write( '{0} {1} {2}\n'.format(i + 1, segments[i][0] + 1, segments[i][1] + 1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    #test_quater_disk()
   # test_rectangle()
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    nPoints = int(argument_list[0])
    h = float(argument_list[1])
    tolerance = float(argument_list[2])
    
    quatercircle2x2(nPoints, h, tolerance)
