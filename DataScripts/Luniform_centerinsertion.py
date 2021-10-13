# Generate a semiuniform triangulation with L-shape with random points inside of
# No use it, I don't remember how it works

import numpy
import pygalmesh
import sys
from shapely.geometry import Point, Polygon, LineString
import random


def square2x2(edge_size, tolerance , number_center, radius_center = 0.2):
    points = numpy.array([[-1.0, 1.0], [-1.0, -1.0], [0.0, -1.0], [0.0, 0.0], [1.0, 0.0], [1.0,1.0]])
    
    
    new_border = points.tolist()
    poly = Polygon(new_border)
    
    s = set()

    #insert in center
    while len(s) + len(new_border) != number_center:
        x = random.uniform(-1.0*radius_center, 1.0*radius_center )
        y = random.uniform(-1.0*radius_center, 1.0*radius_center )
        pt = Point(x,y)
        
        if(poly.contains(pt) and numpy.sqrt(x*x + y*y) <= (radius_center)*(radius_center) ):
            #distance = numpy.round(poly.boundary.distance(pt), len(str(tolerance))-1 )
            if( poly.boundary.distance(pt) < tolerance*radius_center):
                border_pt = poly.exterior.interpolate(poly.boundary.project(pt))
                #border_pt = nearest_points(poly, pt)
                #s.add((border_pt[0].x, border_pt[0].y))
                insert_border_ccw(new_border, border_pt)
            else:
                s.add((x, y))


    border = numpy.array(new_border)
    for i in s:
        border = numpy.append(border, [list(i)], axis=0)
    
    constraints = [[k, k + 1] for k in range(len(new_border)-1)] + [[len(new_border)-1, 0]]

    mesh = pygalmesh.generate_2d(
        border,
        constraints,
        max_edge_size=edge_size,
        num_lloyd_steps=10,
    )

    print(len(mesh.points), len( mesh.get_cells_type("triangle")))
    write_node(mesh.points, edge_size, tolerance, number_center, radius_center)
    write_ele(len(mesh.points), mesh.get_cells_type("triangle"))
    #mesh.write("2x2Luniform_.svg")


def insert_border_ccw(border, pt):
    for i in range(len(border)):
        e1 = border[i % len(border)]
        e2 = border[(i+1) % len(border)]
        line = LineString([Point(e1),Point(e2)])
        #print(e1,e2,line.wkt)
        if (pt.within(line)):
            border.insert((i+1)%len(border),(pt.x, pt.y))
            break


def write_node(vertices, edge_size, tolerance, number_center, radius_center):
    largo =len(vertices)
    f = open('input/2x2LuniformLCI_' + str(largo) + '.node', 'w')
    f.write("# Points {}, Edge_size {},  tolerance {}, number_center {}, radius_center {}\n".format(largo, edge_size,  tolerance, number_center, radius_center))
    f.write("{} 2 0 0\n".format(largo))
    for i in range(0, len(vertices)):
        #print(i +1, vertices[i][0], vertices[i][1])
        f.write('{0} {1} {2}\n'.format(i+1, vertices[i][0],vertices[i][1]))
    f.write('\n')
    f.close()

def write_ele(num_vertices, triangles):
    largo =len(triangles)
    #for i in range(0, largo):
    #    print(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1)
    f = open('input/2x2LuniformLCI_' + str(num_vertices) + '.ele', 'w')
    f.write("{} 3 1\n".format(largo))
    for i in range(0, largo):
        f.write( '{0} {1} {2} {3} 1\n'.format(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    edge_size = float(argument_list[0])
    tolerance = float(argument_list[1])
    number_center = float(argument_list[2])
    radius = float(argument_list[3])
    square2x2(edge_size, tolerance, number_center, radius)
