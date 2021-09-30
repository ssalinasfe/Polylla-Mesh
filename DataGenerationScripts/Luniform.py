# Generate a semi uniform triangulation with L-shape
# Input: edge_len: Max edge-size, the number of points of the final triangulation depend on this parametrer 
# Output: Two files, a .node and a .poly


import numpy
import pygalmesh
import sys

def square2x2(edge_size):
    points = numpy.array([[-1.0, 1.0], [-1.0, -1.0], [0.0, -1.0], [0.0, 0.0], [1.0, 0.0], [1.0,1.0]])
    constraints = [[k, k + 1] for k in range(len(points)-1)] + [[len(points)-1, 0]]
    mesh = pygalmesh.generate_2d(
        points,
        constraints,
        max_edge_size=edge_size,
        num_lloyd_steps=10,
    )

    print(len(mesh.points), len( mesh.get_cells_type("triangle")))
    write_node(mesh.points)
    write_ele(len(mesh.points), mesh.get_cells_type("triangle"))
    #mesh.write("2x2Luniform_.svg")

def write_node(vertices):
    largo =len(vertices)
    f = open('input/2x2Luniform_' + str(largo) + '.node', 'w')
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
    f = open('input/2x2Luniform_' + str(num_vertices) + '.ele', 'w')
    f.write("{} 3 1\n".format(largo))
    for i in range(0, largo):
        f.write( '{0} {1} {2} {3} 1\n'.format(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    edge_size = float(argument_list[0])
    square2x2(edge_size)
