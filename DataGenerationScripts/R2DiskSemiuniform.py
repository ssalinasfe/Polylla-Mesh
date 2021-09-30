# Generate a semiuniform triangulation in a disk of radio 2
# Input: edge_len: Max edge-size, the number of points of the final triangulation depend on this parametrer 
#        h_input: Separation beetwen the points of the disk.
# Output: Two files, a .node and a .poly

import numpy
import pygalmesh
import sys

def generate_disk(h):
    
    n = int(2 * numpy.pi / h)
    
    points = numpy.array(
        [
            [numpy.cos(alpha), numpy.sin(alpha)]
            for alpha in numpy.linspace(0.0, 2 * numpy.pi, n + 1, endpoint=False)
        ]
    )
    constraints = [[k, k + 1] for k in range(n)] + [[n, 0]]
    mesh = pygalmesh.generate_2d(
        points,
        constraints,
        max_edge_size=h,
        num_lloyd_steps=10,
    )
    
    print(len(mesh.points), len( mesh.get_cells_type("triangle")))
    write_node(mesh.points)
    write_ele(len(mesh.points), mesh.get_cells_type("triangle"))
    #mesh.write("diskmesh.svg")

def write_node(vertices):
    largo =len(vertices)
    f = open('input/disk2x2_' + str(largo) + '.node', 'w')
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
    f = open('input/disk2x2_' + str(num_vertices) + '.ele', 'w')
    f.write("{} 3 1\n".format(largo))
    for i in range(0, largo):
        f.write( '{0} {1} {2} {3} 1\n'.format(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    edge_size = float(argument_list[0])
    generate_disk(edge_size)