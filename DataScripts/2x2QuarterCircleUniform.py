# Generate a semi Uniform triangulation inside a square of 2x2 with a quartercircle of radio 0.4 in the left botton corner.
# Input: edge_len: Max edge-size, the number of points of the final triangulation depend on this parametrer 
#        h_input: Separation beetwen the points of the quartercircle.
# Output: Two files, a .node and a .poly


import numpy
import pygalmesh
import getopt, sys

def alejandro(h_input, edge_len):
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
    
    constraints = [[k, k + 1] for k in range(len(points)-1)] + [[len(points)-1, 0]]

    mesh = pygalmesh.generate_2d(
        points,
        constraints,
        max_edge_size=edge_len,
        num_lloyd_steps=10,
    )
    print(len(mesh.points), len( mesh.get_cells_type("triangle")))
    write_node(mesh.points, h_input, edge_len)
    write_ele(len(mesh.points), mesh.get_cells_type("triangle"))
   # mesh.write("alejandro.svg")

def write_node(vertices, h_input, edge_len):
    largo =len(vertices)
    f = open('Uni2x2quartercircle_' + str(largo) + '.node', 'w')
    f.write("#num_vertices {}, h_input: {}, edge_len: {}\n".format(largo, h_input, edge_len))
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
    f = open('Uni2x2quartercircle_' + str(num_vertices) + '.ele', 'w')
    f.write("{} 3 1\n".format(largo))
    for i in range(0, largo):
        f.write( '{0} {1} {2} {3} 1\n'.format(i +1, triangles[i][0] +1 , triangles[i][1] +1, triangles[i][2] +1))
    f.write('\n')
    f.close()

if __name__ == "__main__":
    #test_quater_disk()
   # test_rectangle()
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    h_input = float(argument_list[1])
    edge_len = float(argument_list[0])
    alejandro(h_input, edge_len)
