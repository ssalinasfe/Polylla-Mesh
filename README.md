
# POLYLLA: Polygonal meshing algorithmbased on terminal-edge regions

New algorithm to generate polygonal meshes of arbitrary shape, using any kind of triangulation as input, adaptable to any kind of geometry, no addition of extra points and with a data structure tp an easy implementation to all language programmings.

The input commands are the same as Detri2, the inputs files are

- **[.node](https://www.cs.cmu.edu/~quake/triangle.node.html)**: Point file
- **[.ele](https://www.cs.cmu.edu/~quake/triangle.ele.html)**: Triangle file
- **[.poly](https://www.cs.cmu.edu/~quake/triangle.poly.html)**: PLSG file

Example of use: 

Calls to file with points:

```./DelaunayPolyGenerator -z [name of file].node ```

Calls to PLSG file:

```./DelaunayPolyGenerator -p [name of file].node ```

Calls to triangulation (needs two file in the same folder, .node and .ele)

```./DelaunayPolyGenerator [name of file].ele ```

# Work in progress. TODO

  - [ ] Agregar que acepte cualquier input mientras que tengan un archivo .node. .ele, -neigh
  - [ ] Verificar la estructura de datos del input para que haya una correcta asociación de los triangulos y la adjacencia
  - [ ] Aceptar como input .node and .ele y generar la adjacencia en base a  eso.
  - [ ] Usar safe arrays en vez de array
  - [ ] Cambiar polygon array por list o vector
  - [ ] Eliminar funciones no usadas
  - [ ] Reescribir todo y encapsularlo en una clase para usarlo como librería
  - [ ] Llamar a triangle en vez de a detri2 cuadno solo se tenga de input un .node (detri no genera triangulaciones de más 10^6 puntos)
  - [ ] Generar como output un svg
  - [ ] Construir test


# Input data generation

Following scripts to generate input files were made to test the algorithm using the Virtual Element Method

## Random points files generation (.node)

  - 2x2RandomPoints.py: Generate random points inside a square 2x2

## Random point PLSG (.poly)

  - R2DiskRandom.py: Generate a PLSG of a circle of radio 2 with random points inside of it
  - 2x2Lrandom.py: Generate a PLSG with L-shape with random points inside of it
  - 2x2QuaterCircle.py: Generate a PLSG in a square of 2x2 with a quartercircle of radio 0.4 in the left botton corner and random points inside of it

## Semiuniform triangulations (.node and .ele)
    
  - R2DiskSemiuniform.py: Generate a semiuniform triangulation in a disk of radio 2
  - 2x2QuaterCircleUniform.py: Generate a semi Uniform triangulation inside a square of 2x2 with a quartercircle of radio 0.4 in the left botton corner.
  - Luniform: Generate a semi uniform triangulation with L-shape


 
