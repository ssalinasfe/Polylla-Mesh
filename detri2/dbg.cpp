#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "detri2.h"

using  namespace detri2;

//==============================================================================

bool Triangulation::is_edge(TriEdge &E, int i, int j)
{
  return (((E.org()->idx == i) && (E.dest()->idx == j)) ||
          ((E.org()->idx == j) && (E.dest()->idx == i)));
}

bool Triangulation::is_triangle(Triang *tri, int i, int j, int k)
{
  if (tri->vrt[0]->idx == i) {
    if ((tri->vrt[1]->idx == j) &&
        (tri->vrt[2]->idx == k)) {
      return true;
    }
    if ((tri->vrt[1]->idx == k) &&
        (tri->vrt[2]->idx == j)) {
      return true;
    }
  } else if (tri->vrt[0]->idx == j) {
    if ((tri->vrt[1]->idx == i) &&
        (tri->vrt[2]->idx == k)) {
      return true;
    }
    if ((tri->vrt[1]->idx == k) &&
        (tri->vrt[2]->idx == i)) {
      return true;
    }
  } else if (tri->vrt[0]->idx == k) {
    if ((tri->vrt[1]->idx == i) &&
        (tri->vrt[2]->idx == j)) {
      return true;
    }
    if ((tri->vrt[1]->idx == j) &&
        (tri->vrt[2]->idx == i)) {
      return true;
    }
  }

  return false;
}

Vertex* Triangulation::db_vrt(int v)
{
  for (int i = 0; i < ct_in_vrts; i++) {
    if (in_vrts[i].idx == v) {
      in_vrts[i].print();
      return &(in_vrts[i]);
    }
  }
  if (tr_steiners != NULL) {
    for (int i = 0; i < tr_steiners->used_items; i++) {
      Vertex *vrt = (Vertex *) tr_steiners->get(i);
      if (vrt->is_deleted()) continue;
      if (vrt->idx == v) {
        vrt->print();
        return vrt;
      }
    }
  }
  printf("!! Not exist.\n");
  return NULL;
}

TriEdge Triangulation::db_seg(int v1, int v2)
{
  TriEdge E;
  for (int i = 0; i < tr_tris->used_items; i++) {
    E.tri = (Triang *) tr_tris->get(i);
    if (E.tri->is_deleted()) continue;
    for (E.ver = 0; E.ver < 3; E.ver++) {
      if (((E.org()->idx == v1) && (E.dest()->idx == v2)) ||
          ((E.dest()->idx == v1) && (E.org()->idx == v2))) {
        E.print();
        E.esym().print();
        return E;
      }
    }
  }
  printf("!! Not exist.\n");
  E.tri = NULL;
  return E;
}

TriEdge Triangulation::db_tri(int v1, int v2, int v3)
{
  TriEdge E;
  for (int i = 0; i < tr_tris->used_items; i++) {
    E.tri = (Triang *) tr_tris->get(i);
    if (E.tri->is_deleted()) continue;
    for (E.ver = 0; E.ver < 3; E.ver++) {
      if (((E.org()->idx == v1) && (E.dest()->idx == v2)) ||
          ((E.dest()->idx == v1) && (E.org()->idx == v2))) {
        if (E.apex()->idx == v3) {
          E.tri->print(1);
          return E;
        }
      }
    }
  }
  printf("!! Not exist.\n");
  E.tri = NULL;
  return E;
}

void Triangulation::dump_vertex_star(int v)
{
  Vertex *pt = db_vrt(v);
  if (pt == NULL) return;

  if (pt->typ == UNUSEDVERTEX) {
    return;
  }

  TriEdge E = pt->adj;
  do {
    E.print();
    E = E.eprev_esym();
  } while (E.tri != pt->adj.tri);
}

//==============================================================================
// debug function for get_powercell().

void Triangulation::dump_voronoi_vertices(Vertex *mesh_vertex)
{
  if (mesh_vertex->typ == UNUSEDVERTEX) {
    printf("  It is an UNUSEDVERTEX.\n");
    return;
  }

  TriEdge E = mesh_vertex->adj;
  int vcount = 0; // Count the number of Voronoi vertices.
  do {
    vcount++;
    E = E.eprev_esym(); // ccw
  } while (E.tri != mesh_vertex->adj.tri);

  printf("%d 2 0 0\n", vcount);

  E = mesh_vertex->adj;
  int interiorflag = 1;
  vcount = 0;
  do {
    vcount++;
    if (E.tri->is_dual_on_bdry()) {
      interiorflag = 0;
    } else if (E.tri->is_dual_in_exterior()) {
      interiorflag = -1;
    } else {
      interiorflag = 1;
    }
    printf("%d %g %g # [%d,%d,%d] interior(%d)\n",
           vcount,
           E.tri->cct[0], E.tri->cct[1],
           E.org()->idx, E.dest()->idx, E.apex()->idx,
           interiorflag);
    E = E.eprev_esym(); // ccw
  } while (E.tri != mesh_vertex->adj.tri);
}

void Triangulation::dump_voronoi_vertices(int vertex_index)
{
  Vertex *mesh_vertex;

  int shift = io_firstindex == 0 ? 1 : 0;

  if ((vertex_index + shift) <= ct_in_vrts) {
    mesh_vertex = &(in_vrts[vertex_index - io_firstindex]);
  } else {
    // Find it in the Steiner points.
    printf("  It is a Steiner point.\n");
    return; // to do...
  }

  return dump_voronoi_vertices(mesh_vertex);
}

//==============================================================================
// This function is replaced by construct_voronoi_diagram()

void Triangulation::save_voronoi(int ucd)
{
  // We need a background grid to cut the exterior Voronoi cells.
  bool clean_omt_domain = false;
  if (OMT_domain == NULL) {
    OMT_domain = this; // Use the current triangulation.
    clean_omt_domain = true;
  }

  if (ct_exteriors > 0) {
    // removing exterior triangles before caculating Voronoi cells.
    remove_exteriors();
  }

  // Calculate Voronoi vertices for triangles.
  int i, idx;
  idx = io_firstindex;
  // Calculate circumcenters for hull triangles.
  for (i = 0; i < tr_tris->used_items; i++) {
    Triang* tri = (Triang *) tr_tris->get(i);
    // Ignore exterior triangles.
    if (tri->is_deleted() || tri->is_hulltri()) continue;
    get_tri_orthocenter(tri);
    tri->idx = idx;
    idx++;
  }
  //Tr->op_db_verbose=10; // debug
  // Calculate bisectors for hull triangles.
  for (i = 0; i < tr_tris->used_items; i++) {
    Triang* tri = (Triang *) tr_tris->get(i);
    // Ignore exterior triangles.
    if (tri->is_deleted()) continue;
    if (tri->is_hulltri()) { // A hull triangle.
      get_hulltri_orthocenter(tri);
      tri->idx = idx; // hulltri is also indexed.
      idx++;
    }
  }

  // Store the list of Voronoi vertices and Voronoi edges.
  arraypool *vert_list = new arraypool(sizeof(Vertex), 10);
  arraypool *edge_list = new arraypool(sizeof(Vertex *) * 2, 10);
  arraypool *cell_list = new arraypool(sizeof(Vertex **), 10);
  arraypool *cell_size_list = new arraypool(sizeof(int), 10);

  for (i = 0; i < ct_in_vrts; i++) {
    Vertex *mesh_vertex = &(in_vrts[i]);
    if (mesh_vertex->typ == UNUSEDVERTEX) continue;
    Vertex *ptlist = NULL;
    int ptnum = 0;
    if (get_powercell(mesh_vertex, false, &ptlist, &ptnum)) {
      // Store the cell vertices (pointers).
      int *cell_size = (int *) cell_size_list->alloc();
      *cell_size = ptnum;
      Vertex ***cell = (Vertex ***) cell_list->alloc();
      *cell = new Vertex*[ptnum];
      
      for (int j = 0; j < ptnum; j++) {
        Vertex *pt = (Vertex *) vert_list->alloc();
        pt->init();
        pt->crd[0] = ptlist[j].crd[0];
        pt->crd[1] = ptlist[j].crd[1];
        pt->crd[2] = ptlist[j].crd[2];
        (*cell)[j] = pt;
      }
      for (int j = 0; j < ptnum; j++) {
        Vertex **edge = (Vertex **) edge_list->alloc();
        edge[0] = (*cell)[j];
        edge[1] = (*cell)[(j+1)%ptnum];
      }
    }
    delete [] ptlist;
  }

  if (tr_steiners != NULL) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *mesh_vertex = (Vertex *) tr_steiners->get(i);
      if (mesh_vertex->is_deleted()) continue;
      Vertex *ptlist = NULL;
      int ptnum = 0;
      if (get_powercell(mesh_vertex, false, &ptlist, &ptnum)) {
        // Store the cell vertices (pointers).
        int *cell_size = (int *) cell_size_list->alloc();
        *cell_size = ptnum;
        Vertex ***cell = (Vertex ***) cell_list->alloc();
        *cell = new Vertex*[ptnum];
      
        for (int j = 0; j < ptnum; j++) {
          Vertex *pt = (Vertex *) vert_list->alloc();
          pt->init();
          pt->crd[0] = ptlist[j].crd[0];
          pt->crd[1] = ptlist[j].crd[1];
          pt->crd[2] = ptlist[j].crd[2];
          (*cell)[j] = pt;
        }
        for (int j = 0; j < ptnum; j++) {
          Vertex **edge = (Vertex **) edge_list->alloc();
          edge[0] = (*cell)[j];
          edge[1] = (*cell)[(j+1)%ptnum];
        }
      }
      delete [] ptlist;
    }
  } // if (tr_steiners != NULL)

  // Unifying vertices and counting edges.
  Triangulation *Voro = new Triangulation();

  // debug only.
  strcpy(Voro->io_outfilename, "voro"); // /Users/si/tmp/voro

  Voro->ct_in_vrts = vert_list->objects;
  Voro->in_vrts = new Vertex[vert_list->objects];
  Voro->io_xmin = Voro->io_ymin =  1.e+30;
  Voro->io_xmax = Voro->io_ymin = -1.e+30;

  for (i = 0; i < vert_list->objects; i++) {
    Vertex *p = (Vertex *) vert_list->get(i);
    Vertex *v = & (Voro->in_vrts[i]);
    v->init();
    // Default vertex type is UNUSEDVERTEX (0)
    v->typ = UNUSEDVERTEX;
    Voro->ct_unused_vrts++;
    v->crd[0] = p->crd[0];
    v->crd[1] = p->crd[1];
    v->crd[2] = 0.; // p->crd[0]*p->crd[0] + p->crd[1]*p->crd[1];
    p->Pair = v; // remember this vertex.
    Voro->io_xmin = p->crd[0] < Voro->io_xmin ? p->crd[0] : Voro->io_xmin;
    Voro->io_xmax = p->crd[0] > Voro->io_xmax ? p->crd[0] : Voro->io_xmax;
    Voro->io_ymin = p->crd[1] < Voro->io_ymin ? p->crd[1] : Voro->io_ymin;
    Voro->io_ymax = p->crd[1] > Voro->io_ymax ? p->crd[1] : Voro->io_ymax;
  }
  assert(Voro->ct_unused_vrts == Voro->ct_in_vrts);

  Voro->tr_segs = new arraypool(sizeof(Triang), 10);
  for (i = 0; i < edge_list->objects; i++) {
    Vertex **edge = (Vertex **) edge_list->get(i);
    Triang *seg = (Triang *) Voro->tr_segs->alloc();
    seg->init();
    seg->vrt[0] = edge[0]->Pair;
    seg->vrt[1] = edge[1]->Pair;
    seg->tag = 1; // A non-zero marker
  }
  assert(Voro->ct_segments == 0);

  Voro->incremental_delaunay();
  //Voro->save_triangulation(); // debug only.

  Voro->recover_segments();
  //Voro->save_triangulation(); // debug only.

  if (clean_omt_domain) {
    OMT_domain = NULL;
  }

  char filename[256];
  sprintf(filename, "%s.voro", io_outfilename);
  
  printf("Writing Voronoi diagram to %s\n", filename);
  
  FILE *outfile = fopen(filename, "w");

  int nv = Voro->ct_in_vrts - Voro->ct_unused_vrts;
  int ne = Voro->ct_segments;
  int nc = cell_list->objects;

  fprintf(outfile, "%d %d %d\n", nv, ne, nc);

  printf("  Writing %d vertices.\n", nv);
  idx = 1;  // indexing all vertices.
  for (i = 0; i < Voro->ct_in_vrts; i++) {
    Vertex *v = &(Voro->in_vrts[i]);
    if (v->typ == UNUSEDVERTEX) continue;
    fprintf(outfile, "%g %g %g\n", v->crd[0], v->crd[1], v->crd[2]);
    v->idx = idx;
    idx++;
  }
  //assert(idx == nv);

  printf("  Writing %d edges.\n", ne);
  idx = 1;
  for (i = 0; i < Voro->tr_segs->used_items; i++) {
    Triang *seg = (Triang *) Voro->tr_segs->get(i);
    if (seg->tag == 0) continue;
    Vertex *v1 = seg->vrt[0];
    Vertex *v2 = seg->vrt[1];
    if (seg->vrt[0]->typ == UNUSEDVERTEX) {
      v1 = seg->vrt[0]->Pair;
    }
    if (seg->vrt[1]->typ == UNUSEDVERTEX) {
      v2 = seg->vrt[1]->Pair;
    }
    fprintf(outfile, "%d %d  %d\n", v1->idx, v2->idx, seg->tag);
    idx++;
  }

  printf("  Writing %d cells.\n", nc);
  for (i = 0; i < cell_list->objects; i++) {
    int size = * (int *) cell_size_list->get(i);
    fprintf(outfile, "%d ", size);
    Vertex ***cell = (Vertex ***) cell_list->get(i);
    for (int j = 0; j < size; j++) {
      Vertex *p = (*cell)[j];
      Vertex *v = p->Pair; // This is the mesh vertex in Voro.
      if (v->typ == UNUSEDVERTEX) {
        v = v->Pair;
      }
      fprintf(outfile, " %d", v->idx);
    }
    fprintf(outfile, "\n");
    // release memory.
    delete [] *cell;
  }

  fclose(outfile);

  if (ucd > 0) { // -Iv -Iu
    // Save the skeleton of the Voronoi diagram for Paraview.
    sprintf(filename, "%s_voro.inp", io_outfilename);
    printf("Writing Voronoi diagram to %s\n", filename);
    outfile = fopen(filename, "w");
    
    fprintf(outfile, "%d %d %d 0 0\n", nv, ne, 0);
    
    idx = 1; // UCD index starts from 1.
    for (i = 0; i < Voro->ct_in_vrts; i++) {
      Vertex *v = &(Voro->in_vrts[i]);
      if (v->typ == UNUSEDVERTEX) continue;
      fprintf(outfile, "%d  %g %g 0\n", idx, v->crd[0], v->crd[1]);
      v->idx = idx;
      idx++;
    }
    
    idx = 1;
    for (i = 0; i < Voro->tr_segs->used_items; i++) {
      Triang *seg = (Triang *) Voro->tr_segs->get(i);
      if (seg->tag == 0) continue;
      Vertex *v1 = seg->vrt[0];
      Vertex *v2 = seg->vrt[1];
      if (seg->vrt[0]->typ == UNUSEDVERTEX) {
        v1 = seg->vrt[0]->Pair;
      }
      if (seg->vrt[1]->typ == UNUSEDVERTEX) {
        v2 = seg->vrt[1]->Pair;
      }
      fprintf(outfile, "%d %d line %d %d\n", idx, seg->tag, v1->idx, v2->idx);
      idx++;
    }
    
    fclose(outfile);
  } // if (ucd > 0)

  delete vert_list;
  delete edge_list;
  delete cell_list;
  delete cell_size_list;
  delete Voro;
}

//==============================================================================
/*
// Removing a segment at vertex, update the vertex-to-segment ring.
//   'v' must be an endpoint of a segment, i.e., v->typ is either RIDGEVERTEX or
//   SEGMENTVERTEX. In detri2, v->typ might be STEINERPOINT and v->is_fixed().

int Triangulation::remove_segment_at_vertex(Vertex *v, Triang *seg)
{
    TriEdge nextseg;
    if (seg->vrt[0] == v) {
      nextseg = seg->nei[0];
    } else if (seg->vrt[1] == v) {
      nextseg = seg->nei[1];
    } else {
      assert(0); // seg does not contain v.
    }
    if (nextseg.tri != NULL) {
      assert(nextseg.tri->vrt[nextseg.ver] == v);
    }

    TriEdge prevseg = v->on_bd;
    assert(prevseg.tri->vrt[prevseg.ver] == v);
    if (prevseg.tri == seg) {
      // remove the seg by restting the vertex-to-seg map.
      v->on_bd = nextseg;
      if (nextseg.tri == NULL) {
        // There is no segment at this vertex.
        // v is not a restricted vertex anymore.
        if (v->typ == RIDGEVERTEX) {
          v->typ = FACETVERTEX; // a free input vertex.
        } else if (v->typ == SEGMENTVERTEX) {
          v->typ = FACETVERTEX; // a free input vertex.
        } else if (v->typ == STEINERVERTEX) {
          v->clear_fix();
        } else {
          assert(0); // should not be possible.
        }
      }
      return 1;
    }

    // Search the seg in the vertex-to-segment map.
    TriEdge searchseg = prevseg.tri->nei[prevseg.ver];
    while (searchseg.tri != NULL) {
      assert(searchseg.tri->vrt[searchseg.ver] == v);
      if (searchseg.tri == seg) {
        // Found. Remove seg from the list.
        prevseg.tri->nei[prevseg.ver] = nextseg;
        return 1;
      }
      // Go the next of the next segment.
      TriEdge tmp = searchseg.tri->nei[searchseg.ver];
      searchseg = tmp;
    }

    return 0; // not found.
}
*/

//==============================================================================

REAL Triangulation::get_element_gradient(Triang* tri)
{// get fx, fy, saved in cct[0], [1]
  // The formulas are in Lu-femag2d-2008-09-05.pdf (page 6)
  // The piecewise function on this triangle is equation (8).
  REAL xi = tri->vrt[0]->crd[0];
  REAL yi = tri->vrt[0]->crd[1];
  REAL ui = tri->vrt[0]->fval;

  REAL xj = tri->vrt[1]->crd[0];
  REAL yj = tri->vrt[1]->crd[1];
  REAL uj = tri->vrt[1]->fval;

  REAL xm = tri->vrt[2]->crd[0];
  REAL ym = tri->vrt[2]->crd[1];
  REAL um = tri->vrt[2]->fval;

  //REAL ai = xj*ym - xm*yj;
  //REAL aj = xm*yi - xi*ym;
  //REAL am = xi*yj - xj*yi;

  REAL bi = yj - ym;
  REAL bj = ym - yi;
  REAL bm = yi - yj;
  
  REAL ci = xm - xj;
  REAL cj = xi - xm;
  REAL cm = xj - xi;

  //REAL delta = 0.5 * (ai + aj + am);
  REAL delta = (bi*cj - bj*ci) / 2.;

  /*
  // PostGraph.cpp
  // void MFemPostGraph::CacuEleBxBy(MeshEle& ele, GFloat& Bx, GFloat& By,
  //                                 GFloat& GPx, GFloat& GPy)
  //tai = xj*ym-xm*yj;
  REAL tbi = yj-ym;
  REAL tci = xm-xj;
  //taj = xm*yi-xi*ym;
  REAL tbj = ym-yi;
  REAL tcj = xi-xm;
  //tam = xi*yj-xj*yi;
  REAL tbm = yi-yj;
  REAL tcm = xj-xi;
  REAL delta = (tbi * tcj - tbj * tci) / 2.;
  */
 
  if (delta == 0.0) {
    tri->cct[0] = 0.;
    tri->cct[1] = 0.;
    tri->cct[2] = 0.;
  } else {
    tri->cct[0] = 0.5 * (bi*ui + bj*uj + bm*um) / delta; // grad_x
    tri->cct[1] = 0.5 * (ci*ui + cj*uj + cm*um) / delta; // grad_y
    //tri->cct[0] =  (tbi * ui + tbj * uj + tbm * um) / 2. / delta;
    //tri->cct[1] =  (tci * ui + tcj * uj + tcm * um) / 2. / delta;
    if (fabs(tri->cct[0]) < 1.e-16) tri->cct[0] = 0.0;
    if (fabs(tri->cct[1]) < 1.e-16) tri->cct[1] = 0.0;
    // The Dirichlet energy (the dot product) (= square of H^1 norm) of this element.
    // From Marc Alexa's notation equ (13)
    //   1/2 f^T L_T f, where L_T is the discrete Laplace-Beltrami operator.
    //   It is a 3x3 matrix with components given in Lu's equ (30). (page 10).
    REAL area = fabs(delta);
    // debug only
    //REAL AA = get_tri_area(tri->vrt[0], tri->vrt[1], tri->vrt[2]);
    tri->cct[2] = (tri->cct[0] * tri->cct[0] + tri->cct[1] * tri->cct[1]) * area;
    if (fabs(tri->cct[2]) < 1.e-16) tri->cct[2] = 0.0;
  }

  return tri->cct[2];
}

//==============================================================================

REAL Triangulation::get_Dirichlet_energy() // of the whole triangulation
{
  REAL sum = 0.;

  for (int i = 0; i < tr_tris->used_items; i++) {
    Triang *tri = (Triang *) tr_tris->get(i);
    if (!tri->is_deleted()) {
      if (!tri->is_hulltri() && !tri->is_exterior()) {
        sum += get_element_gradient(tri);
      }
    }
  }

  return sum;
}
/*
REAL Triangulation::dump_Dirichlet_energy() // for visualisation
{
    REAL sum_energy = get_Dirichlet_energy();
    printf(" Energy = %g\n", sum_energy);

    // see file format: detri2/doc/UCD_Format.pdf
    char filename[256];
    sprintf(filename, "%s_energy.inp", io_outfilename);
    FILE *outfile = fopen(filename, "w");

    int ntri = (int) tr_tris->objects - ct_hullsize - ct_exteriors;
    printf("Writing %d triangles to file %s.\n", ntri, filename);
    int nv = ct_in_vrts + (tr_steiners != NULL ? tr_steiners->objects : 0);
    //nv -= ct_unused_vrts;

    //fprintf(outfile, "%d %d %d 0 0\n", nv, ntri, save_val); // nodal data
    fprintf(outfile, "%d %d 1 1 0\n", nv, ntri); // with nodal and cell data

    int i, idx=1; // UCD index starts from 1.
    for (i = 0; i < ct_in_vrts; i++) {
      //if (in_vrts[i].typ == UNUSEDVERTEX) continue;
      fprintf(outfile, "%d %g %g 0\n", idx, in_vrts[i].crd[0], in_vrts[i].crd[1]);
      in_vrts[i].idx = idx;
      idx++;
    }
    if (tr_steiners != NULL) {
      for (i = 0; i < tr_steiners->used_items; i++) {
        Vertex *vrt = (Vertex *) tr_steiners->get(i);
        if (vrt->is_deleted()) continue;
        fprintf(outfile, "%d %g %g 0\n", idx, vrt->crd[0], vrt->crd[1]);
        vrt->idx = idx;
        idx++;
      }
    }

    // UCD assumes vertex index starts from 1.
    //int shift = (io_firstindex == 1 ? 0 : 1);
    idx = 1;
    for (i = 0; i < tr_tris->used_items; i++) {
      Triang* tri = (Triang *) tr_tris->get(i);
      if (!tri->is_deleted()) {
        // ignore a hull triangle.
        if (!tri->is_hulltri() && !tri->is_exterior()) {
          fprintf(outfile, "%d %d tri %d %d %d\n", idx, tri->tag,
                  tri->vrt[0]->idx, // + shift,
                  tri->vrt[1]->idx, // + shift,
                  tri->vrt[2]->idx); // + shift);
          tri->idx = idx;
          idx++;
        }
      }
    }

    // Output nodal data
    fprintf(outfile, "1 1\n");
    fprintf(outfile, "uh, adim\n");

    idx=1;
    for (i = 0; i < ct_in_vrts; i++) {
      //if (in_vrts[i].typ == UNUSEDVERTEX) continue;
      REAL val = in_vrts[i].fval; // in_vrts[i].val;
      if (fabs(val) < 1.e-15) val = 0.0;
      fprintf(outfile, "%d %g\n", idx, val);
      idx++;
    }
    if (tr_steiners != NULL) {
      for (i = 0; i < tr_steiners->used_items; i++) {
        Vertex *vrt = (Vertex *) tr_steiners->get(i);
        if (vrt->is_deleted()) continue;
        REAL val = vrt->fval; // vrt->val;
        if (fabs(val) < 1.e-15) val = 0.0;
        fprintf(outfile, "%d %g\n", idx, val);
        idx++;
      }
    }

    // Output cell data
    fprintf(outfile, "1 1\n");
    fprintf(outfile, "energy, unknown\n");

    idx = 1;
    for (i = 0; i < tr_tris->used_items; i++) {
      Triang* tri = (Triang *) tr_tris->get(i);
      if (!tri->is_deleted()) {
        // ignore a hull triangle.
        if (!tri->is_hulltri() && !tri->is_exterior()) {
          REAL val = tri->cct[2];
          fprintf(outfile, "%d %g\n", idx, val);
          idx++;
        }
      }
    }

    fclose(outfile);

    return sum_energy;
}
*/

//==============================================================================
// see file format: detri2/doc/vtk-legacy-file-formats.pdf
// see examples: detri2/doc/vtk-file-example1.vtk, vtk-file-example2.vtk

void Triangulation::save_sol_to_vtk(int meshidx, bool save_elem_grad)
{
  if (!io_with_sol) {
    printf("Warning:  No solution is available.\n");
    return;
  }

  char filename[256];
  if (meshidx >= 0) {
    sprintf(filename, "%s_%d.vtk", io_outfilename, meshidx);
  } else {
    sprintf(filename, "%s.vtk", io_outfilename);
  }
  FILE *outfile = fopen(filename, "w");
  
  int ntri = (int) tr_tris->objects - ct_hullsize - ct_exteriors;
  int nv = ct_in_vrts + (tr_steiners != NULL ? tr_steiners->objects : 0);
  if (!io_keep_unused) { // no -IJ
    nv -= ct_unused_vrts;
  }
  printf("Writing %d vertices, %d triangles to file %s.\n", nv, ntri, filename);
  
  fprintf(outfile, "# vtk DataFile Version 2.0\n");
  fprintf(outfile, "Mesh with solution\n");
  fprintf(outfile, "ASCII\n");
  fprintf(outfile, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(outfile, "POINTS %d float\n", nv);
  int i, idx=0; // VTK index starts from 0.
  for (i = 0; i < ct_in_vrts; i++) {
    if (!io_keep_unused) { // no -IJ
      if (in_vrts[i].typ == UNUSEDVERTEX) continue;
    }
    fprintf(outfile, "%g %g 0\n", in_vrts[i].crd[0], in_vrts[i].crd[1]);
    in_vrts[i].idx = idx;
    idx++;
  }
  if (tr_steiners != NULL) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *vrt = (Vertex *) tr_steiners->get(i);
      if (vrt->is_deleted()) continue;
      fprintf(outfile, "%g %g 0\n", vrt->crd[0], vrt->crd[1]);
      vrt->idx = idx;
      idx++;
    }
  }

  fprintf(outfile, "CELLS %d %d\n", ntri, ntri*4);
  for (i = 0; i < tr_tris->used_items; i++) {
    Triang* tri = (Triang *) tr_tris->get(i);
    if (!tri->is_deleted()) {
      // ignore a hull triangle.
      if (!tri->is_hulltri() && !tri->is_exterior()) {
        fprintf(outfile, "3 %d %d %d\n",
                tri->vrt[0]->idx,
                tri->vrt[1]->idx,
                tri->vrt[2]->idx);
      }
    }
  }

  fprintf(outfile, "CELL_TYPES %d\n", ntri);
  for (i = 0; i < ntri; i++) {
    fprintf(outfile, "%d\n", 5);
  }

  fprintf(outfile, "POINT_DATA %d\n", nv);
  fprintf(outfile, "SCALARS scalars float 1\n");
  fprintf(outfile, "LOOKUP_TABLE default\n");
  for (i = 0; i < ct_in_vrts; i++) {
    if (!io_keep_unused) { // no -IJ
      if (in_vrts[i].typ == UNUSEDVERTEX) continue;
    }
    fprintf(outfile, "%g\n", in_vrts[i].fval);
  }
  if (tr_steiners != NULL) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *vrt = (Vertex *) tr_steiners->get(i);
      if (vrt->is_deleted()) continue;
      fprintf(outfile, "%g\n", vrt->fval);
    }
  }

  if (save_elem_grad) {
    REAL engery = 0.0;
    fprintf(outfile, "CELL_DATA %d\n", ntri);
    fprintf(outfile, "VECTORS vectors float\n");
    //fprintf(outfile, "LOOKUP_TABLE default\n");
    for (i = 0; i < tr_tris->used_items; i++) {
      Triang* tri = (Triang *) tr_tris->get(i);
      if (!tri->is_deleted()) {
        // ignore a hull triangle.
        if (!tri->is_hulltri() && !tri->is_exterior()) {
          engery += get_element_gradient(tri);
          fprintf(outfile, "%g %g 0.0\n",
                  tri->cct[0],
                  tri->cct[1]);
        }
      }
    }
    printf("Total energy = %.17g\n", engery);
  } // if (save_elem_grad)

  fclose(outfile);
}

//==============================================================================

void Triangulation::save_edge_energy_jump()
{
    REAL sum_energy = get_Dirichlet_energy();
    printf(" Energy = %g\n", sum_energy);

    // see file format: detri2/doc/UCD_Format.pdf
    char filename[256];
    sprintf(filename, "%s_jump.inp", io_outfilename);
    FILE *outfile = fopen(filename, "w");

    //int ntri = (int) tr_tris->objects - ct_hullsize - ct_exteriors;
    //printf("Writing %d triangles to file %s.\n", ntri, filename);
    int nv = ct_in_vrts + (tr_steiners != NULL ? tr_steiners->objects : 0);
    //nv -= ct_unused_vrts;

    int ne = (3 * (tr_tris->objects - ct_hullsize) + ct_hullsize) / 2;
    printf("Writing %d segments to file %s.\n", ne, filename);

    //fprintf(outfile, "%d %d %d 0 0\n", nv, ntri, save_val); // nodal data
    fprintf(outfile, "%d %d 0 1 0\n", nv, ne); // cell data

    int i, idx=1; // UCD index starts from 1.
    for (i = 0; i < ct_in_vrts; i++) {
      //if (in_vrts[i].typ == UNUSEDVERTEX) continue;
      fprintf(outfile, "%d %g %g 0\n", idx, in_vrts[i].crd[0], in_vrts[i].crd[1]);
      in_vrts[i].idx = idx;
      idx++;
    }
    if (tr_steiners != NULL) {
      for (i = 0; i < tr_steiners->used_items; i++) {
        Vertex *vrt = (Vertex *) tr_steiners->get(i);
        if (vrt->is_deleted()) continue;
        fprintf(outfile, "%d %g %g 0\n", idx, vrt->crd[0], vrt->crd[1]);
        vrt->idx = idx;
        idx++;
      }
    }

    // UCD assumes vertex index starts from 1.
    //int shift = (io_firstindex == 1 ? 0 : 1);
    REAL *jump_list = new REAL[ne];
    TriEdge E, N;
    int tag;
    idx = 1;
    for (i = 0; i < tr_tris->used_items; i++) {
      E.tri = (Triang *) tr_tris->get(i);
      if (E.tri->is_deleted()) continue;
      if (!E.tri->is_hulltri()) {
        for (E.ver = 0; E.ver < 3; E.ver++) {
          N = E.esym(); // it neighbor
          if (E.apex()->idx > N.apex()->idx) {
            if(N.is_segment()) {
              Triang* seg = N.get_segment();
              tag = seg->tag;
              // A boundary edge has no jump.
              jump_list[idx-1] = 0.0;
            } else {
              assert(!N.tri->is_hulltri());
              tag = 0;
              jump_list[idx-1] = fabs(sqrt(E.tri->cct[2]) - sqrt(N.tri->cct[2]));
            }
            fprintf(outfile, "%d %d line %d %d\n", idx, tag,
                    E.org()->idx,
                    E.dest()->idx);
            idx++;
          }
        }
      }
    }
    assert(idx = ne + 1);

    // Output metric on nodes or elements
    fprintf(outfile, "1 1\n");
    fprintf(outfile, "jump, unknown\n");

    for (i = 0; i < ne; i++) {
      fprintf(outfile, "%d %g\n", i+1, jump_list[i]);
    }

    fclose(outfile);
}

//==============================================================================
// https://www.geeksforgeeks.org/quick-sort/
/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
int tt_partition (Triang arr[], int low, int high)
{
    double pivot = arr[high].val; // pivot
    int i = (low - 1); // Index of smaller element
  
    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (arr[j].val > pivot)
        {
            i++; // increment index of smaller element
            //swap(&arr[i], &arr[j]);
            {
              Triang tt;
              tt.vrt[0] = arr[i].vrt[0];
              tt.vrt[1] = arr[i].vrt[1];
              tt.val = arr[i].val;
              arr[i].vrt[0] = arr[j].vrt[0];
              arr[i].vrt[1] = arr[j].vrt[1];
              arr[i].val = arr[j].val;
              arr[j].vrt[0] = tt.vrt[0];
              arr[j].vrt[1] = tt.vrt[1];
              arr[j].val = tt.val;
            }
        }
    }
    //swap(&arr[i + 1], &arr[high]);
    {
      Triang tt;
      tt.vrt[0] = arr[i + 1].vrt[0];
      tt.vrt[1] = arr[i + 1].vrt[1];
      tt.val = arr[i + 1].val;
      arr[i + 1].vrt[0] = arr[high].vrt[0];
      arr[i + 1].vrt[1] = arr[high].vrt[1];
      arr[i + 1].val = arr[high].val;
      arr[high].vrt[0] = tt.vrt[0];
      arr[high].vrt[1] = tt.vrt[1];
      arr[high].val = tt.val;
    }
    return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void tt_quickSort(Triang arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
        at right place */
        int pi = tt_partition(arr, low, high);
  
        // Separately sort elements before
        // partition and after partition
        tt_quickSort(arr, low, pi - 1);
        tt_quickSort(arr, pi + 1, high);
    }
}
//==============================================================================

bool Triangulation::get_refine_edge_list(int *p_ne, Triang **ptr_list)
{
  if (!io_with_sol) {
    return false;
  }

  // Calculate the gradients and element energy.
  REAL sum_energy = get_Dirichlet_energy();
  printf("Energy = %g\n", sum_energy);

  int ne = (3 * (tr_tris->objects - ct_hullsize) + ct_hullsize) / 2;
  Triang *edge_jump_list = new Triang[ne];

  TriEdge E, N;
  int i, idx = 1; //int tag;
  for (i = 0; i < tr_tris->used_items; i++) {
    E.tri = (Triang *) tr_tris->get(i);
    if (E.tri->is_deleted()) continue;
    if (!E.tri->is_hulltri()) {
      for (E.ver = 0; E.ver < 3; E.ver++) {
        N = E.esym(); // it neighbor
        if (E.apex()->idx > N.apex()->idx) {
          // get an edge (may be on boundary).
          edge_jump_list[idx-1].init();
          //edge_jump_list[idx-1].omt = E;
          edge_jump_list[idx-1].vrt[0] = E.org();
          edge_jump_list[idx-1].vrt[1] = E.dest();
          if(N.is_segment()) {
            //Triang* seg = N.get_segment();
            //tag = seg->tag;
            // A boundary edge has no jump.
            //jump_list[idx-1] = 0.0;
            edge_jump_list[idx-1].val = 0.0;
          } else {
            assert(!N.tri->is_hulltri());
            //tag = 0;
            edge_jump_list[idx-1].val = fabs(sqrt(E.tri->cct[2]) - sqrt(N.tri->cct[2]));
          }
          //fprintf(outfile, "%d %d line %d %d\n", idx, tag,
          //        E.org()->idx,
          //        E.dest()->idx);
          idx++;
        }
      }
    }
  }
  assert(idx = ne + 1);

  // Sort the list of edges from the highest to the lowest jumps.
  printf("Sort %d edges.\n", ne);
  tt_quickSort(edge_jump_list, 0, ne-1);

  *p_ne = ne;
  *ptr_list = edge_jump_list;

  return true;
}

//==============================================================================

void Triangulation::save_refine_edge_list() // save file .edgem
{
  Triang *edge_jump_list = NULL;
  int ne = 0;

  if (!get_refine_edge_list(&ne, &edge_jump_list)) {
    return;
  }

  /*
  // Calculate the gradients and element energy.
  REAL sum_energy = get_Dirichlet_energy();
  printf("Energy = %g\n", sum_energy);

  int ne = (3 * (tr_tris->objects - ct_hullsize) + ct_hullsize) / 2;
  Triang *edge_jump_list = new Triang[ne];

  TriEdge E, N;
  int i, idx = 1; // .mesh's index start from 1.
  for (i = 0; i < tr_tris->used_items; i++) {
    E.tri = (Triang *) tr_tris->get(i);
    if (E.tri->is_deleted()) continue;
    if (!E.tri->is_hulltri()) {
      for (E.ver = 0; E.ver < 3; E.ver++) {
        N = E.esym(); // it neighbor
        if (E.apex()->idx > N.apex()->idx) {
          // get an edge (may be on boundary).
          edge_jump_list[idx-1].init();
          edge_jump_list[idx-1].omt = E;
          if(N.is_segment()) {
            //Triang* seg = N.get_segment();
            //tag = seg->tag;
            // A boundary edge has no jump.
            //jump_list[idx-1] = 0.0;
            edge_jump_list[idx-1].val = 0.0;
          } else {
            assert(!N.tri->is_hulltri());
            //tag = 0;
            edge_jump_list[idx-1].val = fabs(sqrt(E.tri->cct[2]) - sqrt(N.tri->cct[2]));
          }
          //fprintf(outfile, "%d %d line %d %d\n", idx, tag,
          //        E.org()->idx,
          //        E.dest()->idx);
          idx++;
        }
      }
    }
  }
  assert(idx = ne + 1);

  // Sort the list of edges from the highest to the lowest jumps.
  printf("Sort %d edges.\n", ne);
  tt_quickSort(edge_jump_list, 0, ne-1);
  */

  // Save the list of sorted edges.
  // Index the vertices (to be the same as m.mesh)
  int i, idx = 1; // .mesh's index start from 1.
  for (i = 0; i < ct_in_vrts; i++) {
    if (!io_keep_unused) { // no -IJ
      if (in_vrts[i].typ == UNUSEDVERTEX) continue;
    }
    //fprintf(outfile, "%g %g 0\n", in_vrts[i].crd[0], in_vrts[i].crd[1]);
    in_vrts[i].idx = idx;
    idx++;
  }
  if (tr_steiners != NULL) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *vrt = (Vertex *) tr_steiners->get(i);
      if (vrt->is_deleted()) continue;
      //fprintf(outfile, "%g %g 0\n", vrt->crd[0], vrt->crd[1]);
      vrt->idx = idx;
      idx++;
    }
  }
  
  //int i, idx;
  //int tag;
  idx = 1;
  char filename[256];
  FILE *outfile = NULL;
  strcpy(filename, io_outfilename);
  strcat(filename, ".edgem");
  outfile = fopen(filename, "w");

  printf("Writing %d edges to file %s.\n", ne, filename);

  fprintf(outfile, "%d 0\n", ne);

  for (i = 0; i < ne; i++) {
    fprintf(outfile, "%d  %d %d  # %.17g\n", idx,
            edge_jump_list[i].vrt[0]->idx,
            edge_jump_list[i].vrt[1]->idx, edge_jump_list[i].val);
    idx++;
  }
  
  fclose(outfile);

  delete [] edge_jump_list;
}
