#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "detri2.h"

using  namespace detri2;

//==============================================================================

int Triangulation::remove_skinny_hulltris()
{
  if (op_db_verbose > 1) {
    printf("  Removing skinny hull triangles.\n");
  }
  arraypool *hull_tris = new arraypool(sizeof(TriEdge), 8);

  TriEdge E, N, nn[2], tt[4];
  int count = 0, split_seg_count = 0; // count the removed triangles, segments.
  int fflag, i;

  // Collect all hull triangles.
  for (i = 0; i < tr_tris->used_items; i++) {
    E.tri = (Triang *) tr_tris->get(i);
    if (E.tri->is_deleted()) continue;
    if (E.tri->is_hulltri()) {
      // Get its adjacent triangle.
      for (E.ver = 0; E.ver < 3; E.ver++) {
        if (E.apex() == tr_infvrt) break;
      }
      * (TriEdge *) hull_tris->alloc() = E;
    }
  }

  for (i = 0; i < hull_tris->used_items; i++) {
    E = * (TriEdge *) hull_tris->get(i);
    assert(E.apex() == tr_infvrt); // it is a hull triangle.
    //if (E.is_segment()) continue; // A segment might need to be split.
    N = E.esym();
    if (N.apex() != tr_infvrt) {
      REAL ang = get_angle(N.org(), N.dest(), N.apex());
      if (op_db_verbose > 3) {
        printf("      Check a hull triangle [%d,%d,%d], org_ang = %f degree.\n",
               N.org()->idx, N.dest()->idx, N.apex()->idx, ang / PI * 180.);
      }
      if (ang < op_tol_voronoi_minangle) {
        // Check if N can be removed by a 2-2 flip.
        // We check if N.apex() lies in between N.org() and N.dest().
        REAL L  = get_distance(N.org(), N.dest());
        REAL L1 = get_distance(N.org(), N.apex());
        if (L1 < L) {
          // A 2-2 flip is possible.
          if (op_db_verbose > 3) {
            printf("      Remove a skinny hull triangle.\n");
          }
          // Bakup the two new hull tris.
          nn[0] = N.enext().esym();
          nn[1] = N.eprev().esym();
          
          Segment *seg = NULL;
          int stag = 0; REAL val = 0.;
          
          bool isseg = E.is_segment();
          if (isseg) {
            if (op_db_verbose > 3) {
              printf("      Split a hull segment.\n");
            }
            seg = E.get_segment();
            if (seg == NULL) {
              // report a bug. vertex-to-segment map is wrong.
              clean(); // release memory;
              throw 2;
            }
            stag = seg->tag;
            val = seg->val;
            
            remove_segment(seg); // Remove this segment.
          } // if (isseg)
          
          tt[0] = E;
          fflag = FLIP_22;
          flip(tt, NULL, fflag, NULL);
          
          if (isseg) {
            insert_segment(nn[0], stag, val, NULL);
            insert_segment(nn[1], stag, val, NULL);
            split_seg_count++;
          }
          
          // Add the two new hull edges into list.
          for (int j = 0; j < 2; j++) {
            //assert(nn[j].esym().apex() == tr_infvrt);
            * (TriEdge *) hull_tris->alloc() = nn[j].esym();
          }
          count++;
        } else {
          printf("to debug...\n");
        }
      }
    }
  }

  if (op_db_verbose > 1) {
    printf("  Removed %d skinny hull triangle.\n", count);
    if (split_seg_count > 0) {
      printf("  Split %d hull segments.\n", split_seg_count);
    }
  }

  delete hull_tris;
  return 1;
}

//==============================================================================

int Triangulation::search_point(Vertex *pt, TriEdge &E, int encflag)
{
  int loc = locate_point(pt, E, encflag);
  
  if (!E.tri->is_hulltri()) {
    return loc; //LOC_IN_OUTSIDE;
  }

  // Loop through the set of hull faces, search it from them.
  int i;
  for (i = 0; i < tr_tris->used_items; i++) {
    E.tri = (Triang *) tr_tris->get(i);
    if (!E.tri->is_hulltri()) continue;
    for (E.ver = 0; E.ver < 3; E.ver++) {
      if (E.apex() == tr_infvrt) break;
    }
    assert(E.ver < 3);
    E = E.esym(); // get the triangle inside the domain.
    assert(!E.tri->is_hulltri());
    loc = locate_point(pt, E, encflag);
    if (!E.tri->is_hulltri()) {
      return loc; //LOC_IN_OUTSIDE;
    }
  }

  return LOC_IN_OUTSIDE;
}

//==============================================================================
// We search a boundary (segment) of the OMT_domain which is cut by the
//   line from In_pt to Out_pt. 'S' is an initial boundary (or hull) edge of
//   the OMT_domain. It separates the domain into two halfspaces, such that
//   S.apex() lies in the same halfspace of Out_pt.
// see doc/voronoi-edge-boundary-cut-2018-11-1.pdf

int Triangulation::get_boundary_cut_dualedge(Vertex& In_pt, Vertex& Out_pt, TriEdge& S)
{
  if (op_db_verbose > 3) {
    printf("  S.org(%d) -> S.dest(%d), S.apex(%d)\n",
           S.org()->idx, S.dest()->idx, S.apex()->idx);
    printf("  In (%g,%g) - Out (%g,%g)\n",
           In_pt.crd[0], In_pt.crd[1], Out_pt.crd[0], Out_pt.crd[1]);
  }

  REAL ss = Orient2d(S.org(), S.dest(), &In_pt);

  // assert(ss < 0); // CW orientation
  if (ss >= 0) { // debug this case
    printf("Warning: Wrong orientation of S for finding the cutting edge.\n");
    printf("  S.org()=%d, S.dest()=%d isseg(%d)\n", S.org()->idx, S.dest()->idx, S.is_segment());
    printf("  S.org [%g,%g], S.dest [%g,%g]\n",
           S.org()->crd[0], S.org()->crd[1], S.dest()->crd[0], S.dest()->crd[1]);
    printf("  Inpt [%g,%g], Outpt [%g,%g]\n", In_pt.crd[0], In_pt.crd[1], Out_pt.crd[0], Out_pt.crd[1]);
    // This might be a problem of OMT_domain, e.g., it is too small.
    return 0; //assert(0);
  }

  // maximum vertex degree.
  int max_vert_deg = (OMT_domain != NULL) ? OMT_domain->ct_in_vrts : ct_in_vrts;
  if (max_vert_deg < 1000) max_vert_deg = 1000;
  
  // The maximum possible searching steps.
  int max_iter = (OMT_domain != NULL) ? OMT_domain->ct_hullsize : ct_hullsize;
  if (max_iter < 1000) max_iter = 1000;

  int deg = 0, iter = 0;
  bool bflag = true;

  REAL s1 = Orient2d(&In_pt, &Out_pt, S.org());
  REAL s2 = Orient2d(&In_pt, &Out_pt, S.dest());
  if (op_db_verbose > 3) {
    printf("  [In, out, S.org(%d)]  s1 = %g\n", S.org()->idx, s1);
    printf("  [In, out, S.dest(%d)] s2 = %g\n", S.dest()->idx, s2);
  }

  if (s1 < 0) {
    do {
      iter++; // count a searched boundary edge.
      if (iter >= max_iter) {
        bflag = false; break; // failed to find the cut edge.
      }
      // Must go left of S.org();
      TriEdge searchEdge = S.esym();
      //printf("DBG: adjust S to the left of S.org().\n");
      // Search the next boundary edge at S.org().
      TriEdge workedge = searchEdge.enext();
      searchEdge = workedge;
      while (true) {
        if (op_db_verbose > 3) {
          printf("  workE: (%d,%d,%d)\n", workedge.org()->idx,
                 workedge.dest()->idx, workedge.apex()->idx);
        }
        //assert(workedge.org() == S.org());
        if (!(workedge.org() == S.org())) {
          //quit(QUIT_ON_BUG);
          bflag = false; break;
        }
        if (searchEdge.is_segment()) break;
        if ((searchEdge.esym()).tri->is_hulltri()) break;
        searchEdge = workedge.esym_enext(); // CW
        workedge = searchEdge;
        deg++;
        if (deg >= max_vert_deg) {
          bflag = false; break;
        }
      }
      if (!bflag) {
        break; // failed to find the next boundary edge at S.org().
      }
      S = searchEdge.esym(); // Update S.
      s1 = Orient2d(&In_pt, &Out_pt, S.org());
    } while (s1 < 0);
  } else if (s2 > 0) {
    // Muts go right of S.dest();
    //printf("DBG: adjust S to the right of S.dest().\n");
    do {
      iter++; // count a searched boundary edge.
      if (iter >= max_iter) {
        bflag = false; break; // failed to find the cut edge.
      }
      TriEdge searchEdge = S.esym();
      TriEdge workedge = searchEdge.eprev();
      searchEdge = workedge;
      while (true) {
        if (op_db_verbose > 3) {
          printf("  workE: (%d,%d,%d)\n", workedge.org()->idx,
                 workedge.dest()->idx, workedge.apex()->idx);
        }
        //assert(workedge.dest() == S.dest());
        if (!(workedge.dest() == S.dest())) {
          //quit(QUIT_ON_BUG);
          bflag = false; break;
        }
        if (searchEdge.is_segment()) break;
        if ((searchEdge.esym()).tri->is_hulltri()) break;
        searchEdge = workedge;
        workedge = searchEdge.esym_eprev();
        deg++;
        if (deg >= max_vert_deg) {
          bflag = false; break;
        }
      }
      if (!bflag) {
        break; // failed to find the next boundary edge at S.dest().
      }
      S = searchEdge.esym(); // Update S
      s2 = Orient2d(&In_pt, &Out_pt, S.dest());
    } while (s2 > 0);
  }

  if (!bflag) {
    //quit(QUIT_ON_BUG);
    return 0; // not found.
  }

  if (op_db_verbose > 3) {
    printf("  Found background boundary S: [%d,%d,%d]\n\n",
           S.org()->idx, S.dest()->idx, S.apex()->idx);
  }

  return 1;
}

//==============================================================================
// "tri" must be not a hulltri or an exterior triangle.

int Triangulation::get_tri_orthocenter(Triang* tri)
{
  if (tri->is_hulltri()) {
    return 0; // this is not a triangle.
  }

  REAL Ux, Uy, U_weight;
  REAL Vx, Vy, V_weight;
  REAL Wx, Wy, W_weight;

  Vertex *v1 = tri->vrt[0];
  Vertex *v2 = tri->vrt[1];
  Vertex *v3 = tri->vrt[2];

  Ux = v1->crd[0];
  Uy = v1->crd[1];
  Vx = v2->crd[0];
  Vy = v2->crd[1];
  Wx = v3->crd[0];
  Wy = v3->crd[1];

  // Translate heights to weights.
  U_weight = v1->wei; //Ux*Ux + Uy*Uy - v1->crd[2];
  V_weight = v2->wei; //Vx*Vx + Vy*Vy - v2->crd[2];
  W_weight = v3->wei; //Wx*Wx + Wy*Wy - v3->crd[2];

  if (op_db_verbose > 2) {
    printf(" [%d]: %g,%g,%g\n", v1->idx, Ux, Uy, U_weight);
    printf(" [%d]: %g,%g,%g\n", v2->idx, Vx, Vy, V_weight);
    printf(" [%d]: %g,%g,%g\n", v3->idx, Wx, Wy, W_weight);
  }

  if (op_db_verbose > 3) {
    // debug only
    if (is_triangle(tri, 105, 224, 355)) {
      printf("debug\n");
    }
  }
  
  int success = get_orthocenter(Ux, Uy, U_weight,
                  Vx, Vy, V_weight,
                  Wx, Wy, W_weight,
                  &tri->cct[0], &tri->cct[1], &tri->cct[2],
                  _a11, _a21, _a22);

  if (!success) {
    return 0; // failed to get orthocenter.
  }

  if (OMT_domain != NULL) {
    // Determine if this vertex is inside or outside of the OMT domain (or the subdomain
    //   which contains this triangle -- check segments).
    
    // Must initialise these flags (may be set due to previous iterations).
    tri->clear_dual_in_exterior();
    tri->clear_dual_on_bdry();
    
    Vertex Mass_pt, pt;
    Mass_pt.crd[0] = (Ux + Vx + Wx) / 3.;
    Mass_pt.crd[1] = (Uy + Vy + Wy) / 3.;
    pt.crd[0] = tri->cct[0];
    pt.crd[1] = tri->cct[1];
    // First search the mass center of "tri" in OMT_domain.
    // It must be inside the domain.
    // The background mesh (OMT_domain) may not be convex and it may contain holes,
    //   we must ensure that we can find this point.
    //int loc = OMT_domain->locate_point(&Mass_pt, tri->on, 0, 0); // encflag = 0
    //int loc = OMT_domain->locate_point(&Mass_pt, tri->on, 0); // encflag = 0
    int loc = OMT_domain->search_point(&Mass_pt, tri->omt, 0);
    if (tri->omt.tri->is_hulltri()) {
      // Failed to locate the starting triangle in OMT_domain.
      // One reason is that the given triangle is a very skinny one on the
      // convex hull (see an example: doc, page 8).
      // A possible fix: we try the centers at the neighors of this triangle.
      int j, k;
      for (k = 0; k < 3; k++) {
        if (tri->nei[k].tri->is_hulltri()) continue;
        if (tri->nei[k].is_segment()) continue;
        Vertex *v1 = tri->nei[k].tri->vrt[0];
        Vertex *v2 = tri->nei[k].tri->vrt[1];
        Vertex *v3 = tri->nei[k].tri->vrt[2];
        Mass_pt.crd[0] = (v1->crd[0] + v2->crd[0] + v3->crd[0]) / 3.;
        Mass_pt.crd[1] = (v1->crd[1] + v2->crd[1] + v3->crd[1]) / 3.;
        for (j = 0; j < 3; j++) {
          if (tri->omt.apex() == OMT_domain->tr_infvrt) break;
          tri->omt = tri->omt.enext();
        }
        tri->omt = tri->omt.esym();
        loc = OMT_domain->search_point(&Mass_pt, tri->omt, 0);
        if (!tri->omt.tri->is_hulltri()) {
          break; // Found.
        }
      }
    }
    if (!tri->omt.tri->is_hulltri()) {
      // Locate the orthocenter of "tri" from this triangle (containing its mass center).
      //loc = OMT_domain->locate_point(&pt, tri->on, 0, 1); // encflag = 1 (stop at first segment).
      loc = OMT_domain->locate_point(&pt, tri->omt, 1); // encflag = 1 (stop at first segment).
    } else {
      //loc = LOC_IN_OUTSIDE;
      return 0; // failed to locate the mass point in domain.
    }
    // Set the interior /exterior flag of the dual vertex.
    // Default, the dual (orothocenter) is in the interior.
    if (loc == LOC_IN_OUTSIDE) {
      tri->set_dual_in_exterior();
    } else if (loc == LOC_IN_TRI) {
      if (tri->omt.tri->is_exterior()) {
        tri->set_dual_in_exterior();
      }
    } else if (loc == LOC_ON_EDGE) {
      if (tri->omt.is_segment()) {
        tri->set_dual_on_bdry();
      } else {
        // Check if this is a hull edge.
        if ((tri->omt.tri->is_hulltri()) ||
            (tri->omt.esym().tri->is_hulltri())) {
          tri->set_dual_on_bdry();
        }
      }
    } else if (loc == LOC_ON_VERT) {
      // On a vertex of the background mesh, to be done.
      TriEdge E = tri->omt;
      do {
        if (tri->omt.tri->is_hulltri() || tri->omt.tri->is_exterior()) {
          tri->set_dual_in_exterior(); break;
        } else if (tri->omt.is_segment()) {
          tri->set_dual_on_bdry(); break;
        }
        tri->omt = tri->omt.eprev_esym(); // CCW rotate
        //assert(tri->omt.org() == E.org());
        if (!(tri->omt.org() == E.org())) {
          return 0;
        }
      } while (E.tri != tri->omt.tri);
    } else if (loc == LOC_ENC_SEG) {
      // This dual vertex is behind a segment, which means it lies outside of the
      //   subdomain contains the triangle. We treat it as lying in outside.
      //   This way, a cut point on the segment will be calculated for Voronoi cells.
      //assert(tri->omt.is_segment());
      if (!(tri->omt.is_segment())) {
        return 0;
      }
      tri->set_dual_in_exterior();
    }
  } // if (OMT_domain != NULL)

  return 1;
}

//==============================================================================
// This function assumes that the orthocenter of the adjacent tri to this
//   hulltri has been calculated and located. It must be called after
//   get_tri_orthocenter(tri).

int Triangulation::get_hulltri_orthocenter(Triang* hulltri)
{
  TriEdge E; E.tri = hulltri;
  for (E.ver = 0; E.ver < 3; E.ver++) {
    if (E.apex() == tr_infvrt) break;
  }
  if (!(E.ver < 3)) {
    return 0; //quit(QUIT_ON_BUG);
  }

  REAL Ux, Uy, U_weight;
  REAL Vx, Vy, V_weight;

  Vertex *v1 = E.org();
  Vertex *v2 = E.dest();

  Ux = v1->crd[0];
  Uy = v1->crd[1];
  Vx = v2->crd[0];
  Vy = v2->crd[1];

  // Translate heights to weights.
  U_weight = v1->wei; //Ux*Ux + Uy*Uy - v1->crd[2];
  V_weight = v2->wei; //Vx*Vx + Vy*Vy - v2->crd[2];

  if (op_db_verbose > 2) {
    printf("A Hull edge:\n");
    printf(" [%d]: %g,%g,%g, r(%g)\n", v1->idx, Ux, Uy, U_weight, sqrt(fabs(U_weight)));
    printf(" [%d]: %g,%g,%g, r(%g)\n", v2->idx, Vx, Vy, V_weight, sqrt(fabs(V_weight)));
  }

  if (op_db_verbose > 3) {
    // debug only
    if (((v1->idx == 105) && (v2->idx == 224)) ||
        ((v1->idx == 224) && (v2->idx == 105))) {
      printf("debug\n");
    }
  }
  
  int success = get_bissector(Ux, Uy, U_weight,
                Vx, Vy, V_weight,
                &E.tri->cct[0], &E.tri->cct[1], &E.tri->cct[2],
                _a11, _a21, _a22);

  if (!success) {
    return 0;
  }

  if (OMT_domain) {
    // we must initalise it (it may be set by previous iterations).
    hulltri->clear_dual_in_exterior();
    hulltri->clear_dual_on_bdry();
    
    // The bisector lies exactly in the edge [v1, v2].
    // We move it along the edge normal towards exterior of the triangulaiton.
    TriEdge N = E.esym();
    // The dual vertex of N.tri was already calcualted.
    if (!(N.tri->omt.tri != NULL)) {
      return 0; //quit(QUIT_ON_BUG);
    }
    if (op_db_verbose > 2) {
      printf("  Bisector on the hull edge:  [%g,%g]\n", E.tri->cct[0], E.tri->cct[1]);
      printf("  Ccenter of hull triangle N: [%g,%g]\n", N.tri->cct[0], N.tri->cct[1]);
    }

    // Get the edge normal towards exterior (rotated 90 deg).    
    double outerN[2];
    //outerN[0] =  (E.org()->crd[1] - E.dest()->crd[1]);
    //outerN[1] = -(E.org()->crd[0] - E.dest()->crd[0]);
    //double len = sqrt(outerN[0]*outerN[0] + outerN[1]*outerN[1]);
    //outerN[0] /= len;
    //outerN[1] /= len;

    // Calculate the bisector line (may not be orthogonal if metric is not Euclidean).
    // The line direction vector.
    outerN[0] =   2.*_a22*(Uy-Vy) + 2.*_a21*(Ux-Vx);
    outerN[1] = -(2.*_a11*(Ux-Vx) + 2.*_a21*(Uy-Vy));
    double len = sqrt(outerN[0]*outerN[0] + outerN[1]*outerN[1]);
    outerN[0] /= len;
    outerN[1] /= len;
    // Choose the line direction towards exterior.
    Vertex tmp_pt; // Construct a point on the bisector line.
    tmp_pt.init();
    tmp_pt.crd[0] = E.tri->cct[0] + 2 * len * outerN[0];
    tmp_pt.crd[1] = E.tri->cct[1] + 2 * len * outerN[1];
    REAL check_sign = Orient2d(E.org(), E.dest(), &tmp_pt);
    if (check_sign < 0) {
      // not ccw orientation. reverse the direction.
      outerN[0] = -outerN[0];
      outerN[1] = -outerN[1];
    } else if (check_sign == 0) {
      return 0; // quit(QUIT_ON_BUG); //assert(0); // this is a bug.
    }

    if (op_dt_nearest < 0) { // The farthest-point Voronoi diagram
      // change direction of the outer normal.
      outerN[0] = -outerN[0];
      outerN[1] = -outerN[1];
    }

    len = get_distance(E.org(), E.dest());

    if (!N.tri->is_dual_in_exterior() && !N.tri->is_dual_on_bdry()) {
      // This dual vertex lies inside the domain.
      if (!N.is_segment()) {
        // Move the dual vertex to the boundary of the OMT_domain.
        // For this we need a far point outside the OMT_domain.
        double scale = len;
        /*
        Vertex Mass_pt; // pt;
        REAL Wx = N.apex()->crd[0];
        REAL Wy = N.apex()->crd[1];
        Mass_pt.crd[0] = (Ux + Vx + Wx) / 3.;
        Mass_pt.crd[1] = (Uy + Vy + Wy) / 3.;
        TriEdge S = N.tri->on; // Must be an interior triangle of OMT_domain.
        int loc = OMT_domain->locate_point(&Mass_pt, S, 0, 0);
        */
        Vertex In_pt; // The bisector of this hull edge.
        In_pt.init();
        In_pt.crd[0] = E.tri->cct[0];
        In_pt.crd[1] = E.tri->cct[1];
        TriEdge S = N.tri->omt; // Must be an interior triangle of OMT_domain.
        if (op_db_verbose > 2) {
          printf("  N.tri->on (cct) background tri: [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
        }
        //int loc = OMT_domain->locate_point(&In_pt, S, 0);
        int loc = OMT_domain->search_point(&In_pt, S, 0);
        if (op_db_verbose > 2) {
          printf("  In_pt on background tri: [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
        }
        if (!(loc != LOC_IN_OUTSIDE)) {
          return 0; //quit(QUIT_ON_BUG);
        }

        Vertex Out_pt;
        Out_pt.init();
        while (1) { // we must find a point.
          scale = scale * 2.0;
          Out_pt.crd[0] = E.tri->cct[0] + scale * outerN[0];
          Out_pt.crd[1] = E.tri->cct[1] + scale * outerN[1];
          // Search the Out_pt, stop at the first segment.
          //loc = OMT_domain->locate_point(&Out_pt, S, 0, 1);
          loc = OMT_domain->locate_point(&Out_pt, S, 1);
          if (loc == LOC_IN_OUTSIDE) {
            //assert(S.apex() == OMT_domain->tr_infvrt);
            if (!(S.apex() == OMT_domain->tr_infvrt)) {
              return 0; // quit(QUIT_ON_BUG);
            }
            break;
          } else if (loc == LOC_IN_TRI) {
            if (S.tri->is_exterior()) break;
          } else if (loc == LOC_ENC_SEG) {
            //assert(S.is_segment());
            if (!(S.is_segment())) {
              return 0; // quit(QUIT_ON_BUG);
            }
            break; // hit a segment.
          }
        } // while (1)
        if (op_db_verbose > 3) {
          printf("  Background boundary: [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
        }
        //assert(S.is_segment()); // S is either a segment or a hull edge.
        // It is guaranteed in above.
        //if (!(S.is_segment() || (S.apex() != OMT_domain->tr_infvrt))) {
        //  quit(QUIT_ON_BUG);
        //}
        // Get the boundary edge S in OMT such that [S.org, S.dest] cuts the
        //   line [In_pt, Out_pt].
        if (!get_boundary_cut_dualedge(In_pt, Out_pt, S)) {
          return 0; // failed to get the cut boundary edge.
        }
        E.tri->omt = S; // Remember this triangle in OMT_domain.

        // Calculate the cut point.
        Vertex *e1 = E.tri->omt.org();
        Vertex *e2 = E.tri->omt.dest();
        double X0 = N.tri->cct[0];
        double Y0 = N.tri->cct[1];
        double X1 = Out_pt.crd[0];
        double Y1 = Out_pt.crd[1];
        double X2 = e1->crd[0];
        double Y2 = e1->crd[1];
        double X3 = e2->crd[0];
        double Y3 = e2->crd[1];
        double t1, t2;
        if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
          E.tri->cct[0] = X0 + t1 * (X1 - X0);
          E.tri->cct[1] = Y0 + t1 * (Y1 - Y0);
          // Calulcate rr2. (the weight)
          // The weighted point of this cutpoint must lie on the polar plane
          //   of the weighted point V1 (or V2).
          // The weight of this cut-point can be found as following:
          //   - project this point (E.tri->cct[0], E.tri->cct[1]) vertically to
          //     the tangent plane passing through the lifted point of V1 (or V2),
          //   - then caculate the vertical distance from projection point to the
          //     lifted point of the paraboloid.
          // Let a = e1->crd[0], b = e1->crd[1], wei=e1->wei
          // The tangent plane at e1' is:
          //    z = 2ax + 2by - (a^2+b^2) + wei
          // Let \alpha = E.tri->cct[0], \beta= E.tri->cct[1]
          //    z_proj = 2a(\alpha)+2b(\beta)-(a^2+b^2) + wei
          // Then the weight of the cut-point is:
          //    (\alpha * \alpha) + (\beta * \beta) - z_proj
          double a = Ux; // e1->crd[0];
          double b = Uy; // e1->crd[1];
          double alpha = E.tri->cct[0];
          double beta  = E.tri->cct[1];
          double z_proj = 2.*a*alpha+2.*b*beta - (a*a+b*b) + U_weight;
          E.tri->cct[2] = alpha*alpha+beta*beta - z_proj; // its weight
        } else {
          //printf("!! Warning: Failed at calculating line-line intersection.\n");
          return 0;
        }
      } // if (!N.is_segment()
      else {
        // [2019-07-28] Remember this boundary segment.
        E.tri->omt = N.esym();
      }
      // It is either on a segment or a hull edge of OMT_domain.
      E.tri->set_dual_on_bdry();
    } else {
      // The dual vertex of its neighbor lies outside.
      // Move this dual vertex further than the neighbor dual vertex.
      // This is only for visualization useful.
      E.tri->cct[0] = N.tri->cct[0] + len * outerN[0];
      E.tri->cct[1] = N.tri->cct[1] + len * outerN[1];
      E.tri->cct[2] = 0; // will be calulcated later.
      E.tri->omt = N.tri->omt;
      E.tri->set_dual_in_exterior();
    }
    if (op_db_verbose > 10) {
      printf("p:show_vector(%g,%g,0,%g,%g,0)\n", N.tri->cct[0], N.tri->cct[1],
             E.tri->cct[0], E.tri->cct[1]);
    }
  } // if (OMT_domain)

  return 1;
}

//==============================================================================

int Triangulation::get_cct_location(TriEdge &T, bool hull_flag)
{
  Triang *tri = T.tri;

  // Must initialise these flags (may be set due to previous iterations).
  tri->clear_dual_in_exterior();
  tri->clear_dual_on_bdry();

  Vertex pt;
  pt.crd[0] = tri->cct[0];
  pt.crd[1] = tri->cct[1];

  int loc = LOC_IN_OUTSIDE;

  if (!hull_flag) {
    // tri is an interior triangle.
    tri->omt.tri = tri; // point to itself.
    loc = locate_point(&pt, tri->omt, 1);

    // Set the interior /exterior flag of the dual vertex.
    // Default, the dual (orothocenter) is in the interior.
    if (loc == LOC_IN_OUTSIDE) {
      tri->set_dual_in_exterior();
    } else if (loc == LOC_IN_TRI) {
      if (tri->omt.tri->is_exterior()) {
        tri->set_dual_in_exterior();
      }
    } else if (loc == LOC_ON_EDGE) {
      if (tri->omt.is_segment()) {
        tri->set_dual_on_bdry();
      } else {
        // Check if this is a hull edge.
        if ((tri->omt.tri->is_hulltri()) ||
            (tri->omt.esym().tri->is_hulltri())) {
          tri->set_dual_on_bdry();
        }
      }
    } else if (loc == LOC_ON_VERT) {
      // On a vertex of the mesh -- only possible when this mesh is not a CDT.
      // This circumcenter must be an interior vertex (not crossing any segment).
      if (op_db_verbose) {
        printf("Warning: non-CDT. Circumcenter on vertex.\n");
      }
    } else if (loc == LOC_ENC_SEG) {
      // This dual vertex is behind a segment, which means it lies outside of the
      //   subdomain contains the triangle. We treat it as lying in outside.
      //   This way, a cut point on the segment will be calculated for Voronoi cells.
      //assert(tri->omt.is_segment());
      if (!(tri->omt.is_segment())) {
        return 0; // this should be a bug.
      }
      tri->set_dual_in_exterior();
    }
  } // if (!hull_flag)
  else {
    // The bisector lies exactly in the edge [v1, v2].
    // We move it along the edge normal towards exterior of the triangulaiton.
    TriEdge E = T;
    if (!E.tri->is_hulltri() || !E.is_segment()) {
      return 0; // it is not a hull tri or a segment.
    }
    
    TriEdge N = E.esym();
    if (!N.tri->is_dual_in_exterior() && !N.tri->is_dual_on_bdry()) {
      // It is either on a segment or a hull edge of OMT_domain.
      E.tri->set_dual_on_bdry();
    } else {
      E.tri->set_dual_in_exterior();
    }
  } // if (hull_flag)

  return 1;
}

//==============================================================================
// Both dual vertices of "mesh_edge" are in the exiterior of the domain.
// Check if this dual edge crosses the interior of the domain.
// Return 1 if it is, In_pt returns an interior point of this dual edge.

int Triangulation::get_dual_edge_interior_point(TriEdge mesh_edge, Vertex *In_pt)
{
  Vertex *e1 = mesh_edge.org();
  Vertex *e2 = mesh_edge.dest();
 
  if ((e1 == tr_infvrt) || (e2 == tr_infvrt)) {
    return 0; // a hull edfge.
  }
 
  TriEdge E = mesh_edge;
  TriEdge N = E.esym();

  if (E.tri->omt.tri == N.tri->omt.tri) {
    // These two endpoints of this dual edge lie on the same exterior boundary
    //   of this domain. We think that this dual edge lies completely in exterior.
    return 0;
  }

  Vertex V1, V2;
  V1.init();
  V2.init();

  V1.crd[0] = E.tri->cct[0];
  V1.crd[1] = E.tri->cct[1];
  V2.crd[0] = N.tri->cct[0];
  V2.crd[1] = N.tri->cct[1];
  
  double t1, t2;

  int success = line_line_intersection(
    V1.crd[0], V1.crd[1], V2.crd[0], V2.crd[1],
    e1->crd[0], e1->crd[1], e2->crd[0], e2->crd[1],
    &t1, &t2);

  if (!success) {
    return 0;
  }

  if ((t1 > 0.) && (t1 < 1.)) {
    In_pt->crd[0] = V1.crd[0] + t1 * (V2.crd[0] - V1.crd[0]);
    In_pt->crd[1] = V1.crd[1] + t1 * (V2.crd[1] - V1.crd[1]);

    // Check if this point lies inside of the domain.
    TriEdge S_omt = E.tri->omt.esym(); // Starting from an interior domain triangle.
    int loc = locate_point(In_pt, S_omt, 1); // encflag = 1

    if ((loc != LOC_IN_OUTSIDE) &&
        (loc != LOC_ENC_SEG)) {
      In_pt->on_omt = S_omt;
      return 1;
    }
  }

  return 0;
}

//==============================================================================
// It assumes that a background domain is provided, i.e., OMT_domain != NULL.
//
// In the simplest case, the dual (Voronoi) edge is a line segment, and
//  the two endponts of this line segment are circumcenters (or orthocenters)
//  the two Delaunay triangles sharing at "mesh_edge".
//
// In general, this dual (Voronoi) edge may cut the boundary of the domain.
//   The two endpoints of the returned line segment might be cut points.
//
// If the background domain (OMT_domain) is non-convex, it is possible that
//   a dual edge is cut by the boundary of OMT_domain in many disconnected
//   line segments.
//
// The first line segment of the dual edge is returned in "dual_edge[2]".
// If there are multiple cut line segments, they are returned in the array
//   "pptlist" with length "ptnum", where "pptlist[0]" and "pptlist[1]" is the
//   first line segment, etc.  This array is created inside this function.

int Triangulation::get_dual_edge(TriEdge mesh_edge, Vertex *dual_edge,
                                 Vertex** pptlist, int* ptnum)
{
  TriEdge E = mesh_edge;
  //if (E.is_segment()) return 0; // Skip a segment.
  if (E.org() == tr_infvrt ||
      E.dest() == tr_infvrt) return 0; // Skip a hull edge.

  TriEdge N = E.esym();

  if (op_db_verbose > 3) {
    printf("  Get edge: [%d,%d] - [%d,%d].\n",
           E.org()->idx, E.dest()->idx, E.apex()->idx, N.apex()->idx);
  }
  // Initialize the return values.
  //assert(dual_edge != NULL);
  dual_edge[0].init();
  dual_edge[1].init();
  
  //*pptlist = NULL;
  //*ptnum = 0;

  double a = E.org()->crd[0];
  double b = E.org()->crd[1];
  double wei = E.org()->wei;
  //double alpha, beta, z_proj;

  if (E.tri->is_dual_in_exterior() || E.tri->is_dual_on_bdry()) {
    if (N.tri->is_dual_in_exterior() || N.tri->is_dual_on_bdry()) {
      Vertex In_pt;
      if (get_dual_edge_interior_point(E, &In_pt)) {
        // Get the intersect(s) between this dual edge and OMT_domain.
        // Part of this dual edge lies inside of the domain.
        // We calculate two cut points (see Page 6,7,8).
        //Vertex Cut_pt1, Cut_pt2;
          
        // We're visting the triangles of mesh_vertex in CCW order.
        // The dual edge is directed from N.tri->cct --> E.tri->cct.
        // The first cut point is between N.tri->cct (out) and In_pt.
        Vertex Out_pt;
        Out_pt.init();
        Out_pt.crd[0] = N.tri->cct[0];
        Out_pt.crd[1] = N.tri->cct[1];
          
        TriEdge S_omt = N.tri->omt;
        if (!get_boundary_cut_dualedge(In_pt, Out_pt, S_omt)) { // see Page 3, S = {b1, b2}.
          // Failed to get this boundary edge which cuts the dual edge.
          //assert(0); // should be a bug.
          return 0;
        }

        // Calculate the cut vertex, see Page 3, cut1.
        N.tri->omt = S_omt; // Update E.tri->on
        //printf("  Last Exit [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
          
        Vertex *e1 = N.tri->omt.org();
        Vertex *e2 = N.tri->omt.dest();

        // Sort the vertices to make sure that the cutting point is exactly
        //   the same when we calculate the adjacent cell.
        double X0,Y0,X1,Y1,X2,Y2,X3,Y3;
        if (E.apex()->idx < N.apex()->idx) {
          X0 = E.tri->cct[0];
          Y0 = E.tri->cct[1];
          X1 = N.tri->cct[0];
          Y1 = N.tri->cct[1];
        } else {
          X1 = E.tri->cct[0];
          Y1 = E.tri->cct[1];
          X0 = N.tri->cct[0];
          Y0 = N.tri->cct[1];
        }
        if (e1->idx < e2->idx) {
          X2 = e1->crd[0];
          Y2 = e1->crd[1];
          X3 = e2->crd[0];
          Y3 = e2->crd[1];
        } else {
          X3 = e1->crd[0];
          Y3 = e1->crd[1];
          X2 = e2->crd[0];
          Y2 = e2->crd[1];
        }
        double t1 = 0, t2 = 0;
        if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
          double c0 = X0 + t1 * (X1 - X0);
          double c1 = Y0 + t1 * (Y1 - Y0);
          // Calculate the weight
          double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
          double c2 = c0*c0+c1*c1 - z_proj;

          Vertex *pt = &(dual_edge[0]);
          pt->init();
          pt->crd[0] = c0;
          pt->crd[1] = c1;
          pt->crd[2] = c2; // its weight
          pt->tag = -1; // A cut vertex.
        } else {
          return 0; // quit(QUIT_ON_BUG);
        }

        // The second cut point is between In_pt and E.tri->cct (out).
        Out_pt.crd[0] = E.tri->cct[0];
        Out_pt.crd[1] = E.tri->cct[1];
          
        S_omt = E.tri->omt;
        if (!get_boundary_cut_dualedge(In_pt, Out_pt, S_omt)) { // see Page 3, S = {b1, b2}.
          // Failed to get this boundary edge which cuts the dual edge.
          //assert(0); // should be a bug.
          return 0;
        }

        // Calculate the cut vertex, see Page 3, cut1.
        E.tri->omt = S_omt; // Update E.tri->on
        //printf("  Last Exit [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
          
        e1 = E.tri->omt.org();
        e2 = E.tri->omt.dest();

        if (E.apex()->idx < N.apex()->idx) {
          X0 = E.tri->cct[0];
          Y0 = E.tri->cct[1];
          X1 = N.tri->cct[0];
          Y1 = N.tri->cct[1];
        } else {
          X1 = E.tri->cct[0];
          Y1 = E.tri->cct[1];
          X0 = N.tri->cct[0];
          Y0 = N.tri->cct[1];
        }
        if (e1->idx < e2->idx) {
          X2 = e1->crd[0];
          Y2 = e1->crd[1];
          X3 = e2->crd[0];
          Y3 = e2->crd[1];
        } else {
          X3 = e1->crd[0];
          Y3 = e1->crd[1];
          X2 = e2->crd[0];
          Y2 = e2->crd[1];
        }
        t1 = 0, t2 = 0;

        if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
          double c0 = X0 + t1 * (X1 - X0);
          double c1 = Y0 + t1 * (Y1 - Y0);
          // Calculate the weight
          double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
          double c2 = c0*c0+c1*c1 - z_proj;
    
          Vertex *pt = &(dual_edge[1]);
          pt->init();
          pt->crd[0] = c0;
          pt->crd[1] = c1;
          pt->crd[2] = c2; // its weight
          pt->tag = -1; // A cut vertex.
        } else {
          return 0; // quit(QUIT_ON_BUG);
        }
      } else {
        return 0;
      } // if (!get_dual_edge_interior_point(E, &In_pt))
    } else {
      // N.tri's dual vertex lies inside the domain. (see Page 3, v2->v3)
      // Calculating a cut vertex between this dual edge (from N.tri->cct to E.tri->cct)
      //   and a boundary edge (or a segment) S of the background mesh.
      //   We must ensure that the dual edge and S intersect in S's interior.
      Vertex In_pt, Out_pt;
      In_pt.init();
      Out_pt.init();
      In_pt.crd[0]  = N.tri->cct[0]; // see Page 3, v2
      In_pt.crd[1]  = N.tri->cct[1];
      Out_pt.crd[0] = E.tri->cct[0]; // see Page 3, v3
      Out_pt.crd[1] = E.tri->cct[1];
      
      TriEdge Scut_omt = E.tri->omt;
      if (!get_boundary_cut_dualedge(In_pt, Out_pt, Scut_omt)) { // see Page 3, S = {b1, b2}.
        // Failed to get this boundary edge which cuts the dual edge.
        //assert(0); // should be a bug.
        //delete [] ptlist;
        return 0;
      }
      // Calculate the cut vertex, see Page 3, cut1.
      E.tri->omt = Scut_omt; // Update E.tri->on
      
      Vertex *e1 = E.tri->omt.org();
      Vertex *e2 = E.tri->omt.dest();
      
      // Sort the vertices to make sure that the cutting point is exactly
      //   the same when we calculate the adjacent cell.
      double X0,Y0,X1,Y1,X2,Y2,X3,Y3;
      double t1 = 0, t2 = 0;

      if (E.apex()->idx < N.apex()->idx) {
        X0 = E.tri->cct[0];
        Y0 = E.tri->cct[1];
        X1 = N.tri->cct[0];
        Y1 = N.tri->cct[1];
      } else {
        X1 = E.tri->cct[0];
        Y1 = E.tri->cct[1];
        X0 = N.tri->cct[0];
        Y0 = N.tri->cct[1];
      }
      if (e1->idx < e2->idx) {
        X2 = e1->crd[0];
        Y2 = e1->crd[1];
        X3 = e2->crd[0];
        Y3 = e2->crd[1];
      } else {
        X3 = e1->crd[0];
        Y3 = e1->crd[1];
        X2 = e2->crd[0];
        Y2 = e2->crd[1];
      }
      
      if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
        double c0 = X0 + t1 * (X1 - X0);
        double c1 = Y0 + t1 * (Y1 - Y0);
        // Calculate the weight
        double a = E.org()->crd[0];
        double b = E.org()->crd[1];
        double wei = E.org()->wei;
        double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
        double c2 = c0*c0+c1*c1 - z_proj;
        
        Vertex *pt = &(dual_edge[0]);
        pt->init();
        pt->crd[0] = c0;
        pt->crd[1] = c1;
        pt->crd[2] = c2; // its weight
        pt->tag = -1; // It is a cut vertex.
        pt = &(dual_edge[1]);
        pt->init();
        pt->crd[0] = N.tri->cct[0];
        pt->crd[1] = N.tri->cct[1];
        pt->crd[2] = N.tri->cct[2]; // its weight
        pt->tag = N.tri->idx;
      } else {
        // Failed to get the cut-point (Line-Line intersection).
        return 0; // quit(QUIT_ON_BUG);
      }
    }
  } else {
    // E.tri->cct is an interier point.
    if (N.tri->is_dual_in_exterior() || N.tri->is_dual_on_bdry()) {
      Vertex In_pt, Out_pt;
      In_pt.init();
      Out_pt.init();
      In_pt.crd[0]  = E.tri->cct[0]; // see Page 3, v2
      In_pt.crd[1]  = E.tri->cct[1];
      Out_pt.crd[0] = N.tri->cct[0]; // see Page 3, v3
      Out_pt.crd[1] = N.tri->cct[1];
      
      TriEdge Scut_omt = N.tri->omt;
      if (!get_boundary_cut_dualedge(In_pt, Out_pt, Scut_omt)) { // see Page 3, S = {b1, b2}.
        // Failed to get this boundary edge which cuts the dual edge.
        //assert(0); // should be a bug.
        //delete [] ptlist;
        return 0;
      }
      // Calculate the cut vertex, see Page 3, cut1.
      N.tri->omt = Scut_omt; // Update E.tri->on
      
      Vertex *e1 = N.tri->omt.org();
      Vertex *e2 = N.tri->omt.dest();
      
      // Sort the vertices to make sure that the cutting point is exactly
      //   the same when we calculate the adjacent cell.
      double X0,Y0,X1,Y1,X2,Y2,X3,Y3;
      double t1 = 0, t2 = 0;

      if (E.apex()->idx < N.apex()->idx) {
        X0 = E.tri->cct[0];
        Y0 = E.tri->cct[1];
        X1 = N.tri->cct[0];
        Y1 = N.tri->cct[1];
      } else {
        X1 = E.tri->cct[0];
        Y1 = E.tri->cct[1];
        X0 = N.tri->cct[0];
        Y0 = N.tri->cct[1];
      }
      if (e1->idx < e2->idx) {
        X2 = e1->crd[0];
        Y2 = e1->crd[1];
        X3 = e2->crd[0];
        Y3 = e2->crd[1];
      } else {
        X3 = e1->crd[0];
        Y3 = e1->crd[1];
        X2 = e2->crd[0];
        Y2 = e2->crd[1];
      }
      
      if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
        double c0 = X0 + t1 * (X1 - X0);
        double c1 = Y0 + t1 * (Y1 - Y0);
        // Calculate the weight
        double a = E.org()->crd[0];
        double b = E.org()->crd[1];
        double wei = E.org()->wei;
        double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
        double c2 = c0*c0+c1*c1 - z_proj;
        
        Vertex *pt = &(dual_edge[0]);
        pt->init();
        pt->crd[0] = c0;
        pt->crd[1] = c1;
        pt->crd[2] = c2; // its weight
        pt->tag = -1; // It is a cut vertex.
        pt = &(dual_edge[1]);
        pt->init();
        pt->crd[0] = E.tri->cct[0];
        pt->crd[1] = E.tri->cct[1];
        pt->crd[2] = E.tri->cct[2]; // its weight
        pt->tag = E.tri->idx;
      } else {
        // Failed to get the cut-point (Line-Line intersection).
        return 0; //quit(QUIT_ON_BUG);
      }
    } else {
      // Both circumcenters are in the interior.
      Vertex *pt = &(dual_edge[0]);
      pt->init();
      pt->crd[0] = N.tri->cct[0];
      pt->crd[1] = N.tri->cct[1];
      pt->crd[2] = N.tri->cct[2]; // its weight
      pt->tag = N.tri->idx;
      pt = &(dual_edge[1]);
      pt->init();
      pt->crd[0] = E.tri->cct[0];
      pt->crd[1] = E.tri->cct[1];
      pt->crd[2] = E.tri->cct[2]; // its weight
      pt->tag = E.tri->idx;
    }
  }

  return 1;
}

//==============================================================================
// The power cell of an interior mesh vertex is given by the convex hull of a
//   set of corners, which are Voronoi vertices and cutting vertices (of Voronoi
//   edges with domain boundary edges). The set of corners are returned by
//   an array ``pptlist" ordered in counterclockwise direction.
//
// This functions assumes that a ``OMT_domain" is provided, i.e.,
//   OMT_domain != NULL. It might be the triangulation itself.
//
// If ``calc_cct_flag = true", the orthocenters of (hull) triangles at the
//   "mesh_vertex" will be calculated, and their locations w.r.t. OMT_domain
//   will be determined. Otherwise, these calculations are already done.
//
// see doc/doc-get_powercell.key

int Triangulation::get_powercell(
  Vertex *mesh_vertex,
  bool calc_cct_flag,
  Vertex** pptlist, int* ptnum)
{
  // A background mesh must be provided (maybe itself).
  //assert(OMT_domain != NULL);
  *pptlist = NULL;
  *ptnum = 0;

  TriEdge E = mesh_vertex->adj;
  if (E.tri == NULL) {
    // This vertex does not has an adjacent triangle.
    // It is either an unused vertex or this triangulation is empty.
    return 0;
  }
  if (E.tri->is_deleted()) {
    // This should be a bug.
    //assert(0);
    return 0; //quit(QUIT_ON_BUG);
  }

  int vcount = 0; // Count the number of Voronoi vertices.
  do {
    vcount++;
    //if (E.is_segment()) { //if (E.tri->is_hulltri()) {
    //  if (op_db_verbose > 3) {
    //    printf("  Vertex %d is a boundary (fixed) vertex. Skipped.\n", mesh_vertex->idx);
    //  }
    //  return 0;
    //}
    E = E.eprev_esym(); // ccw
  } while (E.tri != mesh_vertex->adj.tri);

  if (op_db_verbose > 1) {
    printf("  Found %d adjacent triangles at vertex %d.\n", vcount, mesh_vertex->idx);
  }

  if (calc_cct_flag) {
    // Calculate the circumcenters of triangles at this vertex.
    TriEdge E = mesh_vertex->adj;
    do {
      if (!E.tri->is_hulltri()) {
        if (!get_tri_orthocenter(E.tri)) {
          return 0; // failed in calcuation cct.
        }
      }
      E = E.eprev_esym(); // ccw
    } while (E.tri != mesh_vertex->adj.tri);

    E = mesh_vertex->adj;
    do {
      if (E.tri->is_hulltri()) {
        if (!get_hulltri_orthocenter(E.tri)) {
          return 0; // failed in calcuation cct.
        }
      }
      E = E.eprev_esym(); // ccw
    } while (E.tri != mesh_vertex->adj.tri);
  } // if (calc_cct_flag)

  // N_Last_Exit is a boundary edge or a hull edge ([a,b,-1]) of the OMT_domain.
  //   If it is a boundary edge which is not a hull edge, i.e., [a,b,c], c != -1,
  //   it must be an interior segment which separates two subdomains. Then the
  //   vertex 'c' must be in the other subdomain different to this power cell.
  TriEdge N, N_Last_Exit; // [a,b,-1] or [a,b,c]
  N_Last_Exit.tri = NULL;

  // We need to start from an interior dual edge -- need one interior Voronoi vertex.
  // So the initial condition (N_Last_Exit.tri == NULL) is valid.
  // [2019-07-28] We must also check segment.
  //printf("\n DBG start: \n");
  E = mesh_vertex->adj;
  N = E.esym();
  bool is_interor_dual_edge = (!N.tri->is_dual_in_exterior() && !N.tri->is_dual_on_bdry());
  if (!is_interor_dual_edge) {
    do {
      //printf("  E: [%d,%d,%d] %d\n", E.org()->idx, E.dest()->idx, E.apex()->idx, E.triis_dual_in_exterior());
      //printf("  N: [%d,%d,%d] %d\n", N.org()->idx, N.dest()->idx, N.apex()->idx, N.triis_dual_in_exterior());
      E = E.eprev_esym(); // ccw
      N = E.esym();
      is_interor_dual_edge = (!N.tri->is_dual_in_exterior() && !N.tri->is_dual_on_bdry());
      if (is_interor_dual_edge) break;
    } while (E.tri != mesh_vertex->adj.tri);
    mesh_vertex->adj = E; // Must update this.
  }

  // It is possible that all Voronoi vertices are in outside of the domain.
  //if (!is_interor_dual_edge) {
  //  // To fix: the Voronoi cell of this vertex should not be empty
  //  //   (see doc Page 6, vertex 9).
  //  if (op_db_verbose > 3) {
  //    printf("Warning:  All dual edges of this cell lie ouside of the OMT_domain.\n");
  //  }
  //  //return 0;
  //}

  // The weight of this cut-point can be found as following:
  //   - project this point (E.tri->cct[0], E.tri->cct[1]) vertically to
  //     the tangent plane passing through the lifted point of e1 (or e2),
  //   - then caculate the vertical distance from projection point to the
  //     lifted point of the paraboloid.
  // Let a = mesh_vertex->crd[0], b = mesh_vertex->crd[1], wei=mesh_vertex->wei
  // The polar (tangent) plane at of the lifted mesh_vertex' is:
  //    z = 2ax + 2by - (a^2+b^2) + wei
  // Let \alpha = cutpoint->crd[0], \beta= cutpoint->crd[1]
  //    z_proj = 2a(\alpha)+2b(\beta)-(a^2+b^2) + wei
  // Then the weight of the cut-point is:
  //    (\alpha * \alpha) + (\beta * \beta) - z_proj
  double a = mesh_vertex->crd[0];
  double b = mesh_vertex->crd[1];
  double wei = mesh_vertex->wei;
  //double alpha, beta, z_proj;

  // Use an array to save the list of Voronoi or cut vertices (in ccw order).
  // [Remark] Due to the possible cut vertcies, the total number of vertcies
  //   of a power cell may be more than vcount. However, it must be lower than
  //   vcount + hullsize (of OMT)
  Vertex *ptlist = new Vertex[vcount + OMT_domain->ct_hullsize];
  int ptcount = 0;

  do {
    if (op_db_verbose > 3) {
      printf("  Get tri: [%d,%d,%d].\n", E.org()->idx, E.dest()->idx, E.apex()->idx);
    }
    N = E.esym();

    if (E.tri->is_dual_in_exterior()) {
      if (N.tri->is_dual_in_exterior() || N.tri->is_dual_on_bdry()) {
        // We're walking in OUTSIDE (see Page 3, v3->v4).
        Vertex In_pt;
        bool is_cut_domain = get_dual_edge_interior_point(E, &In_pt);
        if (is_cut_domain) {
          // Part of this dual edge lies inside of the domain.
          // We calculate two cut points (see Page 6,7,8).
          Vertex Cut_pt1, Cut_pt2;
          
          // We're visting the triangles of mesh_vertex in CCW order.
          // The dual edge is directed from N.tri->cct --> E.tri->cct.
          // The first cut point is between N.tri->cct (out) and In_pt.
          Vertex Out_pt;
          Out_pt.init();
          Out_pt.crd[0] = N.tri->cct[0];
          Out_pt.crd[1] = N.tri->cct[1];
          
          TriEdge S_omt = N.tri->omt;
          if (!get_boundary_cut_dualedge(In_pt, Out_pt, S_omt)) { // see Page 3, S = {b1, b2}.
            // Failed to get this boundary edge which cuts the dual edge.
            //assert(0); // should be a bug.
            delete [] ptlist;
            return 0;
          }

          // Calculate the cut vertex, see Page 3, cut1.
          N.tri->omt = S_omt; // Update E.tri->on
          //printf("  Last Exit [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
          
          Vertex *e1 = N.tri->omt.org();
          Vertex *e2 = N.tri->omt.dest();

          // Sort the vertices to make sure that the cutting point is exactly
          //   the same when we calculate the adjacent cell.
          double X0,Y0,X1,Y1,X2,Y2,X3,Y3;
          if (E.apex()->idx < N.apex()->idx) {
            X0 = E.tri->cct[0];
            Y0 = E.tri->cct[1];
            X1 = N.tri->cct[0];
            Y1 = N.tri->cct[1];
          } else {
            X1 = E.tri->cct[0];
            Y1 = E.tri->cct[1];
            X0 = N.tri->cct[0];
            Y0 = N.tri->cct[1];
          }
          if (e1->idx < e2->idx) {
            X2 = e1->crd[0];
            Y2 = e1->crd[1];
            X3 = e2->crd[0];
            Y3 = e2->crd[1];
          } else {
            X3 = e1->crd[0];
            Y3 = e1->crd[1];
            X2 = e2->crd[0];
            Y2 = e2->crd[1];
          }
          double t1 = 0, t2 = 0;
          if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
            double c0 = X0 + t1 * (X1 - X0);
            double c1 = Y0 + t1 * (Y1 - Y0);
            // Calculate the weight
            double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
            double c2 = c0*c0+c1*c1 - z_proj;
            /*
            x1 = (int) ( Sxy * N.tri->cct[0] + Cx);
            y1 = (int) (-Sxy * N.tri->cct[1] + Cy);
            x2 = (int) ( Sxy * c0 + Cx);
            y2 = (int) (-Sxy * c1 + Cy);
            painter->drawLine(QPoint(x1, y1), QPoint(x2,y2));
            */
            Vertex *pt = &Cut_pt1;
            pt->init();
            pt->crd[0] = c0;
            pt->crd[1] = c1;
            pt->crd[2] = c2; // its weight
            pt->tag = -1; // A cut vertex.
            //pt->idx = ptcount+1;
            //ptcount++;
          } else {
            quit(QUIT_ON_BUG);
          }

          // The second cut point is between In_pt and E.tri->cct (out).
          Out_pt.crd[0] = E.tri->cct[0];
          Out_pt.crd[1] = E.tri->cct[1];
          
          S_omt = E.tri->omt;
          if (!get_boundary_cut_dualedge(In_pt, Out_pt, S_omt)) { // see Page 3, S = {b1, b2}.
            // Failed to get this boundary edge which cuts the dual edge.
            //assert(0); // should be a bug.
            delete [] ptlist;
            return 0;
          }

          // Calculate the cut vertex, see Page 3, cut1.
          E.tri->omt = S_omt; // Update E.tri->on
          //printf("  Last Exit [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
          
          e1 = E.tri->omt.org();
          e2 = E.tri->omt.dest();

          if (E.apex()->idx < N.apex()->idx) {
            X0 = E.tri->cct[0];
            Y0 = E.tri->cct[1];
            X1 = N.tri->cct[0];
            Y1 = N.tri->cct[1];
          } else {
            X1 = E.tri->cct[0];
            Y1 = E.tri->cct[1];
            X0 = N.tri->cct[0];
            Y0 = N.tri->cct[1];
          }
          if (e1->idx < e2->idx) {
            X2 = e1->crd[0];
            Y2 = e1->crd[1];
            X3 = e2->crd[0];
            Y3 = e2->crd[1];
          } else {
            X3 = e1->crd[0];
            Y3 = e1->crd[1];
            X2 = e2->crd[0];
            Y2 = e2->crd[1];
          }
          t1 = 0, t2 = 0;

          if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
            double c0 = X0 + t1 * (X1 - X0);
            double c1 = Y0 + t1 * (Y1 - Y0);
            // Calculate the weight
            double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
            double c2 = c0*c0+c1*c1 - z_proj;
            /*
            x1 = (int) ( Sxy * N.tri->cct[0] + Cx);
            y1 = (int) (-Sxy * N.tri->cct[1] + Cy);
            x2 = (int) ( Sxy * c0 + Cx);
            y2 = (int) (-Sxy * c1 + Cy);
            painter->drawLine(QPoint(x1, y1), QPoint(x2,y2));
            */
            Vertex *pt = &Cut_pt2;
            pt->init();
            pt->crd[0] = c0;
            pt->crd[1] = c1;
            pt->crd[2] = c2; // its weight
            pt->tag = -1; // A cut vertex.
            //pt->idx = ptcount+1;
            //ptcount++;
          } else {
            quit(QUIT_ON_BUG);
          }

          bool clean_flag = false;
          if (N_Last_Exit.tri == NULL) {
            // This case is due to that all Vronoi vertices are in exterior.
            // This is the only cut line segment of this Voronoi cell  (see Page 8).
            N_Last_Exit = E.tri->omt;
            clean_flag = true;
          }

          if ((N.tri->omt.tri != N_Last_Exit.tri) ||
              ((N.tri->omt.tri == N_Last_Exit.tri) &&
               (N.tri->omt.ver != N_Last_Exit.ver))) {
            // They are on different boundary edges.
            // (see Page 2, the bdry edges containing v2 and v3 are different.)
            // There might be corner vertices to be the vertices of this cell.
            //printf("  !!! Searching background corners (2) DEBUG!!!\n");
            // [Important] The search is by rotating the triangles in the inside of the
            //   power cell. From startEdge to endEdge.
            TriEdge startEdge, endEdge, searchEdge;
            startEdge = N_Last_Exit.esym();
            endEdge = N.tri->omt.esym();
            //===================== Subroutine start =================================
            // maximum vertex degree.
            int max_vert_deg = OMT_domain->ct_in_vrts;
            if (max_vert_deg < 1000) max_vert_deg = 1000;
            // The maximum possible searching steps.
            int max_iter = OMT_domain->ct_hullsize;
            if (max_iter < 1000) max_iter = 1000;
            int deg = 0, iter = 0;
            bool bflag = true;
            // Start searching corner vertices from startEdge towards endEdge.
            // (see Page 2 for an example.)
            searchEdge = startEdge;
            //while ((searchEdge.tri != endEdge.tri)) {
            while ( !((searchEdge.org() == endEdge.org()) &&
                      (searchEdge.dest() == endEdge.dest())) ) {
              iter++; // count a searched boundary edge.
              if (iter >= max_iter) {
                bflag = false; break; // failed to find the cut edge.
              }
              // Found a corner vertex.
              Vertex *pt = &(ptlist[ptcount]);
              pt->init();
              pt->crd[0] = searchEdge.dest()->crd[0];
              pt->crd[1] = searchEdge.dest()->crd[1];
              // Calculate its weight
              double alpha = pt->crd[0];
              double beta  = pt->crd[1];
              double z_proj = 2.*a*alpha+2.*b*beta - (a*a+b*b) + wei;
              pt->crd[2] = alpha*alpha+beta*beta - z_proj; // its weight
              pt->tag = -2; // A background vertex.
              pt->idx = ptcount+1;
              ptcount++;
              // Go to the next boundary edge (around dest()).
              TriEdge workedge = searchEdge.enext();
              if ((workedge.org() == endEdge.org()) &&
                  (workedge.dest() == endEdge.dest())) {
                break; // the walk is finished.
              }
              searchEdge = workedge;
              while (true) {
                if (searchEdge.is_segment()) break;
                if ((searchEdge.esym()).tri->is_hulltri()) break;
                searchEdge = workedge.esym_enext(); // CW
                workedge = searchEdge;
                deg++;
                if (deg >= max_vert_deg) {
                  bflag = false; break;
                }
              }
              //if (searchEdge.is_segment()) {
              //  // hit a segment, we must stop.
              //  // It is possible that this segment is not the endEdge.
              //  break;
              //}
            }
            if (!bflag) {
              // failed to find the next boundary edge at S.org().
              delete [] ptlist;
              return 0;
            }
            //printf("Done!\n");
            //===================== Subroutine end =================================
          }
          // Clear it. We're walking from OUTSIDE to INSIDE already.
          N_Last_Exit.tri = NULL;
          
          // Add the cut line segment of this dual edge.
          Vertex *pt = &(ptlist[ptcount]);
          pt->init();
          pt->crd[0] = Cut_pt1.crd[0];
          pt->crd[1] = Cut_pt1.crd[1];
          pt->crd[2] = Cut_pt1.crd[2]; // its weight
          pt->tag = -1; // A cut vertex.
          pt->idx = ptcount+1;
          ptcount++;

          pt = &(ptlist[ptcount]);
          pt->init();
          pt->crd[0] = Cut_pt2.crd[0];
          pt->crd[1] = Cut_pt2.crd[1];
          pt->crd[2] = Cut_pt2.crd[2]; // its weight
          pt->tag = -1; // A cut vertex.
          pt->idx = ptcount+1;
          ptcount++;

          if (!clean_flag) {
            // Remember this (exit) boundary edge.
            N_Last_Exit = E.tri->omt;
          } else {
            N_Last_Exit.tri = NULL;
          }
        } // if (is_cut_domain)
      } else {
        // N's dual vertex lies inside the domain. (see Page 3, v2->v3)
        // We're walking from INSIDE to OUTSIDE.
        //assert(N_Last_Exit.tri == NULL);
        if (!(N_Last_Exit.tri == NULL)) {
          quit(QUIT_ON_BUG);
        }
        // Calculating a cut vertex between this dual edge (from N.tri->cct to E.tri->cct)
        //    and a boundary edge (or a segment) S of the background mesh.
        //    We must ensure that the dual edge and S intersect in S's interior.
        TriEdge S = E.tri->omt;
        Vertex In_pt, Out_pt;
        In_pt.init();
        Out_pt.init();
        In_pt.crd[0] = N.tri->cct[0]; // see Page 3, v2
        In_pt.crd[1] = N.tri->cct[1];
        Out_pt.crd[0] = E.tri->cct[0]; // see Page 3, v3
        Out_pt.crd[1] = E.tri->cct[1];

        if (!get_boundary_cut_dualedge(In_pt, Out_pt, S)) { // see Page 3, S = {b1, b2}.
          // Failed to get this boundary edge which cuts the dual edge.
          //assert(0); // should be a bug.
          delete [] ptlist;
          return 0;
        }

        // Calculate the cut vertex, see Page 3, cut1.
        E.tri->omt = S; // Update E.tri->on
        //printf("  Last Exit [%d,%d,%d]\n", S.org()->idx, S.dest()->idx, S.apex()->idx);
        
        Vertex *e1 = E.tri->omt.org();
        Vertex *e2 = E.tri->omt.dest();
        /*
        //assert((e1 != NULL) && (e2 != NULL));
        double X0 = E.tri->cct[0];
        double Y0 = E.tri->cct[1];
        double X1 = N.tri->cct[0];
        double Y1 = N.tri->cct[1];
        double X2 = e1->crd[0];
        double Y2 = e1->crd[1];
        double X3 = e2->crd[0];
        double Y3 = e2->crd[1];
        */
        // Sort the vertices to make sure that the cutting point is exactly
        //   the same when we calculate the adjacent cell.
        double X0,Y0,X1,Y1,X2,Y2,X3,Y3;
        if (E.apex()->idx < N.apex()->idx) {
          X0 = E.tri->cct[0];
          Y0 = E.tri->cct[1];
          X1 = N.tri->cct[0];
          Y1 = N.tri->cct[1];
        } else {
          X1 = E.tri->cct[0];
          Y1 = E.tri->cct[1];
          X0 = N.tri->cct[0];
          Y0 = N.tri->cct[1];
        }
        if (e1->idx < e2->idx) {
          X2 = e1->crd[0];
          Y2 = e1->crd[1];
          X3 = e2->crd[0];
          Y3 = e2->crd[1];
        } else {
          X3 = e1->crd[0];
          Y3 = e1->crd[1];
          X2 = e2->crd[0];
          Y2 = e2->crd[1];
        }
        double t1 = 0, t2 = 0;

        if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
          double c0 = X0 + t1 * (X1 - X0);
          double c1 = Y0 + t1 * (Y1 - Y0);
          // Calculate the weight
          double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
          double c2 = c0*c0+c1*c1 - z_proj;
          /*
          x1 = (int) ( Sxy * N.tri->cct[0] + Cx);
          y1 = (int) (-Sxy * N.tri->cct[1] + Cy);
          x2 = (int) ( Sxy * c0 + Cx);
          y2 = (int) (-Sxy * c1 + Cy);
          painter->drawLine(QPoint(x1, y1), QPoint(x2,y2));
          */
          Vertex *pt = &(ptlist[ptcount]);
          pt->init();
          pt->crd[0] = c0;
          pt->crd[1] = c1;
          pt->crd[2] = c2; // its weight
          pt->tag = -1; // A cut vertex.
          pt->idx = ptcount+1;
          ptcount++;
        } else {
          quit(QUIT_ON_BUG);
        }

        // Remember this (exit) boundary edge.
        N_Last_Exit = E.tri->omt;
        //printf("  2 Last Exit [%d,%d,%d]\n", N_Last_Exit.org()->idx,
        //       N_Last_Exit.dest()->idx, N_Last_Exit.apex()->idx);
      }
    } else if (E.tri->is_dual_on_bdry()) {
      if (N.tri->is_dual_in_exterior() || N.tri->is_dual_on_bdry()) {
        // We're walking from OUTSIDE to INSIDE.
        // (see Page 2, E.tri (p,-1,p2)'s cct is v3, N.tri (p,p1,-1)'s cct is v2.)
        //assert(N_Last_Exit.tri != NULL);
        if (!(N_Last_Exit.tri != NULL)) {
          quit(QUIT_ON_BUG);
        }
        // Since E's dual vertex lies exactly on boundary, then E.tri->on must contain this dual vertex.
        //   This is guaranteed by the function get_hulltri_orthocenter().
        //   No need to call get_boundary_cut_dualedge().
        //printf("\n DBG start: \n");
        //printf("  N_Last_Exit: [%d,%d,%d]\n", N_Last_Exit.org()->idx, N_Last_Exit.dest()->idx, N_Last_Exit.apex()->idx);
        //printf("  E.tri->on:   [%d,%d,%d]\n", E.tri->on.org()->idx, E.tri->on.dest()->idx, E.tri->on.apex()->idx);
        //printf("\n DBG end: \n");
        if ((E.tri->omt.tri != N_Last_Exit.tri) ||
            ((E.tri->omt.tri == N_Last_Exit.tri) && (E.tri->omt.ver != N_Last_Exit.ver))) {
          // They are on different boundary edges.
          // (see Page 2, the bdry edges containing v2 and v3 are different.)
          // There might be corner vertices to be the vertices of this cell.
          //printf("  !!! Searching background corners (2) DEBUG!!!\n");
          // [Important] The search is by rotating the triangles in the inside of the
          //   power cell. From startEdge to endEdge.
          TriEdge startEdge, endEdge, searchEdge;
          startEdge = N_Last_Exit.esym();
          endEdge = E.tri->omt.esym();
          //===================== Subroutine start =================================
          // Start searching corner vertices from startEdge towards endEdge.
          // (see Page 2 for an example.)
          // maximum vertex degree.
          int max_vert_deg = OMT_domain->ct_in_vrts;
          if (max_vert_deg < 1000) max_vert_deg = 1000;
          // The maximum possible searching steps.
          int max_iter = OMT_domain->ct_hullsize;
          if (max_iter < 1000) max_iter = 1000;
          int deg = 0, iter = 0;
          bool bflag = true;

          searchEdge = startEdge;
          //while ((searchEdge.tri != endEdge.tri)) {
          while ( !((searchEdge.org() == endEdge.org()) &&
                    (searchEdge.dest() == endEdge.dest())) ) {
            iter++; // count a searched boundary edge.
            if (iter >= max_iter) {
              bflag = false; break; // failed to find the cut edge.
            }
            // Found a corner vertex.
            Vertex *pt = &(ptlist[ptcount]);
            pt->init();
            pt->crd[0] = searchEdge.dest()->crd[0];
            pt->crd[1] = searchEdge.dest()->crd[1];
            // Calculate its weight
            double alpha = pt->crd[0];
            double beta  = pt->crd[1];
            double z_proj = 2.*a*alpha+2.*b*beta - (a*a+b*b) + wei;
            pt->crd[2] = alpha*alpha+beta*beta - z_proj; // its weight
            pt->tag = -2; // A background vertex.
            pt->idx = ptcount+1;
            ptcount++;
            // Go to the next boundary edge (around dest()).
            TriEdge workedge = searchEdge.enext();
            if ((workedge.org() == endEdge.org()) &&
                (workedge.dest() == endEdge.dest())) {
              break; // the walk is finished.
            }
            searchEdge = workedge;
            while (true) {
              if (searchEdge.is_segment()) break;
              if ((searchEdge.esym()).tri->is_hulltri()) break;
              searchEdge = workedge.esym_enext(); // CW
              workedge = searchEdge;
              deg++;
              if (deg >= max_vert_deg) {
                bflag = false; break;
              }
            }
            //if (searchEdge.is_segment()) {
            //  // hit a segment, we must stop.
            //  // It is possible that this segment is not the endEdge.
            //  break;
            //}
          }
          if (!bflag) {
            // failed to find the next boundary edge at S.org().
            delete [] ptlist;
            return 0;
          }
          //printf("Done!\n");
          //===================== Subroutine end =================================
        }
        // Clear it. We're walking from OUTSIDE to INSIDE already.
        N_Last_Exit.tri = NULL;
      } else {
        // N's dual vertex lies inside the domain.
        // We're walking to the direction from INSIDE to OUTSIDE.
        // (see Page 2, E.tri (p,p1,-1)'s cct is v2, N.tri (p,p3,p1)'s cct is v1.)
        // Remeber this boundary edge.
        //assert(N_Last_Exit.tri == NULL);
        if (!(N_Last_Exit.tri == NULL)) {
          quit(QUIT_ON_BUG);
        }
        N_Last_Exit = E.tri->omt; // (see Page 2, N_Last_Exit = startEdge.)
        //printf("  2 Last Exit [%d,%d,%d]\n", N_Last_Exit.org()->idx,
        //       N_Last_Exit.dest()->idx, N_Last_Exit.apex()->idx);
        /*
        // They are in the same subdomain.
        x1 = (int) ( Sxy * E.tri->cct[0] + Cx);
        y1 = (int) (-Sxy * E.tri->cct[1] + Cy);
        x2 = (int) ( Sxy * N.tri->cct[0] + Cx);
        y2 = (int) (-Sxy * N.tri->cct[1] + Cy);
        painter->drawLine(QPoint(x1, y1), QPoint(x2,y2));
        */
      }

      Vertex *pt = &(ptlist[ptcount]);
      pt->init();
      pt->crd[0] = E.tri->cct[0]; // see Page 2, v2
      pt->crd[1] = E.tri->cct[1];
      // Calculate its weight
      double alpha = pt->crd[0];
      double beta  = pt->crd[1];
      double z_proj = 2.*a*alpha+2.*b*beta - (a*a+b*b) + wei;
      pt->crd[2] = alpha*alpha+beta*beta - z_proj; // its weight
      //pt->crd[2] = E.tri->cct[2]; // its weight
      pt->tag = E.tri->idx; // A Voronoi vertex.
      pt->idx = ptcount+1;
      ptcount++;
    } else {
      // E's dual vertex lies inside the domain (or subdomain).
      // There must exist a (or part of) dual edge.
      if (!N.tri->is_dual_in_exterior() || N.tri->is_dual_on_bdry()) {
        // We're walking inside of the OMT_domain.
        //assert(N_Last_Exit.tri == NULL);
        if (!(N_Last_Exit.tri == NULL)) {
          quit(QUIT_ON_BUG);
        }
        /*
        // An interior voronoi edge.
        x1 = (int) ( Sxy * E.tri->cct[0] + Cx);
        y1 = (int) (-Sxy * E.tri->cct[1] + Cy);
        x2 = (int) ( Sxy * N.tri->cct[0] + Cx);
        y2 = (int) (-Sxy * N.tri->cct[1] + Cy);
        painter->drawLine(QPoint(x1, y1), QPoint(x2,y2));
        */
      } else {
        // The adjacent voronoi vertex is in exterior of the OMT.
        // We're walking from OUTSIDE to INSIDE.
        // (see Page 3, E,tri (p,p5,p6)'cct is v5, N.tri (p,p4,p5)'cct is v4,
        //    the walk is from v4 to v5.)
        //assert(N_Last_Exit.tri != NULL);
        if (!(N_Last_Exit.tri != NULL)) {
          quit(QUIT_ON_BUG);
        }
        // Calculating a cut vertex between this dual edge (from N.tri->cct to E.tri->cct)
        //    and a boundary edge (or a segment) S of the background mesh.
        //    We must ensure that the dual edge and S intersect in S's interior.
        TriEdge S = N.tri->omt;
        Vertex In_pt, Out_pt;
        In_pt.init();
        Out_pt.init();
        In_pt.crd[0] = E.tri->cct[0]; // see Page 3, v5
        In_pt.crd[1] = E.tri->cct[1];
        Out_pt.crd[0] = N.tri->cct[0]; // see Page 3, v4
        Out_pt.crd[1] = N.tri->cct[1];

        if (!get_boundary_cut_dualedge(In_pt, Out_pt, S)) {
          delete [] ptlist;
          return 0;
        }

        N.tri->omt = S; // Update N.tri->on

        //printf("\n DBG start: \n");
        //printf("  N_Last_Exit: [%d,%d,%d]\n", N_Last_Exit.org()->idx, N_Last_Exit.dest()->idx, N_Last_Exit.apex()->idx);
        //printf("  N.tri->on:   [%d,%d,%d]\n", N.tri->on.org()->idx, N.tri->on.dest()->idx, N.tri->on.apex()->idx);
        //printf("\n DBG end: \n");
        if ((N.tri->omt.tri != N_Last_Exit.tri) ||
            ((N.tri->omt.tri == N_Last_Exit.tri) && (N.tri->omt.ver != N_Last_Exit.ver))) {
          // They are on different boundary.
          // There might be corner vertices to be the vertices of this cell.
          //printf("  !!! Searching background corners (2) DEBUG!!!\n");
          // (see Page 2 for an example.)
          TriEdge startEdge, endEdge, searchEdge;
          startEdge = N_Last_Exit.esym();
          endEdge = N.tri->omt.esym();
          //===================== Subroutine start =================================
          // Start searching corner vertices from startEdge towards endEdge.
          // maximum vertex degree.
          int max_vert_deg = OMT_domain->ct_in_vrts;
          if (max_vert_deg < 1000) max_vert_deg = 1000;
          // The maximum possible searching steps.
          int max_iter = OMT_domain->ct_hullsize;
          if (max_iter < 1000) max_iter = 1000;
          int deg = 0, iter = 0;
          bool bflag = true;
          
          searchEdge = startEdge;
          //while (searchEdge.tri != endEdge.tri) {
          while ( !((searchEdge.org() == endEdge.org()) &&
                    (searchEdge.dest() == endEdge.dest())) ) {
            iter++; // count a searched boundary edge.
            if (iter >= max_iter) {
              bflag = false; break; // failed to find the cut edge.
            }
            // Found a corner vertex.
            Vertex *pt = &(ptlist[ptcount]);
            pt->init();
            pt->crd[0] = searchEdge.dest()->crd[0];
            pt->crd[1] = searchEdge.dest()->crd[1];
            // Calculate its weight
            double alpha = pt->crd[0];
            double beta  = pt->crd[1];
            double z_proj = 2.*a*alpha+2.*b*beta - (a*a+b*b) + wei;
            pt->crd[2] = alpha*alpha+beta*beta - z_proj; // its weight
            pt->tag = -2; // A background vertex.
            pt->idx = ptcount+1;
            ptcount++;
            // Go to the next boundary edge (around dest()).
            TriEdge workedge = searchEdge.enext();
            if ((workedge.org() == endEdge.org()) &&
                (workedge.dest() == endEdge.dest())) {
              break; // the walk is finished.
            }
            searchEdge = workedge;
            while (true) {
              if (searchEdge.is_segment()) break;
              if ((searchEdge.esym()).tri->is_hulltri()) break;
              searchEdge = workedge.esym_enext(); // CW
              workedge = searchEdge;
              deg++;
              if (deg >= max_vert_deg) {
                bflag = false; break;
              }
            }
            //if (searchEdge.is_segment()) {
            //  // to do ...
            //}
          }
          if (!bflag) {
            // failed to find the next boundary edge at S.org().
            delete [] ptlist;
          }
          //printf("Done!\n");
          //===================== Subroutine end =================================
        }

        // Calculate the cut2 vertex (see Page 3, cut2)
        Vertex *e1 = N.tri->omt.org();
        Vertex *e2 = N.tri->omt.dest();
        /*
        assert((e1 != NULL) && (e2 != NULL));
        double X0 = E.tri->cct[0];
        double Y0 = E.tri->cct[1];
        double X1 = N.tri->cct[0];
        double Y1 = N.tri->cct[1];
        double X2 = e1->crd[0];
        double Y2 = e1->crd[1];
        double X3 = e2->crd[0];
        double Y3 = e2->crd[1];
        */
        // Sort the vertices to make sure that the cutting point is exactly
        //   the same when we calculate the adjacent cell.
        double X0,Y0,X1,Y1,X2,Y2,X3,Y3;
        if (E.apex()->idx < N.apex()->idx) {
          X0 = E.tri->cct[0];
          Y0 = E.tri->cct[1];
          X1 = N.tri->cct[0];
          Y1 = N.tri->cct[1];
        } else {
          X1 = E.tri->cct[0];
          Y1 = E.tri->cct[1];
          X0 = N.tri->cct[0];
          Y0 = N.tri->cct[1];
        }
        if (e1->idx < e2->idx) {
          X2 = e1->crd[0];
          Y2 = e1->crd[1];
          X3 = e2->crd[0];
          Y3 = e2->crd[1];
        } else {
          X3 = e1->crd[0];
          Y3 = e1->crd[1];
          X2 = e2->crd[0];
          Y2 = e2->crd[1];
        }
        double t1 = 0, t2 = 0;

        if (line_line_intersection(X0, Y0, X1, Y1, X2, Y2, X3, Y3, &t1, &t2)) {
          double c0 = X0 + t1 * (X1 - X0);
          double c1 = Y0 + t1 * (Y1 - Y0);
          // Calculate weight
          double z_proj = 2.*a*c0+2.*b*c1 - (a*a+b*b) + wei;
          double c2 = c0*c0+c1*c1 - z_proj;
          /*
          x1 = (int) ( Sxy * E.tri->cct[0] + Cx);
          y1 = (int) (-Sxy * E.tri->cct[1] + Cy);
          x2 = (int) ( Sxy * c0 + Cx);
          y2 = (int) (-Sxy * c1 + Cy);
          painter->drawLine(QPoint(x1, y1), QPoint(x2,y2));
          */
          Vertex *pt = &(ptlist[ptcount]);
          pt->init();
          pt->crd[0] = c0;
          pt->crd[1] = c1;
          pt->crd[2] = c2; // its weight
          pt->tag = -1; // A cut vertex.
          pt->idx = ptcount+1;
          ptcount++;
        } else {
          quit(QUIT_ON_BUG);
        }
      }

      Vertex *pt = &(ptlist[ptcount]);
      pt->init();
      pt->crd[0] = E.tri->cct[0]; // see Page 3, v5
      pt->crd[1] = E.tri->cct[1];
      // Calculate its weight
      double alpha = pt->crd[0];
      double beta  = pt->crd[1];
      double z_proj = 2.*a*alpha+2.*b*beta - (a*a+b*b) + wei;
      pt->crd[2] = alpha*alpha+beta*beta - z_proj;
      //pt->crd[2] = E.tri->cct[1]; // its weight
      pt->tag = E.tri->idx; // A Voronoi vertex.
      pt->idx = ptcount+1;
      ptcount++;

      // Clear it, we're inside the OMT_domain again.
      N_Last_Exit.tri = NULL;
    }
  
    // Go to the next edge.
    E = E.eprev_esym(); // ccw
  } while (E.tri != mesh_vertex->adj.tri);

  *pptlist = ptlist;
  *ptnum = ptcount;

  return 1;
}

//==============================================================================

int Triangulation::get_constrained_power_cell(
  Vertex*  mesh_vertex,
  TriEdge& S_start,
  TriEdge& S_end,
  Vertex** pptlist, int* ptnum)
{
  // This function uses the triangulation itself as domain.
  if (OMT_domain != NULL) {
    quit(QUIT_ON_BUG);
  }

  return 1;
}

//==============================================================================

int Triangulation::get_mass_center(Vertex* ptlist, int ptnum, REAL mc[2])
{
  // Calculate cell mass center
  mc[0] = mc[1] = 0.;
  
  // The formula from Wikipedia.
  double Area = 0.; // signed area.
  for (int i = 0; i < ptnum; i++) {
    Vertex *v1 = &(ptlist[i]);
    Vertex *v2 = &(ptlist[(i+1)%ptnum]);
    double x0 = v1->crd[0];
    double y0 = v1->crd[1];
    double x1 = v2->crd[0];
    double y1 = v2->crd[1];
    Area += 0.5 * (x0 * y1 - x1 * y0);
  }

  double mcx = 0., mcy = 0.;
  for (int i = 0; i < ptnum; i++) {
    Vertex *v1 = &(ptlist[i]);
    Vertex *v2 = &(ptlist[(i+1)%ptnum]);
    double x0 = v1->crd[0];
    double y0 = v1->crd[1];
    double x1 = v2->crd[0];
    double y1 = v2->crd[1];
    mcx += (x0 + x1)*(x0 * y1 - x1 * y0);
    mcy += (y0 + y1)*(x0 * y1 - x1 * y0);
  }
  mcx /= (6 * Area);
  mcy /= (6 * Area);

  mc[0] = mcx;
  mc[1] = mcy;
  return 1;
}

//==============================================================================
/*
int Triangulation::get_voronoi_diagram_from_dual_edges()
{
  // The set of cut points and corner points will be generated and stored.
  Voronoi->tr_steiners = new arraypool(sizeof(Vertex), 10);
  // The set of Voronoi edges are stored as segments in Voronoi.
  Voronoi->tr_segs = new arraypool(sizeof(Segment), 10);

  arraypool *fqueue = new arraypool(sizeof(TriEdge), 8);
  TriEdge Evoro, tt[4];
  Vertex *edge[2], masspt;
  int loc, dir, i, j;

  // Get the list of mesh edges (for finding their dual edges).
  int ne = (3 * (tr_tris->objects - ct_hullsize) + ct_hullsize) / 2;
  TriEdge *mesh_edge_list = new TriEdge[ne];
  
  int idx = 0; // io_firstindex; see io.cpp, save_edges()
  for (i = 0; i < tr_tris->used_items; i++) {
    Triang* tri = (Triang *) tr_tris->get(i);
    if (tri->is_deleted() || tri->is_hulltri()) continue;
    tt[0].tri = tri;
    for (tt[0].ver = 0; tt[0].ver < 3; tt[0].ver++) {
      tt[1] = tt[0].esym(); // tt[1] might be a hull edge (apex is -1).
      if (tt[0].apex()->idx > tt[1].apex()->idx) { // must use '>'
        mesh_edge_list[idx++] = tt[0];
      }
    }
  } // i
  assert(idx == ne);

  // Get the dual (Voronoi) edge of each mesh edge one by one.
  Vertex voro_edge[2];
  Vertex *ptlist = NULL;
  int ptcount = 0;
  bool release_ptlist_flag;
  TriEdge E, N;
  
  for (i = 0; i < ne; i++) {
    
    Vertex **vrtlist = new Vertex*[ptcount];

    for (j = 0; j < ptcount; j++) {
      Evoro.tri = NULL;
      loc = Voronoi->locate_point(&(ptlist[j]), Evoro, 0); // encflag = 0
      
      if (loc == LOC_ON_VERT) {
        vrtlist[j] = Evoro.org();
      } else if ((loc == LOC_IN_OUTSIDE) ||
                 (loc == LOC_IN_TRI)) {
        // Add a Steiner points (a cut point or a corner point).
        vrtlist[j] = (Vertex *) Voronoi->tr_steiners->alloc();
        vrtlist[j]->init();
        vrtlist[j]->crd[0] = ptlist[j].crd[0];
        vrtlist[j]->crd[1] = ptlist[j].crd[1];
        vrtlist[j]->wei = ptlist[j].crd[2];
        vrtlist[j]->idx = Voronoi->ct_in_vrts + Voronoi->tr_steiners->objects;
        vrtlist[j]->typ = FREEVERTEX;
        //vrtlist[j]->tag = in_vrts[i].idx; // remember it.
        tt[0] = Evoro;
        int fflag = FLIP_13;
        Voronoi->flip(tt, &(vrtlist[j]), fflag, fqueue);
        Voronoi->lawson_flip(vrtlist[j], 1, fqueue);
      } else if (loc == LOC_ON_EDGE) {
        // Add a Steiner points (a cut point or a corner point).
        vrtlist[j] = (Vertex *) Voronoi->tr_steiners->alloc();
        vrtlist[j]->init();
        vrtlist[j]->crd[0] = ptlist[j].crd[0];
        vrtlist[j]->crd[1] = ptlist[j].crd[1];
        vrtlist[j]->wei = ptlist[j].crd[2];
        vrtlist[j]->idx = Voronoi->ct_in_vrts + Voronoi->tr_steiners->objects;
        vrtlist[j]->typ = FREEVERTEX;
        //vrtlist[j]->tag = in_vrts[i].idx; // remember it.
        tt[0] = Evoro;
        int fflag = FLIP_24;
        Voronoi->flip(tt, &(vrtlist[j]), fflag, fqueue);
        Voronoi->lawson_flip(vrtlist[j], 1, fqueue);
      } else {
        assert(0); // report a bug.
      }
    } // j

    int edge_num = ptcount / 2;
    
    for (j = 0; j < edge_num; j++) {
      edge[0] = vrtlist[j*2];
      edge[1] = vrtlist[j*2+1];
      dir = Voronoi->recover_edge(edge[0], edge[1], Evoro, fqueue);
      if (dir != INTERSECT_SHARE_EDGE) {
        assert(0); // report a bug.
      }
      if (!Evoro.is_segment()) {
        Voronoi->insert_segment(Evoro, -1, 0, NULL);
      }
    }

    delete [] vrtlist;

    if (release_ptlist_flag) {
      delete [] ptlist;
    }
  } // i
  
  delete [] mesh_edge_list;
  delete fqueue;
  
  return 1;
}
*/

//==============================================================================
// Get the Voronoi diagram by computing the Voronoi cells.

int Triangulation::get_voronoi_diagram_from_powrcells()
{
  // The set of cut points and corner points will be generated and stored.
  Voronoi->tr_steiners = new arraypool(sizeof(Vertex), 10);
  // The set of Voronoi edges are stored as segments in Voronoi.
  Voronoi->tr_segs = new arraypool(sizeof(Segment), 10);

  arraypool *fqueue = new arraypool(sizeof(TriEdge), 8);
  TriEdge Evoro, tt[4];
  Vertex *edge[2], masspt;
  int loc, dir, i, j;
  
  // Count the number of failures.
  int failed_count = 0;

  // Get the list of mesh vertices (for find their power cells).
  int nv = ct_in_vrts;
  if (tr_steiners != NULL) {
    nv += tr_steiners->objects;
  }
  Vertex **mesh_vertex_list = new Vertex*[nv];
  
  for (i = 0; i < ct_in_vrts; i++) {
    mesh_vertex_list[i] = &(in_vrts[i]);
  }
  if (tr_steiners) {
    for (i = 0; i < tr_steiners->used_items; i++) {
      Vertex *mesh_vertex = (Vertex *) tr_steiners->get(i);
      if (mesh_vertex->is_deleted()) continue;
      mesh_vertex_list[i+ct_in_vrts] = mesh_vertex;
    }
  }

  // Every mesh vertex has a unique index. All triangles in Voronoi will have
  //   a marker which belong to its Voronoi cell.
  //   The index of vertex is the same as save_nodes().
  int vertex_idx = io_firstindex;

  for (i = 0; i < nv; i++) {
    Vertex *mesh_vertex = mesh_vertex_list[i];

    if (mesh_vertex->typ == UNUSEDVERTEX) {
      if (io_keep_unused) { // with -IJ
        // Index it even it is an UNUSEDVERTEX vertex.
        mesh_vertex->idx = vertex_idx;
        vertex_idx++;
      }
      continue;
    }
    
    mesh_vertex->idx = vertex_idx;
    vertex_idx++;
 
    Vertex *ptlist = NULL;
    int ptnum = 0;
    bool success = true;

    if (get_powercell(mesh_vertex, false, &ptlist, &ptnum)) {
      // Step (1), Get the vertices of this power cell.
      Vertex **vrtlist = new Vertex*[ptnum];
      masspt.crd[0] = masspt.crd[1] = 0.;
      for (j = 0; j < ptnum; j++) {
        Evoro.tri = NULL;
        loc = Voronoi->locate_point(&(ptlist[j]), Evoro, 0); // encflag = 0
        
        if (loc == LOC_ON_VERT) {
          vrtlist[j] = Evoro.org();
        } else if ((loc == LOC_IN_OUTSIDE) ||
                   (loc == LOC_IN_TRI)) {
          // Add a Steiner points (a cut point or a corner point).
          vrtlist[j] = (Vertex *) Voronoi->tr_steiners->alloc();
          vrtlist[j]->init();
          vrtlist[j]->crd[0] = ptlist[j].crd[0];
          vrtlist[j]->crd[1] = ptlist[j].crd[1];
          vrtlist[j]->wei = ptlist[j].crd[2];
          vrtlist[j]->idx = Voronoi->ct_in_vrts + Voronoi->tr_steiners->objects;
          vrtlist[j]->typ = FREEVERTEX;
          //vrtlist[j]->tag = in_vrts[i].idx; // remember it.
          tt[0] = Evoro;
          int fflag = FLIP_13;
          Voronoi->flip(tt, &(vrtlist[j]), fflag, fqueue);
          Voronoi->lawson_flip(vrtlist[j], 1, fqueue);
        } else if (loc == LOC_ON_EDGE) {
          // Add a Steiner points (a cut point or a corner point).
          vrtlist[j] = (Vertex *) Voronoi->tr_steiners->alloc();
          vrtlist[j]->init();
          vrtlist[j]->crd[0] = ptlist[j].crd[0];
          vrtlist[j]->crd[1] = ptlist[j].crd[1];
          vrtlist[j]->wei = ptlist[j].crd[2];
          vrtlist[j]->idx = Voronoi->ct_in_vrts + Voronoi->tr_steiners->objects;
          vrtlist[j]->typ = FREEVERTEX;
          //vrtlist[j]->tag = in_vrts[i].idx; // remember it.
          tt[0] = Evoro;
          int fflag = FLIP_24;
          Voronoi->flip(tt, &(vrtlist[j]), fflag, fqueue);
          Voronoi->lawson_flip(vrtlist[j], 1, fqueue);
        } else {
          //assert(0); // report a bug.
          success = false;
          break;
        }
        masspt.crd[0] += vrtlist[j]->crd[0];
        masspt.crd[1] += vrtlist[j]->crd[1];
      }
      masspt.crd[0] /= ptnum;
      masspt.crd[1] /= ptnum;

      if (!success) {
        delete [] vrtlist;
        failed_count++;
        continue;
      }

      // Step (2), Insert the Voronoi edges (as line segments) in Voronoi.
      for (j = 0; j < ptnum; j++) {
        edge[0] = vrtlist[j];
        edge[1] = vrtlist[(j+1)%ptnum];
        if (edge[0] != edge[1]) {
          dir = Voronoi->recover_edge(edge[0], edge[1], Evoro, fqueue);
          if (dir != INTERSECT_SHARE_EDGE) {
            //assert(0); // report a bug.
            success = false;
            break;
          }
          if (!Evoro.is_segment()) {
            Voronoi->insert_segment(Evoro, -1, 0, NULL);
          }
        }
      }

      if (!success) {
        delete [] vrtlist;
        failed_count++;
        continue;
      }

      // Step (3), mark all triangles in this Voronoi cell.
      // use masspt to locate a triangle.
      loc = Voronoi->locate_point(&masspt, Evoro, 0); // encflag = 0

      if (loc != LOC_IN_TRI) {
        // assert(0); // report a bug.
        success = false;
        delete [] vrtlist;
        failed_count++;
        continue;
      }

      // Let the mesh vertex remember this omt triangle.
      mesh_vertex->on_omt = Evoro;
      
      // Mark all triangles in this cell (re-use fqueue).
      Evoro.tri->set_infect();
      * (TriEdge *) fqueue->alloc() = Evoro;
      for (j = 0; j < fqueue->objects; j++) {
        Evoro = * (TriEdge *) fqueue->get(j);
        for (Evoro.ver = 0; Evoro.ver < 3; Evoro.ver++) {
          if (!Evoro.is_segment()) {
            tt[0] = Evoro.esym();
            if (!tt[0].tri->is_infected()) {
              tt[0].tri->set_infect();
              * (TriEdge *) fqueue->alloc() = tt[0];
            }
          }
        }
      }
      for (j = 0; j < fqueue->objects; j++) {
        Evoro = * (TriEdge *) fqueue->get(j);
        Evoro.tri->tag = mesh_vertex->idx;
        Evoro.tri->clear_infect();
      }
      fqueue->clean();

      delete [] vrtlist;
    } else {
      // get_powercell(...) failed.
      failed_count++;
      continue;
    }
    delete [] ptlist;
  } // i

  delete fqueue;
  delete [] mesh_vertex_list;

  if (op_db_verbose) {
    if (failed_count > 0) {
      printf("  !! There are %d failures. The Voronoi diagram is not complete.", failed_count);
    }
  }

  return failed_count == 0;
}

int Triangulation::construct_voronoi_diagram(int option)
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

  remove_skinny_hulltris();

  // Calculate Voronoi vertices for triangles.
  int i, idx;
  idx = 0; // io_firstindex; // index the triangles (0-based).
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

  if (Voronoi != NULL) {
    delete Voronoi;
  }

  Voronoi = new Triangulation();
  tr_free_voronoi = true;

  sprintf(Voronoi->io_outfilename, "%s.voro", io_outfilename);

  // The set of vertices of the Voronoi diagram is the set of circumcenters
  //   (of triangles) and bisectors (of boundary edges).
  Voronoi->ct_in_vrts = tr_tris->objects; // include hull triangles.
  Voronoi->in_vrts = new Vertex[Voronoi->ct_in_vrts];

  for (i = 0; i < tr_tris->used_items; i++) {
    Triang* tri = (Triang *) tr_tris->get(i);
    if (tri->is_deleted()) continue;
    // A Voronoi vertex.
    Vertex *vrt = &(Voronoi->in_vrts[tri->idx]);
    vrt->init();
    vrt->idx = tri->idx;
    vrt->crd[0] = tri->cct[0];
    vrt->crd[1] = tri->cct[1];
    vrt->crd[2] = 0.;
    // tri->cct[2] is the weight of this circumcenter (might be negative).
    // see function get_orthocenter().
    vrt->wei = tri->cct[2];
    if (tri->is_dual_in_exterior()) {
      // Do not insert this vertex into the triangulation of Voronoi.
      vrt->typ = UNUSEDVERTEX;
      Voronoi->ct_unused_vrts++;
    } else {
      // It will be a vertex in the triangulation of Voronoi.
      vrt->typ = FREEVERTEX; // flip13() and flip24().
    }
  }

  Voronoi->io_xmin = OMT_domain->io_xmin;
  Voronoi->io_xmax = OMT_domain->io_xmax;
  Voronoi->io_ymin = OMT_domain->io_ymin;
  Voronoi->io_ymax = OMT_domain->io_ymax;

  Voronoi->io_diagonal = OMT_domain->io_diagonal;
  Voronoi->io_diagonal2 = OMT_domain->io_diagonal2;

  // For constructing the DT. (ignore the weights of Voronoi vertices).
  Voronoi->io_with_wei = 0; // needed by incremental_delaunay().
  Voronoi->op_metric = METRIC_Euclidean_no_weight;

  // Construct the DT of Voronoi vertices (skip exterior vertices).
  int arysize = Voronoi->ct_in_vrts - Voronoi->ct_unused_vrts;
  Vertex **vrtarray = new Vertex*[arysize];

  idx = 0;
  for (i = 0; i < Voronoi->ct_in_vrts; i++) {
    if (Voronoi->in_vrts[i].typ != UNUSEDVERTEX) {
      vrtarray[idx++] = &(Voronoi->in_vrts[i]);
    }
  }
  assert(idx == arysize);

  Voronoi->insert_vertices(vrtarray, arysize, 1); // hullflag = 1

  delete [] vrtarray;

  if (option == 0) { // default
    get_voronoi_diagram_from_powrcells();
  } else {
    // to do...
    //get_voronoi_diagram_from_dual_edges();
  }

  if (clean_omt_domain) {
    OMT_domain = NULL;
  }
  return 1;
}
