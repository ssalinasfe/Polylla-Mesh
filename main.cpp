#include <iostream>
#include <stdio.h>   
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <vector> 
#include <set>
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include <list>
#include <string>
#include <math.h>  



#include "delaunay.h"
#include "io.h"
#include "consts.h"
#include "triang.h"
#include "polygon.h"
#include "mesh.h"
#include "BET_elimitation.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)

//make &&  python3 rp.py 100 0 0 && ./DelaunayPolyGenerator autodata.node && geomview output/.off
//make &&  python3 rp.py 506 0 0 && ./DelaunayPolyGenerator -p -q unicorn.poly && geomview output/.off

//detri funca
// python3 2x2quatercircle.py 63000 && make &&  ./DelaunayPolyGenerator -p input/RP2x2quartercircle_63000.poly 
// make && ./DelaunayPolyGenerator input/2x2_10.node && geomview output/2x2_10.off 
//algp da vacios
// make &&  ./DelaunayPolyGenerator -z input/unisquare2x2_41857.node && geomview output/unisquare2x2_41857.off

int main(int argc, char* argv[]){

	char* ppath;
	int print_triangles = 0;

	// int nparam = 3;
    //char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("test.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506randompoints.node")};
	//char* params[] = {const_cast<char*> ("./detri2"), const_cast<char*> ("-z"), const_cast<char*> ("506equilateral.node")};
    //char* ppath = const_cast<char*> ("test");

	//number elements
	int tnumber, pnumber, i,j, num_border;
	//arrays inicilization
	double *r;
	int *triangles;
	int *adj;
    int *seed;
    int *max;
	int *mesh;
	int *trivertex;
	
	auto tb_delaunay = std::chrono::high_resolution_clock::now();
	generate_delaunay_from_random_points(argc, argv, pnumber,tnumber);
	auto te_delaunay = std::chrono::high_resolution_clock::now();

	r = (double *)malloc(2*pnumber*sizeof(double)); // cambiar por pnumber
    triangles = (int *)malloc(3*tnumber*sizeof(int));
	adj = (int *)malloc(3*tnumber*sizeof(int));
	seed = (int *)malloc(tnumber*sizeof(int));
    max = (int *)malloc(tnumber*sizeof(int));
	mesh = (int *)malloc(3*tnumber*sizeof(int));
	trivertex = (int *)malloc(pnumber*sizeof(int));

	int *border = (int *)malloc(2*tnumber*sizeof(int));

	copy_delaunay_arrays(tnumber, r, triangles, adj);

	//Asociate each vertex to an adjacent  triangle
	for(i = 0; i < pnumber; i++){
		for (j = 0; j < tnumber; j++)
		{
			if(i == triangles[3*j + 0] ||  i == triangles[3*j + 1] || i == triangles[3*j + 2]){
				trivertex[i] = j;
				break;
			}
		}
		//std::cout<<"trivertex["<<i<<"] "<<trivertex[i]<<std::endl;
	}

	num_border = get_border_points(pnumber,tnumber, border,triangles, adj, r);
	
	//stats
	int i_mesh = 0;	
	int num_BE = 0;
	int est_total_be = 0;
	int est_min_triangles_be = 2147483641;
	int est_max_triangles_be = 0;
	int est_poly_with_be = 0;
	double est_ratio_be = 0;
	int tcost_be = 0;
	

	//initialize array
	for(i = 0; i < tnumber; i++){
		seed[i] = FALSE;
	}
	
	auto t1 = std::chrono::high_resolution_clock::now();
	auto tb_label =std::chrono::high_resolution_clock::now();
	/* Etapa 1: Encontrar aristas máximas. */
	debug_msg("Etapa 1: Encontrar aristas máximas. \n");
	for(i = 0; i < tnumber; i++)
	{
		//max[i] marca la arista entre los puntos P0-P1 como 0, P1-P2 como 1 y P2-P0 como 2 del arreglo triangles
		max[i] = max_edge_index(i, r, triangles); 
	}	
	
	//Marcar triangulos semilla
	for(i = 0; i < tnumber; i++){
		for(j = 0; j < 3; j++){
			if(adj[3*i +j] != -1 && is_max_max(i, adj[3*i + j], triangles, max) == TRUE ){
				if(adj[3*i + j] < i){
					seed[i] = TRUE;
					break;
				}
			}
			//(j + 1)%3 es debido a que la arista P0-P1 de triangle es 1 del arreglo de adjacencia, y así con el resto
			if(adj[3*i +j] == -1  && max[i] == (j +1)%3){ 
				seed[i] = TRUE;
				break;
			}
		}
	}

	/* Etapa 2: Desconectar arcos asociados a aristas nomáx-nomáx. */
	debug_msg("Etapa 2: Desconectar arcos asociados a aristas nomáx-nomáx. \n");
	
	int num_terminal_edges =0;
	int num_terminal_border_edges=0;
	int num_frontier_edges=0;
	int num_frontier_border_edges=0;
	int num_interior_edges=0;
		
	for(i = 0; i < tnumber; i++)
	{
		for(j = 0; j < 3; j++)
		{
			
			if(adj[3*i + j] < 0){
				if ((j + 1)%3 == max[i])
					num_terminal_border_edges++;
				else
					num_frontier_border_edges++;
			}else
			{
				//If has frontieredge
				if(is_nomax_nomax(i, adj[3*i + j], triangles, max))
					num_frontier_edges++;
				//If has terminal_edge
				else if(is_max_max(i, adj[3*i + j], triangles, max))
					num_terminal_edges++;
				else //if is interioredge
					num_interior_edges++;
				//if(is_max_nomax(i, adj[3*i + j], triangles, max))
			}
			

			//Marcación real
			if(adj[3*i +j] < 0 || is_nomax_nomax(i, adj[3*i + j], triangles, max))
			{
				adj[3*i + j] = NO_ADJ;
			}
			
		}
	}
	auto te_label =std::chrono::high_resolution_clock::now();
	
	free(max);
	
	
    int length_poly;
	std::list <int> seed_bet;  
	std::vector <int> seed_bet_mark(3*tnumber, 0);
	
	debug_msg("Etapa 5: Generar poligonos\n");
	int poly[1000];	
	//std::bitset<tnumber> tuhermana;
	auto tb_travel = std::chrono::high_resolution_clock::now();
	for(i = 0; i < tnumber; i++)
	{
		if(seed[i] == TRUE){			

			length_poly = generate_polygon(i, poly, triangles, adj, r);
			num_BE = count_BarrierEdges(poly, length_poly);
	
					
			if(num_BE>0){
				//i_mesh = save_to_mesh(mesh, poly, i_mesh, length_poly);
				//printf("%d %d\n", num_BE, length_poly);
				est_total_be += num_BE;
				est_min_triangles_be = (num_BE < est_min_triangles_be) ? num_BE : est_min_triangles_be;
				est_max_triangles_be = (num_BE > est_max_triangles_be) ? num_BE : est_max_triangles_be;
				est_poly_with_be++;
				est_ratio_be += (float)num_BE/length_poly;
			}
			
		
			debug_msg("Poly: "); debug_block(print_poly(poly, length_poly););
			if( num_BE > 0){
				seed[i] = FALSE;
				//printf("Se dectecto %d BE\n", num_BE);
				auto tb_be = std::chrono::high_resolution_clock::now();
				i_mesh = Remove_BE3(1,poly, length_poly, num_BE, triangles, adj, r, tnumber, mesh, i_mesh, trivertex, seed_bet, seed_bet_mark);
				//i_mesh = Remove_BE2(1,poly, length_poly, num_BE, triangles, adj, r, tnumber, mesh, i_mesh, trivertex, seed_bet);
				//i_mesh = Remove_BE(1,poly, length_poly, num_BE, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
				auto te_be = std::chrono::high_resolution_clock::now();
				tcost_be += std::chrono::duration_cast<std::chrono::milliseconds>( te_be - tb_be ).count();
				//i_mesh = save_to_mesh(mesh, poly, i_mesh, length_poly, r);	
			}else{
				debug_msg("Guardando poly\n");
				i_mesh = save_to_mesh(mesh, poly, i_mesh, length_poly, r);	
				
			}
		}
	}
	auto te_travel = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();

	int num_region = count_regions(mesh,i_mesh);

	std::string name(argv[argc-1]);
	//name.erase(0,6);
	name.erase(name.end()-5,name.end());

	write_geomview(name, r, triangles, pnumber, tnumber, i_mesh, mesh, seed, num_region, print_triangles);
	//write_alejandro(name, r, triangles, pnumber, tnumber, i_mesh, mesh, num_region);
	//write_alejandro_quater_circle(name, r, triangles, pnumber, tnumber, i_mesh, mesh, num_region);
	//write_alejandro_custom(name, r, triangles, pnumber, tnumber, i_mesh, mesh, num_region, border, num_border);
	//write_VEM(name, r, triangles, pnumber, tnumber, i_mesh, mesh, seed, num_region, print_triangles);
	//write_VEM_triangles(name, r, triangles, adj, pnumber, tnumber, i_mesh, mesh, seed, num_region, seed_bet);
	//write_GID(name, r, triangles, adj, pnumber, tnumber);
	//write_triangulation(name, r, triangles, adj, pnumber, tnumber);
	int t_delaunay = std::chrono::duration_cast<std::chrono::milliseconds>(te_delaunay - tb_delaunay).count();
	int t_label = std::chrono::duration_cast<std::chrono::milliseconds>(te_label - tb_label).count();
	int t_total = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1 ).count();
	int t_travel_and_opt = std::chrono::duration_cast<std::chrono::milliseconds>(te_travel - tb_travel).count();
	int t_travel = std::chrono::duration_cast<std::chrono::milliseconds>(te_travel - tb_travel).count() - tcost_be;
	
	write_metrics(name,r, triangles, pnumber, tnumber,i_mesh,  mesh,  num_region,  num_border,  num_terminal_edges,  num_terminal_border_edges,  num_frontier_edges,  num_frontier_border_edges,  num_interior_edges,  t_delaunay,  t_label,  t_total,  t_travel_and_opt,  t_travel, tcost_be, num_BE,  est_total_be,  est_min_triangles_be,  est_max_triangles_be,  est_poly_with_be, est_ratio_be);
	
	free(trivertex);
	free(r);
	free(triangles);
	free(adj);
	free(seed);
	free(mesh);    
	free(border);
	return EXIT_SUCCESS;
}
    

