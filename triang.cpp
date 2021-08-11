/* Funciones para manejar triangulación. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "triang.h"
#include "consts.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)

int Equality(double a, double b, double epsilon)
{
  return fabs(a - b) < epsilon;
}

int GreaterEqualthan(double a, double b, double epsilon){
	return Equality(a,b,epsilon) || a > b;
}

/* dist
 * 
 * Retorna la distancia que hay entre el punto (x0,y0)
 * y (x1,y1).
 * */

double dist(double x0, double y0, double x1, double y1)
{
	return sqrt(pow(x0 - x1, 2.0) + pow(y0 - y1, 2.0));
}



/* max_edge_index
 * 
 * Retorna el índice k de la arista máxima de un triángulo i, 
 * descrito por los puntos p0p1p2. Será 0 si p0p1 es máxima.
 * Será 1 si p1p2 lo es. Será 2 si p2p0 lo es.
 * */

int max_edge_index(int i, double *r, int *p)
{
	double l0;
	double l1;
	double l2;
	
	int p0;
	int p1;
	int p2;
	
	p0 = p[3*i + 0];
	p1 = p[3*i + 1];
	p2 = p[3*i + 2];
	
	l0 = dist(r[2*p0 + 0], r[2*p0 + 1], r[2*p1 + 0], r[2*p1 + 1]);
	l1 = dist(r[2*p1 + 0], r[2*p1 + 1], r[2*p2 + 0], r[2*p2 + 1]);
	l2 = dist(r[2*p2 + 0], r[2*p2 + 1], r[2*p0 + 0], r[2*p0 + 1]);
	//se aumento la precisión para las mallas de 25k en 2x2
	double epsion = 0.00000001f;

	//if((l0 >= l1 && l1 >= l2) || (l0 >= l2 && l2 >= l1))
	if( (GreaterEqualthan(l0,l1,epsion) && GreaterEqualthan(l1,l2,epsion)) || ( GreaterEqualthan(l0,l2,epsion) && GreaterEqualthan(l2,l1,epsion)))
	{
		return 0;
	}
	//else if((l1 >= l0 && l0 >= l2) || (l1 >= l2 && l2 >= l0))
	else if((GreaterEqualthan(l1,l0,epsion) && GreaterEqualthan(l0,l2,epsion)) || ( GreaterEqualthan(l1,l2,epsion) && GreaterEqualthan(l2,l0,epsion)))
	{
		return 1;
	}
	else
	{
		return 2;
	}



	//if((l2 >= l0 && l0 >= l1) || (l2 >= l1 && l1 >= l0))
	/*
	if( (GreaterEqualthan(l2,l0,epsion) && GreaterEqualthan(l0,l1,epsion)) || ( GreaterEqualthan(l2,l1,epsion) && GreaterEqualthan(l1,l0,epsion)))
	{
		return 2;
	}
	//else if((l0 >= l1 && l1 >= l2) || (l0 >= l2 && l2 >= l1))
	else if((GreaterEqualthan(l0,l1,epsion) && GreaterEqualthan(l1,l2,epsion)) || ( GreaterEqualthan(l0,l2,epsion) && GreaterEqualthan(l2,l1,epsion)))
	{
		return 0;
	}
	else
	{
		return 1;
	}
	*/

	/*	
	if((l1 >= l0 && l0 >= l2) || (l1 >= l2 && l2 >= l0))
	if( (GreaterEqualthan(l1,l0,epsion) && GreaterEqualthan(l0,l2,epsion)) || ( GreaterEqualthan(l1,l2,epsion) && GreaterEqualthan(l2,l0,epsion)))
	{
		return 1;
	}
	//else if((l0 >= l1 && l1 >= l2) || (l0 >= l2 && l2 >= l1))
	else if((GreaterEqualthan(l0,l1,epsion) && GreaterEqualthan(l1,l2,epsion)) || ( GreaterEqualthan(l0,l2,epsion) && GreaterEqualthan(l2,l1,epsion)))
	{
		return 0;
	}
	else
	{
		return 2;
	}
	*/
}



/* same_edge
 * 
 * Indica para las aristas {u,v} y {w,x} si son iguales o no.
 * */
 
int same_edge(int u, int v, int w, int x)
{
	return (u == w && v == x) || (u == x && v == w);
}



/* get_edge_index
 * 
 * Entrega el índice de la arista {u,v} respecto del triángulo i.
 * */

static int get_edge_index(int u, int v, int i, int *p)
{
	int p0;
	int p1;
	int p2;
	
	p0 = p[3*i + 0];
	p1 = p[3*i + 1];
	p2 = p[3*i + 2];
	
	if(same_edge(u, v, p0, p1))
	{
		return 0;
	}
	else if(same_edge(u, v, p1, p2))
	{
		return 1;
	}
	else if(same_edge(u, v, p2, p0))
	{
		return 2;
	}
	else
	{
		fprintf(stderr, "%s:%d:%s() ** ERROR ** get_edge_index: Arista {%d,%d} no pertenece al triángulo %d.\n", __FILE__,  __LINE__, __func__, u, v, i);
		exit(EXIT_FAILURE);
	}
}


/* is_max_nomax
 * 
 * Indica si la arista compartida entre los triángulos i y j
 * es internal-edge.
 * */

int is_max_nomax(int i, int j, int *p, int *max)
{
	int p0i;
	int p1i;
	int p2i;
	int p0j;
	int p1j;
	int p2j;
	
	p0i = p[3*i + 0];
	p1i = p[3*i + 1];
	p2i = p[3*i + 2];
	
	p0j = p[3*j + 0];
	p1j = p[3*j + 1];
	p2j = p[3*j + 2];
	
	int ij;
	int ii;
	
	if(same_edge(p0i, p1i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 2;
	}
	else
	{
		fprintf(stderr, "** ERROR ** is_nomax_nomax: Problema insperado para triángulos %d y %d.\n", i, j);
		exit(EXIT_FAILURE);
	}
	
	return ((ij == max[j]) && (ii != max[i]) || (ij != max[j]) && (ii == max[i]));
}


/* is_nomax_nomax
 * 
 * Indica si la arista compartida entre los triángulos i y j
 * es nomáx-nomáx.
 * */

int is_nomax_nomax(int i, int j, int *p, int *max)
{
	int p0i;
	int p1i;
	int p2i;
	int p0j;
	int p1j;
	int p2j;
	
	p0i = p[3*i + 0];
	p1i = p[3*i + 1];
	p2i = p[3*i + 2];
	
	p0j = p[3*j + 0];
	p1j = p[3*j + 1];
	p2j = p[3*j + 2];
	
	int ij;
	int ii;
	
	if(same_edge(p0i, p1i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 2;
	}
	else
	{
		fprintf(stderr, "** ERROR ** is_nomax_nomax: Problema insperado para triángulos %d y %d.\n", i, j);
		exit(EXIT_FAILURE);
	}
	
	return (ij != max[j]) && (ii != max[i]);
}



/* is_max_max
 * 
 * Indica si la arista compartida entre los triángulos i y j
 * es máx-máx.
 * */

int is_max_max(int i, int j, int *p, int *max)
{
	int p0i;
	int p1i;
	int p2i;
	
	int p0j;
	int p1j;
	int p2j;
	
	p0i = p[3*i + 0];
	p1i = p[3*i + 1];
	p2i = p[3*i + 2];
	
	p0j = p[3*j + 0];
	p1j = p[3*j + 1];
	p2j = p[3*j + 2];
	
	int ij;
	int ii;
	
	if(same_edge(p0i, p1i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p0j, p1j))
	{
		ij = get_edge_index(p0j, p1j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p1j, p2j))
	{
		ij = get_edge_index(p1j, p2j, j, p);
		ii = 2;
	}
	else if(same_edge(p0i, p1i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 0;
	}
	else if(same_edge(p1i, p2i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 1;
	}
	else if(same_edge(p2i, p0i, p2j, p0j))
	{
		ij = get_edge_index(p2j, p0j, j, p);
		ii = 2;
	}
	else
	{
		fprintf(stderr, "** ERROR ** is_max_max: Problema insperado para triángulos %d y %d.\n", i, j);
		exit(EXIT_FAILURE);
	}
	
	return (ij == max[j]) && (ii == max[i]);
}



/* edge_belongs_to
 * 
 * Indica si arista {k,l} pertenece al triángulo i.
 * */

static int edge_belongs_to(int k, int l, int i, int *p)
{
	return same_edge(k, l, p[3*i + 0], p[3*i + 1])
					|| same_edge(k, l, p[3*i + 1], p[3*i + 2])
					|| same_edge(k, l, p[3*i + 2], p[3*i + 0]);
}

/* Given one triangle i, return the edge index that containts u and v*/
int get_shared_edge(int i, int u, int v, int *p){
	int j, ind1,ind2;
	for(j = 0; j < 3; j++){
		ind1 = 3*i + j;
		ind2 = 3*i + (j+1)%3;
		//debug_print("%d %d %d %d %d\n", ind1, ind2, ind3, (p[ind1] == u || p[ind2] == u), (p[ind1] == v || p[ind2] == v));
		if( (p[ind1] == u || p[ind2] == u) && (p[ind1] == v || p[ind2] == v))
			return (j+2)%3;
	}
	fprintf(stderr, "ERROR get_edge: No se encontro el edge %d - %d del triangulo %d", u,v,i);
	exit(0);
	return EXIT_FAILURE;

}


/* get_adjacent_triangle
 * 
 * Retorna el identificador del triángulo que es adyacente al
 * triángulo i, mediante la arista {k,l}.
 * 
 * Si no hay triángulo, retorna NO_ADY (aún si es porque {k,l}
 * es de borde de triangulación).
 * */

int get_adjacent_triangle(int i, int k, int l, int *p, int *adj){
	return adj[3*i + get_shared_edge(i, k, l, p)];
}




/*Indica si un triangulo contiene al punto endpoint*/
int is_continuous(int i, int endpoint, int *p ){
	int p0, p1, p2;
	if (i != -1){
		p0 = p[3*i + 0];
		p1 = p[3*i + 1];
		p2 = p[3*i + 2];
				
		if(endpoint == p0){
			return  0; /* indica que está en p0*/
		}else if (endpoint == p1){
			return  1;  /* indica que está en p1*/
		}else if(endpoint == p2){
			return 2;  /* indica que está en p2*/
		}
	}
	return -1;
}


// advance i triangles arround vertex endpoint
int advance_i_adjacents_triangles_share_endpoint(int adv, int t, int origen, int endpoint, int *p, int *adj){
	int aux;
	while(adv > 0){
		printf("%d %d\n", t, origen) ;
		aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, endpoint, p, adj);
        origen = aux;
		adv--;
	}
	printf("%d %d\n", t, origen) ;
	return t;
}

// advance i triangles arround vertex endpoint
int get_next_and_prev_triangles(int adv, int t, int origen, int endpoint, int *p, int *adj){
	int aux;
	while(adv > 0){
		printf("%d %d\n", t, origen) ;
		aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, endpoint, p, adj);
        origen = aux;
		adv--;
	}
	printf("%d %d\n", t, origen) ;
	return t;
}
/* 
	Busca un triangulo adjacente que comparte el mismo endpoint.
	Origen es el triangulo de donde se viene, -1 si se quiere que se pueda devolver a triangulo anterior.
*/
int get_adjacent_triangle_share_endpoint(int i, int origen, int endpoint, int *p, int *adj){
	int p0 = p[3*i + 0];
	int p1 = p[3*i + 1];
	int p2 = p[3*i + 2];
	
	/* consigue los triangulos adyacentes */
	//int i0 = get_adjacent_triangle(i, p0, p1, p, adj);
	//int i1 = get_adjacent_triangle(i, p1, p2, p, adj);
	//int i2 = get_adjacent_triangle(i, p2, p0, p, adj);
	
	int i0 = adj[3*i + 2];
	int i1 = adj[3*i + 0];
	int i2 = adj[3*i + 1];

	/*verifica si los triangulos son continuos al endpoint */
	int ic0 = is_continuous(i0 ,endpoint, p);
	int ic1 = is_continuous(i1 ,endpoint, p);
	int ic2 = is_continuous(i2 ,endpoint, p);
	
	debug_print("FUNCTION i0 ic0 %d %d   || i1 ic1 %d %d || i2 ic2 %d %d  \n", i0, ic0, i1,ic1,  i2,ic2);
	debug_print("T %d endpoint %d | Triangles %d %d %d | ADJ  %d %d %d\n", i, endpoint, p[3*i + 0], p[3*i + 1], p[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2] );
	if(ic0 != -1 &&  i0 != origen && i0 != -1){ /*Si hay contuinidad y no retrocede al origen */
		return i0;
	}else if(ic1 != -1 && i1 != origen  && i1 != -1){
		return i1;
	}else if(ic2 != -1 && i2 != origen  && i2 != -1){
		return i2;
	}
	return -2;
}


int count_FrontierEdges(int triangle, int *adj){
    int adj_counter = 0;
    int j;
    for(j = 0; j < 3; j++){ 
        if(adj[3*triangle + j] == NO_ADJ){
            adj_counter++;
        }
    }
    return adj_counter;
}


//Search a triangle asociated to a barrier-edge that contains vertex v
//This doesnt use trivertex to find a triangle asociated to v, instead of, search in the list of the triangles one that containts v as frontier-edge
//Input:  index vertex v, array of triangles, array of neigh, number of triangles, array that asociated each triangle to a vertex
//output: index of triangle that have one frontier-edge that contains v
int search_triangle_by_vertex_with_FrontierEdge(int v, int *triangles, int *adj, int tnumber){
	int i,j;
	for (i = 0; i < tnumber; i++)
		for (j = 0; j < 3; j++){
			//If the triangle contains v and has 2 fronter-edges
			if(triangles[3*i +j] == v  && ( adj[3*i + ((j + 1)%3)] == -1 || adj[3*i + ((j + 2)%3)] == -1 )){
				debug_print("v %d |t %d | Triangles %d %d %d | ADJ  %d %d %d\n", v, i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);
				return i;
			}
		}
				
	fprintf(stderr,"%s:%d:%s() No se encontro triangulo que contiene el vertice %d \n",__FILE__,  __LINE__, __func__, v);
    exit(0);
	return -1;
}

//Search a triangle asociated to a barrier-edge that contains vertex v
//This use trivertex to find a triangle asociated to v and travel through  adjacents triangles until find one with frontie-edges
//Input:  index vertex v, array of triangles, array of neigh, number of triangles, array that asociated each triangle to a vertex
//output: index of triangle that have one frontier-edge that contains v
int search_triangle_by_vertex_with_FrontierEdge_from_trivertex(int v, int *triangles, int *adj, int tnumber, int* trivertex){
	int t = trivertex[v];
	int origen = -1;
	int j, aux;
	while (1)
	{
		for (j = 0; j < 3; j++){
			//If the triangle contains v and has 2 fronter-edges
			if(triangles[3*t +j] == v  && ( adj[3*t + ((j + 1)%3)] == -1 || adj[3*t + ((j + 2)%3)] == -1 ))
				return t;
		}
		//avanza al siguiente triangulo
		aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v, triangles, adj);
        origen = aux;

	}	
	fprintf(stderr,"%s:%d:%s() No se encontro triangulo que contiene el vertice %d \n",__FILE__,  __LINE__, __func__, v);
    exit(0);
	return -1;
}

//Given a triangle i, return the triangle adjacent to the triangle origen that containts the vertex v
int search_prev_vertex_to_split(int i, int v, int origen, int *triangles, int *adj){
	int t0, t1,t2;
	int a0, a1, a2;

	t0 = triangles[3*i + 0];
	t1 = triangles[3*i + 1];
	t2 = triangles[3*i + 2];

	a0 = adj[3*i + 0];
	a1 = adj[3*i + 1];
	a2 = adj[3*i + 2];

	debug_print("origen %d, actual %d  | Triangles %d %d %d | ADJ  %d %d %d\n", origen,i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2] );

	if(t1 == v && origen == a0)
			return t2;
	else  if(t2 == v && origen == a0)
			return t1;
	else  if(t0 == v && origen == a1)
			return t2;
	else  if(t2 == v && origen == a1)
			return t0;
	else  if(t0 == v && origen == a2)
			return t1;
	else  if(t1 == v && origen == a2)
			return t0;	
	
	fprintf(stderr,"%s:%d:%s()No se pudo encontrar el vertice anterior para partición \n",__FILE__,  __LINE__, __func__);
    exit(0);
	return -1;
}


//Given a triangle i, return the triangle no adjacent to the triangle origen that containts the vertex v
int search_next_vertex_to_split(int i, int v, int origen, int *triangles, int *adj){
	int t0, t1,t2;
	int a0, a1, a2;

	t0 = triangles[3*i + 0];
	t1 = triangles[3*i + 1];
	t2 = triangles[3*i + 2];

	a0 = adj[3*i + 0];
	a1 = adj[3*i + 1];
	a2 = adj[3*i + 2];

	debug_print("v %d origen %d, actual %d  | Triangles %d %d %d | ADJ  %d %d %d\n", v, origen,i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);

	if(a0 != NO_ADJ && t1 == v && origen != a0)
			return t2;
	else if(a0 != NO_ADJ && t2 == v && origen != a0)
			return t1;		
	else if(a1 != NO_ADJ && t0 == v && origen != a1)
			return t2;
	else if(a1 != NO_ADJ && t2 == v && origen != a1)
			return t0;
	else if(a2 != NO_ADJ && t0 == v && origen != a2)
			return t1;
	else if(a2 != NO_ADJ && t1 == v && origen != a2)
			return t0;

	/*caso particular en poligonos grandes, ya no hay más triangulos para avanzar */
	if(get_adjacent_triangle_share_endpoint(i, origen, v, triangles, adj) == -2){
		debug_msg("No se encuentran más triangulos para avanzar\n");
		return -2;
	}

	fprintf(stderr,"%s:%d:%s()No se pudo encontrar el vertice siguiente para partición \n",__FILE__,  __LINE__, __func__);
    exit(0);
	return -1;
}



/*
error en 147
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 7 | Triangles 41 123 37 | ADJ  -1 31 6
triang.c:440:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 7 | Triangles 41 123 37 | ADJ  -1 31 6
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 8 | Triangles 95 91 37 | ADJ  197 33 42
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 31 | Triangles 31 41 37 | ADJ  7 -1 -1
triang.c:440:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 31 | Triangles 31 41 37 | ADJ  7 -1 -1
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 33 | Triangles 31 95 37 | ADJ  8 -1 -1
triang.c:440:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 33 | Triangles 31 95 37 | ADJ  8 -1 -1
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 196 | Triangles 66 123 37 | ADJ  -1 197 195
triang.c:440:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 196 | Triangles 66 123 37 | ADJ  -1 197 195
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 37 |t 197 | Triangles 66 91 37 | ADJ  8 196 -1

triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 123 |t 2 | Triangles 123 48 74 | ADJ  -1 6 3
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 123 |t 3 | Triangles 123 48 32 | ADJ  -1 195 2
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 123 |t 6 | Triangles 41 123 74 | ADJ  2 -1 7
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 123 |t 7 | Triangles 41 123 37 | ADJ  -1 31 
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 123 |t 195 | Triangles 66 123 32 | ADJ  3 -1 196
triang.c:438:search_triangle_by_vertex_with_FrontierEdge(): v 123 |t 196 | Triangles 66 123 37 | ADJ  -1 197 195

Recorrido, se salta los triangulos 7 y 6

polygon.c:277:generate_polygon(): T 31 Tiene 2 Barrier edge, es oreja, se usa como semilla para generar el poly
polygon.c:325:generate_polygon(): T_inicial 31 | Triangles 31 41 37 | ADJ  7 -1 -1
polygon.c:326:generate_polygon(): initial_point 37 endpoint 41
polygon.c:337:generate_polygon(): origen 31| t 7 | Triangles 41 123 37 | ADJ  -1 31 6
polygon.c:337:generate_polygon(): origen 7| t 6 | Triangles 41 123 74 | ADJ  2 -1 7
polygon.c:337:generate_polygon(): origen 6| t 2 | Triangles 123 48 74 | ADJ  -1 6 3
polygon.c:337:generate_polygon(): origen 2| t 3 | Triangles 123 48 32 | ADJ  -1 195 2
polygon.c:337:generate_polygon(): origen 3| t 195 | Triangles 66 123 32 | ADJ  3 -1 196
polygon.c:337:generate_polygon(): origen 195| t 196 | Triangles 66 123 37 | ADJ  -1 197 195
polygon.c:337:generate_polygon(): origen 196| t 197 | Triangles 66 91 37 | ADJ  8 196 -1
polygon.c:337:generate_polygon(): origen 197| t 8 | Triangles 95 91 37 | ADJ  197 33 42
polygon.c:337:generate_polygon(): origen 8| t 42 | Triangles 63 95 91 | ADJ  8 -1 -1
polygon.c:337:generate_polygon(): origen 42| t 8 | Triangles 95 91 37 | ADJ  197 33 42
polygon.c:337:generate_polygon(): origen 8| t 33 | Triangles 31 95 37 | ADJ  8 -1 -1

*/