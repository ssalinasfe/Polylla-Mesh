#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>

#include "consts.h"
#include "triang.h"
#include "polygon.h"
#include "mesh.h"
#include "hashtable.h"
#include "BET_elimitation.h"



#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)


int Remove_BE2(int option, int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int i_mesh, int* trivertex, std::list <int> &seed_bet){
    
    int i,j,k,x,y, auxind_poly, ind_poly, ind_poly_after;
    int v_be, v_other;
    int t1, t2, triangle;
    node* hashtable_seed[MAX_HASH];
    //MAX_HASH = num_BE +1;
    //node* hashtable_seed = (node *)malloc(MAX_HASH*sizeof(int));
    node* aux;
    int index;
    //Hash table is initialize
    generate_hash_table(hashtable_seed);
    debug_print("Removiendo %d barrier-edges de ", num_BE); debug_block(print_poly(poly, length_poly); );
    //search by all barrier edge tips and insert edge in the middle    
    for (i = 0; i < length_poly; i++)
    {
        x = i;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y]){
            
            v_be= poly[(i+1) %length_poly];
            debug_print("Encontrado v_be %d %d %d\n", poly[x], v_be, poly[y]);

            
            t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, triangles, adj, tnumber, trivertex);
            //t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
            v_other = optimice2_middle_edge(&t1, v_be, triangles, adj);
            //v_other = optimice2_middle_edge_no_memory(&t1, v_be, triangles, adj);
            if(v_other == -2){
                fprintf(stderr, "Caso critico especial, no encuentra vertices para avanzar en la busqueda de eliminación de barries edge, pero es la primera iteración\n");
                debug_print("v_be - v_other: %d - %d | t1 - t2: %d \n", v_be, v_other, t1);
            }
            t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);
            debug_print("v_be - v_other: %d - %d | t1 - t2: %d - %d\n", v_be, v_other, t1, t2);
            //"<<std::endl;
            // Agregar arista
            if(t2 >= 0){
                adj[3*t1 + get_shared_edge(t1, v_be, v_other, triangles)] = NO_ADJ;
                adj[3*t2 + get_shared_edge(t2, v_be, v_other, triangles)] = NO_ADJ;

                //Guardar triangulos adjacentes a hash t able de largo k;
                index = hashy(t1);
                if(!search(hashtable_seed[index], t1))
                    hashtable_seed[index] = prepend(hashtable_seed[index],t1);
                index = hashy(t2);
                if(!search(hashtable_seed[index], t2))
                    hashtable_seed[index] = prepend(hashtable_seed[index],t2);
            }else{
                index = hashy(t1);
                if(!search(hashtable_seed[index], t1))
                    hashtable_seed[index] = prepend(hashtable_seed[index],t1);
            }
        }
    }
    //Use triangles of hash table to generate polygons
    for (i = 0, ind_poly = 0; i < MAX_HASH; i++)
    {
        while(hashtable_seed[i] != NULL){ //Select the i nth linked list
            
            triangle = hashtable_seed[i]->triangle; //get the head node of the list
            aux =  hashtable_seed[i]; 
            hashtable_seed[i] = hashtable_seed[i]->next; //change the head for the next element
            free(aux); //delete the original head
            //generate polygon with index ind_poly +1 
            ind_poly_after = generate_polygon_from_BE(triangle, poly,triangles,adj,r,ind_poly+1, hashtable_seed); 
            poly[ind_poly] = ind_poly_after - ind_poly - 1; // calculate lenght poly and save it before their vertex
            ind_poly = ind_poly_after;
            seed_bet.push_front(triangle);
        } 
        
    } 
    //save new polygons in mesh
    for(int i = 0; i <ind_poly; i++){
        mesh[i_mesh + i] = poly[i];
    }
    return i_mesh + ind_poly;
    
}

int generate_polygon_from_BE(int i, int * poly, int * triangles, int * adj, double *r, int ind_poly, node** hashtable_seed){
//    int ind_poly = 0;
	
	int initial_point = 0;
	int end_point = 0;
	
	int t0;
	int t1;	
	int t2;
    int ind0;
    int ind1;
    int ind2;
	int continuous;
	int k, j, aux;
	int origen;

    int num_FrontierEdges = count_FrontierEdges(i, adj);
    debug_print("Generando polinomio con triangulo %d FE %d\n", i, num_FrontierEdges);
    /*si tiene 3 se agregan y se corta el ciclo*/
    if (num_FrontierEdges == 3) {
        debug_print("T %d Tiene 3 Frontier edge, se guardan así\n", i);
        poly[ind_poly] = triangles[3 * i + 0];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 1];
        ind_poly++;
        poly[ind_poly] = triangles[3 * i + 2];
        ind_poly++;

        //visited[i] = TRUE;
        return ind_poly;
    } else if(num_FrontierEdges == 2) {
        debug_print("T %d Tiene 2 Frontier edge, es oreja, se usa como semilla para generar el poly\n", i);
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            ind0 = 3*i + j;
            ind1 = 3*i + (j+1)%3;
            ind2 = 3*i + (j+2)%3;
            if(adj[ind0] == NO_ADJ && adj[ind1] == NO_ADJ){
                poly[ind_poly] = triangles[ind1];
                ind_poly++;
                poly[ind_poly] = triangles[ind2];
                ind_poly++;

                initial_point = triangles[ind1];
                end_point = triangles[ind0];  
            }
        }
    }else if (num_FrontierEdges == 1){
        debug_print("T %d Tiene 1 Frontier edge,se usa como FE initial\n", i);
        /*si tiene dos FE se agregan y se empieza el ciclo*/
        for(j = 0; j<3; j++){
            if(adj[3*i + j] == NO_ADJ){
                poly[ind_poly] = triangles[3*i + (j+1)%3];
                ind_poly++;
                initial_point = triangles[3*i + (j+1)%3];

                end_point = triangles[3*i + (j+2)%3];  
            }
        }
    }else {
        end_point = triangles[3*i + 0];
        initial_point = triangles[3*i + 0];
    }
    
    
    /*se marca como visitado */
    //visited[i] = TRUE;
    num_FrontierEdges = 0;
    k = i;
    aux = k;
    k = get_adjacent_triangle_share_endpoint(k, k, end_point, triangles, adj); /* cambia el indice */
    continuous = is_continuous(k, end_point, triangles);
    origen = aux;
//        debug_print("k %d origen %d, conti %d\n", k, origen, continuous);
    debug_print("T_inicial %d | Triangles %d %d %d | ADJ  %d %d %d\n", i, triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2], adj[3*i + 0], adj[3*i + 1], adj[3*i + 2]);
    debug_print("initial_point %d endpoint %d | T_sig %d\n", initial_point, end_point, k);

    int triangugulo_initial = i;
    while (initial_point != end_point || triangugulo_initial != k) {

        /*se marca el triangulo visto como visitado y se suma al area del poligono */
        searchandremove(hashtable_seed[hashy(k)], k);
      //  visited[k] = TRUE;
        t0 = adj[3 * k + 0];
        t1 = adj[3 * k + 1];
        t2 = adj[3 * k + 2];

        num_FrontierEdges = count_FrontierEdges(k, adj);
        debug_print("FE %d | origen %d t %d | Triangles %d %d %d | ADJ  %d %d %d\n", num_FrontierEdges, origen, k, triangles[3*k + 0], triangles[3*k + 1], triangles[3*k + 2], adj[3*k + 0], adj[3*k + 1], adj[3*k + 2]);
        if(origen == -2)
            exit(0);
        if (num_FrontierEdges == 2 && continuous != -1) {
            /* ///////////////////si tiene 2 frontier edge se agregan a poly //////////////////////////////////// */

            if (t0 == NO_ADJ && t1 == NO_ADJ) {
                /*si endpoint es continua a t0  y t0-t1 son fe*/
                if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];
                }
            } else if (t2 == NO_ADJ && t0 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2-t0 son fe*/
                if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ && t2 == NO_ADJ) {
                /*si endpoint es continua a t1 y t1-t2 son fe*/
                if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                }
            } else {
                fprintf(stderr, "** ERROR ** Adding 2 fronter edges\n");
                fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d \n", k, t0, t1, t2, initial_point, end_point);
            }

            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, -1, end_point, triangles, adj); /* se le permite volver al triangulo anterior */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;

        } else if (num_FrontierEdges == 1 && continuous != -1) {
            /* ///////////////////si solo se tiene 1 frontier edge //////////////////////////////////// */
            if (t0 == NO_ADJ) {
                /*si endpoint es continua a t0  y t0 es fe*/
                if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                }
            } else if (t2 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2 es fe*/
                if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ) {
                /*si endpoint es continua a t1  y t1 es fe*/
                if (continuous == 2) {
                    poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                }
            } else {
                fprintf(stderr, "** ERROR ** Adding 1 fronter edges\n");
                fprintf(stderr, "** ERROR ** k %d t %d %d %d ini %d end %d conti %d\n", k, t0, t1, t2, initial_point, end_point, continuous);
            }
            /*si es continuo y tiene 1 fe no puede volver, ind si se guarda  o no*/
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        } else {
            /*si no es continuo no puede regresar de donde venía*/
            aux = k;
            k = get_adjacent_triangle_share_endpoint(k, origen, end_point, triangles, adj); /* cambia el indice */
            continuous = is_continuous(k, end_point, triangles);
            origen = aux;
        }

    }
    
    return ind_poly;
}

//Remove_BE, recursive with exponencial new memory
int Remove_BE(int option, int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int i_mesh, int* trivertex){

    debug_msg("Removiendo barrier edge de "); debug_block(print_poly(poly, length_poly); );

    int *poly1 = (int *)malloc(length_poly*sizeof(int));
	int *poly2 = (int *)malloc(length_poly*sizeof(int));
	int length_poly1;
	int length_poly2;
    int num_BE_poly1;
    int num_BE_poly2;
    int v_be, v_other;
    int t1, t2;


    v_be = get_vertex_BarrierEdge(poly, length_poly);
    
    if(option == 0){
        //t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
        t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, triangles, adj, tnumber, trivertex);
        v_other = optimice1_max_area_criteria(&t1, v_be, poly, length_poly, poly1, &length_poly1, poly2, &length_poly2, num_BE, triangles, adj, r, tnumber);
        t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);
    }else if(option == 1){
        t1 = search_triangle_by_vertex_with_FrontierEdge(v_be, triangles, adj, tnumber);
        //t1 = search_triangle_by_vertex_with_FrontierEdge_from_trivertex(v_be, triangles, adj, tnumber, trivertex);
        v_other = optimice2_middle_edge(&t1, v_be, triangles, adj);
        //v_other = optimice2_middle_edge_no_memory(&t1, v_be, triangles, adj);
        //printf("v %d, t %d  | Triangles %d %d %d | ADJ  %d %d %d\n", v_be, t1, triangles[3*t1 + 0], triangles[3*t1 + 1], triangles[3*t1 + 2], adj[3*t1 + 0], adj[3*t1 + 1], adj[3*t1 + 2]);
        t2 = get_adjacent_triangle(t1, v_other, v_be, triangles, adj);
    }

    debug_print("Eliminando arista de %d - %d de los triangulos  %d y %d", v_be, v_other, t1, t2);
    adj[3*t1 + get_shared_edge(t1, v_be, v_other, triangles)] = NO_ADJ;
    adj[3*t2 + get_shared_edge(t2, v_be, v_other, triangles)] = NO_ADJ;

    //posible bug, si se parte un BE puede omitir un vertice
    //solución, verificar si el triangulo que particiona es válido, sino se cambia.


    length_poly1 = generate_polygon(t1, poly1, triangles, adj, r);
    length_poly2 = generate_polygon(t2, poly2, triangles, adj, r);

    
    num_BE_poly1 = count_BarrierEdges(poly1, length_poly1);
    num_BE_poly2 = count_BarrierEdges(poly2, length_poly2);

    debug_print("num_BE_poly1 %d, num_BE_poly2 %d\n", num_BE_poly1, num_BE_poly2);

    if(num_BE_poly1 > 0 && num_BE_poly2 == 0){						
        debug_msg("Guardando poly2 y enviando recursivamente poly1\n");
        i_mesh = save_to_mesh(mesh, poly2, i_mesh, length_poly2,r);	
        i_mesh = Remove_BE(option, poly1, length_poly1, num_BE_poly1, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
    }else if(num_BE_poly2 > 0 && num_BE_poly1 == 0){
        debug_msg("Guardando poly1 y enviando recursivamente poly2\n");
        i_mesh = save_to_mesh(mesh, poly1, i_mesh, length_poly1,r);	
        i_mesh = Remove_BE(option, poly2, length_poly2, num_BE_poly2, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
    }else if(num_BE_poly1 > 0 && num_BE_poly2 > 0){
        debug_msg("Enviando recursivamente poly1 y poly2\n");
        i_mesh = Remove_BE(option, poly1, length_poly1, num_BE_poly1, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
        i_mesh = Remove_BE(option, poly2, length_poly2, num_BE_poly2, triangles, adj, r, tnumber, mesh, i_mesh, trivertex);
    }else{
        debug_msg("Guardando poly1 y poly2\n");
        i_mesh = save_to_mesh(mesh, poly1, i_mesh, length_poly1,r);	
        i_mesh = save_to_mesh(mesh, poly2, i_mesh, length_poly2,r);	
    }
    free(poly1);
    free(poly2);
    return i_mesh;
}



/* Method 1 to split polygon, max area criteria */

/* Dado un poly con barrier edges
Optimiza la división de este y devuelve poly1 y poly2*/
int optimice1_max_area_criteria(int *t_original, int v_be, int *poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int num_BE, int *triangles, int *adj, double *r, int tnumber){
    double A_poly, A1, A2, opt, r_prev, r_act;
    int v_other, aux, origen, t;
    t = *t_original;

    /*se calculca el valor optimo para el poligono */
    A_poly = get_signed_area_poly(poly, length_poly,r);
    opt = fabs(A_poly/(num_BE+1));

    debug_print("Area poly: %.2lf, opt = %.2lf\n", A_poly, opt);

    /*se calcula el otro vertice para partir poly*/
    v_other = search_next_vertex_to_split(t, v_be, -2, triangles, adj);

    debug_print("Agregar edge %d - %d del Triangulo %d \n", v_be, v_other,t);
    
    if(v_other == -1 || v_other == -2){
        debug_print("No se puede agregar el edge %d - %d del Triangulo %d\n", v_be, v_other, t);
        if(v_other == -2)
            fprintf(stderr, "Caso critico especial, no encuentra vertices para avanzar en la busqueda de eliminación de barries edge, pero es la primera iteración");
        exit(0);
    }

    debug_msg("Dividiendo poligono\n");
    /* Se divide el polygono en dos */
    split_poly(poly, length_poly, poly1, &(*length_poly1), poly2, &(*length_poly2), v_be, v_other);

    debug_msg("poly1: "); debug_block(print_poly(poly1, *length_poly1););
	debug_msg("poly2: "); debug_block( print_poly(poly2, *length_poly2););


    A1 = get_signed_area_poly(poly1, *length_poly1,r);
    A2 = get_signed_area_poly(poly2, *length_poly2,r);

    /* se calcula el r */
    r_prev = fabs(fmin(fabs(A1), fabs(A2)) - opt);
    r_act = 0.0;

    debug_print("A: %.2lf, A1: %.2lf, A2:  %.2lf, A1/A = %.2lf, A2/A = %.2lf, r_prev = %.2lf, r_act = %.2lf\n", A_poly, A1 , A2,  A1/A_poly, A2/A_poly, r_prev, r_act);
    origen = t;
    while (1){
        
        aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v_be, triangles, adj);
        origen = aux;
        v_other = search_next_vertex_to_split(t, v_be, origen, triangles, adj);
        

        debug_print("Agregar edge %d - %d del nuevo triangulo %d | origen = %d \n", v_be, v_other,t, origen); 
        
        if(v_other != -2){
            debug_msg("Dividiendo poligono de nuevo\n");

            split_poly(poly, length_poly, poly1, &(*length_poly1), poly2, &(*length_poly2), v_be, v_other);
            A1 =get_signed_area_poly(poly1, *length_poly1,r);
            A2 = get_signed_area_poly(poly2,*length_poly2,r);
            debug_msg("poly1: "); debug_block(print_poly(poly1, *length_poly1););
            debug_msg("poly2: "); debug_block( print_poly(poly2, *length_poly2););
            
            r_act = fabs(fmin(fabs(A1), fabs(A2)) - opt);
            debug_print("A: %.2lf, A1: %.2lf, A2:  %.2lf, A1/A = %.2lf, A2/A = %.2lf, r_prev = %.2lf, r_act = %.2lf\n", A_poly, A1 , A2,  A1/A_poly, A2/A_poly, r_prev, r_act);
        }
        if (r_act <= r_prev && v_other != -2){
            r_prev = r_act;
            debug_msg("Solución optima no encontrada, repitiendo\n");
        }else{
            debug_print("Se encontro la optimización con r_act %.2lf\n", r_act);
            v_other = search_prev_vertex_to_split(t, v_be, origen, triangles, adj);
            *t_original = t;
            return v_other;
            /*
            debug_print("Agregar edge %d - %d del nuevo triangulo %d | origen = %d \n", v_be, v_other,t, origen);
            split_poly(poly, length_poly, poly1, &(*length_poly1), poly2, &(*length_poly2), v_be, v_other);
            debug_msg("Poligonos generados\n");
			debug_msg("poly1: "); debug_block(print_poly(poly1, *length_poly1););
			debug_msg("poly2: "); debug_block( print_poly(poly2, *length_poly2); );
            debug_block(
            A1 =get_signed_area_poly(poly1, *length_poly1,r);
            A2 = get_signed_area_poly(poly2,*length_poly2,r););
            debug_print("A: %.2lf, A1: %.2lf, A2:  %.2lf, A1/A = %.2lf, A2/A = %.2lf, r_prev = %.2lf, r_act = %.2lf\n", A_poly, A1 , A2,  A1/A_poly, A2/A_poly, r_prev, r_act);
            debug_msg("División optima terminada\n");
            */
            
        }
        
    }
    exit(0);
    return EXIT_FAILURE;
}

int optimice2_middle_edge(int *t_original, int v_be, int *triangles, int *adj){
    
    int aux, origen,i;
    int t_incident[10000];
    int t = *t_original;
    i =0;
    origen = -1;

    while(1){
        debug_print("%d %d %d\n", *t_original, t, origen );

        t_incident[i] = t;
        i++;
        aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v_be, triangles, adj);
        origen = aux;
        if (t<0)
            break;
    }
    debug_print("%d %d %d\n", *t_original, t, origen );
    debug_msg("t_incident"); debug_block(print_poly(t_incident, i); );
    if(i == 1){
        *t_original = t_incident[0];
        //return search_prev_vertex_to_split(t_incident[0], v_be, -1, triangles, adj);
        for(int j = 0; j < 3; j++){
            if(triangles[3*t_incident[0] + j] == v_be)
                return triangles[3*t_incident[0] + (j+1)%3];
        }
    }
    if(i%2 == 0){ //if the triangles surrounding the BET are even 
        i = floor(i/2-1);
        *t_original = t_incident[i];
        //Choose the edge in common with the two middle triangles
        debug_print("search_prev i %d t_incident[i+1] %d v_be %d t_incident[i] %d \n", i, t_incident[i+1], v_be, t_incident[i]);
        return search_prev_vertex_to_split(t_incident[i+1], v_be, t_incident[i], triangles, adj);
    }else{   
        //if the triangles surrounding the BET are odd, edges are even 
        i = floor(i/2);
        *t_original = t_incident[i];
        //Choose any edge of the triangle in the middle; prov is choose due to this always exists
        debug_print("search_next i %d t_incident[i] %d v_be %d t_incident[i-1] %d \n", i, t_incident[i], v_be, t_incident[i-1]);
        return search_next_vertex_to_split(t_incident[i], v_be, t_incident[i-1], triangles, adj);
    }   
}

int optimice2_middle_edge_no_memory(int *t_original, int v_be, int *triangles, int *adj){
    
    int aux, origen,adv;

    int t_incident;
    int t = *t_original;
    t_incident = t;
    adv = 0;
    origen = -1; 
    int t_next = -1, t_prev;
    while (1)
    {
        debug_print("%d %d %d\n", *t_original, t, origen) ;
        adv++;
        aux = t;
        t = get_adjacent_triangle_share_endpoint(t, origen, v_be, triangles, adj);
        origen = aux;
        if (t<0)
            break;
    }
    debug_print("%d %d %d\n", *t_original, t, origen );
    //print_poly(t_incident, i);
    if(adv == 1){
        *t_original = t_incident;
        //return search_prev_vertex_to_split(t_incident[0], v_be, -1, triangles, adj);
        for(int j = 0; j < 3; j++){
            if(triangles[3*t_incident + j] == v_be)
                return triangles[3*t_incident + (j+1)%3];
        }
    }
    if(adv%2 == 0){ //if the triangles surrounding the BET are even 
        adv = adv/2 - 1;
        origen = -1;
        t_prev = advance_i_adjacents_triangles_share_endpoint(adv,t_incident, origen, v_be, triangles, adj);
        *t_original = t_prev;
        //Choose the edge in common with the two middle triangles   
        t_next = get_adjacent_triangle_share_endpoint(t_prev, origen, v_be, triangles, adj);
        debug_print("search_prev adv %d t_next %d v_be %d t_prev %d \n", adv, t_next, v_be, t_prev);    
        return search_prev_vertex_to_split(t_next, v_be, t_prev, triangles, adj);
    }else{   
        //if the triangles surrounding the BET are odd, edges are even
        //Choose any edge of the triangle in the middle; prov is choose due to this always exists
        adv = adv/2;
        t_next = advance_i_adjacents_triangles_share_endpoint(adv,t_incident, t_prev, v_be, triangles, adj);
        *t_original = t_next;
        //t_next = get_adjacent_triangle_share_endpoint(t_prev, -1, v_be, triangles, adj);
        debug_print("search_next adv %d t_next %d v_be %d t_prev %d \n", adv, t_next, v_be, t_prev);
        return search_next_vertex_to_split(t_next, v_be, t_prev, triangles, adj);
    }   
}

/*
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:35:Remove_BE2(): Removiendo barrier edge de (17) 2 66 39 61 90 54 82 87 74 45 21 45 74 87 24 97 19
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 129 -1
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 142 129
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 132 142
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 123 132
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 149 123
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 153 149
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:491:optimice2_middle_edge(): 129 147 153
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:501:optimice2_middle_edge(): 129 -2 147
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:502:optimice2_middle_edge(): t_incident(7) 129 142 132 123 149 153 147
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:515:optimice2_middle_edge(): search_next i 3 t_incident[i] 123 v_be 21 t_incident[i-1] 132 
/home/cuyi/Dropbox/Doctorado/paper_examen/algo_with_detri2/BET_elimitation.cpp:50:Remove_BE2(): v_be - v_other: 21 - 66 | t1 - t2: 123 - 149
*/
