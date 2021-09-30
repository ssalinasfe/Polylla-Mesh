/* Funciones para manejar poligonos. */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "consts.h"
#include "triang.h"
#include "polygon.h"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)


// Function to reverse elements of an array
void reverse(int arr[], int n)
{
    for (int low = 0, high = n - 1; low < high; low++, high--)
    {
        int temp = arr[low];
        arr[low] = arr[high];
        arr[high] = temp;
    }
}

//Verifica si es válido como triangulo semilla del poligono
int is_BarrierEdge(int i, int *adj, int *adj_copy, int *root_id){
    //debug_print("i = %d, root_id[i] = %d, root_id[t0_adj] = %d, root_id[t1_adj] = %d, root_id[t2_adj] = %d\n", i, root_id[i], root_id[t0_adj], root_id[t1_adj], root_id[t2_adj]  );
    //if( root_id[i] != root_id[t0_adj] && root_id[i] != root_id[t1_adj] && root_id[i] != root_id[t2_adj] ){
    int j;
    int num_fe = count_FrontierEdges(i, adj);
    debug_print("FE %d, Triangulo %d\n", num_fe, i);
    if(num_fe == 3){
        return 1;
    }else if (num_fe == 1)
    {
        for (j = 0; j < 3; j++)
        {
            if (adj_copy[3*i + j] != TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[3*i + j]])
                    return 1;
            }else
            {
                return 1;
            }
            
                
        }   
    }else if(num_fe == 2)
    {
        for (j = 0; j < 3; j++)
        {
            int ind = 3*i + j;
            int ind2 = 3*i + (j + 1)%3;
            if (adj_copy[ind] != TRIANG_BORDER && adj_copy[ind2] != TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[ind]] && root_id[i] !=  root_id[adj_copy[ind2]] )
                    return 1;
            }else if(adj_copy[ind] == TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[ind2]] )
                    return 1;
            }else if(adj_copy[ind2] == TRIANG_BORDER){
                if (root_id[i] !=  root_id[adj_copy[ind]] )
                    return 1;
            }
        }   
    }
    return 0;
}




/* Divide un poly dado un vertice e1-e2
    resultados poly1 y poly*/
void split_poly(int *original_poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int e1, int e2){
    int pos1, pos2,i;
    pos1= -1;
    pos2= -1;
    for(i =0; i< length_poly; i++)
        if(original_poly[i] == e1 || original_poly[i] == e2){
            pos1 = i;
            break;
        }
    for(i =pos1 + 1; i< length_poly; i++)
        if(original_poly[i] == e1  || original_poly[i] == e2){
            pos2 = i;
            break;
        }
    debug_print("Divide pos1: %d, pos2: %d \n", pos1, pos2);
    if(pos1 == -1 || pos2 == -1){
        fprintf(stderr, "%s:%d:%s(): No se encontro pos1 o pos2\n",__FILE__,  __LINE__, __func__);
        exit(0);
    }
    
    *length_poly1 = abs(pos1-pos2) +1;
    *length_poly2 = length_poly - *length_poly1 +2;

    for (i = 0; i < *length_poly1 ; i++)
        poly1[i] = original_poly[(pos1 + i) %length_poly];

    for (i = 0; i < *length_poly2 ; i++)
        poly2[i] = original_poly[(pos2 + i) %length_poly];
}

int copy_poly(int *in, int *out, int len){
    int i;
    for (i = 0; i < len; i++)
        out[i] = in[i];
    return len;
}

void print_poly(int *poly, int length_poly){
    int i;
    fprintf(stderr,"(%d)", length_poly);
    for (i = 0; i < length_poly; i++)
        fprintf(stderr," %d", poly[i]); 
    fprintf(stderr,"\n");
}

double get_signed_area_poly(int *poly, int length_poly, double *r){
    double area = 0.0;
    double x1,y1,x2,y2;
    int i,j;
    for (i = 0; i < length_poly; i++)
    {
        j = (i+1) %length_poly;
        x1=r[2*poly[i] + 0];
        y1=r[2*poly[i] + 1];
        x2=r[2*poly[j] + 0];
        y2=r[2*poly[j] + 1];
        area += (x1 + x2)*(y2 - y1);
    }
    return area/2;
}




int count_BarrierEdges(int *poly, int length_poly){
    int count = 0;
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            count++;
    }
    return count;
}

int has_BarrierEdgeTip(int *poly, int length_poly){
    int count = 0;
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            return 1;
    }
    return 0;
}


int get_vertex_BarrierEdge(int *poly, int length_poly){
    int x, y,i;
    for (i = 0; i < length_poly; i++)
    {
        x = i % length_poly;
        y = (i+2) % length_poly;
        if (poly[x] == poly[y])
            return poly[(i+1) %length_poly];
    }
    fprintf(stderr,"num_BE %d\n", count_BarrierEdges(poly, length_poly));
    fprintf(stderr,"%s:%d:%s(): No se encontro vertice BarrierEdge\n",__FILE__,  __LINE__, __func__);
    exit(0);
    return -1;
}


int generate_polygon(int i, int * poly, int * triangles, int * adj, double *r) {
    int ind_poly = 0;
	
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
