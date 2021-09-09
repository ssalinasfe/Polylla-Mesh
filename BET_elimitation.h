typedef struct node node;

//Method 1, MAX area criteria
int Remove_BE(int option, int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int i_mesh, int* trivertex);
int Remove_BE2(int option, int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int i_mesh, int* trivertex, std::list <int> &seed_bet);
int optimice1_max_area_criteria(int *t_original, int v_be, int *poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int num_BE, int *triangles, int *adj, double *r, int tnumber);
int optimice2_middle_edge(int *t_original, int v_be, int *triangles, int *adj);
int optimice2_middle_edge_no_memory(int *t_original, int v_be, int *triangles, int *adj);
int generate_polygon_from_BE(int i, int * poly, int * triangles, int * adj, double *r, int ind_poly, node** hashtable_seed);
int Remove_BE3(int option, int *poly, int length_poly, int num_BE, int *triangles, int *adj, double *r, int tnumber, int *mesh, int i_mesh, int* trivertex, std::list <int> &seed_bet, std::vector <int> &seed_bet_mark);
int generate_polygon_from_BE_with_vector(int i, int * poly, int * triangles, int * adj, double *r, int ind_poly, std::vector <int> &seed_bet_mark);