/* Prototipos de funciones para manejar triangulación. */

double dist(double x0, double y0, double x1, double y1);
int max_edge_index(int i, double *r, int *p);
int is_nomax_nomax(int i, int j, int *p, int *max);
int is_max_max(int i, int j, int *p, int *max);
int get_adjacent_triangle(int i, int k, int l, int *p, int *ady);
int same_edge(int u, int v, int w, int x);
int is_continuous(int i, int endpoint, int *p );
int get_adjacent_triangle_share_endpoint(int i, int origen, int endpoint, int *p, int *adj);
int count_FrontierEdges(int triangle, int *adj);
int search_triangle_by_vertex_with_FrontierEdge(int v, int *triangles, int *adj, int tnumber);
int search_next_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
int search_prev_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
int  get_shared_edge(int i, int u, int v, int *p);
int Equality(float a, float b, float epsilon);
int GreaterEqualthan(float a, float b, float epsilon);
int is_max_nomax(int i, int j, int *p, int *max);
int advance_i_adjacents_triangles_share_endpoint(int adv, int t, int origen, int endpoint, int *p, int *adj);
int search_triangle_by_vertex_with_FrontierEdge_from_trivertex(int v, int *triangles, int *adj, int tnumber, int* trivertex);