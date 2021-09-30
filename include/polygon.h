void reverse(int arr[], int n);
int generate_polygon(int i, int * poly, int * triangles, int * adj, double *r);
double get_signed_area_poly(int *poly, int length_poly, double *r);

void print_poly(int *poly, int length_poly);
int copy_poly(int *in, int *out, int len);

void split_poly(int *original_poly, int length_poly, int *poly1, int *length_poly1, int *poly2, int *length_poly2, int e1, int e2);

int has_BarrierEdgeTip(int *poly, int length_poly);
int count_BarrierEdges(int *poly, int length_poly);
int get_vertex_BarrierEdge(int *poly, int length_poly);
int is_BarrierEdge(int i, int *adj, int *adj_copy, int *root_id);
