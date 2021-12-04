/* Prototipos de funciones para manejar entrada y salida de datos. */

void write_geomview(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles);
void write_svg(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles);
void write_VEM(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles);
void write_VEM_triangles(std::string name, double *r, int *triangles, int *adj, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, std::list <int> &seed_bet);
void write_triangulation(std::string name, double *r, int *triangles, int *adj, int pnumber, int tnumber);
int look_triangles(int i, std::set<int> &s, int * triangles, int * adj, double *r);
void write_alejandro(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region);
void write_alejandro_custom(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region, int *border, int num_boder);
void write_alejandro_quater_circle(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region);

void write_GID(std::string name, double *r, int *triangles, int *adj, int pnumber, int tnumber);


void write_metrics(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region, int num_border, int num_terminal_edges, int num_terminal_border_edges, int num_frontier_edges, int num_frontier_border_edges, int num_interior_edges, int t_delaunay, int t_label, int t_total, int t_travel_and_opt, int t_travel, unsigned int tcost_be, int num_BE, int est_total_be, int est_min_triangles_be, int est_max_triangles_be, int est_poly_with_be, double est_ratio_be);

void read_from_triangle(std::string node_file, std::string ele_file, std::string neigh_file, int &pnumber, int &tnumber, double *&points, int *&triangles, int *&neigh, int *&trivertex);