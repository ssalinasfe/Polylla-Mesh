void generate_delaunay_from_random_points(int argc, char* argv[], int &pnumber, int &tnumber);
void copy_delaunay_arrays(int tnumber, double *r, int* triangles, int* adj);
int get_border_points(int pnumber, int tnumber, int *border, int * triangles, int * adj, double *r);
void generate_constrained_vonoronoi();
void free_detri2();
/*
Por hacer 

int is_nomax_nomax(int i, int j, int *p, int *max);
int is_max_max(int i, int j, int *p, int *max);
int get_adjacent_triangle(int i, int k, int l, int *p, int *ady);
int same_edge(int u, int v, int w, int x);

------------------------------------------------------
int is_continuous(int i, int endpoint, int *p );
int get_adjacent_triangle_share_endpoint(int i, int origen, int endpoint, int *p, int *adj);
int count_FrontierEdges(int triangle, int *adj);
int search_triangle_by_vertex_with_FrontierEdge(int v, int *triangles, int *adj, int tnumber);
int search_next_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
int search_prev_vertex_to_split(int i, int v, int origen, int *triangles, int *adj);
int  get_shared_edge(int i, int u, int v, int *p);
*/

//#include <math.h>
//#define NO_ADJ -1
/*
double dist(double x0, double y0, double x1, double y1);

double dist(double x0, double y0, double x1, double y1)
{
	return sqrt(pow(x0 - x1, 2.0) + pow(y0 - y1, 2.0));
}

//usar composición clases

class TMesh : public detri2::Triangulation{

public:
    detri2::Triangulation *trimesh;
    
    int tnumber;
    int pnumber;
    
    TMesh(int nparam, char* params[]){
        trimesh = new detri2::Triangulation();
        trimesh->parse_commands(nparam, params);
        trimesh->read_nodes();
        trimesh->incremental_delaunay();

        this->tnumber =  trimesh->tr_tris->objects - trimesh->ct_hullsize;
        this->pnumber = trimesh->ct_in_vrts + (trimesh->tr_steiners != NULL ? trimesh->tr_steiners->objects : 0);          
        if (!trimesh->io_keep_unused) { // no -IJ
            this->pnumber -= trimesh->ct_unused_vrts;
        }
    
    }

    ~TMesh(){
        delete trimesh;
        delete max;
    }

    void print(); 

    detri2::Triang *getTriangle(int i) const{
        return (detri2::Triang *) trimesh->tr_tris->get(i);
    }

    int max_edge_index(int i);
    int is_nomax_nomax(int i, int j);

    int same_triangle(detri2::Triang *i, detri2::Triang *j)
    {
        if (i->is_deleted() || i->is_hulltri() || j->is_deleted() || j->is_hulltri()) 
            return -1;
        return i->vrt[0] == j->vrt[0] && i->vrt[1] == j->vrt[1] && i->vrt[2] == j->vrt[2];
    }

};
*/
/*
int TMesh::is_nomax_nomax(int i, int j)
{
	detri2::Triang* const tri1 = getTriangle(i);
    detri2::Triang* const tri2 = getTriangle(j);

    int p0,p1,p2,p3;
    int indi,indj;
    for(i = 0; i < 3; i++)
    {
        if( same_triangle(tri1->nei[indi].esym().tri, tri2) )
            break;
    }
    for(indj = 0; i < 3; i++)
    {
        if( same_triangle(tri2->nei[indj].esym().tri, tri1) )
            break;
    }

    return (indj != max[j]) && (indi != max[i]);
    
}

//Devuelve el indice del máximo edge
//Input; ind triangulo
//outout; ind edge más largo
int TMesh::max_edge_index(int i){
	double l0;
	double l1;
	double l2;
	
	detri2::Vertex *p0;
	detri2::Vertex *p1;
	detri2::Vertex *p2;
	detri2::Triang* const tri = getTriangle(i);
	
    p0 = tri->vrt[0];
    p1 = tri->vrt[1];
    p2 = tri->vrt[2];
	
	l0 = dist(p0->crd[0], p0->crd[1], p1->crd[0], p1->crd[1]);
    l0 = dist(p1->crd[0], p1->crd[1], p2->crd[0], p2->crd[1]);
	l0 = dist(p2->crd[0], p2->crd[1], p0->crd[0], p0->crd[1]);
	
	if((l0 >= l1 && l1 >= l2) || (l0 >= l2 && l2 >= l1))
	{
		return 0;
	}
	else if((l1 >= l0 && l0 >= l2) || (l1 >= l2 && l2 >= l0))
	{
		return 1;
	}
	else
	{
		return 2;
	}
}
*/

/*
void TMesh::print(){
    int idx = trimesh->io_firstindex;

    std::cout<<"vertices:"<<trimesh->ct_in_vrts<<std::endl;
    for (int i = 0; i < trimesh->ct_in_vrts; i++)
    {   
        detri2::Vertex *v = &(trimesh->in_vrts[i]);
        std::cout<<"v "<<v->idx<<" ("<<v->crd[0]<<", "<<v->crd[1]<<")"<<std::endl;
        idx++;
    }

    std::cout<<"Triangles: "<< trimesh->tr_tris->objects - trimesh->ct_hullsize<<std::endl;

    // Index all triangles (hull triangles all have index -1).
    idx = trimesh->io_firstindex;
    for (int i = 0; i < trimesh->tr_tris->used_items; i++) {
        detri2::Triang* tri = (detri2::Triang *) trimesh->tr_tris->get(i);
        if (tri->is_deleted()) continue;
        if (tri->is_hulltri()) {
            tri->idx = -1;
        } else {
            tri->idx = idx;
            idx++;
        }
    }
    idx = trimesh->io_firstindex;
    for (int i = 0; i < trimesh->tr_tris->used_items; i++)
    {
        
        detri2::Triang* tri = (detri2::Triang *) trimesh->tr_tris->get(i);
        if (tri->is_deleted() || tri->is_hulltri()) continue;
        std::cout<<"t "<<idx<<" | "<< tri->vrt[0]->idx <<" "<< tri->vrt[1]->idx<<" "<<tri->vrt[2]->idx<<" | "<<
            tri->nei[0].tri->idx<<" "<<tri->nei[1].tri->idx<<" "<<tri->nei[2].tri->idx<<std::endl;
        idx ++;
    }

    trimesh->save_neighbors();
    trimesh->save_triangulation();   
}
*/