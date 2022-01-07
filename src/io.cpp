#include <iostream>
#include <vector> 
#include <algorithm>
#include <set>
#include <list>
#include <string>
#include <math.h>  
#include <float.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "consts.h"
#include "triang.h"
#include "io.h"
#include "metrics.h"


#define filespath "input/"
#define filespathoutput "output/"

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif

#define debug_block(fmt) do { if (DEBUG_TEST){ fmt }} while (0)
#define debug_print(fmt, ...) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__); } while (0)
#define debug_msg(fmt) do { if (DEBUG_TEST) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__,  __LINE__, __func__); } while (0)


void read_from_triangle(std::string node_file, std::string ele_file, std::string neigh_file, int &pnumber, int &tnumber, double *&points, int *&triangles, int *&neigh, int *&trivertex){
    std::string line;
    std::ifstream nodefile(node_file);
    double a1, a2, a3, a4;
    int i = 0;
    
    //std::cout<<"Node file"<<std::endl;
    if (nodefile.is_open())
    {
        nodefile >> pnumber ;
        //std::cout<<pnumber<<std::endl;

        std::getline(nodefile, line); 
        points = (double *)malloc(2*pnumber*sizeof(double));
        while (std::getline(nodefile, line))
        {
            std::stringstream(line) >> std::ws >> a1 >> std::ws >> a2 >> std::ws >> a3;   
            points[2*i + 0] = a2;
            points[2*i + 1] = a3;
            //std::cout<<points[2*i + 0]<<" "<<points[2*i + 1]<<std::endl;
            i++;
        }        
    }
    else 
        std::cout << "Unable to open node file"; 

    nodefile.close();


    trivertex =(int *)malloc(pnumber*sizeof(int));
    //std::cout<<"Ele file"<<std::endl;
    std::ifstream elefile(ele_file);
    int t1, t2, t3, t4;
    i = 0;
    if(elefile.is_open()){
        elefile >> tnumber ;
        triangles = (int *)malloc(3*tnumber*sizeof(int));
        std::getline(elefile, line); 
        while (elefile >> t1 >> t2 >> t3 >> t4)
        {
            //std::cout<<t2<<" "<<t3<<" "<<t4<<std::endl;
            triangles[3*i + 0] = t2;
            triangles[3*i + 1] = t3;
            triangles[3*i + 2] = t4;
            trivertex[t2] = i;
            trivertex[t3] = i;
            trivertex[t4] = i;
            //std::cout<<triangles[3*i + 0]<<" "<<triangles[3*i + 1]<<" "<<triangles[3*i + 2]<<std::endl;
            i++;
            std::getline(elefile, line); 
        }
    }else std::cout << "Unable to open ele file";

    elefile.close();

    //std::cout<<"Neigh file"<<std::endl;
    std::ifstream neighfile(neigh_file);
    i = 0;
    if(neighfile.is_open()){
        std::getline(neighfile, line); 
        neigh =(int *)malloc(3*tnumber*sizeof(int));
        while (neighfile >> t1 >> t2 >> t3 >> t4 )
        {
            neigh[3*i + 0] = t2;
            neigh[3*i + 1] = t3;
            neigh[3*i + 2] = t4;
            //std::cout<<t2<<" "<<t3<<" "<<t4<<std::endl;
            i++;
        }
    }else std::cout << "Unable to open neigh file";
    neighfile.close();

    ////std::cout<<"Neigh file"<<std::endl;
    //std::ifstream trivertexfile(name + ".trivertex");
    //i = 0;
    //if(trivertexfile.is_open()){
    //    std::getline(trivertexfile, line); 
    //    trivertex =(int *)malloc(pnumber*sizeof(int));
    //    while (trivertexfile >> t1 >> t2)
    //    {
    //        trivertex[i] = t2;
    //        i++;
    //    }
    //}else std::cout << "Unable to open neigh file";
    //trivertexfile.close();
}

/*geomview output*/
void write_geomview(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles){

    int i,j;
    char cmd[1024] = "\0";
    strcat(cmd, name.c_str());
    strcat(cmd,".off");
    //printf("Wiriting off file in %s\n", cmd);
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr, "{ appearance  {+edge +face linewidth 2} LIST\n");
    fprintf(fptr, "OFF\n");
    fprintf(fptr,"%d %d 0\n", pnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.5f %.5f 0\n", r[2*i + 0], r[2*i + 1]);

  //imprimir polginos
    i = 0;
    while(i < i_mesh){
        int length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i]);
            i++;
        }
        fprintf(fptr, "\n");
    }

    if(print_triangles){
        
        fprintf(fptr, "{ appearance  {+edge -face linewidth 2} LIST\n");
        int p0, p1,p2;
        for(i = 0; i < tnumber; i++){
            p0 = 3*i + 0;
            p1 = 3*i + 1;
            p2 = 3*i + 2;
            fprintf(fptr,"# %d %d\n", p0, p1);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.5f %.5f 0\n%.5f %.5f 0\n",	r[2*triangles[p0]+0], r[2*triangles[p0]+1], 
                                r[2*triangles[p1]+0], r[2*triangles[p1]+1]);
            fprintf(fptr,"0 1 1 1\n");

            fprintf(fptr,"# %d %d\n", p1, p2);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.5f %.5f 0\n%.5f %.5f 0\n",	r[2*triangles[p1]+0], r[2*triangles[p1]+1], 
                                r[2*triangles[p2]+0], r[2*triangles[p2]+1]);
            fprintf(fptr,"0 1 1 1\n");

            fprintf(fptr,"# %d %d\n", p0, p2);
            fprintf(fptr,"VECT 1 2 1 2 1\n");
            fprintf(fptr,"%.5f %.5f 0\n%.5f %.5f 0\n",	r[2*triangles[p0]+0], r[2*triangles[p0]+1], 
                                r[2*triangles[p2]+0], r[2*triangles[p2]+1]);
            fprintf(fptr,"0 1 1 1\n");
        }
    }
    
    fprintf(fptr," }\n");
    if(print_triangles){
        fprintf(fptr," }\n");
        fprintf(fptr," }\n");
    }
    /*
    std::sort( border_point.begin(), border_point.end() );
    border_point.erase( std::unique( border_point.begin(), border_point.end() ), border_point.end() );
    fprintf(fptr,"#Border vertices\n#");
    for ( i=0; i<border_point.size(); i++)
    {
        fprintf(fptr,"%d ", border_point[i]);
    }
    */
    fprintf(fptr,"\n");
    fclose(fptr);
   // std::cout<<"Output saved in: "<<cmd<<std::endl;
}

int r8_to_i4 ( double xmin, double xmax, double x, int ixmin, int ixmax )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
//
//  Discussion:
//
//    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
//    IX := min ( IX, max ( IXMIN, IXMAX ) )
//    IX := max ( IX, min ( IXMIN, IXMAX ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XMIN, XMAX, the real range.  XMAX and XMIN must not be
//    equal.  It is not necessary that XMIN be less than XMAX.
//
//    Input, double X, the real number to be converted.
//
//    Input, int IXMIN, IXMAX, the allowed range of the output
//    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
//    It is not necessary that IXMIN be less than IXMAX.
//
//    Output, int R8_TO_I4, the value in the range [IXMIN,IXMAX] that
//    corresponds to X.
//
{
    int ix;
    double temp;

    /*
    if ( xmax == xmin )
    {
    cerr << "\n";
    cerr << "R8_TO_I4 - Fatal error!\n";
    cerr << "  XMAX = XMIN, making a zero divisor.\n";
    cerr << "  XMAX = " << xmax << "\n";
    cerr << "  XMIN = " << xmin << "\n";
    exit ( 1 );
    }*/

    temp =
        ( ( xmax - x        ) * ( double ) ixmin
        + (        x - xmin ) * ( double ) ixmax )
        / ( xmax     - xmin );

    if ( 0.0 <= temp )
    {
        temp = temp + 0.5;
    }
    else
    {
        temp = temp - 0.5;
    }

    ix = ( int ) temp;

return ix;
}

//basado en: https://people.math.sc.edu/Burkardt/cpp_src/triangulation_svg/triangulation_svg.cpp
void write_svg(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles){
    int i,j;
    int j4, j4_min, j4_max;
    int i4, i4_max, i4_min;
    double x_max, y_max, x_min, y_min, x_scale, y_scale;    

    //Determinar escala
    //busca máximos y min
    for(i = 0; i < pnumber; i++){
        //search range x
        if(r[2*i + 0] > x_max )
            x_max = r[2*i + 0];
        if(r[2*i + 0] < x_min )
            x_min = r[2*i + 0];
        //search range y
        if(r[2*i + 1] > y_max )
            y_max = r[2*i + 1];
        if(r[2*i + 1] < y_min )
            y_min = r[2*i + 1];
    }
    
    x_scale = x_max - x_min;
    x_max = x_max + 0.05 * x_scale;
    x_min = x_min - 0.05 * x_scale;
    x_scale = x_max - x_min;
    y_scale = y_max - y_min;
    y_max = y_max + 0.05 * y_scale;
    y_min = y_min - 0.05 * y_scale;
    y_scale = y_max - y_min;

    i4_min = 1;
    j4_min = 1;

    if ( x_scale < y_scale )
    {
        i4_max = ( int ) ( 0.5 + 500.0 * x_scale / y_scale );
        j4_max = 500;
    }
    else
    {
        i4_max = 500;
        j4_max = ( int ) ( 0.5 + 500.0 * y_scale / x_scale );
    }
    char cmd[1024] = "\0";
    strcat(cmd, name.c_str());
    strcat(cmd,".svg");
    printf("Writting svg in %s\n", cmd);
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr, "<svg width=\"%d\" height=\"%d\" viewbox=\"%d %d %d %d \" fill=\"white\">\n", i4_max, j4_max,  i4_min, j4_min, i4_max, j4_max);
  //imprimir polginos
    i = 0;
    while(i < i_mesh){
        int length_poly = mesh[i];
        i++;
        fprintf(fptr, "<polygon points=\"");
        for(j=0; j < length_poly;j++){
            i4 = r8_to_i4( x_min, x_max, r[0+mesh[i]*2], i4_min, i4_max);
            j4 = r8_to_i4( y_max, y_min, r[1+mesh[i]*2], j4_min, j4_max);
            fprintf(fptr, " %d %d", i4, j4);
            i++;
        }
        fprintf(fptr, "\" stroke=\"black\" stroke-width=\"2\" />\n");
    }
    
    fprintf(fptr,"</svg>\n");
    fclose(fptr);
}

void write_metrics(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region, int num_border, int num_terminal_edges, int num_terminal_border_edges, int num_frontier_edges, int num_frontier_border_edges, int num_interior_edges, int t_delaunay, int t_label, int t_total, int t_travel_and_opt, int t_travel, unsigned int tcost_be, int num_BE, int est_total_be, int est_min_triangles_be, int est_max_triangles_be, int est_poly_with_be, double est_ratio_be)
{
    
    int i,j,length_poly = 0,k;
    int poly[100];

    double avg_edge_size = 0;
    int edges, est_max_edges, est_min_edges,est_total_edges;
	edges = mesh[0];
	est_max_edges = edges;
	est_min_edges = edges;
	est_total_edges = edges;
	i = 0;
	while(i < i_mesh){
        edges = mesh[i];
		est_max_edges = edges > est_max_edges ? edges : est_max_edges;
		est_min_edges = edges < est_min_edges ? edges : est_min_edges;
		est_total_edges += edges;
        i++;
        for(j=0; j < edges;j++){
            poly[j] = mesh[i];              
			i++;
        }
        for(j=0; j < edges;j++){
            k = (j + 1)%edges;
            avg_edge_size += dist(r[2*poly[j]+0],r[2*poly[j]+1],r[2*poly[k]+0],r[2*poly[k]+1]);
        }
    }
    avg_edge_size = avg_edge_size/est_total_edges;

    double metrics[11];
	double mesh_sparcing = 0.0;
	i = 0;

    double max_area = 0, min_area = DBL_MAX, avg_area = 0;
    double max_perimeter = 0, min_perimeter = DBL_MAX, avg_perimeter = 0;
    double max_apr = 0, min_apr = DBL_MAX, avg_apr = 0;
    double min_angle = 360.0, max_angle = 0;
    double min_edge = DBL_MAX, max_edge = 0;
    double min_edgeratio = DBL_MAX, max_edgeratio = 0, avg_edgeratio = 0;
    double min_pointDistance = DBL_MAX, max_pointDistance = 0;
    double min_radius = DBL_MAX, max_radius = 0, avg_radius = 0;
    double area, perimeter, apr,angle,edge, edgeratio, pointDistance, radius;

	while(i < i_mesh){
        length_poly = mesh[i];
        i++;
        for(j=0; j < length_poly; j++){
			poly[j] = mesh[i];            
			i++;
        }

        get_metrics_of_polygon(metrics, poly, length_poly, r);


        //area metrics
        //std::cout<<"area: "<<metrics[0]<<std::endl;
        area = metrics[0];
        avg_area += abs(area);
        if(area < min_area) min_area = area;
        if(area > max_area) max_area = area;

        //std::cout<<"perimeter: "<<metrics[1]<<std::endl;
        perimeter = metrics[1];
        avg_perimeter += perimeter;
        if(perimeter < min_perimeter) min_perimeter = perimeter;
        if(perimeter > max_perimeter) max_perimeter = perimeter;

        //std::cout<<"area_perimeter_ratio: "<<metrics[2]<<std::endl;
        apr = metrics[2];
        avg_apr += apr;
        if(apr < min_apr) min_apr = apr;
        if(apr > max_apr) max_apr = apr;

        //std::cout<<"min_angle: "<<metrics[3]<<std::endl;
        angle = metrics[3];
        if(angle < min_angle) min_angle = angle;

        //std::cout<<"max_angle: "<<metrics[4]<<std::endl;
        angle = metrics[4];
        if(angle > max_angle) max_angle = angle;

        //std::cout<<"min_edge: "<<metrics[5]<<std::endl;
        edge = metrics[5];
        if(edge < min_edge) min_edge = edge;

        //std::cout<<"max_edge: "<<metrics[6]<<std::endl;
        edge = metrics[6];
        if(edge > max_edge) max_edge = edge;

        //std::cout<<"edge_rati: "<<metrics[7]<<std::endl;
        edgeratio = metrics[7];
        avg_edgeratio += edgeratio;
        if(edgeratio < min_edgeratio) min_edgeratio = edgeratio;
        if(edgeratio > max_edgeratio) max_edgeratio = edgeratio;

        //std::cout<<"min_pd: "<<metrics[8]<<std::endl;
        pointDistance = metrics[8];
        if(pointDistance < min_pointDistance) min_pointDistance = pointDistance;
        //std::cout<<"max_pd: "<<metrics[9]<<std::endl;
        pointDistance = metrics[9];
        if(pointDistance > max_pointDistance) max_pointDistance = pointDistance;        

        //std::cout<<"circumscribed:circle: "<<metrics[10]<<std::endl;   
        radius = metrics[1];
        avg_radius += radius;
        if(radius < min_radius) min_radius = radius;
        if(radius > max_radius) max_radius = radius;

    }

    mesh_sparcing = std::sqrt(avg_area/num_region);
    avg_area = avg_area/num_region;
    avg_perimeter = avg_perimeter/num_region;
    avg_edgeratio = avg_edgeratio/num_region;
    avg_radius = avg_radius/num_region;
	
    char cmd[1024] = "\0";
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_metrics.json");
    printf("Writting metrics in %s\n",cmd);
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }      
    fprintf(fptr, "{\n");
    fprintf(fptr, "\"pnumber\": %d,\n", pnumber);
    fprintf(fptr, "\"tnumber\": %d,\n", tnumber);
    fprintf(fptr, "\"num_regions\": %d,\n", num_region);
    fprintf(fptr, "\"num_border\": %d,\n", num_border);
    fprintf(fptr, "\"avg_area\": %.8f,\n", avg_area);
    fprintf(fptr, "\"mesh_sparcing\": %.12f,\n", mesh_sparcing);
    fprintf(fptr, "\"min_area\": %.8f,\n", min_area);
    fprintf(fptr, "\"max_area\": %.8f,\n", max_area);
    fprintf(fptr, "\"avg_perimeter\": %.9f,\n", avg_perimeter);
    fprintf(fptr, "\"min_perimeter\": %.9f,\n", min_perimeter);
    fprintf(fptr, "\"max_perimeter\": %.9f,\n", max_perimeter);
    fprintf(fptr, "\"avg_radius\": %.9f,\n", avg_radius);
    fprintf(fptr, "\"min_radius\": %.9f,\n", min_radius);
    fprintf(fptr, "\"max_radius\": %.9f,\n", max_radius);
    fprintf(fptr, "\"avg_area_perimeter_ratio\": %.9f,\n", avg_apr);
    fprintf(fptr, "\"min_area_perimeter_ratio\": %.9f,\n", min_apr);
    fprintf(fptr, "\"max_area_perimeter_ratio\": %.9f,\n", max_apr);
    fprintf(fptr, "\"min_angle\": %.9f,\n", min_angle);
    fprintf(fptr, "\"max_angle\": %.9f,\n", max_angle);
    fprintf(fptr, "\"avg_edge_size\": %.9f,\n", avg_edge_size);
    fprintf(fptr, "\"min_edge\": %.9f,\n", min_edge);
    fprintf(fptr, "\"max_edge\": %.9f,\n", max_edge);
    fprintf(fptr, "\"avg_edgeratio\": %.9f,\n", avg_edgeratio);
    fprintf(fptr, "\"min_edgeratio\": %.9f,\n", min_edgeratio);
    fprintf(fptr, "\"max_edgeratio\": %.9f,\n", max_edgeratio);
    fprintf(fptr, "\"min_pointDistance\": %.9f,\n", min_pointDistance);
    fprintf(fptr, "\"max_pointDistance\": %.9f,\n", max_pointDistance);
    fprintf(fptr, "\"poly_with_be\": %d,\n", est_poly_with_be);
    fprintf(fptr, "\"total_be\": %d,\n", est_total_be);
    fprintf(fptr, "\"min_poly_be\": %d,\n", est_min_triangles_be);
    fprintf(fptr, "\"max_poly_be\": %d,\n", est_max_triangles_be);
    fprintf(fptr, "\"ratio_be_per_poly\": %.9f,\n", (est_poly_with_be > 0 ? est_ratio_be/est_poly_with_be : 0.0));
    fprintf(fptr, "\"num_total_edges\": %d,\n", est_total_edges);
    fprintf(fptr, "\"num_max_edges\": %d,\n", est_max_edges);
    fprintf(fptr, "\"num_min_edges\": %d,\n", est_min_edges);
    fprintf(fptr, "\"num_edges_by_poly\": %.9f,\n", (float)est_total_edges/num_region);
    fprintf(fptr, "\"num_terminal_edges_no_border\": %d,\n", num_terminal_edges/2);
    fprintf(fptr, "\"num_terminal_border_edges\": %d,\n", num_frontier_border_edges);
    fprintf(fptr, "\"num_terminal_edges_total\": %d,\n", num_frontier_border_edges + num_terminal_edges/2);
    fprintf(fptr, "\"num_frontier_edges_no_border\": %d,\n", num_frontier_edges/2);
    fprintf(fptr, "\"num_frontier_border_edges\": %d,\n", num_frontier_border_edges);
    fprintf(fptr, "\"num_frontier_edges\": %d,\n", num_frontier_edges/2 + num_frontier_border_edges);
    fprintf(fptr, "\"num_interior_edges\": %d,\n", num_interior_edges/2);
    fprintf(fptr, "\"t_delaunay\": %d,\n",t_delaunay);
    fprintf(fptr, "\"t_total\": %d,\n", t_total);
    fprintf(fptr, "\"t_label\": %d,\n", t_label);
    fprintf(fptr, "\"t_travel_and_opt\": %d,\n", t_travel_and_opt);
    fprintf(fptr, "\"t_travel\": %d,\n", t_travel);
    fprintf(fptr, "\"t_opt\": %d\n", tcost_be);
    fprintf(fptr, "}\n");
    fclose(fptr);

	//std::cout << std::setprecision(4) << std::fixed;
    //std::cout<<"pnumber tnumber num_reg num_border mesh_sparcing min_radius avg_radius min_angle max_angle min_edge max_edge min_pointDistance min_edgeratio poly_with_be total_be talgorithm funca"<<std::endl;
	//std::cout<<pnumber<<" "<<tnumber<<" "<<num_region;
	//std::cout<<" "<<num_border;
	//printf(" %.12f", mesh_sparcing);
    //std::cout<<" "<<min_radius<<" "<<avg_radius;
    //std::cout<<" "<<min_angle;
    //std::cout<<" "<<max_angle;
    //std::cout<<" "<<min_edge;
    //std::cout<<" "<<max_edge;
    //std::cout<<" "<<min_pointDistance;
    //std::cout<<" "<<min_edgeratio;
    //std::cout<<" "<<est_poly_with_be;
    //std::cout<<" "<<est_total_be; 
    //std::cout<<" "<<t_total; 
    //std::cout<<"\n";
    //std::cout<<"Voronoi_Sites tnumber num_reg num_terminal_edge_regions Sites_per_poly Triangles_per_poly Edge_per_poly max_poly_be Total_BE"<<std::endl;
    std::cout<<pnumber;
    std::cout<<" "<<tnumber;
    std::cout<<" "<<num_region;
    std::cout<<" "<<num_terminal_border_edges + num_terminal_edges/2;
    std::cout<<" "<<t_delaunay;
    std::cout<<" "<<t_total;
    std::cout<<" "<<t_label;
    std::cout<<" "<<t_travel;
    std::cout<<" "<<tcost_be;
    //std::cout<<" "<<(float) pnumber/num_region;
    //std::cout<<" "<<(float) tnumber/num_region;
    //std::cout<<" "<<(float) est_total_edges/num_region;
    //std::cout<<" "<<est_max_triangles_be;
    //std::cout<<" "<<est_total_be;
    std::cout<<"\n";
	//std::cout<<" "<<est_poly_with_be<<" "<<est_total_be<<" "<<est_min_triangles_be<<" "<<est_max_triangles_be;
	//std::cout<<" "<<(est_poly_with_be > 0 ? est_ratio_be/est_poly_with_be : 0.0);
	//std::cout<<" "<<est_total_edges<<" "<<est_max_edges<<" "<<est_min_edges<<" "<<(float)est_total_edges/num_region;
	//std::cout<<" "<<num_terminal_edges/2<<" "<<num_terminal_border_edges<<" "<<num_frontier_edges/2<<" //"<<num_frontier_border_edges<<" "<<num_interior_edges/2;

	//std::cout<<num_terminal_border_edges<<" "<<3*pnumber - 3 - (num_terminal_border_edges + num_frontier_border_edges) <<" = "<<num_terminal_edges/2 + num_terminal_border_edges + num_frontier_edges/2 + num_frontier_border_edges + num_interior_edges/2<<" "<<(3*pnumber - 3 - (num_terminal_border_edges + num_frontier_border_edges) == num_terminal_edges/2 + num_terminal_border_edges + num_frontier_edges/2 + num_frontier_border_edges + num_interior_edges/2)<<std::endl;
     //std::cout<<"Metric saved in: "<<cmd<<std::endl;
}

/*geomview output*/
void write_VEM(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, int print_triangles){

    int i,j;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_vertice.txt");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr,"%d %d\n", pnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.5f %.5f\n", r[2*i + 0], r[2*i + 1]);

  //imprimir polginos
    i = 0;
    while(i < i_mesh){
        int length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i] + 1);
            i++;
        }
        fprintf(fptr, "\n");
    }


    fprintf(fptr,"\n");
    fclose(fptr);
}

/*geomview output*/
void write_VEM_triangles(std::string name, double *r, int *triangles, int *adj, int pnumber, int tnumber, int i_mesh, int *mesh, int *seed, int num_region, std::list <int> &seed_bet){

    int i,j, lenght_poly;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_triangulos.txt");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    
    //imprimir puntos
    fprintf(fptr,"%d %d %d 0\n", pnumber, tnumber, num_region);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.5f %.5f\n", r[2*i + 0], r[2*i + 1]);

        //imprimir triangulos
    for(i = 0; i < tnumber; i++)
        fprintf(fptr,"%d %d %d\n", triangles[3*i + 0] +1, triangles[3*i + 1] +1, triangles[3*i + 2] +1);

  //imprimir polginos
    std::set<int> s;
    j = 1;
    for( i = 0; i < tnumber; i++)
	{
		if(seed[i] == TRUE ){
			lenght_poly = look_triangles(i, s, triangles, adj, r);
            fprintf(fptr, "%d %ld ", j, s.size());
            
            for ( auto it = s.begin(); it != s.end(); it++ )
                fprintf(fptr, "%d ", *it + 1);
            fprintf(fptr, "\n");
            j++;
            s.clear();
        }
    }
    std::list <int> :: iterator ite;
    //std::cout<<"iterando "<<seed_bet.size()<<" de "<<j<<std::endl;
    for(ite = seed_bet.begin(); ite != seed_bet.end(); ++ite){
            //std::cout<<"seed: "<<*ite<<std::endl;
			lenght_poly = look_triangles(*ite, s, triangles, adj, r);
            fprintf(fptr, "%d %ld ", j, s.size());
            
            for ( auto it = s.begin(); it != s.end(); it++ )
                fprintf(fptr, "%d ", *it + 1);
            fprintf(fptr, "\n");
            j++;
            s.clear();
    }

    fprintf(fptr,"\n");
    fclose(fptr);
}

/*geomview output*/
void write_triangulation(std::string name, double *r, int *triangles, int *adj, int pnumber, int tnumber){
    int i,j, lenght_poly;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_tri.off");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    fprintf(fptr, "{ appearance  {+edge +face linewidth 2} LIST\n");

    fprintf(fptr, "OFF\n");
    fprintf(fptr,"%d %d 0\n", pnumber, tnumber);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.5f %.5f 0\n", r[2*i + 0], r[2*i + 1]);

        //imprimir triangulos
    for(i = 0; i < tnumber; i++)
        fprintf(fptr,"3 %d %d %d\n", triangles[3*i + 0], triangles[3*i + 1], triangles[3*i + 2]);
    fprintf(fptr," }\n");
    fclose(fptr);
}

void write_GID(std::string name, double *r, int *triangles, int *adj, int pnumber, int tnumber){
    int i,j, lenght_poly;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_gid.msh");
    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    fprintf(fptr, "MESH    dimension 2 ElemType Triangle  Nnode 3\n");

    fprintf(fptr, "Coordinates\n");
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%d %.16f %.16f\n", i + 1, r[2*i + 0], r[2*i + 1]);
     fprintf(fptr, "end coordinates\n      \nElements\n");
    for(i = 0; i < tnumber; i++)
        fprintf(fptr,"%d %d %d %d 1\n",i+1, triangles[3*i + 0]+1, triangles[3*i + 1]+1, triangles[3*i + 2]+1);
    fprintf(fptr,"end elements\n");
    fclose(fptr);
}

void write_alejandro_custom(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region, int *border, int num_boder)
{
    int i;
    double xmax,ymax, xmin,ymin;
    xmax = r[0];
    xmin = r[0];
    ymax = r[1];
    ymin = r[1];
    std::list<int> bottom;
    std::list<int> top;
    std::list<int> left;
    std::list<int> right;

    //busca máximos y min
    for(i = 0; i < pnumber; i++){
        //search range x
        if(r[2*i + 0] > xmax )
            xmax = r[2*i + 0];
        if(r[2*i + 0] < xmin )
            xmin = r[2*i + 0];
        //search range y
        if(r[2*i + 1] > ymax )
            ymax = r[2*i + 1];
        if(r[2*i + 1] < ymin )
            ymin = r[2*i + 1];
    }

    //agrega elementos de borde
    for(i = 0; i < pnumber; i++){
        if(r[2*i + 1] == ymin)
            bottom.push_back(i);
        if(r[2*i + 1] == ymax)
            top.push_back(i);
        if(r[2*i + 0] == xmin)
            left.push_back(i);
        if(r[2*i + 0] == xmax)
            right.push_back(i);
    }

    
    int j,length_poly = 0;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_alejandro.txt");

    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    fprintf(fptr, "# domain type\n");
    fprintf(fptr, "Custom\n");
    fprintf(fptr, "# nodal coordinates: number of nodes followed by the coordinates \n");
    fprintf(fptr,"%d\n", pnumber);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.16f %.16f\n", r[2*i + 0], r[2*i + 1]);
    fprintf(fptr, "# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel) \n");
    //imprimir polginos
    fprintf(fptr,"%d\n", num_region);
    i = 0;
    while(i < i_mesh){
        length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i]+1);
            i++;
        }
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "# indices of nodes located on the Dirichlet boundary\n");
    for (i = 0; i< num_boder; i++)
        fprintf(fptr, "%d ", border[i] + 1);
    fprintf(fptr, "\n");
    fprintf(fptr, "# indices of nodes located on the Neumann boundary\n0\n");
    fprintf(fptr, "# xmin, xmax, ymin, ymax of the bounding box\n");
    fprintf(fptr, "%.16f %.16f %.16f %.16f\n", xmin, xmax, ymin, ymax);
    fclose(fptr);
}

void write_alejandro(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region){
    int i;
    double xmax,ymax, xmin,ymin;
    xmax = r[0];
    xmin = r[0];
    ymax = r[1];
    ymin = r[1];
    std::list<int> bottom;
    std::list<int> top;
    std::list<int> left;
    std::list<int> right;

    //busca máximos y min
    for(i = 0; i < pnumber; i++){
        //search range x
        if(r[2*i + 0] > xmax )
            xmax = r[2*i + 0];
        if(r[2*i + 0] < xmin )
            xmin = r[2*i + 0];
        //search range y
        if(r[2*i + 1] > ymax )
            ymax = r[2*i + 1];
        if(r[2*i + 1] < ymin )
            ymin = r[2*i + 1];
    }

    //agrega elementos de borde
    for(i = 0; i < pnumber; i++){
        if(r[2*i + 1] == ymin)
            bottom.push_back(i);
        if(r[2*i + 1] == ymax)
            top.push_back(i);
        if(r[2*i + 0] == xmin)
            left.push_back(i);
        if(r[2*i + 0] == xmax)
            right.push_back(i);
    }

    //ordernar
    //bottom menor a mayor
    bottom.sort([&r](int a, int b) {return r[2*a + 0] < r[2*b + 0]; });
    //top mayor a menor
    top.sort([&r](int a, int b) {return r[2*a + 0] > r[2*b + 0]; });
    ////left mayor a menor
    left.sort([&r](int a, int b) {return r[2*a + 1] > r[2*b + 1]; });
    ////right menor a mayor
    right.sort([&r](int a, int b) {return r[2*a + 1] < r[2*b + 1]; });
    
    int j,length_poly = 0;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_alejandro_");
    strcat(cmd, std::to_string(num_region).c_str());
    strcat(cmd, "_polygons.txt");

    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    fprintf(fptr, "# domain type\n");
    fprintf(fptr, "RectangularDomain\n");
    fprintf(fptr, "# nodal coordinates: number of nodes followed by the coordinates \n");
    fprintf(fptr,"%d\n", pnumber);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.16f %.16f\n", r[2*i + 0], r[2*i + 1]);
    fprintf(fptr, "# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel) \n");
    //imprimir polginos
    fprintf(fptr,"%d\n", num_region);
    i = 0;
    while(i < i_mesh){
        length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i]+1);
            i++;
        }
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "# indices of nodes located on the bottom boundary\n");
    for (auto v : bottom)
        fprintf(fptr, "%d ", v+1);

    std::cout<<"\n";
    fprintf(fptr, "\n");

    fprintf(fptr, "# indices of nodes located on the top boundary\n");
    for (auto v : top)
        fprintf(fptr, "%d ", v+1);
    fprintf(fptr, "\n");    

    fprintf(fptr, "# indices of nodes located on the left boundary\n");
    for (auto v : left)
        fprintf(fptr, "%d ", v+1);

    fprintf(fptr, "\n");    
    fprintf(fptr, "# indices of nodes located on the right boundary\n");
    for (auto v : right)
        fprintf(fptr, "%d ", v+1);
    fprintf(fptr, "\n");    

    fprintf(fptr, "# xmin, xmax, ymin, ymax of the bounding box\n");
    fprintf(fptr, "%.16f %.16f %.16f %.16f\n", xmin, xmax, ymin, ymax);
    fclose(fptr);
}

void write_alejandro_quater_circle(std::string name, double *r, int *triangles, int pnumber, int tnumber, int i_mesh, int *mesh, int num_region){
    int i;
    double xmax,ymax, xmin,ymin;
    xmax = r[0];
    xmin = r[0];
    ymax = r[1];
    ymin = r[1];
    std::list<int> bottom;
    std::list<int> top;
    std::list<int> left;
    std::list<int> right;

    //busca máximos y min
    for(i = 0; i < pnumber; i++){
        //search range x
        if(r[2*i + 0] > xmax )
            xmax = r[2*i + 0];
        if(r[2*i + 0] < xmin )
            xmin = r[2*i + 0];
        //search range y
        if(r[2*i + 1] > ymax )
            ymax = r[2*i + 1];
        if(r[2*i + 1] < ymin )
            ymin = r[2*i + 1];
    }

    //agrega elementos de borde
    for(i = 0; i < pnumber; i++){
        if(r[2*i + 1] == ymin)
            bottom.push_back(i);
        if(r[2*i + 1] == ymax)
            top.push_back(i);
        if(r[2*i + 0] == xmin)
            left.push_back(i);
        if(r[2*i + 0] == xmax)
            right.push_back(i);
    }

    //ordernar
    //bottom menor a mayor
    bottom.sort([&r](int a, int b) {return r[2*a + 0] < r[2*b + 0]; });
    //top mayor a menor
    top.sort([&r](int a, int b) {return r[2*a + 0] > r[2*b + 0]; });
    ////left mayor a menor
    left.sort([&r](int a, int b) {return r[2*a + 1] > r[2*b + 1]; });
    ////right menor a mayor
    right.sort([&r](int a, int b) {return r[2*a + 1] < r[2*b + 1]; });
    
    int j,length_poly = 0;
    char cmd[1024] = "\0";
    strcat(cmd, filespathoutput);
    strcat(cmd, name.c_str());
    //strcat(cmd, ppath);
    strcat(cmd,"_quarter_circle_alejandro.txt");

    FILE *fptr;
    fptr = fopen(cmd, "w");
    if(fptr == NULL){
        printf("Hmnito, no abrio el archivo :C\n");
        perror("fopen");
        exit(0);
    }   
    fprintf(fptr, "# domain type\n");
    fprintf(fptr, "PlateWithHoleQuarterDormain\n");
    fprintf(fptr, "# nodal coordinates: number of nodes followed by the coordinates \n");
    fprintf(fptr,"%d\n", pnumber);
    for(i = 0; i < pnumber; i++)
        fprintf(fptr,"%.16f %.16f\n", r[2*i + 0], r[2*i + 1]);
    fprintf(fptr, "# element connectivity: number of elements followed by the connectivity (each line: nodes_per_element(nel) node1 node 2 ... node_nel) \n");
    //imprimir polginos
    fprintf(fptr,"%d\n", num_region);
    i = 0;
    while(i < i_mesh){
        length_poly = mesh[i];
        i++;
        fprintf(fptr, "%d ", length_poly);
        for(j=0; j < length_poly;j++){
            fprintf(fptr, "%d ", mesh[i]+1);
            i++;
        }
        fprintf(fptr, "\n");
    }

    fprintf(fptr, "# indices of nodes located on the Dirichlet left boundary\n");
    for (auto v : left)
        fprintf(fptr, "%d ", v+1);
    fprintf(fptr, "\n");

    fprintf(fptr, "# indices of nodes located on the Dirichlet bottom boundary\n");
    for (auto v : bottom)
        fprintf(fptr, "%d ", v+1);
    fprintf(fptr, "\n");

    fprintf(fptr, "# indices of nodes located on the Neumann top boundary\n");
    for (auto v : top)
        fprintf(fptr, "%d ", v+1);
    fprintf(fptr, "\n"); 

    fprintf(fptr, "# indices of nodes located on the Neumann right boundary\n");
    for (auto v : right)
        fprintf(fptr, "%d ", v+1);
    fprintf(fptr, "\n");    


    fprintf(fptr, "# xmin, xmax, ymin, ymax of the bounding box\n");
    fprintf(fptr, "%.16f %.16f %.16f %.16f\n", xmin, xmax, ymin, ymax);
    fclose(fptr);
}

int look_triangles(int i, std::set<int> &s, int * triangles, int * adj, double *r) {


    s.insert(i);

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
        //poly[ind_poly] = triangles[3 * i + 0];
        ind_poly++;
        ////poly[ind_poly] = triangles[3 * i + 1];
        ind_poly++;
        ////poly[ind_poly] = triangles[3 * i + 2];
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
                ////poly[ind_poly] = triangles[ind1];
                ind_poly++;
                ////poly[ind_poly] = triangles[ind2];
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
                ////poly[ind_poly] = triangles[3*i + (j+1)%3];
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
        s.insert(k);

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
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    ////poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];
                }
            } else if (t2 == NO_ADJ && t0 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2-t0 son fe*/
                if (continuous == 0) {
                    ////poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ && t2 == NO_ADJ) {
                /*si endpoint es continua a t1 y t1-t2 son fe*/
                if (continuous == 2) {
                    ////poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;
                    ////poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    ////poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;
                    //poly[ind_poly] = triangles[3 * k + 0];
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
                    //poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 2];

                } else if (continuous == 2) {
                    //poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                }
            } else if (t2 == NO_ADJ) {
                /*si endpoint es continua a t2  y t2 es fe*/
                if (continuous == 0) {
                    //poly[ind_poly] = triangles[3 * k + 0];
                    ind_poly++;

                    end_point = triangles[3 * k + 1];

                } else if (continuous == 1) {
                    //poly[ind_poly] = triangles[3 * k + 1];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                }
            } else if (t1 == NO_ADJ) {
                /*si endpoint es continua a t1  y t1 es fe*/
                if (continuous == 2) {
                    //poly[ind_poly] = triangles[3 * k + 2];
                    ind_poly++;

                    end_point = triangles[3 * k + 0];

                } else if (continuous == 0) {
                    //poly[ind_poly] = triangles[3 * k + 0];
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