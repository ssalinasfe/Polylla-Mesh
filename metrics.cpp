#include <iostream>
#include <math.h>  
#include <assert.h>     /* assert */
#include <float.h>

#include "SmallestEnclosingCircle.hpp"
#include "metrics.h"
#include "polygon.h"
#include "triang.h"



/*
int main(int argc, char* argv[]){


    //Min angle tests
    //test puntos colineales
    double r_col[] = {0,0,0.5,0,1,0};
    int col[] = {0,1,2};
    assert( get_min_angle_non_zero(col, 3, r_col ) == 360.0);

    //test square
    double r_square[] = {0,0,0.5,0,1,0,1,1,0,1};
    int square[] = {0,1,2,3,4};
    assert( get_min_angle_non_zero(square, 5, r_square) - 90 < 0.001);
    //    std::cout<<"min angle: "<<get_min_angle_non_zero(square, len_poly, r_square)<<std::endl;

    //test reflex angle
    double r_L[] = {-1,-1,0,-1,0,0,1,0,1,1,0,1,-1,1,-1,0};
    int L[] = {0,1,2,3,4,5,6,7};
    assert(std::fabs(get_min_angle_non_zero(L, 8, r_L) - 90) < 0.001);

    //angle reflex test
    assert(std::fabs(get_angle_three_points(r_L, 1,2,3) - 270) < 0.001);

    //std::cout<<get_max_edge(square,5,r_square)<<std::endl;
    //std::cout<<get_min_edge(square,5,r_square)<<std::endl;

    double metrics[10];

    get_metrics_of_polygon(metrics,square, 5,r_square);

    std::cout<<"\nsquare\narea: "<<metrics[0]<<std::endl;
    std::cout<<"perimeter: "<<metrics[1]<<std::endl;
    std::cout<<"area_perimeter_ratio: "<<metrics[2]<<std::endl;
    std::cout<<"min_angle: "<<metrics[3]<<std::endl;
    std::cout<<"max_angle: "<<metrics[4]<<std::endl;
    std::cout<<"min_edge: "<<metrics[5]<<std::endl;
    std::cout<<"max_edge: "<<metrics[6]<<std::endl;
    std::cout<<"edge_rati: "<<metrics[7]<<std::endl;
    std::cout<<"min_pd: "<<metrics[8]<<std::endl;
    std::cout<<"max_pd: "<<metrics[9]<<std::endl;
    
    get_metrics_of_polygon(metrics,L, 8,r_L);
    std::cout<<"\nL\narea: "<<metrics[0]<<std::endl;
    std::cout<<"perimeter: "<<metrics[1]<<std::endl;
    std::cout<<"area_perimeter_ratio: "<<metrics[2]<<std::endl;
    std::cout<<"min_angle: "<<metrics[3]<<std::endl;
    std::cout<<"max_angle: "<<metrics[4]<<std::endl;
    std::cout<<"min_edge: "<<metrics[5]<<std::endl;
    std::cout<<"max_edge: "<<metrics[6]<<std::endl;
    std::cout<<"edge_rati: "<<metrics[7]<<std::endl;
    std::cout<<"min_pd: "<<metrics[8]<<std::endl;
    std::cout<<"max_pd: "<<metrics[9]<<std::endl;
    
    get_metrics_of_polygon(metrics,col, 3,r_col);
    std::cout<<"\ncolinear\narea: "<<metrics[0]<<std::endl;
    std::cout<<"perimeter: "<<metrics[1]<<std::endl;
    std::cout<<"area_perimeter_ratio: "<<metrics[2]<<std::endl;
    std::cout<<"min_angle: "<<metrics[3]<<std::endl;
    std::cout<<"max_angle: "<<metrics[4]<<std::endl;
    std::cout<<"min_edge: "<<metrics[5]<<std::endl;
    std::cout<<"max_edge: "<<metrics[6]<<std::endl;
    std::cout<<"edge_rati: "<<metrics[7]<<std::endl;
    std::cout<<"min_pd: "<<metrics[8]<<std::endl;
    std::cout<<"max_pd: "<<metrics[9]<<std::endl;
    return 0;
}
*/


void get_metrics_of_polygon(double *metrics, int *poly, int len_poly, double *r){
    double ax,ay,bx,by;
    int j,k, a,b,c;
    double d;
    double max_edge, min_edge, edge_rati;
    double perimeter = 0.0, area, area_perimeter_ratio;
    double  min_angle = 360.0, max_angle = 0, angle;

    area = get_signed_area_poly(poly, len_poly, r);
    std::vector<Point> circumscribed_circle_points;

    for (size_t i = 0; i < len_poly; i++)
    {
        j = (i + 1)%len_poly;
        k = (i + 2)%len_poly;
        a = poly[i];
        b = poly[j];
        c = poly[k];

        ax = r[2*a+0];
        ay = r[2*a+1];

        bx = r[2*b+0];
        by = r[2*b+1]; 

        circumscribed_circle_points.push_back(Point{ax,ay});
        d = dist(ax,ay,bx,by);
        perimeter += d;

        if(i == 0 || d > max_edge)
            max_edge = d;
        if(i == 0 || d < min_edge)
            min_edge = d;

        angle = get_angle_three_points(r,a,b,c);

        if( angle != 0 && angle < min_angle )
            min_angle = angle;
        if(angle > max_angle )
            max_angle = angle;
    }

    double min_pd = DBL_MAX, max_pd = 0;

    for (size_t i = 0; i < len_poly; i++)
    {
        for (size_t j = 0; j < len_poly; j++)
        {
            if(i!=j){
                a = poly[i];
                b = poly[j];
                ax = r[2*a+0];
                ay = r[2*a+1];

                bx = r[2*b+0];
                by = r[2*b+1];

                d = dist(ax,ay,bx,by);

                if(d > max_pd)
                    max_pd = d;
                if(d < min_pd)
                    min_pd = d;
            }
        }
        
    }


    edge_rati = min_edge/max_edge;
    area_perimeter_ratio =  (2*M_PI*area)/pow(perimeter,2);
    Circle cir = makeSmallestEnclosingCircle(circumscribed_circle_points);

    metrics[0] = area;
    metrics[1] = perimeter;
    metrics[2] = area_perimeter_ratio;
    metrics[3] = min_angle;
    metrics[4] = max_angle;
    metrics[5] = min_edge;
    metrics[6] = max_edge;
    metrics[7] = edge_rati;
    metrics[8] = min_pd;
    metrics[9] = max_pd;
    metrics[10] = cir.r;

}


double get_angle_three_points(double *r, int i, int j, int k){
    double ax,ay,bx,by,cx,cy;    
    ax = r[2*i+0];
    ay = r[2*i+1];

    bx = r[2*j+0];
    by = r[2*j+1];

    cx = r[2*k+0];
    cy = r[2*k+1];

    //int a = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    if(ax * (by - cy) + bx * (cy - ay) + cx * (ay - by) == 0)
        return 0.0;


    double angba = atan2((bx - ax), (by - ay));
    double angbc = atan2((bx - cx), (by - cy));
    double rslt =  angbc - angba;
    double rs = (rslt * 180) / 3.141592;
    rs = rs < 0 ? rs + 360.0 : rs;
    return rs;
}
/*
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
*/
double get_max_edge(int *poly, int len_poly, double *r){
    double ax,ay,bx,by;
    int j;
    double d,max = 0;
    for (size_t i = 0; i < len_poly; i++)
    {
        j = (i + 1)%len_poly;

        ax = r[2*i+0];
        ay = r[2*i+1];

        bx = r[2*j+0];
        by = r[2*j+1]; 

        d = dist(ax,ay,bx,by);
        //std::cout<<i<<"-"<<j<<": "<<d<<std::endl;
        if(d > max)
            max = d;
    }
    return max;
}

double get_min_edge(int *poly, int len_poly, double *r){
    double ax,ay,bx,by;
    int j, a,b;
    double d,min;
    for (size_t i = 0; i < len_poly; i++)
    {
        j = (i + 1)%len_poly;

        a = poly[i];
        b = poly[j];

        ax = r[2*a+0];
        ay = r[2*a+1];

        bx = r[2*b+0];
        by = r[2*b+1]; 

        d = dist(ax,ay,bx,by);
        //std::cout<<i<<"-"<<j<<": "<<d<<std::endl;
        if(i == 0 || d < min)
            min = d;
    }
    return min;
}

double get_min_angle_non_zero(int *poly, int len_poly, double *r){
    int j,k,a,b,c;
    double min = 360.0;
    double angle;
    for (size_t i = 0; i < len_poly; i++)
    {   
        
        j = (i + 1)%len_poly;
        k = (i + 2)%len_poly;
        
        a = poly[i];
        b = poly[j];
        c = poly[k];

        angle = get_angle_three_points(r,a,b,c);
        if( angle != 0 && angle < min )
            min = angle;
    }
    return min;
}




/*
int Equality(double a, double b, double epsilon)
{
  return fabs(a - b) < epsilon;
}


int GreaterEqualthan(double a, double b, double epsilon){
	return Equality(a,b,epsilon) || a > b;
}

double dist(double x0, double y0, double x1, double y1)
{
	return sqrt(pow(x0 - x1, 2.0) + pow(y0 - y1, 2.0));
}
*/