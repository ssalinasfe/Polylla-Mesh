#include "polygon.h"
#include <iostream>

int count_regions(int *mesh, int i_mesh){
    int i = 0;
	int num_region = 0;
    while(i < i_mesh){
		num_region++;
		i = mesh[i] + i + 1;
    }
    return num_region;
}

int save_to_mesh(int *mesh, int *poly, int i_mesh, int length_poly, double *r){    

    mesh[i_mesh] = length_poly;
    //std::cout<<get_signed_area_poly(poly,length_poly,r)<<std::endl;
    if(get_signed_area_poly(poly,length_poly,r) < 0.0) {
        reverse(poly, length_poly);
        std::cout<<"hay que voltear esta"<<std::endl;
    }

    i_mesh++;
    for(int i = 0; i <length_poly; i++){
        mesh[i_mesh + i] = poly[i];
    }
    
    return i_mesh + length_poly;
}
