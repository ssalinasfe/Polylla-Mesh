// reading a text file
#include <iostream>
#include <fstream>
#include <string>

void read_from_triangle(std::string name, int &pnumber, int &tnumber, double *&points, int *&triangles, int *&neigh){
    std::string line;
    std::ifstream nodefile(name + ".node");
    double a1, a2, a3, a4;
    int i = 0;
    
    //std::cout<<"Node file"<<std::endl;
    if (nodefile.is_open())
    {
        nodefile >> pnumber ;
        //std::cout<<pnumber<<std::endl;

        std::getline(nodefile, line); 
        points = (double *)malloc(2*pnumber*sizeof(double));
        while (nodefile >> a1 >> a2 >> a3 >> a4)
        {
            points[2*i + 0] = a2;
            points[2*i + 1] = a3;
            //std::cout<<points[2*i + 0]<<" "<<points[2*i + 1]<<std::endl;
            i++;
            //std::cout<<a2<<" "<<a3<<std::endl;
            
        }
        
    }
    else 
        std::cout << "Unable to open node file"; 

    nodefile.close();


    //std::cout<<"Ele file"<<std::endl;
    std::ifstream elefile(name + ".ele");
    int t1, t2, t3, t4;
    i = 0;
    if(elefile.is_open()){
        elefile >> tnumber ;
        triangles = (int *)malloc(3*tnumber*sizeof(int));
        std::getline(elefile, line); 
        while (elefile >> t1 >> t2 >> t3 >> t4 )
        {
            //std::cout<<t2<<" "<<t3<<" "<<t4<<std::endl;
            triangles[3*i + 0] = t2;
            triangles[3*i + 1] = t3;
            triangles[3*i + 2] = t4;
            //std::cout<<triangles[3*i + 0]<<" "<<triangles[3*i + 1]<<" "<<triangles[3*i + 2]<<std::endl;
            i++;
        }
    }else std::cout << "Unable to open ele file";

    elefile.close();

    //std::cout<<"Neigh file"<<std::endl;
    std::ifstream neighfile(name + ".neigh");
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
}

__global__ void initialize_memory(int* cu_trivertex, int* cu_triangles, int pnumber, int tnumber){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int j;
    if(i < pnumber){  
        for (j = 0; j < tnumber; j++){
            if(i == cu_triangles[3*j + 0] ||  i == cu_triangles[3*j + 1] || i == cu_triangles[3*j + 2]){
                cu_trivertex[i] = j;
                break;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    double *points;
    int *triangles;
    int *neigh;
    int *trivertex;
    int *cu_triangles;
    int *cu_trivertex;

    int pnumber, tnumber;
    std::string name(argv[1]);
	std::cout<<name<<std::endl;
	read_from_triangle(name, pnumber, tnumber, points, triangles, neigh);
	std::cout << " " << tnumber << " " << pnumber << "\n";

    trivertex = (int *)malloc(pnumber*sizeof(int));
	cudaMalloc((void**) &cu_triangles, 3*tnumber*sizeof(int));
	cudaMemcpy(cu_triangles, triangles, 3*tnumber*sizeof(int), cudaMemcpyHostToDevice);
	cudaMalloc((void**) &cu_trivertex, pnumber*sizeof(int));
	int numThreads = 128;  // max register per block is 65536, 65536/512
    std::cout<<"generating trivertex"<<std::endl;
	initialize_memory<<<(pnumber + (numThreads-1))/numThreads, numThreads>>>( cu_trivertex, cu_triangles, pnumber, tnumber);
    
	cudaMemcpy(trivertex, cu_trivertex, pnumber*sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(cu_trivertex);
	cudaFree(cu_triangles);
    
    std::cout<<"storing .trivertex"<<std::endl;
    std::ofstream myfile;
    myfile.open (name + ".trivertex");
    myfile << pnumber<<"\n";
    for (size_t i = 0; i < pnumber; i++)
    {
        myfile << i<<" "<<trivertex[i]<<"\n";
    }
    
    myfile.close();

    return 0;
}
