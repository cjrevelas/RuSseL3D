#include <iostream>
#include "/home/cjrevelas/temp/include/fem.h"
#include "/home/cjrevelas/temp/include/mesh.h"

int main(int argc, char* argv[]){

    int gid = atoi(argv[1]);

    Mesh mesh(268491, 48400);
    Mesh* ptr1 = &mesh;

    ptr1->elementsContainingNode(gid);
    std::cout << ptr1->computeMeshVolume() << std::endl;

    return 0;
}
