#include <iostream>
#include "/Users/sergeypanov/bin/mpi/include/mpi.h"
#include <fstream>
#include <vector>

using namespace std;

typedef struct Me {
    int my_id;
    char from_node;
    char to_node;
    bool is_back;
} Me;


vector<Me> constructEdges(string tree){

    int nodes = tree.size();

    vector<Me> edges;

    int procid = 0;

    tree = " " + tree;

    for (int i = 1; i < tree.size(); ++i) {
        char parent_node = tree[i];

        char left_son;
        char right_son;


        left_son = tree[2*i];
        right_son = tree[2*i + 1];

        Me forward_left_edge;
        forward_left_edge.my_id = procid;
        forward_left_edge.from_node = parent_node;
        forward_left_edge.to_node = left_son;
        forward_left_edge.is_back = false;

        if (2*nodes-2 > edges.size())
            edges.push_back(forward_left_edge);

        ++procid;

        Me back_left_edge;
        back_left_edge.my_id = procid;
        back_left_edge.from_node = left_son;
        back_left_edge.to_node = parent_node;
        back_left_edge.is_back = true;

        if (2*nodes-2 > edges.size())
            edges.push_back(back_left_edge);

        ++procid;

        Me forward_right_edge;
        forward_right_edge.my_id = procid;
        forward_right_edge.from_node = parent_node;
        forward_right_edge.to_node = right_son;
        forward_right_edge.is_back = false;

        if (2*nodes-2 > edges.size())
            edges.push_back(forward_right_edge);

        ++procid;

        Me back_right_edge;
        back_right_edge.my_id = procid;
        back_right_edge.from_node = right_son;
        back_right_edge.to_node = parent_node;
        back_right_edge.is_back = true;

        if (2*nodes-2 > edges.size())
            edges.push_back(back_right_edge);

        ++procid;

    }
    return edges;
}
int main(int argc, char** argv) {

    int numprocs;
    int myid;


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0){
        // Read tree
        string file = "tree";
        ifstream infile;
        string input;
        infile.open(file);
        getline(infile, input);
        infile.close();


        vector<Me> edges = constructEdges(input);

        for (int i = 0; i < edges.size(); ++i) {
            cout << "id: " << edges[i].my_id << endl;
            cout << "from: " << (char) edges[i].from_node << endl;
            cout << "to: " << (char) edges[i].to_node << endl;
            cout << "is_back: " << edges[i].is_back << endl;
            cout << "------------" << endl;
        }
    }


    MPI_Finalize();
    return 0;
}