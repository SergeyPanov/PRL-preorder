#include <iostream>
#include "/Users/sergeypanov/bin/mpi/include/mpi.h"
#include <fstream>
#include <vector>
#include <map>

using namespace std;

typedef struct Edge {
    int my_id;
    char from_node;
    char to_node;
    bool is_back;
} Edge;

void display_edge(pair< Edge, Edge > p){
    cout << "first: "<< endl;
    cout << p.first.my_id << " " << p.first.from_node << " " << p.first.to_node << " " << p.first.is_back << endl;

    cout << "second: "<< endl;
    cout << p.second.my_id << " " << p.second.from_node << " " << p.second.to_node << " " << p.second.is_back << endl;
}

Edge make_edge(int id, char from_node, char to_node, bool is_back){
    Edge edge;
    edge.my_id = id;
    edge.from_node = from_node;
    edge.to_node = to_node;
    edge.is_back = is_back;

    return edge;
}

pair< Edge, Edge > make_pair(Edge from, Edge to){
    pair< Edge, Edge > pr;
    pr.first = from;
    pr.second = to;
    return pr;
};

void show_map(map< char, vector< pair< Edge, Edge > > > adjacencies){

    for ( auto it = adjacencies.begin(); it != adjacencies.end(); ++it  ) {
        std::cout << it->first <<  std::endl;

        for (int i = 0; i < it->second.size(); ++i) {
            display_edge(it->second[i]);
        }
    }

}

vector<Edge> constructEdges(string tree){

    int nodes = tree.size();

    vector<Edge> edges;

    int procid = 0;

    map< char, vector< pair<Edge, Edge> > >  adjacency_list;

    tree = " " + tree;

    for (int i = 1; i < tree.size(); ++i) {

        vector< pair< Edge, Edge > > list_for_parent;
        vector< pair< Edge, Edge > > list_for_left;
        vector< pair< Edge, Edge > > list_for_right;


        char parent_node = tree[i];

        char left_son;
        char right_son;


        left_son = tree[2*i];
        right_son = tree[2*i + 1];

        Edge forward_left_edge = make_edge(procid, parent_node, left_son, false);

        if (2*nodes-2 > edges.size())
            edges.push_back(forward_left_edge);

        ++procid;

        Edge back_left_edge = make_edge(procid, left_son, parent_node, true);


        if (2*nodes-2 > edges.size())
            edges.push_back(back_left_edge);

        ++procid;

        pair< Edge, Edge > left_pair = make_pair(forward_left_edge, back_left_edge);
        list_for_parent.push_back(left_pair);

        // Insert pair for the left son
        pair< Edge, Edge > left_sons_pair = make_pair(back_left_edge, forward_left_edge);
        list_for_left.push_back(left_sons_pair);


        Edge forward_right_edge = make_edge(procid, parent_node, right_son, false);

        if (2*nodes-2 > edges.size())
            edges.push_back(forward_right_edge);

        ++procid;

        Edge back_right_edge = make_edge(procid, right_son, parent_node, true);


        if (2*nodes-2 > edges.size())
            edges.push_back(back_right_edge);


        pair< Edge, Edge > right_pair = make_pair(forward_right_edge, back_right_edge);
        list_for_parent.push_back(right_pair);

        // Insert pair for the right son
        pair< Edge, Edge > right_sons_pair = make_pair(back_right_edge, forward_right_edge);
        list_for_right.push_back(right_sons_pair);



        map< char, vector< pair<Edge, Edge> > >::iterator it = adjacency_list.find(parent_node);

        if (it != adjacency_list.end()){

            for (int j = 0; j < list_for_parent.size(); ++j) {

                if (list_for_parent[j].first.from_node != 0 && list_for_parent[j].first.to_node != 0)
                    it->second.push_back(list_for_parent[j]);
            }


        } else{
            // Insert parent
            adjacency_list.insert(pair< char, vector< pair<Edge, Edge> > >(parent_node, list_for_parent));

            // Insert left son
            adjacency_list.insert(pair< char, vector< pair< Edge, Edge > > >(left_son, list_for_left));

            // Insert right son
            adjacency_list.insert(pair< char, vector< pair< Edge, Edge > > >(right_son, list_for_right));

        }

        ++procid;
    }

    show_map(adjacency_list);
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


        vector<Edge> edges = constructEdges(input);

//        for (int i = 0; i < edges.size(); ++i) {
//            cout << "id: " << edges[i].my_id << endl;
//            cout << "from: " << (char) edges[i].from_node << endl;
//            cout << "to: " << (char) edges[i].to_node << endl;
//            cout << "is_back: " << edges[i].is_back << endl;
//            cout << "------------" << endl;
//        }
    }


    MPI_Finalize();
    return 0;
}