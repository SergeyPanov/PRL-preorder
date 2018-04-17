#include <iostream>
#include "/Users/sergeypanov/bin/mpi/include/mpi.h"
#include <fstream>
#include <vector>
#include <map>

#define TAG 0

using namespace std;

typedef struct Edge {
    int my_id;
    char from_node;
    char to_node;
    bool is_back;
} Edge;

void display_one_direction(Edge e){
    cout << e.my_id << " " << e.from_node << " " << e.to_node << " " << e.is_back << endl;
}

void display_edge(pair< Edge, Edge > p){
    cout << "first: "<< endl;
    display_one_direction(p.first);

    cout << "second: "<< endl;
    display_one_direction(p.second);
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
}

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

    for (int i = 1; i <= tree.size(); ++i) {

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


        Edge forward_right_edge = make_edge(procid, parent_node, right_son, false);

        if (2*nodes-2 > edges.size())
            edges.push_back(forward_right_edge);

        ++procid;

        Edge back_right_edge = make_edge(procid, right_son, parent_node, true);

        if (2*nodes-2 > edges.size())
            edges.push_back(back_right_edge);

        ++procid;
    }
    return edges;
}


pair< Edge, Edge > find_pair(vector< Edge > edges, char from, char to){

    Edge forward_edge = {-1};
    Edge back_edge = {-1};

    for (int i = 0; i < edges.size(); ++i) {

        if (edges[i].from_node == from && edges[i].to_node == to){
            forward_edge = edges[i];
        }


        if (edges[i].from_node == to && edges[i].to_node == from){
            back_edge = edges[i];
        }
    }

    return pair< Edge, Edge >(forward_edge, back_edge);

}

pair< Edge, Edge > find_parent(vector< Edge > edges, char node){
    Edge parent_forward = {-1};
    Edge parent_back = {-1};

    for (int i = 0; i < edges.size(); ++i) {

        if (edges[i].is_back && edges[i].from_node == node){
            parent_back = edges[i];
        }

        if (!edges[i].is_back && edges[i].to_node == node){
            parent_forward = edges[i];
        }
    }

    return pair< Edge, Edge > (parent_back, parent_forward);
}

void display_etour(vector< pair<Edge, Edge> > tour){

    for (int i = 0; i < tour.size(); ++i) {
        cout << "After: ";
        display_one_direction(tour[i].first);
        cout << "Is: ";
        display_one_direction(tour[i].second);
    }


}


map< char, vector< pair< Edge, Edge > > > construct_adacency_list(const vector< Edge > &edges, string nodes){

    int n = static_cast<int>(nodes.size());
    nodes = " " + nodes;


    map< char, vector< pair< Edge, Edge > > > adjacency_list;

    int inned_iterator = 0;

    for (int i = 1; i < nodes.size(); ++i) {

        char node = nodes[i];
        char left_son = nodes[2*i];
        char right_son = nodes[2*i + 1];
        inned_iterator += 3;

        vector< pair< Edge, Edge > > adacency;


        pair< Edge, Edge > left = find_pair(edges, node, left_son);
        pair< Edge, Edge > right = find_pair(edges, node, right_son);
        pair< Edge, Edge > parent = find_parent(edges, node);


        if (parent.first.my_id >= 0 && parent.second.my_id >= 0){
            adacency.push_back(parent);
        }


        if (left.first.my_id >= 0 && left.second.my_id >= 0 ){
            adacency.push_back(left);
        }

        if (right.first.my_id >= 0 && right.second.my_id >= 0){
            adacency.push_back(right);
        }

        adjacency_list.insert(pair< char, vector< pair< Edge, Edge > > >(node, adacency));
    }

    return adjacency_list;

}

void display_vector(vector< pair< Edge, Edge > > vec){
    for (int i = 0; i < vec.size(); ++i) {

        cout << "[\nFrom edge: ";
        display_one_direction(vec[i].first);
        cout << "To edge: ";
        display_one_direction(vec[i].second);
        cout << "]" << endl;

    }
}

Edge next(Edge e, map< char, vector< pair< Edge, Edge > > > adjacency_list){

    Edge next_e = {-1};

    for (auto &it : adjacency_list) {
//        cout << it.first <<  endl;
        int j = 0;
//        cout << "Node: " << it.first << endl;
//        display_vector(it.second);
//        cout << "Size of vector is: " << it.second.size() << endl;

        for (int i = 0; i < it.second.size(); ++i) {
//            cout << it.second[i].first.from_node;

            if (e.my_id == it.second[i].first.my_id &&
                    e.from_node == it.second[i].first.from_node){

//                cout << "Edge: " << e.my_id << " is from node " << it.first << endl;
                if (i + 1 >= it.second.size()){
//                    cout << "No next for " << e.my_id << endl;

                }else{
//                    cout << "Next for: " << e.my_id << " is " << it.second[i+1].first.my_id << endl;
                    next_e = it.second[i+1].first;
                }
            }
        }
//        cout << "" << endl;
    }
    return next_e;

}


Edge get_reversed(Edge e, vector< Edge > edges){
    Edge reverse = {-1};

    for (int i = 0; i < edges.size(); ++i) {
        if ( edges[i].is_back == !e.is_back
                && edges[i].from_node == e.to_node
             && edges[i].to_node == e.from_node){
            reverse = edges[i];
        }
    }
    return reverse;

}

Edge get_first_from_list(map< char, vector< pair< Edge, Edge > > > adj_list, Edge e){

    Edge edge = {-1};

    for (auto &it : adj_list) {
//        cout << it.first <<  endl;
        int j = 0;
//        cout << "Node: " << it.first << endl;
//        display_vector(it.second);
//        cout << "Size of vector is: " << it.second.size() << endl;

        for (int i = 0; i < it.second.size(); ++i) {
//            cout << it.second[i].first.from_node;

            if (e.my_id == it.second[i].first.my_id &&
                e.from_node == it.second[i].first.from_node){

                return it.second[0].first;

            }
        }
//        cout << "" << endl;
    }

    return e;
}


vector< pair< Edge, Edge > > construct_etour(const map< char, vector< pair< Edge, Edge > > > &adj_list, vector< Edge > edges, string nodes) {

    nodes = " " + nodes;
//    show_map(adj_list);

    vector< pair<Edge, Edge> > etour;
    for (int i = 0; i < edges.size(); ++i) {

        Edge rev = get_reversed(edges[i], edges);


//        cout << "For edge: " << edges[i].my_id << " reversed is " << rev.my_id << endl;

        Edge next_edge = next(rev, adj_list);

        if (next_edge.my_id == -1){
            next_edge = get_first_from_list(adj_list, rev);
        }

        etour.push_back(pair< Edge, Edge >(edges[i], next_edge));

    }

//    display_etour(etour);

    return etour;
}

vector< pair< Edge, int > > calculate_positions(vector< pair< Edge, Edge > > etour, char root){


    int edge_index = 0;
    vector< pair< Edge, int > > positiones_edges;

    Edge next_edge;

    // Find first edge
    for (int j = 0; j < etour.size(); ++j) {
        if (etour[j].first.my_id == 0 && etour[j].first.from_node == root){
            next_edge = etour[j].first;
        }
    }
//    display_one_direction(next_edge);

    positiones_edges.push_back(pair< Edge, int > (next_edge, edge_index) );

    ++edge_index;

    while (edge_index < etour.size()){

        for (int i = 0; i < etour.size(); ++i) {

            if (next_edge.my_id == etour[i].first.my_id){
                next_edge = etour[i].second;
                positiones_edges.push_back(pair< Edge, int > (next_edge, edge_index) );
                ++edge_index;
                break;
            }

        }
    }

//    for (int k = 0; k < positiones_edges.size(); ++k) {
//        cout << "Edge: " << positiones_edges[k].first.my_id << " on position: " << positiones_edges[k].second << endl;
//    }

    return positiones_edges;

}


int main(int argc, char** argv) {

    int numprocs;
    int myid;

    MPI_Request  request;
    MPI_Status stat;


    Edge myedge;



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
        map< char, vector< pair< Edge, Edge > > > list = construct_adacency_list(edges, input);
        vector< pair< Edge, Edge > > tour = construct_etour(list, edges, input);
        calculate_positions(tour, input[0]);



        // Broadcast profiles to all processes
        for (int i = 0; i < numprocs; ++i) {
            MPI_Send(&edges[i], sizeof(Edge), MPI_UNSIGNED, i, TAG, MPI_COMM_WORLD);
        }

//        cout << "------------" << endl;

//        for (int i = 0; i < edges.size(); ++i) {
//            cout << "id: " << edges[i].my_id << endl;
//            cout << "from: " << (char) edges[i].from_node << endl;
//            cout << "to: " << (char) edges[i].to_node << endl;
//            cout << "is_back: " << edges[i].is_back << endl;
//            cout << "------------" << endl;
//        }
    }


    // Receive information about myself
    MPI_Recv(&myedge, sizeof(Edge), MPI_UNSIGNED, 0, TAG ,MPI_COMM_WORLD, &stat);

    cout << "I'm: " << myid << " my edge is: "<< endl;
    cout << myedge.my_id << " " << myedge.from_node << " " << myedge.to_node << " " << myedge.is_back << endl;

    MPI_Finalize();
    return 0;
}