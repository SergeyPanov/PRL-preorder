#include <iostream>
#include "/Users/sergeypanov/bin/mpi/include/mpi.h"
#include <fstream>
#include <vector>
#include <map>

#define TAG 0

using namespace std;


template<class T>
std::vector<unsigned char> serialize_object(T object)
{
    T moved_object = std::move(object);
    std::vector<unsigned char> byte_v;
    auto *data_ptr = reinterpret_cast<unsigned char *>(&moved_object);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        byte_v.push_back(data_ptr[i]);
    }
    return byte_v;
}

template<class T>
T deserialize_to_object(std::vector<unsigned char>& byte_v)
{
    return *reinterpret_cast<T*>(byte_v.data());
}

template<class T>
std::vector<unsigned char> serialize_pair(std::pair<T, T> p)
{
    std::vector<unsigned char> serialized_pair;
    std::vector<unsigned char> first_serialized = serialize_object(p.first);
    std::vector<unsigned char> second_serialized = serialize_object(p.second);
    serialized_pair.insert(serialized_pair.end(), first_serialized.begin(), first_serialized.end());
    serialized_pair.insert(serialized_pair.end(), second_serialized.begin(), second_serialized.end());
    return serialized_pair;
}

template<class T>
std::pair<T, T> deserialize_to_pair_of(std::vector<unsigned char> &byte_v)
{
    std::vector<unsigned char> first_serialized(byte_v.begin(), byte_v.begin() + byte_v.size()/2);
    std::vector<unsigned char> second_serialized(byte_v.begin() + byte_v.size()/2, byte_v.end());
    return std::pair<T, T>(deserialize_to_object<T>(first_serialized), deserialize_to_object<T>(second_serialized));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

void display_adacency_list(map< char, vector< pair< Edge, Edge > > > adjacencies){

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
        int j = 0;

        for (int i = 0; i < it.second.size(); ++i) {

            if (e.my_id == it.second[i].first.my_id &&
                    e.from_node == it.second[i].first.from_node){
                if (i + 1 >= it.second.size()){

                }else{
                    next_e = it.second[i+1].first;
                }
            }
        }
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
        int j = 0;

        for (int i = 0; i < it.second.size(); ++i) {

            if (e.my_id == it.second[i].first.my_id &&
                e.from_node == it.second[i].first.from_node){

                return it.second[0].first;

            }
        }
    }

    return e;
}


Edge get_my_next(const map< char, vector< pair< Edge, Edge > > > &adj_list, vector< Edge > edges, Edge me) {
    Edge my_next = {-1};

    Edge rev = get_reversed(me, edges);

    Edge next_edge = next(rev, adj_list);

    if (next_edge.my_id == -1){
        next_edge = get_first_from_list(adj_list, rev);
    }

    return next_edge;
}


vector< pair< Edge, Edge > > construct_etour(const map< char, vector< pair< Edge, Edge > > > &adj_list, vector< Edge > edges, string nodes) {

    nodes = " " + nodes;

    vector< pair<Edge, Edge> > etour;
    for (int i = 0; i < edges.size(); ++i) {
        Edge rev = get_reversed(edges[i], edges);

        Edge next_edge = next(rev, adj_list);

        if (next_edge.my_id == -1){
            next_edge = get_first_from_list(adj_list, rev);
        }
        etour.push_back(pair< Edge, Edge >(edges[i], next_edge));
    }


    return etour;
}

map<int, int> calculate_positions(vector< pair< Edge, Edge > > etour, char root){


    int edge_index = 0;
    vector< pair< Edge, int > > positiones_edges;


    map< int, int > positions;

    Edge next_edge;

    // Find first edge
    for (int j = 0; j < etour.size(); ++j) {
        if (etour[j].first.my_id == 0 && etour[j].first.from_node == root){
            next_edge = etour[j].first;
        }
    }

    positiones_edges.push_back(pair< Edge, int > (next_edge, edge_index) );

    positions.insert(pair< int, int > (edge_index, next_edge.my_id));

    ++edge_index;

    while (edge_index < etour.size()){

        for (int i = 0; i < etour.size(); ++i) {

            if (next_edge.my_id == etour[i].first.my_id){
                next_edge = etour[i].second;
                positiones_edges.push_back(pair< Edge, int > (next_edge, edge_index) );
                positions.insert(pair< int, int > (edge_index, next_edge.my_id));
                ++edge_index;
                break;
            }

        }
    }

    return positions;

}

int get_id_for_pos(int pos,  map< int, int > id_pos){

    for (auto const& x : id_pos) {
        if (x.second == pos){
            return x.first;
        }
    }
    return -1;

}


int main(int argc, char** argv) {

    int numprocs;
    int myid;

    MPI_Request  request;
    MPI_Status stat;

    int tour_size;

    Edge myedge;

    vector< pair< Edge, Edge > > my_tour;
    vector< pair< Edge, int > > my_positioned;



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

        // Broadcast profiles to all processes
        for (int i = 0; i < numprocs; ++i) {
            MPI_Send(&edges[i], sizeof(Edge), MPI_UNSIGNED, i, TAG, MPI_COMM_WORLD);
        }

        // Send broadcast input
        int input_size = static_cast<int>(input.size());
        for (int j = 0; j < numprocs; ++j) {
            MPI_Send(&input_size, 1, MPI_INT, j, TAG, MPI_COMM_WORLD);
            MPI_Send(input.c_str(), static_cast<int>(input.size()), MPI_CHAR, j, TAG, MPI_COMM_WORLD);
        }

    }


    // Receive information about myself
    MPI_Recv(&myedge, sizeof(Edge), MPI_UNSIGNED, 0, TAG ,MPI_COMM_WORLD, &stat);

    int input_size;
    MPI_Recv(&input_size, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat);



    string input;
    input.resize(input_size);
    MPI_Recv((void *) input.data(), input_size, MPI_CHAR, 0, TAG, MPI_COMM_WORLD, &stat);



    vector<Edge> edges = constructEdges(input); // Create vector of edges
    map< char, vector< pair< Edge, Edge > > > list = construct_adacency_list(edges, input); // Create adjacency list


    Edge my_next = get_my_next(list, edges, myedge);

    pair< Edge, Edge > me_mynext = pair<Edge, Edge>(myedge, my_next);

//    cout << "ME: " << myid << endl;


    vector< pair< Edge, Edge > > complete_etour;

    // Broadcast information about my next to all processors. Create complete etour
    for (int k = 0; k < numprocs; ++k) {
        vector< unsigned char > vec = serialize_pair(me_mynext);
        int local_size = static_cast<int>(vec.size());

        MPI_Send(&local_size, 1, MPI_INT, k, TAG, MPI_COMM_WORLD);

        MPI_Send(vec.data(), static_cast<int>(vec.size()), MPI_UNSIGNED_CHAR, k, TAG, MPI_COMM_WORLD);
    }

    // Receive broadcasted edges
    for (int l = 0; l < numprocs; ++l) {
        int loc_size;
        MPI_Recv(&loc_size, 1, MPI_INT, l, TAG, MPI_COMM_WORLD, &stat);

        vector<unsigned char > vec;
        vec.resize(sizeof(pair< Edge, Edge >));

        MPI_Recv(vec.data(), vec.size(), MPI_UNSIGNED_CHAR, l, TAG , MPI_COMM_WORLD, &stat);
        pair< Edge, Edge > pr = deserialize_to_pair_of< Edge >(vec);

        complete_etour.push_back(pr);
    }


//    vector< pair< Edge, int > > positiones_edges;

    map<int, int> positiones_edges;
    positiones_edges = calculate_positions(complete_etour, input[0]);   // Calculate positions of each edge



    ////////////////////////////////// Pre order algorithm //////////////////////////////////


    /////////////////////////// Stage 1 //////////////////////////////////////////////////////

    ///////// Broadcast info about direction(0 is backward, 1 is forward) /////////
    for (int k = 0; k < numprocs; ++k) {
        int direction = (!myedge.is_back);
        MPI_Send(&direction, 1, MPI_INT, k, TAG, MPI_COMM_WORLD);
    }


    ///////// Receive directions(0 is backward, 1 is forward) /////////
    map< int, int > directions;
//    vector< pair< int, int > > directions;
    for (int l = 0; l < numprocs; ++l) {
        int direction;
        MPI_Recv(&direction, 1, MPI_INT, l, TAG, MPI_COMM_WORLD, &stat);

        pair< int, int > pr = pair< int, int > (l, direction);
        directions.insert(pr);
    }



    /////////////////////////// Stage 2 //////////////////////////////////////////////////////
    ///////// Calculate suffix sum based on positions and direction /////////

    int index = positiones_edges.size() - 1;    // Start from last
    map< int, int > suffix_sum;  // Map with suffix sums

    int e_id = positiones_edges.at(index);

    suffix_sum.insert(pair< int, int >(index, 0));  // Suffix sum for last edge is 0

    --index;

    int total_suffix_sum = 0;

    while (index >= 0){

        int edge_id = positiones_edges.at(index);

        total_suffix_sum += directions.at(edge_id);

        suffix_sum.insert( pair<int, int>(edge_id, total_suffix_sum));
        --index;
    }



    /////////////////////////// Stage 3 //////////////////////////////////////////////////////
    ///////// Corection /////////

    int correction_value = input.size() - 1;
    for(auto& suffix : suffix_sum) {
        suffix.second = correction_value - suffix.second;
    }



//    cout << "I'm: " << myid << endl;
//    for (auto const& x : suffix_sum) {
//        std::cout << x.first
//                  << ':'
//                  << x.second
//                  << std::endl ;
//    }


    map< int, int > forward_edges;

    // Get only forward edges in structure [position, id]
    for (auto const& x : directions) {
        if (x.second == 1){
            int pos = suffix_sum.at(x.first);
            forward_edges.insert(pair< int, int >(x.first, pos));
        }
    }

//    cout << "I'm: " << myid << endl;
//    for (auto const& x : forward_edges) {
//        std::cout << x.first
//                  << ':'
//                  << x.second
//                  << std::endl ;
//    }




    if (!myedge.is_back){

        int my_pos = forward_edges.at(myid);

        int next_id = get_id_for_pos(my_pos + 1, forward_edges);

        if (myid == 0){
//            cout << "I'm: " << myid << " my next: " << next_id << endl;

            cout << myedge.from_node << myedge.to_node;

            MPI_Send(&next_id, 1, MPI_INT, next_id, TAG, MPI_COMM_WORLD);

        } else{
            int dummy;
            int prev_id = get_id_for_pos(my_pos - 1, forward_edges);

//            cout << "I'm: " << myid << " my next: " << next_id << endl;

            MPI_Recv(&dummy, 1, MPI_INT, prev_id, TAG, MPI_COMM_WORLD, &stat);
            cout << myedge.to_node;

            if (next_id != -1){
                MPI_Send(&next_id, 1, MPI_INT, next_id, TAG, MPI_COMM_WORLD);
            }

        }

    }



//    cout << "------------------------------------------------------" << endl;
//    for (int m = 0; m < directions.size(); ++m) {
//        cout << "Proc: " << directions[m].first << " has direction: " << directions[m].second << endl;
//    }



//    cout << "I'm: " << myid << endl;
//
//    for (int l = 0; l < positiones_edges.size(); ++l) {
//        cout << "Edge: " << positiones_edges[l].first.my_id << " on position: " << positiones_edges[l].second << endl;
//    }



    MPI_Finalize();
    return 0;
}