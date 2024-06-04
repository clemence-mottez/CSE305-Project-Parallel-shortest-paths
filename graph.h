
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <climits>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <ctime>
#include <chrono>
#include <thread>
#include <mutex>
#include <set>
#include <atomic>
#include <functional>
#include <type_traits>
#include <cmath>
#include <future>
#include <condition_variable>
#include <functional>
#include <cstdlib> 



std::random_device rd;
std::mt19937 gen(rd());
const int INF = INT_MAX;


template <typename T>
struct Edge {
    int dest;
    T weight;
    Edge(int d, T w) : dest(d), weight(w) {}
};


// template <typename T>
// class Graph {
// private:
//     int vertices;
//     std::vector<std::list<Edge<T>>> adj_list;
//     int avgDegree;
//     void computeAverageDegree();

// public:
//     Graph(int n);
//     int size() const;
//     int nb_buckets(int delta) const;
//     const std::list<Edge<T>>& get_adjacent(int u) const;
//     void add_edge(int u, int v, T weight);
//     void print_graph() const;
//     void gen_small_graph_int();
//     void gen_small_graph_real();
//     void gen_random_graph(int type_weight, int num_vertices, int num_edges, int min_weight, int max_weight);
//     void gen_graph_from_txt(const std::string& filename);
//     int suggestOptimalNumberOfThreads() const;
//     int findDelta();
// };

// Custom hash function for pairs
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2> &pair) const {
        auto h1 = std::hash<T1>{}(pair.first);
        auto h2 = std::hash<T2>{}(pair.second);
        return h1 ^ h2;
    }
};

template <typename T>
struct CompareEdge {
    bool operator()(const Edge<T>& a, const Edge<T>& b) const {
        return a.dest > b.dest; // Prioritize lower destination values
    }
};


template <typename T>
class Graph {
private:
    int vertices; // nb of vertices in the graph
    std::vector<std::list<Edge<T>>> adj_list; // adjacency list
    int avgDegree; // Average degree of the graph (to compute)
    
    void computeAverageDegree() { // average number of edges connected to each vertex
        int totalEdges = 0;
        for (const auto& list : adj_list) {
            totalEdges += list.size();
        }
        if (vertices != 0){
            avgDegree = std::round(static_cast<double>(totalEdges) / vertices);
        } else{
            avgDegree = 0;
        }
    }


public:

    Graph(int n) : vertices(n), adj_list(n), avgDegree(0) {
        if (vertices > 0) {
            computeAverageDegree();
        }
    }

    // Get the size of the graph (nb of vertices)
    int size() const {
        return vertices;
    }

    // computes the sufficient nb of buckets to use
    int nb_buckets(int delta) const {   
        int max_bucket = 0;
        for (const auto& list : adj_list) {
            for (const auto& edge : list) {
                int bucket = static_cast<int>(std::ceil(edge.weight / delta));
                if (bucket > max_bucket) {
                    max_bucket = bucket;
                }
            }
        }
        return max_bucket + 1; 
    }


    // Access the adjacency list of a vertex
    const std::list<Edge<T>>& get_adjacent(int u) const {
        return adj_list[u];
    }

    std::vector<std::list<Edge<T>>>  get_adj_list() const
    {
        return adj_list;
    }

    // Add an edge between u and v with weight
    void add_edge(int u, int v, T weight) {
        if (u >= 0 && u < size() && v >= 0 && v < size()) {
            Edge<T> edge(v, weight);
            adj_list[u].push_back(edge);
        }
    }


    // Find delta function
    int findDelta() {
        // Find the minimum edge weight
        T delta0 = std::numeric_limits<T>::max();
        for (const auto& list : adj_list) {
            for (const auto& edge : list) {
                if (edge.weight < delta0) delta0 = edge.weight;
            }
        }

        if (delta0 == std::numeric_limits<T>::max()) {
            return 0; // No edges in the graph
        }

        T deltaCur = delta0;
        std::unordered_map<std::pair<int, int>, T, pair_hash> found;
        bool change = true;

        while (change) {
            change = false;
            std::unordered_map<int, std::list<std::pair<int, T>>> T_list;

            for (int u = 0; u < vertices; ++u) {
                for (const auto& edge : adj_list[u]) {
                    int v = edge.dest;
                    T weight = edge.weight;
                    if (weight <= deltaCur) {
                        int j = static_cast<int>(std::floor(weight / delta0));
                        T_list[j].emplace_back(u, v);
                        if (found.find({u, v}) == found.end() || found[{u, v}] > weight) {
                            found[{u, v}] = weight;
                            change = true;
                        }
                    }
                }
            }
            deltaCur *= 2;
        }
        deltaCur /= 2;
        int delta = static_cast<int>(deltaCur);
        return delta;
    }


    int suggestOptimalNumberOfThreads() const {
        int numPhysicalCores = std::thread::hardware_concurrency();
        int suggestedThreads = numPhysicalCores; 

        // std::cout<<numPhysicalCores<<std::endl;
        // std::cout<<vertices<<std::endl;
        // increase threads if the graph is large and dense
        if (vertices >= 1000) {
            suggestedThreads = 2 * numPhysicalCores;
        }
        // Of reduce number of thread if the graph is sparse
        if (vertices < 10) {
            suggestedThreads = numPhysicalCores / 2; 
        }

        return suggestedThreads;
    }


    // Generate a small and basic graph
    void gen_small_graph_int() {
        add_edge(0, 1, 10);
        add_edge(0, 3, 3);
        add_edge(0, 4, 10);
        add_edge(1, 2, 5);
        add_edge(2, 4, 1);
        add_edge(3, 2, 2);
        add_edge(3, 1, 6);
        add_edge(4, 2, 3);
        add_edge(4, 3, 3);
        add_edge(4, 1, 1);
        add_edge(5, 4, 10);
        add_edge(5, 2, 2);
        add_edge(5, 0, 6);
    }

    void gen_small_graph_real() {
        add_edge(0, 1, 10.1);
        add_edge(0, 3, 3.2);
        add_edge(0, 4, 10.);
        add_edge(1, 2, 5.5);
        add_edge(2, 4, 1.1);
        add_edge(3, 2, 2.4);
        add_edge(3, 1, 6.6);
        add_edge(4, 2, 3.7);
        add_edge(4, 3, 3.);
        add_edge(4, 1, 1.8);
        add_edge(5, 4, 10.9);
        add_edge(5, 2, 2.3);
        add_edge(5, 0, 6.);
    }

    
    void gen_random_graph(int type_weight , int num_vertices, int num_edges, int min_weight, int max_weight) {
        if (num_edges > (num_vertices * (num_vertices - 1))) {
        // Too many edges for acyclic graph
            std::cout << "Number of edges exceeds limit for acyclic graph" << std::endl;
        }
        std::set<int> visited; //set to keep track of visited nodes to avoid cycles

        for (int i = 0; i < num_edges; i++) {
        int source = std::uniform_int_distribution<int>(0, num_vertices - 1)(gen); //pick a random source
        int target;

        // finding a valid target 
        while (visited.count(target) > 0 && source != target) {
            target = std::uniform_int_distribution<int>(0, num_vertices - 1)(gen);
        }

        T weight; 
        if (type_weight == 1){
            weight = std::uniform_real_distribution<double>(min_weight, max_weight)(gen); //random real positive weight
        } 
        else {
            weight = std::uniform_int_distribution<int>(min_weight, max_weight)(gen); //random integer weight
        }
        add_edge(source, target, weight); // adds edge to the graph
        visited.insert(target); // marks the target node as visited
        }
    }

    void gen_graph_from_txt(std::string filename) {
        std::ifstream file(filename);
        std::string line;

        if (!file.is_open()) {
            throw std::runtime_error("Could not open file");
        }

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int u, v;
            T weight;
            if (!(iss >> u >> v >> weight)) {
                throw std::invalid_argument("Invalid line format");
            }
            add_edge(u, v, weight);
        }
        
        file.close();
    }

    void gen_graph_from_txt_without_weight(std::string filename) {
        std::ifstream file(filename);
        std::string line;

        if (!file.is_open()) {
            throw std::runtime_error("Could not open file");
        }

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int u, v;
            T weight = 1;
            if (!(iss >> u >> v)) {
                throw std::invalid_argument("Invalid line format");
            }
            add_edge(u, v, 1);
        }
        
        file.close();
    }

    

    // Print the graph
    void print_graph() const {
        for (int i = 0; i < vertices; i++) {
            std::cout << "Vertex " << i << " has edges to:\n";
            for (const auto& edge : adj_list[i]) {
                std::cout << " -> " << edge.dest << " with weight " << edge.weight << "\n";
            }
        }
    }


    
};

#endif
