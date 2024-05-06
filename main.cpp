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

using namespace std;


// Struct Edge holds destination vertex and weight of the edge
struct Edge {
    int to;
    double weight;
};

class Graph {
private:
    int vertices; // nb of vertices in the graph
    std::vector<std::list<Edge>> adj_list; // adjacency list

public:
    Graph(int n) : vertices(n), adj_list(n) {}

    // Get the size of the graph (nb of vertices)
    int size() const {
        return vertices;
    }

    // Access the adjacency list of a vertex
    const std::list<Edge>& get_adjacent(int u) const {
        return adj_list[u];
    }

    // Add an edge between u and v with weight
    void add_edge(int u, int v, double weight) {
        if (u >= 0 && u < vertices && v >= 0 && v < vertices) {
            adj_list[u].push_back({v, weight});
        }
    }

    // Generate a small and basic graph
    void gen_small_graph() {
        add_edge(0, 1, 10);
        add_edge(0, 3, 30);
        add_edge(0, 4, 100);
        add_edge(1, 2, 50);
        add_edge(2, 4, 10);
        add_edge(3, 2, 20);
        add_edge(3, 4, 60);
    }



    void gen_graph_from_txt(std::string filename) {
        std::ifstream file(filename);
        std::string line;

        if (!file.is_open()) {
            throw std::runtime_error("Could not open file");
        }

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int u, v, weight;
            if (!(iss >> u >> v >> weight)) {
                throw std::invalid_argument("Invalid line format");
            }
            add_edge(u, v, weight);
        }

        file.close();
    }


    // Print the graph
    void print_graph() const {
        for (int i = 0; i < vertices; i++) {
            std::cout << "Vertex " << i << " has edges to:\n";
            for (const auto& edge : adj_list[i]) {
                std::cout << " -> " << edge.to << " with weight " << edge.weight << "\n";
            }
        }
    }
    
};



const int INF = INT_MAX;

std::vector<int> delta_stepping(int source, Graph& graph, int delta) {
    int n = graph.size();
    std::vector<int> dist(n, INF);
    std::vector<std::list<int>> buckets((INF / delta) + 1);

    dist[source] = 0;
    buckets[0].push_back(source);

    int idx = 0;
    while (idx < buckets.size()) {
        while (!buckets[idx].empty()) {
            int u = buckets[idx].front();
            buckets[idx].pop_front();

            for (const auto& e : graph.get_adjacent(u)) {
                int v = e.to;
                int weight = e.weight;
                int newDist = dist[u] + weight;

                if (newDist < dist[v]) {
                    if (dist[v] != INF) {
                        int oldBucketIdx = dist[v] / delta;
                        buckets[oldBucketIdx].remove(v);
                    }
                    dist[v] = newDist;
                    int newBucketIdx = newDist / delta;
                    buckets[newBucketIdx].push_back(v);
                }
            }
        }
        idx++;
    }

    // Print the distances
    for (int i = 0; i < n; ++i) {
        std::cout << "Distance from " << source << " to " << i << " is ";
        if (dist[i] == INF)
            std::cout << "infinity" << std::endl;
        else
            std::cout << dist[i] << std::endl;
    }

    return dist;
}





int main() {
    int type_graph = 1;  // 0 for small graph, 1 for txt graph, 2 for random graph 
    int run_both_algo = 0; // 0 run both, 1 run dijkstra, 2 run delta stepping
    int nb_vertices = 6;

    Graph g(6);

    if (type_graph == 0){
        std::cout << "\nGenerating a small graph:\n";
        g.gen_small_graph(); 
        g.print_graph(); 
    }
    else if (type_graph == 1){
        std::cout << "\nGenerating a small graph via text file:\n";
        //Graph g(6); 
        g.gen_graph_from_txt("small_graph_txt.txt");
        g.print_graph(); 
    }
    
    std::chrono::steady_clock::time_point begin_delta_stepping = std::chrono::steady_clock::now();
    std::vector<int> dist_delta_stepping;
    if (run_both_algo!=1){
        // Run delta stepping algo
        std::cout << "\nResults with delta stepping algo";
        int delta = 10; 
        std::cout << "\nDelta = " << delta << std::endl;
        dist_delta_stepping = delta_stepping(0, g, delta);
    }
    std::chrono::steady_clock::time_point end_delta_stepping = std::chrono::steady_clock::now();
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping - begin_delta_stepping).count() << " ms" << std::endl;
    
    return 0;
}