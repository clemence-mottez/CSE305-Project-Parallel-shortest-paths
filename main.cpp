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


// To declare an edge : 
// Edge<T> edge(v, weight);  

// To define a function using an Edge : 

// template <typename T>
// void some_function(const Edge<T>& edge) {
//   int dest = edge.dest;
//   T weight = edge.weight;



template <typename T>
class Graph {
private:
    int vertices; // nb of vertices in the graph
    std::vector<std::list<Edge<T>>> adj_list; // adjacency list
    int avgDegree; // Average degree of the graph (to compute)
    
    void computeAverageDegree() {
        int totalEdges = 0;
        for (const auto& list : adj_list) {
            totalEdges += list.size();
        }
        avgDegree = vertices ? std::round(static_cast<double>(totalEdges) / vertices) : 0;
    }


public:
    // Graph(int n) : vertices(n), adj_list(n) {}

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

        std::cout<<numPhysicalCores<<std::endl;
        std::cout<<vertices<<std::endl;
        // increase threads if the graph is large and dense
        if (vertices >= 1000) {
            suggestedThreads = std::min(2 * numPhysicalCores, vertices / 50);
        }
        // Of reduce number of thread if the graph is sparse
        if (vertices < 10) {
            suggestedThreads = std::max(1, numPhysicalCores / 2); 
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

    
    void gen_random_graph(int num_vertices, int num_edges, int min_weight, int max_weight) {
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
        if (std::is_same<T, double>::value){
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


//_________________________________________________________________________________________________________________________



// called by relaxRequests
// updates the tentative distance of node v (of neighbor u) if a shorter path is discovered
template <typename T>
void relax(int u, int v, T weight, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    T newDist = dist[u] + weight; //distance through neighbor
    if (newDist < dist[v]) { 
        if (dist[v] != INT_MAX) {  //new shortest path was found
            int oldBucketIdx = static_cast<int>(dist[v] / delta);
            buckets[oldBucketIdx].remove(v); // remove v from its current bucket
        }
        dist[v] = newDist;  // update new distance
        int newBucketIdx = static_cast<int>(newDist / delta);
        buckets[newBucketIdx].push_back(v); // assign v to a new bucket
    }
}


// creates requests for the edges of type isLight and of the nodes in R
// returns a set of edges
template <typename T>
std::set<int> findRequests(const std::list<int>& R, const Graph<T>& graph, int delta, const std::vector<T>& dist, bool isLight) {
    std::set<int> requests;
    for (int u : R) {
        for (const auto& e : graph.get_adjacent(u)) {
            int v = e.dest; 
            T weight = e.weight;
            if ((isLight && weight <= delta) || (!isLight && weight > delta)) { // checks if e has the correct type : light or heavy
                if (dist[u] + weight < dist[v]) { // found a new dist
                    requests.insert(v);
                }
            }
        }
    }
    return requests;
}

// do relaxations, may move nodes between buckets
template <typename T>
void relaxRequests(const std::set<int>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    for (int v : requests) {
        for (const auto& e : graph.get_adjacent(v)) {
            relax(v, e.dest, e.weight, dist, buckets, delta);
        }
    }
}

template <typename T>
std::vector<T> delta_stepping(int source, const Graph<T>& graph, int delta, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INT_MAX);
    int b = graph.nb_buckets(delta);
    //int b = int(ceil(INT_MAX / delta)) + 1;
    std::vector<std::list<int>> buckets(b); 

    dist[source] = 0;
    buckets[0].push_back(source); // Insert source node with distance 0

    for (int i = 0; i < buckets.size(); ++i) {
        while (!buckets[i].empty()) {   
            std::list<int> R = move(buckets[i]);
            for (int u : R) {
                for (const auto& e : graph.get_adjacent(u)) {
                    relax(u, e.dest, e.weight, dist, buckets, delta);
                }
            }

            // Relax light edges
            // Handle reinsertion here using findRequests and relaxRequests
            auto lightRequests = findRequests(R, graph, delta, dist, true);
            relaxRequests(lightRequests, graph, dist, buckets, delta);
            
            // Relax heavy edges
            auto heavyRequests = findRequests(R, graph, delta, dist, false);
            relaxRequests(heavyRequests, graph, dist, buckets, delta);
        }
    }

    // Print the distances
    if (print_dist){
        for (int i = 0; i < n; ++i) {
            std::cout << "Distance from " << source << " to " << i << " is ";
            if (dist[i] == INF){
                std::cout << "infinity" << std::endl;
            }
            else{
                std::cout << dist[i] << std::endl;
            }
        }
    }

    return dist;
}




// Deletion and edge relaxation for an entire bucket can be done in parallel and in arbitrary order 
// as long as an individual relaxation is atomic, i.e., the relaxations for a particular node are done sequentially

// relaxations for a particular node are done sequentially
// we use a mutex to ensure that updates to the distance vector and bucket lists are atomic
template <typename T>
void relaxPar(int u, int v, T weight, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta, std::mutex& m) {
    std::unique_lock<std::mutex> lock(m);
    relax(u, v, weight, dist, buckets, delta);
    lock.unlock();
}



// edge relaxation for an entire bucket can be done in parallel so we use threads
// we use multiple threads to process the requests and a mutex to manage shared resources safely
template <typename T>
void relaxRequestsPar(
    const std::set<int>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta, int nb_threads){
    std::vector<std::thread> threads; 
    std::mutex m; // mutex for synchronizing access to shared resources (dist and buckets)
    
    int batchSize = requests.size() / nb_threads; // compute size of each batch of requests to be processed by each thread

    auto it = requests.begin(); // iterator to access elements in the set of requests

    for (int i = 0; i < nb_threads; ++i) {
        // start iterator for the current thread
        auto start = std::next(it, i * batchSize);
        // end iterator
        auto end = (i == nb_threads - 1) ? requests.end() : std::next(start, batchSize);

        threads.emplace_back([&graph, start, end, &dist, &buckets, delta, &m](){
            for (auto it = start; it != end; ++it) {
                // iterate over all adjacent edges of the current node
                for (const auto& e : graph.get_adjacent(*it)) {
                    // perform relaxation step for current edge
                    relaxPar(*it, e.dest, e.weight, dist, buckets, delta, m);
                }
            }
        });
    }
    // wait for threads
    for (auto& t : threads) {
        t.join();
    }
}


// Previous version for relaxRequestsPar that might work better?
// // edge relaxation for an entire bucket can be done in parallel so we use threads
// template <typename T>
// void relaxRequestsPar(const std::set<int>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
//     std::vector<std::thread> threads;
//     std::mutex mutex;
//     for (int v : requests) {
//         threads.push_back(std::thread([&, v]() {
//             for (const auto& e : graph.get_adjacent(v)) {
//                 relaxPar(v, e.dest, e.weight, dist, buckets, delta, mutex);
//             }
//         }));
//     }
    
//     // wait for threads
//     for (std::thread& t : threads) {
//         if (t.joinable()) {
//             t.join();
//         }
//     }
// }



// Previous version that parallelize less but seems to run in less time
template <typename T>
std::vector<T> delta_stepping_Par(int source, const Graph<T>& graph, int delta, int nb_threads, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INT_MAX);
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b); 
    std::mutex mutex;

    dist[source] = 0;
    buckets[0].push_back(source); // insert source node with distance 0

    for (int i = 0; i < buckets.size(); ++i) {
        while (!buckets[i].empty()) {   
            
            std::list<int> R = move(buckets[i]);
            
            for (int u : R) {
                for (const auto& e : graph.get_adjacent(u)) {
                    relaxPar(u, e.dest, e.weight, dist, buckets, delta, mutex);
                }
            }

            // Relax light edges in parallel
            // Handle reinsertion here using findRequests and relaxRequests
            auto lightRequests = findRequests(R, graph, delta, dist, true);
            relaxRequestsPar(lightRequests, graph, dist, buckets, delta, nb_threads);
            
            // Relax heavy edges in parallel
            auto heavyRequests = findRequests(R, graph, delta, dist, false);
            relaxRequestsPar(heavyRequests, graph, dist, buckets, delta, nb_threads);
        }
    }

    // Print the distances
    if (print_dist){
        for (int i = 0; i < n; ++i) {
            std::cout << "Distance from " << source << " to " << i << " is ";
            if (dist[i] == INF){
                std::cout << "infinity" << std::endl;
            }
            else{
                std::cout << dist[i] << std::endl;
            }
        }
    }

    return dist;
}

// // More parallelize but takes more time to work
// template <typename T>
// std::vector<T> delta_stepping_Par(int source, const Graph<T>& graph, int delta, int nb_threads, bool print_dist) {
//     int n = graph.size();
//     std::vector<T> dist(n, INT_MAX);
//     int b = graph.nb_buckets(delta);
//     std::vector<std::list<int>> buckets(b); 
//     std::mutex mutex;

//     dist[source] = 0;
//     buckets[0].push_back(source); // insert source node with distance 0

//     for (int i = 0; i < buckets.size(); ++i) {
//         while (!buckets[i].empty()) {
            
//             std::list<int> R = move(buckets[i]);

//             std::vector<std::thread> threads;
//             auto start_iter = R.begin();

//             int num_nodes_per_thread = R.size() / nb_threads;
//             int remainder = R.size() % nb_threads;

//             for (int t = 0; t < nb_threads; ++t) {
//                 int num_nodes_this_thread = num_nodes_per_thread + (t < remainder ? 1 : 0);
//                 auto end_iter = std::next(start_iter, num_nodes_this_thread);

//                 threads.emplace_back([&graph, &dist, &buckets, &mutex, delta, start_iter, end_iter]() {
//                     for (auto it = start_iter; it != end_iter; ++it) {
//                         int u = *it;
//                         for (const auto& e : graph.get_adjacent(u)) {
//                             relaxPar(u, e.dest, e.weight, dist, buckets, delta, mutex);
//                         }
//                     }
//                 });

//                 start_iter = end_iter;
//             }

//             for (auto& thread : threads) {
//                 thread.join();
//             }

//             // Relax light edges in parallel
//             auto lightRequests = findRequests(R, graph, delta, dist, true);
//             relaxRequestsPar(lightRequests, graph, dist, buckets, delta, nb_threads);
            
//             // Relax heavy edges in parallel
//             auto heavyRequests = findRequests(R, graph, delta, dist, false);
//             relaxRequestsPar(heavyRequests, graph, dist, buckets, delta, nb_threads);
//         }
//     }

//     // Print the distances
//     if (print_dist){
//         for (int i = 0; i < n; ++i) {
//             std::cout << "Distance from " << source << " to " << i << " is ";
//             if (dist[i] == INT_MAX) {
//                 std::cout << "infinity" << std::endl;
//             } else {
//                 std::cout << dist[i] << std::endl;
//             }
//         }
//     }

//     return dist;
// }

//______________________________________________________________________________________________________________________________



template <typename T>
std::vector<T> dijkstra(int source, Graph<T>& graph, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INT_MAX); 
    std::vector<bool> visited(n, false); 
    //std::priority_queue<std::pair<int, T>, std::vector<std::pair<int, T>>, std::greater<std::pair<int, T>>> pq;
    std::priority_queue< Edge<T>, std::vector<Edge<T>>, CompareEdge<T>> pq;

    // Initialize priority queue with the source node
    Edge<T> e1(0, source);
    pq.push(e1); // (distance, vertex)
    dist[source] = 0;

    while (!pq.empty()) {
        T u = pq.top().weight;
        pq.pop();

        // If vertex already been visited, continue to the next
        if (visited[u]) continue;
        visited[u] = true;

        // Check each adjacent vertex of u
        for (const auto& edge : graph.get_adjacent(u)) {
            int v = edge.dest;
            T weight = edge.weight;

            // Only consider this vertex if it has not been visited
            if (!visited[v] && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                Edge<T> e2(dist[v],v);
                pq.push(e2); // Push updated distance and vertex
            }
        }
    }


    if (print_dist){
        for (int i = 0; i < n; ++i) {
            std::cout << "Distance from " << source << " to " << i << " is ";
            if (dist[i] == INF)
                std::cout << "infinity" << std::endl;
            else
                std::cout << dist[i] << std::endl;
        }
    }

    return dist;
}

// Compare distances from the 2 algorithms to check is they match
template <typename T>
void compare_distances(const std::vector<T>& dist1, const std::vector<T>& dist2) {
    if (dist1.size() != dist2.size()) {
        std::cerr << "Error: Distance vectors are of different sizes." << std::endl;
        return;
    }

    //double mse = 0;
    int diff_count = 0;
    for (size_t i = 0; i < dist1.size(); i++) {
        if (dist1[i] != dist2[i]) {
            diff_count++;
            // if (dist1[i] != INT_MAX && dist2[i] != INT_MAX) { 
            //     mse += pow(dist1[i] - dist2[i], 2);
            // }
        }
    }

    // if (diff_count > 0) {
    //     mse /= diff_count; 
    // }

    // Print MSE and nb of different values
    // std::cout << "Mean Squared Error: " << mse << std::endl;
    std::cout << "Number of different values: " << diff_count << std::endl;
}




template <typename T>
int continue_main(Graph<T> g, int run_algo, int delta, int print_graph, int print_dist, int num_threads){

    if (print_graph){
        g.print_graph(); 
    }

    double t1 = 0;
    double t2 = 0;
    double t3 = 0;
    std::vector<T> dist_dijkstra;
    std::vector<T> dist_delta_stepping;
    std::vector<T> dist_delta_stepping_threads;

    // running and chronometring the algorithms: 

    if (run_algo == 1 || run_algo == 4 || run_algo == 5 || run_algo == 0){// Run dijkstra algo
        std::cout << "\nResults with Dijkstra algo\n";
        std::chrono::steady_clock::time_point begin_dijkstra = std::chrono::steady_clock::now();
        dist_dijkstra = dijkstra(0, g, print_dist);
        std::chrono::steady_clock::time_point end_dijkstra = std::chrono::steady_clock::now();
        t1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_dijkstra - begin_dijkstra).count() ;
        std::cout << "Total time with Dijkstra: " << t1 << " ms" << std::endl;

    }

    if (run_algo == 2 || run_algo == 4 || run_algo == 6 || run_algo == 0){ // Run delta-stapping
        std::cout << "\nResults with delta stepping algo, delta = " << delta << std::endl;
        std::chrono::steady_clock::time_point begin_delta_stepping = std::chrono::steady_clock::now();
        dist_delta_stepping = delta_stepping(0, g, delta, print_dist);
        std::chrono::steady_clock::time_point end_delta_stepping = std::chrono::steady_clock::now();
        t2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping - begin_delta_stepping).count();
        std::cout << "Total time with Delta stepping: " << t2 << " ms" << std::endl;

    }

    if (run_algo == 3 || run_algo == 5 || run_algo == 6 || run_algo == 0){ // Run delta-stepping threads
        std::cout << "\nResults with delta stepping threads algo, delta = " << delta << " , nb of threads = " << num_threads << std::endl;
        std::chrono::steady_clock::time_point begin_delta_stepping_threads = std::chrono::steady_clock::now();
        // dist_delta_stepping_threads = delta_stepping_parallel(0, g, delta, num_threads, print_dist);
        dist_delta_stepping_threads = delta_stepping_Par(0, g, delta, num_threads, print_dist);
        std::chrono::steady_clock::time_point end_delta_stepping_threads = std::chrono::steady_clock::now();
        t3 = std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping_threads - begin_delta_stepping_threads).count() ;
        std::cout << "Total time with Delta stepping threads: " << t3 << " ms" << std::endl;

    }

    //comparison of algorithms (results & runtimes): 

    if (run_algo == 4 || run_algo == 0){// compare dijkstra & delta-stepping: 
        std::cout << "\nComparing Dijkstra / delta-stepping "<< std::endl;
        double speed_up = t1/t2; 
        std::cout << "Speed up: " << speed_up << std::endl;
        compare_distances(dist_dijkstra, dist_delta_stepping);
    }

    // if (run_algo == 5 || run_algo == 0){// compare dijkstra & delta-stepping threads:
    //     std::cout << "\nComparing Dijkstra / delta-stepping threads "<< std::endl; 
    //     double speed_up = t1/t3; 
    //     std::cout << "Speed up: " << speed_up << std::endl;
    //     compare_distances(dist_dijkstra, dist_delta_stepping_threads);
    // }

    if (run_algo == 6 || run_algo == 0){// compare delta-stepping & delta-stepping threads: 
        std::cout << "\nComparing delta-stepping / delta-stepping threads "<< std::endl;
        double speed_up = t2/t3; 
        std::cout << "Speed up: " << speed_up << std::endl;
        compare_distances(dist_delta_stepping, dist_delta_stepping_threads);
    }

    return 0;
}





// For all type of graphs
// int type_graph = 0 for small graph, 1 for txt graph, 2 for random graph 
// bool print_dist = if want to print the resulting distances or not, it affects the running time so put 0 preferably
// bool print_graph = whether or not want to print the graph
// int run_algo = 1 : dijkstra ; 2 : delta-stepping ; 3 : DS-threads ; 4 : compare dijkstra & DS ; 5 : compare dijkstra & DS threads ; 6 : compare DS & DS threads ; 0 : compare all 
// int delta = 0 if want to use computed value, or = value if want a specific value
// int num_threads = 0 if want to use computed value, or = value if want a specific value 

// For txt graph
// std::string name_of_txt = "txt_graph_1000.txt" for example
// int num_vertices = nb vertices

// For random graph
// int num_vertices = nb vertices
// int num_edges =  nb edges
// int min_weight = min weight (>0 positive weights)
// int max_weight = max weight



// Examples to run:

// First run: g++ main.cpp -o test

// Then:
// For testing an existing small graph
// ./test 0 [run_algo] [delta] [num_threads] [print_dist] [print_graph]
// ./test 0 0 5 6 1 1

// For testing a graph from a txt file
// ./test 1 name_of_txt_file [num_vertices] [run_algo] [delta] [num_threads] [print_dist] [print_graph]
// ./test 1 txt_graph_1000.txt 1000 0 0 0 0 0

// For testing a random graph
// ./test 2 [num_vertices] [num_edges] [min_weight] [max_weight] [run_algo] [delta] [num_threads] [print_dist] [print_graph]
// ./test 2 1000 10000 1 50 0 10 10 0 0

int main(int argc, char* argv[]) {

    int type_graph = std::atoi(argv[1]);

    std::string name_of_txt;
    int num_vertices, num_edges, min_weight, max_weight;
    int run_algo, delta, num_threads;
    bool print_dist, print_graph;

    switch (type_graph) {
    case 0: // Small graph
        break;

    case 1: // txt graph
        name_of_txt = argv[2];
        num_vertices = std::atoi(argv[3]);
        break;

    case 2: // Random graph
        num_vertices = std::atoi(argv[2]);
        num_edges = std::atoi(argv[3]);
        min_weight = std::atoi(argv[4]);
        max_weight = std::atoi(argv[5]);
        break;

    default:
        std::cerr << "Invalid type_graph, valid options are 0, 1, or 2." << std::endl;
        return 1;
    }

    run_algo = std::atoi(argv[argc-5]);
    delta = std::atoi(argv[argc-4]);
    num_threads = std::atoi(argv[argc-3]);
    print_dist = std::atoi(argv[argc-2]);
    print_graph = std::atoi(argv[argc-1]);


    if (type_graph == 0){
        std::cout << "\nGenerating a small graph\n";
        if (run_algo == 2 || run_algo == 3 || run_algo ==  6 ){ // creating a graph with real-valued positive edge weights (run_algo does not call dijkstra)
            Graph<double> g(6); //change with the nb of vertices in the fixed graph
            g.gen_small_graph_real();
            if (delta==0){
                delta = g.findDelta();
            }
            if (num_threads==0){
                num_threads = g.suggestOptimalNumberOfThreads();
            }
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        } 
        else {
            Graph<int> g(6);
            g.gen_small_graph_int();
            if (delta==0){
                delta = g.findDelta();
            }
            if (num_threads==0){
                num_threads = g.suggestOptimalNumberOfThreads();
            }
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads); 
        }
    }

    else if (type_graph == 1){
        std::cout << "\nGenerating a graph via text file\n";
        if (run_algo == 2 || run_algo == 3 || run_algo ==  6 ){ // creating a graph with real-valued positive edge weights (run_algo does not call dijkstra)
            Graph<double> g(num_vertices);
            g.gen_graph_from_txt(name_of_txt);
            if (delta==0){
                delta = g.findDelta();
            }
            if (num_threads==0){
                num_threads = g.suggestOptimalNumberOfThreads();
            }
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        } 
        else {
            Graph<int> g(num_vertices); //change with the nb of vertices in the txt graph
            g.gen_graph_from_txt(name_of_txt);
            if (delta==0){
                delta = g.findDelta();
            }
            if (num_threads==0){
                num_threads = g.suggestOptimalNumberOfThreads();
            }
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        }
    }

    else if (type_graph == 2){
        std::cout << "\nGenerating a random graph\n";
        if (run_algo == 2 || run_algo == 3 || run_algo ==  6 ){ // creating a graph with real-valued positive edge weights (run_algo does not call dijkstra)
            Graph<double> g(num_vertices);
            g.gen_random_graph(num_vertices, num_edges, min_weight, max_weight);
            if (delta==0){
                delta = g.findDelta();
            }
            if (num_threads==0){
                num_threads = g.suggestOptimalNumberOfThreads();
            }
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        } 
        else {
            Graph<int> g(num_vertices);
            g.gen_random_graph(num_vertices, num_edges, min_weight, max_weight);
            if (delta==0){
                delta = g.findDelta();
            }
            if (num_threads==0){
                num_threads = g.suggestOptimalNumberOfThreads();
            }
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        }
    }

}

