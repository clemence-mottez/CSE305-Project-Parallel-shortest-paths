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


std::random_device rd;
std::mt19937 gen(rd());
const int INF = INT_MAX;

template <typename T>
struct Edge {
    int dest;
    T weight;
    Edge(int d, T w) : dest(d), weight(w) {}
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

public:
    Graph(int n) : vertices(n), adj_list(n) {}

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


// int calculate_b(const Graph& graph, double delta) {
//     int max_bucket = 0;
//     for (const auto& list : graph.get_adj_list()) { // Access each adjacency list
//         for (const auto& edge : list) { // Access each edge in the list
//             double weight = edge.second;
//             int bucket = static_cast<int>(std::ceil(weight / delta));
//             if (bucket > max_bucket) {
//                 max_bucket = bucket;
//             }
//         }
//     }
//     return max_bucket + 1; // +1 as per the formula provided
// }

//_________________________________________________________________________________________________________________________


//called by relaxRequests
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
            // Handle reinsertion here using findRequests and relaxRequests
            auto lightRequests = findRequests(R, graph, delta, dist, true);
            relaxRequests(lightRequests, graph, dist, buckets, delta);

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

/*
// Function to perform relaxations on light edges
void relaxLightEdges(const std::list<int>& R, const Graph& graph, std::vector<double>& dist, std::vector<std::list<int>>& buckets, int delta, std::mutex& mtx) {
    std::set<int> lightRequests = findRequests(R, graph, delta, dist, true);
    std::unique_lock<std::mutex> lock(mtx);
    relaxRequests(lightRequests, graph, dist, buckets, delta);
    lock.unlock();
}

// Function to perform relaxations on heavy edges
void relaxHeavyEdges(const std::list<int>& R, const Graph& graph, std::vector<double>& dist, std::vector<std::list<int>>& buckets, int delta, std::mutex& mtx) {
    std::set<int> heavyRequests = findRequests(R, graph, delta, dist, false);
    std::unique_lock<std::mutex> lock(mtx);
    relaxRequests(heavyRequests, graph, dist, buckets, delta);
    lock.unlock();
}

// Function to parallelize delta stepping
std::vector<double> delta_stepping_parallel(int source, const Graph& graph, int delta, bool print_dist) {
    int n = graph.size();
    std::vector<double> dist(n, INT_MAX);
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b); //int(ceil(INT_MAX / delta)) + 1);

    dist[source] = 0;
    buckets[0].push_back(source); // Insert source node with distance 0

    std::mutex mtx;

    std::vector<std::thread> threads;

    for (int i = 0; i < buckets.size(); ++i) {
        while (!buckets[i].empty()) {   
            std::list<int> R = move(buckets[i]);
            for (int u : R) {
                for (const auto& e : graph.get_adjacent(u)) {
                    relax(u, e.first, e.second, dist, buckets, delta);
                }
            }
            // Handle reinsertion here using findRequests and relaxRequests
            auto lightRequests = findRequests(R, graph, delta, dist, true);
            relaxRequests(lightRequests, graph, dist, buckets, delta);

            auto heavyRequests = findRequests(R, graph, delta, dist, false);
            relaxRequests(heavyRequests, graph, dist, buckets, delta);
        }
    }

        threads.emplace_back(std::thread(relaxLightEdges, std::ref(R), std::ref(graph), std::ref(dist), std::ref(buckets), delta, std::ref(mtx)));
        threads.emplace_back(std::thread(relaxHeavyEdges, std::ref(R), std::ref(graph), std::ref(dist), std::ref(buckets), delta, std::ref(mtx)));

        for (auto& t : threads) {
            t.join();
        }
        threads.clear();
    
  
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
*/


// Utility function to perform relaxation in parallel
template <typename T>
void parallel_relax(const std::list<int>& R, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta, std::mutex& mutex) {
    
    for (int u : R) {
        for (const auto& e : graph.get_adjacent(u)) {
            T newDist = dist[u] + e.weight; //distance through neighbor
            if (newDist < dist[e.dest]) {
                std::lock_guard<std::mutex> lock(mutex); // Lock for critical section
                if (newDist < dist[e.dest]) { // Double-checked locking
                    if (dist[e.dest] != INT_MAX) { //new shortest path was found
                        int oldBucketIdx = static_cast<int>(dist[e.dest] / delta);
                        buckets[oldBucketIdx].remove(e.dest); // remove v from its current bucket
                    }
                    dist[e.dest] = newDist;  // update new distance
                    int newBucketIdx = static_cast<int>(newDist / delta);
                    buckets[newBucketIdx].push_back(e.dest); // assign v to a new bucket
                }
            }
        }
    }
    return;
}


template <typename T>
std::vector<T> delta_stepping_parallel(int source, const Graph<T>& graph, int delta, int num_threads, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INF);
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b);
    std::mutex mutex;
    
    dist[source] = 0;
    buckets[0].push_back(source);  

    std::vector<std::thread> threads(num_threads);
    for (int i = 0; i < buckets.size(); ++i) {
        while (!buckets[i].empty()) {
            std::list<int> R = std::move(buckets[i]);

            // Split R into almost equal parts for each thread
            size_t part_size = R.size() / num_threads;
            auto it = R.begin();
            for (int t = 0; t < num_threads; ++t) {
                auto start = it;
                std::advance(it, t == num_threads - 1 ? std::distance(it, R.end()) : part_size);
                std::list<int> sub_R(start, it);
                //threads[t] = std::thread(parallel_relax, std::cref(sub_R), std::cref(graph), std::ref(dist), std::ref(buckets), delta, std::ref(mutex));
                threads[t] = std::thread([&]() {
                    parallel_relax(sub_R, graph, dist, buckets, delta, mutex);
                }); //use lambda function to handle template 
            }

            for (auto& thread : threads) {
                if (thread.joinable()) {
                    thread.join();
                }
            }
        }
    }

    if (print_dist) {
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


//__________________________________________________________________________________________________________________________________________________________________

// ATTEMPT 1 : ton code Clémence

// void relaxRequestsAux(const std::set<int>& requests, std::mutex &tentMutex, const Graph& graph, std::vector<double>& dist, std::vector<std::list<int>>& buckets, int delta){
//     tentMutex.lock();
//     for (int v : requests) {
//         for (const auto& e : graph.get_adjacent(v)) {
//             relax(v, e.first, e.second, dist, buckets, delta);
//         }
//     }
//     tentMutex.unlock();
// }



// std::vector<double> parDeltaStepping(int source, const Graph& graph, double delta, int num_threads, bool print_dist) {
//     int n = graph.size();
//     std::vector<std::atomic<double>> dist(n);
//     int b = graph.nb_buckets(delta);
//     std::vector<std::priority_queue<int>> buckets(b);
//     std::vector<std::thread> threads(num_threads);

//     // Initialize distances
//     for (int i = 0; i < n; ++i) {
//         double dist[i];
//         if (i==source){
//             dist[i] = 0.;
//         } else{
//             dist[i] =  std::numeric_limits<double>::max();
//         }
//     }
//     buckets[0].push(source);

//     int current_bucket = 0;
//     while (current_bucket < b) {
//         if (buckets[current_bucket].empty()) {
//             ++current_bucket;
//             continue;
//         }

//         std::vector<std::pair<int, double>> requests;
//         // Collect all requests
//         while (!buckets[current_bucket].empty()) {
//             int u = buckets[current_bucket].top();
//             buckets[current_bucket].pop();
//             for (const auto& e : graph.get_adjacent(u)) {
//                 requests.emplace_back(u, e.first);
//             }
//         }

//         // Divide requests among threads
//         int per_thread = (requests.size() + num_threads - 1) / num_threads;
//         for (int i = 0; i < num_threads && !requests.empty(); ++i) {
//             int start_idx = i * per_thread;
//             int end_idx = std::min(start_idx + per_thread, static_cast<int>(requests.size()));
//             threads[i] = std::thread(relaxRequestsAux, std::vector<std::pair<int, double>>(requests.begin() + start_idx, requests.begin() + end_idx), dist.data(), buckets.data(), std::ref(delta), std::ref(current_bucket));
//         }

//         // Wait for threads to finish
//         for (auto& t : threads) {
//             if (t.joinable()) {
//                 t.join();
//             }
//         }
//     }

//     std::vector<double> final_dists(n);
//     for (int i = 0; i < n; ++i) {
//         final_dists[i] = dist[i];
//     }

//     // Print the distances
//     if (print_dist){
//         for (int i = 0; i < n; ++i) {
//             std::cout << "Distance from " << source << " to " << i << " is ";
//             if (dist[i] == INF){
//                 std::cout << "infinity" << std::endl;
//             }
//             else{
//                 std::cout << dist[i] << std::endl;
//             }
//         }
//     }
//     return final_dists;
    
// }




// // PREVIOUS CODE THAT WORKED BEFORE (WITH OLD DELTA-STEPPING)
// void delta_stepping_worker(int source, const Graph& graph, int delta, int start_idx, int end_idx, std::vector<int>& dist, std::vector<std::list<int>>& buckets, std::mutex& mtx) {
//     while (start_idx < end_idx) {
//         while (!buckets[start_idx].empty()) {
//             int u = buckets[start_idx].front();
//             buckets[start_idx].pop_front();

//             for (const auto& e : graph.get_adjacent(u)) {
//                 int v = e.first;
//                 int weight = e.second;
//                 int newDist = dist[u] + weight;

//                 // Thread-safe update using mutex
//                 mtx.lock();
//                 if (newDist < dist[v]) {
//                     if (dist[v] != INF) {
//                         buckets[dist[v] / delta].remove(v);
//                     }
//                     dist[v] = newDist;
//                     buckets[newDist / delta].push_back(v);
//                 }
//                 mtx.unlock();
//             }
//         }
//         start_idx++;
//     }
// }

// std::vector<double> delta_stepping_threads(int source, const Graph& graph, int delta, int num_threads, bool print_dist) {
//     int n = graph.size();
//     std::vector<double> dist(n, INF);
//     // Divide vertices into buckets based on distance ranges
//     std::vector<std::list<int>> buckets((INF / delta) + 1);
//     dist[source] = 0;
//     buckets[0].push_back(source);

//     // Mutex for thread safety when accessing shared data
//     std::mutex mtx;

//     // Divide work among threads
//     int bucket_per_thread = (buckets.size() + num_threads - 1) / num_threads;
//     std::vector<std::thread> workers(num_threads);

//     int start_idx = 0;
//     for (int i = 0; i < num_threads; ++i) {
//         int end_idx = std::min(start_idx + bucket_per_thread, (int)buckets.size());
//         workers[i] = std::thread(delta_stepping_worker, source, std::cref(graph), delta, start_idx, end_idx, std::cref(dist), std::cref(buckets), std::ref(mtx));
//         start_idx = end_idx;
//     }


//     // Wait for all threads to finish
//     for (int i = 0 ; i<num_threads; i++) {
//         workers[i].join();
//     }

//     // Print the distances
//     if (print_dist){
//     for (int i = 0; i < n; ++i) {
//         std::cout << "Distance from " << source << " to " << i << " is ";
//         if (dist[i] == INF){
//             std::cout << "infinity" << std::endl;
//         }
//         else{
//             std::cout << dist[i] << std::endl;
//         }
//     }}
// 
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

    double mse = 0;
    int diff_count = 0;
    for (size_t i = 0; i < dist1.size(); i++) {
        if (dist1[i] != dist2[i]) {
            diff_count++;
            if (dist1[i] != INT_MAX && dist2[i] != INT_MAX) { 
                mse += pow(dist1[i] - dist2[i], 2);
            }
        }
    }

    if (diff_count > 0) {
        mse /= diff_count; 
    }

    // Print MSE and nb of different values
    std::cout << "Mean Squared Error: " << mse << std::endl;
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
        dist_delta_stepping_threads = delta_stepping_parallel(0, g, delta, num_threads, print_dist);
        std::chrono::steady_clock::time_point end_delta_stepping_threads = std::chrono::steady_clock::now();
        t3 = std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping_threads - begin_delta_stepping_threads).count() ;
        std::cout << "Total time with Delta stepping threads: " << t3 << " ms" << std::endl;

    }

    //comparison of algorithms (results & runtimes): 

    if (run_algo == 4 || run_algo == 0){// compare dijkstra & delta-stepping: 
        std::cout << "\n\nComparing Dijkstra / delta-stepping "<< std::endl;
        double speed_up = t1/t2; 
        std::cout << "\n\nSpeed up Dijkstra / delta-stepping: " << speed_up << std::endl;
        compare_distances(dist_dijkstra, dist_delta_stepping);
    }

    if (run_algo == 5 || run_algo == 0){// compare dijkstra & delta-stepping threads:
        std::cout << "\n\nComparing Dijkstra / delta-stepping threads "<< std::endl; 
        double speed_up = t1/t3; 
        std::cout << "\n\nSpeed up Dijkstra / delta-stepping threads: " << speed_up << std::endl;
        compare_distances(dist_dijkstra, dist_delta_stepping_threads);
    }

    if (run_algo == 6 || run_algo == 0){// compare delta-stepping & delta-stepping threads: 
        std::cout << "\n\nComparing delta-stepping / delta-stepping threads "<< std::endl;
        double speed_up = t2/t3; 
        std::cout << "\n\nSpeed up delta-stepping / delta-stepping threads: " << speed_up << std::endl;
        compare_distances(dist_delta_stepping, dist_delta_stepping_threads);
    }

    return 0;
}





int main() {
    int type_graph = 2;  // 0 for small graph, 1 for txt graph, 2 for random graph 
    std::string name_of_txt = "txt_graph.txt";


    int run_algo = 4; // 1 : dijkstra ; 2 : delta-stepping ; 3 : DS-threads ; 4 : compare dijkstra & DS ; 5 : compare dijkstra & DS threads ; 6 : compare DS & DS threads ; 0 : compare all 
    
    int delta = 5; 
    int num_threads = 10;

    // for generating a random graph, type_graph = 2 
    int num_vertices = 20;
    int num_edges =  (num_vertices * (num_vertices - 1)) ;
    int min_weight = 1; //positive weights
    int max_weight = 10;

    bool print_dist = 1; // if want to print the resulting distances or not, it affects the running time so put 0 preferably
    bool print_graph = 0; // Whether or not want to print the graph
    
    

    if (type_graph == 0){
        std::cout << "\nGenerating a small graph\n";
        if (run_algo == 2 || run_algo == 3 || run_algo ==  6 ){ // creating a graph with real-valued positive edge weights (run_algo does not call dijkstra)
            Graph<double> g(6); //change with the nb of vertices in the fixed graph
            g.gen_small_graph_real();
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        } 
        else {
            Graph<int> g(6);
            g.gen_small_graph_int();
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads); 
        }
    }

    else if (type_graph == 1){
        std::cout << "\nGenerating a small graph via text file\n";
        if (run_algo == 2 || run_algo == 3 || run_algo ==  6 ){ // creating a graph with real-valued positive edge weights (run_algo does not call dijkstra)
            Graph<double> g(999);
            g.gen_graph_from_txt(name_of_txt);
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        } 
        else {
            Graph<int> g(999); //change with the nb of vertices in the txt graph
            g.gen_graph_from_txt(name_of_txt);
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        }
    }

    else if (type_graph == 2){
        std::cout << "\nGenerating a random graph\n";
        if (run_algo == 2 || run_algo == 3 || run_algo ==  6 ){ // creating a graph with real-valued positive edge weights (run_algo does not call dijkstra)
            Graph<double> g(num_vertices);
            g.gen_random_graph(num_vertices, num_edges, min_weight, max_weight);
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        } 
        else {
            Graph<int> g(num_vertices);
            g.gen_random_graph(num_vertices, num_edges, min_weight, max_weight);
            return continue_main(g, run_algo, delta, print_graph, print_dist, num_threads);
        }
    }

}

