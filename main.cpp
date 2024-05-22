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


const int INF = INT_MAX;

typedef std::pair<int, double> Edge; // an edge from a node (destination, weight)


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

    // computes the sufficient nb of buckets to use
    int nb_buckets(int delta) const {   
        int max_bucket = 0;
        for (const auto& list : adj_list) {
            for (const auto& edge : list) {
                int bucket = static_cast<int>(std::ceil(edge.second / delta));
                if (bucket > max_bucket) {
                    max_bucket = bucket;
                }
            }
        }
        return max_bucket;
    }

    // Access the adjacency list of a vertex
    const std::list<Edge>& get_adjacent(int u) const {
        return adj_list[u];
    }

    std::vector<std::list<Edge>>  get_adj_list() const
    {
        return adj_list;
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

    // Generate a random graph with num_vertices and num_edges and weights between min and max_weight
    void gen_random_graph(int num_vertices, int num_edges, int min_weight, int max_weight) {
        std::mt19937 rng(static_cast<unsigned int>(time(nullptr))); 
        std::uniform_int_distribution<> vert_dist(0, num_vertices - 1);
        std::uniform_int_distribution<> weight_dist(min_weight, max_weight);


        for (int i = 0; i < num_edges; i++) {
            int u = vert_dist(rng);
            int v = vert_dist(rng);
            double weight = weight_dist(rng);
            add_edge(u, v, weight);
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
                std::cout << " -> " << edge.first << " with weight " << edge.second << "\n";
            }
        }
    }


    
};

//_________________________________________________________________________________________________________________________


//called by relaxRequests
// updates the tentative distance of node v (of neighbor u) if a shorter path is discovered
void relax(int u, int v, double weight, std::vector<double>& dist, std::vector<std::list<int>>& buckets, int delta) {
    double newDist = dist[u] + weight; //distance through neighbor
    if (newDist < dist[v]) { 
        if (dist[v] != INT_MAX) {  //new shortest path was found
            int oldBucketIdx = dist[v] / delta;
            buckets[oldBucketIdx].remove(v); // remove v from its current bucket
        }
        dist[v] = newDist;  // update new distance
        int newBucketIdx = newDist / delta;
        buckets[newBucketIdx].push_back(v); // assign v to a new bucket
    }
}


// creates requests for the edges of type isLight and of the nodes in R
// returns a set of edges
std::set<int> findRequests(const std::list<int>& R, const Graph& graph, int delta, const std::vector<double>& dist, bool isLight) {
    std::set<int> requests;
    for (int u : R) {
        for (const auto& e : graph.get_adjacent(u)) {
            int v = e.first; 
            double weight = e.second;
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
void relaxRequests(const std::set<int>& requests, const Graph& graph, std::vector<double>& dist, std::vector<std::list<int>>& buckets, int delta) {
    for (int v : requests) {
        for (const auto& e : graph.get_adjacent(v)) {
            relax(v, e.first, e.second, dist, buckets, delta);
        }
    }
}


std::vector<double> delta_stepping(int source, const Graph& graph, int delta, bool print_dist) {
    int n = graph.size();
    std::vector<double> dist(n, INT_MAX);
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b); //int(ceil(INT_MAX / delta)) + 1);

    dist[source] = 0;
    buckets[0].push_back(source); // Insert source node with distance 0

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





// ATTEMPT 2 : 



// void relax_thread_function(int thread_id, const Graph& graph, std::vector<double>& dist, std::vector<std::list<int>>& buckets, int delta, std::atomic<bool>& stop) {
//     while (!stop.load()) {
//         // Get a bucket (round-robin style)
//         int bucket_idx = (thread_id + buckets.size()) % buckets.size();

//         // Process elements in the bucket while it's not empty and stop signal not received
//         while (!buckets[bucket_idx].empty() && !stop.load()) {
//             std::list<int> R = move(buckets[bucket_idx]);

//             // Find light and heavy requests (can be parallelized further)
//             auto lightRequests = findRequests(R, graph, delta, dist, true);
//             auto heavyRequests = findRequests(R, graph, delta, dist, false);

//             // Relax edges within the requests (can be parallelized further)
//             relaxRequests(lightRequests, graph, dist, buckets, delta);
//             relaxRequests(lightRequests, graph, dist, buckets, delta);

//         }
//     }
// }

// std::vector<double> parDeltaStepping(int source, const Graph& graph, int delta, int num_threads, bool print_dist) {
//     int n = graph.size();
//     std::vector<double> dist(n, INT_MAX);
//     int b = graph.nb_buckets(delta);
//     std::vector<std::list<int>> buckets(b);

//     dist[source] = 0;
//     buckets[0].push_back(source); // Insert source node with distance 0

//     std::vector<std::thread> threads(num_threads);
//     std::atomic<bool> stop(false); // Flag to signal thread termination

//     // Spawn threads
//     for (int i = 0; i < num_threads; ++i) {
//         threads[i] = std::thread(relax_thread_function, i, std::ref(graph), std::ref(dist), std::ref(buckets), delta, std::ref(stop));
//     }

//     // Main thread processing (can be used for heavy nodes)
//     bool all_buckets_empty = false;
//     while (!all_buckets_empty && !stop.load()) {
//         all_buckets_empty = true;
//         for (int i = 0; i < buckets.size(); ++i) {
//             if (!buckets[i].empty()) {
//                 all_buckets_empty = false;
//                 std::list<int> R = move(buckets[i]);
//                 for (int u : R) {
//                     for (const auto& e : graph.get_adjacent(u)) {
//                         relax(u, e.first, e.second, dist, buckets, delta);
//                     }
//                 }
//             }
//         }
//     }

//     // Signal termination to threads
//     stop.store(true);

//     for (auto& thread : threads) {
//     thread.join();
//     }

//   // Print the distances
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

//     return dist;
// }




//ATTEMPT 3 : 



// std::mutex mutex;
// void relax_atomic(int u, int v, double weight, std::vector<std::atomic<double>>& dist, std::vector<std::list<int>>& buckets, int delta) {
//     double newDist = dist[u].load(std::memory_order_relaxed) + weight;
//     if (newDist < dist[v].load(std::memory_order_relaxed)) {
//         dist[v].store(newDist, std::memory_order_release);
//         int oldBucketIdx = static_cast<int>(dist[v].load(std::memory_order_relaxed) / delta);
//         buckets[oldBucketIdx].remove_if([v](int node) { return node == v; }); 
//     }
// }

// void relaxRequestsAux(const std::vector<std::pair<int, double>>& requests, std::vector<std::atomic<double>>& dist, const Graph& graph, std::vector<std::list<int>>& buckets, int delta) {
//     mutex.lock();

//     for (const auto& [u, v] : requests) {
//         double weight = std::numeric_limits<double>::infinity(); // Initialize weight to positive infinity
//         for (const auto& edge : graph.get_adjacent(u)) {
//             if (edge.first == v) { // Found the edge between u and v
//                 weight = edge.second;
//                 break; // No need to iterate further once the edge is found
//             }
//         }
//         relax_atomic(u, v, weight, dist, buckets, delta);  // Use weight from loop or infinity
//     }

//     mutex.unlock();
// }

// std::vector<double> parDeltaStepping(int source, const Graph& graph, double delta, int num_threads, bool print_dist) {
//     int n = graph.size();
//     std::vector<std::atomic<double>> dist(n, std::numeric_limits<double>::max());
//     int b = graph.nb_buckets(delta);
//     std::vector<std::list<int>> buckets(b); // Use lists for buckets (consistent with sequential version)

//     // Initialize distances
//     dist[source] = 0.;
//     buckets[0].push_back(source);

//     int current_bucket = 0;
//     while (current_bucket < b) {
//         if (buckets[current_bucket].empty()) {
//             ++current_bucket;
//             continue;
//         }

//         std::vector<std::pair<int, double>> requests;
//         // Collect all requests
//         while (!buckets[current_bucket].empty()) {
//             int u = buckets[current_bucket].front();
//             buckets[current_bucket].pop_front();
//             for (const auto& e : graph.get_adjacent(u)) {
//                 requests.emplace_back(u, e.first);
//             }
//         }

//         // Divide requests among threads
//         int per_thread = (requests.size() + num_threads - 1) / num_threads;
//         std::vector<std::thread> threads(num_threads);
//         for (int i = 0; i < num_threads && !requests.empty(); ++i) {
//             int start_idx = i * per_thread;
//             int end_idx = std::min(start_idx + per_thread, static_cast<int>(requests.size()));
//             threads[i] = std::thread(relaxRequestsAux, std::vector<std::pair<int, double>>(requests.begin() + start_idx, requests.begin() + end_idx), std::ref(dist),std::ref(graph), std::ref(buckets), std::ref(delta));
//         }

//         // Wait for threads to finish
//         for (auto& t : threads) {
//             if (t.joinable()) {
//                 t.join();
//             }
//         }

//         // Check if all buckets are empty
//         bool all_empty = true;
//         for (int i = 0; i < b; ++i) {
//                 if (!buckets[i].empty()) {
//                 all_empty = false;
//                 break;
//             }
//         }
//         if (all_empty) break;

//         ++current_bucket;
//     }

//     std::vector<double> final_dists(n);
//     for (int i = 0; i < n; ++i) {
//         final_dists[i] = dist[i];
//     }

//     // Print the distances
//     if (print_dist) {
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






//______________________________________________________________________________________________________________________________


std::vector<double> dijkstra(int source, Graph& graph, bool print_dist) {
    int n = graph.size();
    std::vector<double> dist(n, INT_MAX); 
    std::vector<bool> visited(n, false); 
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    // Initialize priority queue with the source node
    pq.push({0, source}); // (distance, vertex)
    dist[source] = 0;

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        // If vertex already been visited, continue to the next
        if (visited[u]) continue;
        visited[u] = true;

        // Check each adjacent vertex of u
        for (const auto& edge : graph.get_adjacent(u)) {
            int v = edge.first;
            int weight = edge.second;

            // Only consider this vertex if it has not been visited
            if (!visited[v] && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pq.push({dist[v], v}); // Push updated distance and vertex
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
void compare_distances(const std::vector<double>& dist1, const std::vector<double>& dist2) {
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





int main() {
    int type_graph = 1;  // 0 for small graph, 1 for txt graph, 2 for random graph 
    std::string name_of_txt = "txt_graph.txt";
    int run_all_algo = 0; // 0 run both 1 and 2, 1 run dijkstra, 2 run delta stepping, 3 run delta-stepping w/ threads
    int nb_vertices = 1000; // To change
    int delta = 10; 
    int num_threads = 1;

    bool print_dist = 0; // if want to print the resulting distances or not, it affects the running time so put 0 preferably
    bool print_graph = 0; // Whether or not want to print the graph
    
    Graph g(nb_vertices);

    if (type_graph == 0){
        std::cout << "\nGenerating a small graph:\n";
        g.gen_small_graph(); 
    }
    else if (type_graph == 1){
        std::cout << "\nGenerating a small graph via text file:\n";
        //Graph g(6); 
        g.gen_graph_from_txt(name_of_txt);
    }
    else if (type_graph == 2){
        std::cout << "\nGenerating a random graph:\n";
        g.gen_random_graph(10000, 100000, 1, 100); // Generate a random graph with 5 vertices and 10 edges
        //g.gen_random_graph(25, 50, 1, 100);
    }

    if (print_graph){
        g.print_graph(); 
    }
    std::chrono::steady_clock::time_point begin_dijkstra = std::chrono::steady_clock::now();
    std::vector<double> dist_dijkstra;
    if (run_all_algo!=2 && run_all_algo!=3){
        // Run dijkstra algo
        std::cout << "\nResults with dijkstra algo\n";
        dist_dijkstra = dijkstra(0, g, print_dist);
    }
    std::chrono::steady_clock::time_point end_dijkstra = std::chrono::steady_clock::now();
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_dijkstra - begin_dijkstra).count() << " ms" << std::endl;

    std::chrono::steady_clock::time_point begin_delta_stepping = std::chrono::steady_clock::now();
    std::vector<double> dist_delta_stepping;
    if (run_all_algo!=1 && run_all_algo!=3){
        // Run delta stepping algo
        std::cout << "\nResults with delta stepping algo";
        std::cout << "\nDelta = " << delta << std::endl;
        dist_delta_stepping = delta_stepping(0, g, delta, print_dist);
    }
    std::chrono::steady_clock::time_point end_delta_stepping = std::chrono::steady_clock::now();
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping - begin_delta_stepping).count() << " ms" << std::endl;
    
    // std::chrono::steady_clock::time_point begin_delta_stepping_threads = std::chrono::steady_clock::now();
    // std::vector<double> dist_delta_stepping_threads;
    // if (run_all_algo!=1 && run_all_algo!=2){
    //     // Run delta stepping algo
    //     std::cout << "\nResults with delta stepping w/ threads algo";
    //     std::cout << "\nDelta = " << delta << std::endl;
    //     std::cout << "\nNb of threads = " << num_threads << std::endl;
    //     dist_delta_stepping_threads = parDeltaStepping(0, g, delta, num_threads, print_dist);
    // }
    // std::chrono::steady_clock::time_point end_delta_stepping_threads = std::chrono::steady_clock::now();
    // std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping_threads - begin_delta_stepping_threads).count() << " ms" << std::endl;


    if (run_all_algo == 0){
        std::cout << "\n\nComparison 3 algorithms";
        double t1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_dijkstra - begin_dijkstra).count();
        double t2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping - begin_delta_stepping).count();
        // double t3 = std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping_threads - begin_delta_stepping_threads).count();
        double speed_up_1 = t1/t2; 
        // double speed_up_2 = t1/t3; 
        // double speed_up_3 = t2/t3; 

        std::cout << "\n\nTime Dijkstra: "<<t1 << "ms" << std::endl;
        std::cout << "Time delta stepping: "<<t2 << "ms" << std::endl;
        // std::cout << "Time delta stepping threads: "<<t3 << "ms" << std::endl;

        std::cout << "\n\nSpeed up Dijkstra/delta-stepping: " << speed_up_1 << std::endl;
        // std::cout << "Speed up Dijkstra/delta-stepping-threads: " << speed_up_2 << std::endl;
        // std::cout << "Speed up delta-stepping/delta-stepping-threads: " << speed_up_3 << std::endl;
        
        std::cout << "\n\nCompare results Dijkstra/delta-stepping "<< std::endl;
        compare_distances(dist_dijkstra, dist_delta_stepping);
        // std::cout << "\nCompare results delta-stepping/delta-stepping-threads "<< std::endl;
        // compare_distances(dist_delta_stepping, dist_delta_stepping_threads);
    }

    return 0;
}