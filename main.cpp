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

using namespace std;

typedef std::pair<int, double> Pii; // (vertex, weight)


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
                std::cout << " -> " << edge.to << " with weight " << edge.weight << "\n";
            }
        }
    }

    
    int nb_buckets(int delta) {
        int max_bucket = 0;
        for (const auto& list : adj_list) {
            for (const auto& edge : list) {
                int bucket = static_cast<int>(std::ceil(edge.weight / delta));
                if (bucket > max_bucket) {
                    max_bucket = bucket;
                }
            }
        }
        return max_bucket;
    }
    
};



const int INF = INT_MAX;

// std::vector<int> delta_stepping(int source, Graph& graph, int delta, bool print_dist) {
//     int n = graph.size();
//     std::vector<int> dist(n, INF);
//     std::vector<std::list<int>> buckets((INF / delta) + 1);

//     dist[source] = 0;
//     buckets[0].push_back(source);

//     int idx = 0;
//     while (idx < buckets.size()) {
//         while (!buckets[idx].empty()) {
//             int u = buckets[idx].front();
//             buckets[idx].pop_front();

//             for (const auto& e : graph.get_adjacent(u)) {
//                 int v = e.to;
//                 int weight = e.weight;
//                 int newDist = dist[u] + weight;

//                 if (newDist < dist[v]) {
//                     if (dist[v] != INF) {
//                         int oldBucketIdx = dist[v] / delta;
//                         buckets[oldBucketIdx].remove(v);
//                     }
//                     dist[v] = newDist;
//                     int newBucketIdx = newDist / delta;
//                     buckets[newBucketIdx].push_back(v);
//                 }
//             }
//         }
//         idx++;
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

//     return dist;
// }



// void delta_stepping_worker(int source, Graph& graph, int delta, int start_idx, int end_idx, std::vector<int>& dist, std::vector<std::list<int>>& buckets, std::mutex& mtx) {
//     while (start_idx < end_idx) {
//         while (!buckets[start_idx].empty()) {
//             int u = buckets[start_idx].front();
//             buckets[start_idx].pop_front();

//             for (const auto& e : graph.get_adjacent(u)) {
//                 int v = e.to;
//                 int weight = e.weight;
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

// std::vector<double> delta_stepping_threads(int source, Graph& graph, int delta, int num_threads, bool print_dist) {
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
//     std::vector<std::thread> threads(num_threads);

//     int start_idx = 0;
//     for (int i = 0; i < num_threads; ++i) {
//         int end_idx = std::min(start_idx + bucket_per_thread, (int)buckets.size());
//         threads[i] = std::thread(delta_stepping_worker, source, std::ref(graph), delta, start_idx, end_idx, std::ref(dist), std::ref(buckets), std::ref(mtx));
//         start_idx = end_idx;
//     }

//     // Wait for all threads to finish
//     for (int i = 0 ; i<num_threads; i++) {
//         threads[i].join();
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

//     return dist;
// }




// New attempt delta stepping



void findRequests(std::vector<std::pair<int, double>> &requests, std::priority_queue<int> bucket, bool light_edge, std::vector<std::vector<std::pair<int, double>>> adj_list, std::vector<double> tent, int delta)
{
    while (!bucket.empty())
    {
        int i = bucket.top();
        bucket.pop();

        for (auto edge : adj_list[i])
        {
            int j = edge.first;
            double weight = edge.second;

            if (weight > 0 && ((light_edge && weight <= delta) || (!light_edge && weight > delta)))
            {
                requests.push_back(std::make_pair(j, tent[i] + weight));
            }
        }
    }
}


void relax(int w, int x, std::vector<double> &tent, std::vector<std::priority_queue<int> > &B, int delta)
{
    if (x < tent[w])
    {
        if (tent[w] != INF) // If there is something to pop
        {
            B[std::floor(tent[w] / delta)].pop(); // remove from old bucket
        }
        B[std::floor(x / delta)].push(w); // insert in new bucket
        tent[w] = x;
    }
}

void relaxRequests(std::vector<std::pair<int, double>> requests, std::vector<double> &tent, std::vector<std::priority_queue<int> > &B, int delta)
{
    for (int i = 0; i < requests.size(); ++i)
    {
        relax(requests[i].first, requests[i].second, tent, B, delta);
    }
}


bool isNotEmpty(std::vector<std::priority_queue<int> > bucket_vec)
{
    for (size_t i = 0; i < bucket_vec.size(); ++i)
    {
        if (!bucket_vec[i].empty())
        {
            return true;
        }
    }
    return false;
}



std::vector<double> delta_stepping(int source, Graph& graph, int delta, bool print_dist) {
    int n = graph.size();
    int b = 100; //graph.nb_buckets(delta);
    std::vector<std::priority_queue<int> > B(b); // b buckets (vectors) stored in B, priority queues
    std::vector<double> tent(n, INF);                   // tentative distances
    std::vector<std::vector<std::pair<int, double>>> adj_list; //= graph.get_adjacent();

    tent[source] = 0;
    B[0].push(source);
    int i = b;
    while (isNotEmpty(B))
    {
        for (size_t j = 0; j < b; ++j)
        {
            if (!(B[j].empty()))
            {
                i = j;
                break;
            }
        }
        std::priority_queue<int> R;
        while (!(B[i].empty()))
        {
            std::vector<std::pair<int, double>> Req_light;
            findRequests(Req_light, B[i], true, adj_list, tent, delta);
            R = B[i];
            B[i] = std::priority_queue<int>();
            relaxRequests(Req_light, tent, B, delta);
        }
        std::vector<std::pair<int, double>> Req_heavy;
        findRequests(Req_heavy, R, false, adj_list, tent, delta);
        relaxRequests(Req_heavy, tent, B, delta);
    }
    return tent;
}




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
            int v = edge.to;
            int weight = edge.weight;

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
    }}

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
    int type_graph = 2;  // 0 for small graph, 1 for txt graph, 2 for random graph 
    string name_of_txt = "./graphs/rmat_medium.txt";
    int run_both_algo = 0; // 0 run both 1 and 2, 1 run dijkstra, 2 run delta stepping, 3 run delta-stepping w/ threads
    int nb_vertices = 1000; // To change
    int delta = 2; 
    int num_threads = 10;

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
        g.gen_random_graph(1000, 10000, 1, 100); // Generate a random graph with 5 vertices and 10 edges
        //g.gen_random_graph(25, 50, 1, 100);
    }

    if (print_graph){
        g.print_graph(); 
    }
    std::chrono::steady_clock::time_point begin_dijkstra = std::chrono::steady_clock::now();
    std::vector<double> dist_dijkstra;
    if (run_both_algo!=2 && run_both_algo!=3){
        // Run dijkstra algo
        std::cout << "\nResults with dijkstra algo\n";
        dist_dijkstra = dijkstra(0, g, print_dist);
    }
    std::chrono::steady_clock::time_point end_dijkstra = std::chrono::steady_clock::now();
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_dijkstra - begin_dijkstra).count() << " ms" << std::endl;

    std::chrono::steady_clock::time_point begin_delta_stepping = std::chrono::steady_clock::now();
    std::vector<double> dist_delta_stepping;
    if (run_both_algo!=1 && run_both_algo!=3){
        // Run delta stepping algo
        std::cout << "\nResults with delta stepping algo";
        std::cout << "\nDelta = " << delta << std::endl;
        dist_delta_stepping = delta_stepping(0, g, delta, print_dist);
    }
    std::chrono::steady_clock::time_point end_delta_stepping = std::chrono::steady_clock::now();
    std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping - begin_delta_stepping).count() << " ms" << std::endl;
    
    // std::chrono::steady_clock::time_point begin_delta_stepping_threads = std::chrono::steady_clock::now();
    // std::vector<double> dist_delta_stepping_threads;
    // if (run_both_algo!=1 && run_both_algo!=2){
    //     // Run delta stepping algo
    //     std::cout << "\nResults with delta stepping w/ threads algo";
    //     std::cout << "\nDelta = " << delta << std::endl;
    //     std::cout << "\nNb of threads = " << num_threads << std::endl;
    //     dist_delta_stepping_threads = delta_stepping_threads(0, g, delta, num_threads, print_dist);
    // }
    // std::chrono::steady_clock::time_point end_delta_stepping_threads = std::chrono::steady_clock::now();
    // std::cout << "Total time : " << std::chrono::duration_cast<std::chrono::milliseconds>(end_delta_stepping_threads - begin_delta_stepping_threads).count() << " ms" << std::endl;


    if (run_both_algo == 0){
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
        std::cout << "\nCompare results delta-stepping/delta-stepping-threads "<< std::endl;
        // compare_distances(dist_delta_stepping, dist_delta_stepping_threads);
    }

    return 0;
}