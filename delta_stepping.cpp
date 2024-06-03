
#include "graph.h"

// called by relaxRequests
// updates the tentative distance of node v (of neighbor u) if a shorter path is discovered
template <typename T>
void relax(int u, int v, T weight, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    T newDist = dist[u] + weight; //distance through neighbor
    if (newDist < dist[v]) { 
        if (dist[v] != INT_MAX) {  //new shortest path was found
            int oldBucketIdx = static_cast<int>(std::floor(dist[v] / delta));
            buckets[oldBucketIdx].remove(v); // remove v from its current bucket
        }
        dist[v] = newDist;  // update new distance
        int newBucketIdx = static_cast<int>(std::floor(newDist / delta));
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

