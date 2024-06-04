
#include "graph.h"

// relax node u with new dist x
template <typename T>
void relax(int w, T x, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    if (x < dist[w]) { 
        if (dist[w] != INT_MAX) {  //new shortest path was found
            int oldBucketIdx = static_cast<int>(std::floor(dist[w] / delta));
            buckets[oldBucketIdx].remove(w); // remove v from its current bucket
        }
        int newBucketIdx = static_cast<int>(std::floor(x / delta));
        buckets[newBucketIdx].push_back(w); // assign v to a new bucket
        dist[w] = x;  // update new distance
    }
}


// creates requests for the edges of type isLight and of the nodes in R
// returns a set of edges
template <typename T>
std::vector<Edge<T>> findRequests(const std::list<int>& Vprime, const Graph<T>& graph, int delta, const std::vector<T>& dist, bool isLight) {
    std::vector<Edge<T>> requests;
    for (int v : Vprime) {
        for (const auto& e : graph.get_adjacent(v)) { 
            int w = e.dest; 
            T weight = e.weight;  // weight of edge from v to w 
            if ((isLight && weight <= delta) || (!isLight && weight > delta)) { // checks if e has the correct type : light or heavy
                Edge<T> e(w, dist[v] + weight);
                requests.push_back(e);
            }
        }
    }
    return requests;
}

// do relaxations, may move nodes between buckets
template <typename T>
void relaxRequests(const std::vector<Edge<T>>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    for (auto req : requests) {
        int w = req.dest;
        T x = req.weight;
        relax(w, x, dist, buckets, delta);
    }
}


int find_smallest_non_empty_bucket(const std::vector<std::list<int>>& buckets) {
  int min_size = INT_MAX; 
  int min_bucket = -1;

  for (size_t i = 0; i < buckets.size(); ++i) {
    if (!buckets[i].empty()) {
      int size = buckets[i].size();
      if (size < min_size) {
        min_size = size;
        min_bucket = i;
      }
    }
  }
  return min_bucket;
}


template <typename T>                                                        //  for reference : fig 1 p. 123 of paper 
std::vector<T> delta_stepping(int source, const Graph<T>& graph, int delta, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INT_MAX);                                         // line 1 
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b); 

    dist[source] = 0;                                                        // line 2
    buckets[0].push_back(source); // Insert source node with distance 0

    while (true){                                                           // line 3
        int i = find_smallest_non_empty_bucket(buckets);                    // line 4
        if (i == -1) break;                                                 
        std::list<int> R ;                                                  // line 5
        while (!buckets[i].empty()) {                                       // line 6
            std::vector<Edge<T>> lightRequests = findRequests(buckets[i], graph, delta, dist, true); // line 7
            // move elements of buckets[i] to the end of R, empties buckets[i] : 
            R.splice(R.end(), buckets[i])  ;                                // line 8 & 9
            relaxRequests(lightRequests, graph, dist, buckets, delta);      // line 10
        }        
        // Relax heavy edges
        std::vector<Edge<T>> heavyRequests = findRequests(R, graph, delta, dist, false);  // line 11
        relaxRequests(heavyRequests, graph, dist, buckets, delta);        // line 12
    }

    // Print the distances
    if (print_dist){
        for (int i = 0; i < n; ++i) {
            std::cout << "Distance from " << source << " to " << i << " is ";
            if (dist[i] == INT_MAX){
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
void relaxPar(int u, T x, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta, std::mutex& m) {
    std::unique_lock<std::mutex> lock(m);
    relax(u, x, dist, buckets, delta);
    lock.unlock();
}



// edge relaxation for an entire bucket can be done in parallel so we use threads
// we use multiple threads to process the requests and a mutex to manage shared resources safely

template <typename T>
void relaxRequestsPar(const std::vector<Edge<T>>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta, int nb_threads) {
    // mutex to protect concurrent access to buckets
    std::mutex mutex;

    // divide requests into chunks for each thread
    int chunk_size = (requests.size() + nb_threads - 1) / nb_threads;
    std::vector<std::thread> threads(nb_threads);

    // Thread function to process a chunk of requests
    auto relax_chunk = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
        auto req = requests[i];
        int w = req.dest;
        T x = req.weight;
        relaxPar(w, x, dist, buckets, delta, mutex); 
        }
    };

    for (int i = 0; i < nb_threads; ++i) {
        int start = i * chunk_size;
        int end = std::min(start + chunk_size, (int)requests.size());
        threads[i] = std::thread(relax_chunk, start, end);
    }

    for (auto& thread : threads) {
        thread.join();
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



// Version that parallelize less but seems to run in less time
template <typename T>
std::vector<T> delta_stepping_Par(int source, const Graph<T>& graph, int delta, int nb_threads, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INT_MAX);                                         // line 1 
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b); 

    dist[source] = 0;                                                        // line 2
    buckets[0].push_back(source); // Insert source node with distance 0

    while (true){                                                           // line 3
        int i = find_smallest_non_empty_bucket(buckets);                    // line 4
        if (i == -1) break;                                                 
        std::list<int> R ;                                                  // line 5
        while (!buckets[i].empty()) {                                       // line 6
            std::vector<Edge<T>> lightRequests = findRequests(buckets[i], graph, delta, dist, true); // line 7
            // move elements of buckets[i] to the end of R, empties buckets[i] : 
            R.splice(R.end(), buckets[i])  ;                                // line 8 & 9
            relaxRequestsPar(lightRequests, graph, dist, buckets, delta, nb_threads);      // line 10
        }        
        // Relax heavy edges
        std::vector<Edge<T>> heavyRequests = findRequests(R, graph, delta, dist, false);  // line 11
        relaxRequestsPar(heavyRequests, graph, dist, buckets, delta, nb_threads);        // line 12
    }

    // Print the distances
    if (print_dist){
        for (int i = 0; i < n; ++i) {
            std::cout << "Distance from " << source << " to " << i << " is ";
            if (dist[i] == INT_MAX){
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

