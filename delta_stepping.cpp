
#include "graph.h"

// relax node u with new dist x
template <typename T>
void relax(int w, T x, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    if (x < dist[w]) {                                                        // if x < tent(w) then
        if (dist[w] != INT_MAX) {                                             // if in : 
            int oldBucketIdx = static_cast<int>(std::floor(dist[w] / delta)); // index of old bucket
            buckets[oldBucketIdx].remove(w);                                  // remove w from its current bucket  
        }
        int newBucketIdx = static_cast<int>(std::floor(x / delta));           // index of new bucket
        buckets[newBucketIdx].push_back(w);                                   // insert w into new bucket
        dist[w] = x;                                                          // update new distance
    }
}


// creates requests for the edges of type isLight and of the nodes in R
// returns a set of edges
template <typename T>
std::vector<Edge<T>> findRequests(const std::list<int>& Vprime, const Graph<T>& graph, int delta, const std::vector<T>& dist, bool isLight) {
    std::vector<Edge<T>> requests; 
    for (int v : Vprime) {                                   // v in V'
        for (const auto& e : graph.get_adjacent(v)) {        // get (v, w)
            int w = e.dest; 
            T weight = e.weight;                             // c(v, w) which is different from c(w,v), that is what caused problems earlier !!
            if ((isLight && weight <= delta) || (!isLight && weight > delta)) { // (v,w) in E_kind (edge of the good type, light or heavy)
                Edge<T> e(w, dist[v] + weight);              // tent(v) + c(v,w)
                requests.push_back(e);                       // put (w, newDist) in set of requests
            }
        }
    }
    return requests;                                        // return {(w, tent(v)+c(v,w)) | v in V' ∧ (v,w) in E_kind}
}

// do relaxations, may move nodes between buckets
template <typename T>
void relaxRequests(const std::vector<Edge<T>>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta) {
    for (auto req : requests) {            // for each (w ,x) in Req 
        int w = req.dest;
        T x = req.weight;
        relax(w, x, dist, buckets, delta);  // do relax(w, x)
    }
}

// helper function to find the smallest non-empty bucket for delta-stepping
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
    std::vector<T> dist(n, INT_MAX);                                        // line 1: for each v in V do tent(v) = inf 
    int b = graph.nb_buckets(delta);
    std::vector<std::list<int>> buckets(b); 

    dist[source] = 0;                                                       // line 2
    buckets[0].push_back(source); // Insert source node with distance 0     // line 2 : relax(s,0)

    while (true){                                                           // line 3 : while !isEmpty(B) do
        int i = find_smallest_non_empty_bucket(buckets);                    // line 4 : i := min {j>= 0 | B[j] != ∅}
        if (i == -1) break;                                                 // line 3 
        std::list<int> R ;                                                  // line 5 : R := ∅
        while (!buckets[i].empty()) {                                       // line 6 : while B[i] != ∅ do
            std::vector<Edge<T>> lightRequests = findRequests(buckets[i], graph, delta, dist, true); // line 7 : Req := findRequests(B[i], light)
            // move elements of buckets[i] to the end of R, empties buckets[i] : 
            R.splice(R.end(), buckets[i])  ;                                // line 8 & 9 : R := R ∪ B[i], B[i] := ∅
            relaxRequests(lightRequests, graph, dist, buckets, delta);      // line 10 : relaxRequests(Req)
        }        
        // Relax heavy edges
        std::vector<Edge<T>> heavyRequests = findRequests(R, graph, delta, dist, false);  // line 11 : Req := findRequests(R, heavy)
        relaxRequests(heavyRequests, graph, dist, buckets, delta);          // line 12 : relaxRequests(Req)
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



// Previous version for relaxRequestsPar
// ONE THREAD PER REQUEST TO RELAX : OK FOR SMALL GRAPHS BUT OTHERWIZE TOO LONG
// edge relaxation for an entire bucket can be done in parallel so we use threads
// template <typename T>
// void relaxRequestsPar(const std::vector<Edge<T>>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::list<int>>& buckets, int delta, int nb_threads) {
//     std::vector<std::thread> threads;
//     std::mutex mutex;
//     for (auto req : requests) {
//         threads.push_back(std::thread([&, req]() {
//             relaxPar(req.dest, req.weight, dist, buckets, delta, mutex);
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

