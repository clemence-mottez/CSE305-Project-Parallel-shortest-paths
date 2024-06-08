
#include "graph.h"


// relax node u with new dist x
template <typename T>
void relax(int w, T x, std::vector<T>& dist, std::vector<std::deque<int>>& buckets, int delta) {
    if (x >= dist[w]) { // early exit if no update needed
        return;
    }
    // if x < tent(w) then
    if (dist[w] != INT_MAX) { 
        int oldBucketIdx = static_cast<int>(dist[w] / delta); // index of old bucket
        //std::cout<<"remove"<<std::endl;
        auto& oldBucket = buckets[oldBucketIdx];
        oldBucket.erase(std::remove(oldBucket.begin(), oldBucket.end(), w), oldBucket.end()); // remove w from its current bucket  
    }
    int newBucketIdx = static_cast<int>(x / delta);  // index of new bucket
    
    if (newBucketIdx >= buckets.size()) {
        buckets.resize(newBucketIdx + 1);
    }

    buckets[newBucketIdx].push_back(w);                                   // insert w into new bucket

    dist[w] = x;                                                          // update new distance
    
}


// creates requests for the edges of type isLight and of the nodes in R
// returns a set of edges
template <typename T>
std::vector<Edge<T>> findRequests(const std::deque<int>& Vprime, const Graph<T>& graph, int delta, const std::vector<T>& dist, bool isLight) {
    //std::cout<<"findres"<<std::endl;
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
void relaxRequests(const std::vector<Edge<T>>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::deque<int>>& buckets, int delta) {
    //std::cout<<"relaxreq"<<std::endl;
    for (auto req : requests) {            // for each (w ,x) in Req 
        int w = req.dest;
        T x = req.weight;
        relax(w, x, dist, buckets, delta);  // do relax(w, x)
    }
}

// helper function to find the smallest non-empty bucket for delta-stepping
int find_smallest_non_empty_bucket(const std::vector<std::deque<int>>& buckets) {
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
    std::vector<std::deque<int>> buckets(b); 

    dist[source] = 0;                                                       // line 2
    buckets[0].push_back(source); // Insert source node with distance 0     // line 2 : relax(s,0)

    while (true){                                                           // line 3 : while !isEmpty(B) do
        int i = find_smallest_non_empty_bucket(buckets);                    // line 4 : i := min {j>= 0 | B[j] != ∅}
        if (i == -1) break;                                                 // line 3 
        std::deque<int> R;                                                  // line 5 : R := ∅
        while (!buckets[i].empty()) {                                       // line 6 : while B[i] != ∅ do
            std::vector<Edge<T>> lightRequests = findRequests(buckets[i], graph, delta, dist, true); // line 7 : Req := findRequests(B[i], light)
            // move elements of buckets[i] to the end of R, empties buckets[i] : 
            R.insert(R.end(), std::make_move_iterator(buckets[i].begin()), std::make_move_iterator(buckets[i].end()));
            buckets[i].clear();                             // line 8 & 9 : R := R ∪ B[i], B[i] := ∅
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
void relaxPar2(int u, T x, std::vector<T>& dist, std::vector<std::deque<int>>& buckets, int delta, std::mutex& m) {
    std::lock_guard<std::mutex> lock(m);
    relax(u, x, dist, buckets, delta);
}




// edge relaxation for an entire bucket can be done in parallel so we use threads
// we use multiple threads to process the requests and a mutex to manage shared resources safely
template <typename T>
void relaxRequestsPar(const std::vector<Edge<T>>& requests, const Graph<T>& graph, std::vector<T>& dist, std::vector<std::deque<int>>& buckets, int delta, int nb_threads) {
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
            relaxPar2(w, x, dist, buckets, delta, mutex); 
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


template <typename T>
std::vector<Edge<T>> findRequestsPar(const std::deque<int>& Vprime, const Graph<T>& graph, int delta, const std::vector<T>& dist, bool isLight, int nb_threads) {
    
    std::vector<Edge<T>> requests;
    // Use a mutex to protect the shared requests vector
    std::mutex mutex;
    size_t chunk_size = (Vprime.size() + nb_threads - 1) / nb_threads;
    std::vector<std::thread> threads(nb_threads);

    auto process_chunk = [&](int start, int end) {
        // Create a local vector to store requests for each thread
        std::vector<Edge<T>> local_requests;
        for (int i = start; i < end; ++i) {
            
            auto it = std::next(Vprime.begin(), i);
            int v = *it;
            for (const auto& e : graph.get_adjacent(v)) {        
                int w = e.dest; 
                T weight = e.weight;                             
                if ((isLight && weight <= delta) || (!isLight && weight > delta)) { 
                    Edge<T> e(w, dist[v] + weight); 
                    local_requests.push_back(e);                     
                }
            }
        }
        // Acquire the mutex and add the local requests to the shared vector
        std::lock_guard<std::mutex> lock(mutex);
        requests.insert(requests.end(), local_requests.begin(), local_requests.end());
    };

    for (int i = 0; i < nb_threads; ++i) {
        int start = i * chunk_size;
        int end = std::min((int)(start + chunk_size), (int)Vprime.size());
        threads[i] = std::thread(process_chunk, start, end);
    } 

    for (auto& thread : threads) {
        thread.join();
    }

    return requests;
}


template <typename T>
std::vector<T> delta_stepping_Par(int source, const Graph<T>& graph, int delta, int nb_threads, bool print_dist) {
    int n = graph.size();
    std::vector<T> dist(n, INT_MAX);                                         // line 1 
    int b = graph.nb_buckets(delta);
    std::vector<std::deque<int>> buckets(b); 

    dist[source] = 0;                                                       // line 2
    buckets[0].push_back(source); // Insert source node with distance 0

    while (true){                                                           // line 3
        int i = find_smallest_non_empty_bucket(buckets);                    // line 4
        if (i == -1) break;                                                 
        std::deque<int> R;                                             // line 5
        while (!buckets[i].empty()) {                                       // line 6
            // std::vector<Edge<T>> lightRequests = findRequests(buckets[i], graph, delta, dist, true); 
            std::vector<Edge<T>> lightRequests = findRequestsPar(buckets[i], graph, delta, dist, true, nb_threads); // line 7
            
            // move elements of buckets[i] to the end of R, empties buckets[i] :   
            R.insert(R.end(), std::make_move_iterator(buckets[i].begin()), std::make_move_iterator(buckets[i].end()));
            buckets[i].clear();                             // line 8 & 9
            relaxRequestsPar(lightRequests, graph, dist, buckets, delta, nb_threads);      // line 10
        }        
        // Relax heavy edges
        // std::vector<Edge<T>> heavyRequests = findRequests(R, graph, delta, dist, false); 
        std::vector<Edge<T>> heavyRequests = findRequestsPar(R, graph, delta, dist, false, nb_threads);  // line 11
        
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
//______________________________________________________________________________________________________________________________

