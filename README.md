# CSE305-Project-Parallel-shortest-paths

- Authors:
Clémence Mottez
Garance Perrot

- Introduction
Project as part of our Concurrent and Distributed Computing class at Ecole Polytechnique.
Finding shortest distances in a graph is one the fundamental problems in computer science with
numerous applications (route planning by CityMapper and Google/Yandex-maps, for example).
We have already seen classical algorithms (BFS, Dijkstra, etc) for this task in our algorithms course.
However, these algorithms are inherently sequential and hard to parallelize (although parallele versions exist).

- Aim
In this project, we implement and benchmark one of the most standard shortest path algorithms, ∆-stepping algorithm.
We implement the algorithm, run it on a set of benchmarks, and compare it with non-parallel algorithms including single-thread ∆-stepping and Dijkstra.

- How to run the project
