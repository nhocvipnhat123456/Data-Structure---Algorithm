#ifndef GRAPH_ALGS
#define GRAPH_ALGS

#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <algorithm>
#include "weighted_graph.hpp"
#include "easy_weighted_graph_algorithms.cpp"

//Returns true if the graph is connected, false otherwise.
template <typename vertex>
bool is_connected(const weighted_graph<vertex>& g){
	int size = g.num_vertices();
	std::unordered_set<vertex> vertices = g.get_vertices();
	
	if (size == 0){	//Check if g is a graph or not.
		return true;
	} else {
		//Create a visited with an unordered map type to store both vertex and index of the vertex.
		std::unordered_map<vertex, int> visited;
		//Accessing to first vertex.
		auto begin_vertex = vertices.begin();
		visited[*begin_vertex] = 0;
		
		while (visited.size() < size){
			int count = 0;
			for( const auto& iter : visited ) {
				if (iter.second == 0) {
					//Create an iterator to store neighbours of the vertex that being checked.
					for (auto neighbour = g.neighbours_begin(iter.first); neighbour != g.neighbours_end(iter.first); neighbour++){
						auto got = visited.find(neighbour->first);
						if (got == visited.end()){
							visited[neighbour->first] = 0;
							count++;
						}
					}
				}
				//Prevent being checked again.
				visited.at(iter.first) = 1;
			}
			//If the for loop does not add anything in visited, it will break the while loop and return false.
			if (!count) break;
		}
		if (visited.size() < size) return false;
		else return true;
		
	}
	return false;
}

//Returns a vector of weighted graphs, where each weighted graph is a connected
//component of the input graph.
template <typename vertex>
std::vector<weighted_graph<vertex>> connected_components(const weighted_graph<vertex>& g){

	std::unordered_set<vertex> visited;
	std::unordered_set<vertex> unvisited = g.get_vertices();
	std::vector<weighted_graph<vertex>> result;
	
	for (auto it = unvisited.begin(); it != unvisited.end(); it++){
		//Check if the vertex is visited.
		if (visited.count(*it) == 0){
			std::vector<vertex> connected = depth_first(g, *it);
			// create connected component graph here
			weighted_graph<vertex> sub_graph;
			//Input vertex
			for (int i = 0; i < connected.size(); i++){
				sub_graph.add_vertex(connected[i]);
			}
			//Input edge
			for (int i = 0; i < connected.size(); i++){
				for (int j = 1; j < connected.size(); j++){
					if (g.are_adjacent(connected[i], connected[j])){
						sub_graph.add_edge(connected[i], connected[j], g.get_edge_weight(connected[i], connected[j]));
					}
				}
			}
			//Add sub graph into the result.
			result.push_back(sub_graph);
			
			for (auto it_connected = connected.begin(); it_connected != connected.end(); it_connected++){
				visited.insert(*it_connected);
			}
		}
	}
	
	return result;
}

//Returns a map of the vertices of the weighted graph g and their distances from
//the given starting vertex v.
template <typename vertex> 
std::map<vertex, int> dijkstras(const weighted_graph<vertex>& g, const vertex& v){
	
	std::unordered_set<vertex> labelled;
	std::map<vertex, int> length;
	std::unordered_set<vertex> vertices = g.get_vertices();

	//Insert all of the vertices into the map
	for (auto it : vertices) {
		length.insert(std::pair<vertex, int>(it, g.total_weight() + 1));
	}
	//Check if the graph is empty or not
	if (g.num_vertices() != 0){
		//Set v as the start vertex
		length[v] = 0;
		//If the unvisited is not empty then do stuffs.
		while (labelled.size() != g.num_vertices()) {
			int min = -1;
			auto selected_vertex = vertex();
			//look for the minimum distance of length
			for (auto it : vertices) {
				if (labelled.count(it) == 0 && length[it] != g.total_weight() + 1 && (length[it] < min || min == -1)) {
					min = length[it];
					selected_vertex = it;
				}
			}
			//Label the selected vertex to not run it again
			labelled.insert(selected_vertex);
			//Get all neighbours of the selected vertex
			for (auto nit = g.cneighbours_begin(selected_vertex); nit != g.cneighbours_end(selected_vertex); nit++) {
				//If the vertex is not labelled and there is an edge between the selected vertex and its neighbour vertex
				if (labelled.count(nit->first) == 0 && g.get_edge_weight(selected_vertex, nit->first) > 0) {
					auto temp_dist = length[selected_vertex] + g.get_edge_weight(selected_vertex, nit->first);
					//Look for the minimum distance of the selected vertex with its neighbours
					if (length[nit->first] == g.total_weight() + 1 || length[nit->first] > temp_dist) {
						length[nit->first] = temp_dist;
					}
				}
			}
		}
	}
	
	return length;
}

//Returns a vector containing all the articulation points of the
//input weighted graph g.
template <typename vertex>
std::vector<vertex> articulation_points(const weighted_graph<vertex>& g){
	std::vector<vertex> result;
	std::vector<weighted_graph<vertex>> g_components = connected_components(g);
	size_t size_g = g_components.size();
	std::unordered_set<vertex> vertices = g.get_vertices();
	
	for (auto it = vertices.begin(); it != vertices.end(); it++){
		weighted_graph<vertex> test;
		//Copy g to test graph.
		test.copy_graph(g);
		//Test if removing the vertex causes to increase the number of connected components
		test.remove_vertex(*it);
		std::vector<weighted_graph<vertex>> test_components = connected_components(test);
		if (test_components.size() > size_g){
			result.push_back(*it);
		}
	}
	
	return result;
}

#endif
