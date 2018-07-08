#include <iostream>
#include <vector>
#include <string>
#include "weighted_graph.hpp"

int main(){
	weighted_graph<int> g;
	int r = std::rand()%20 + 1;
	std::cout << r << std::endl;
	int re[r][r];

	for (int i = 0; i < r; ++i){
		for (int  j = 0; j < r; ++j){
			re[i][j] = 0;
		}
	}

	for (int i = 1; i < r; ++i){
		for (int j = i+1; j < r;  ++j){
			if (std::rand()%2 == 1){
				int weight = (std::rand()%10) + 1;
				re[i][j] = weight;
				re[j][i] = weight;
			}
		}
	}

	for (int i = 0; i < r; ++i){
		g.add_vertex(i);
	}
	
	for (int i = 0; i < r; ++i) {
		std::cout << g.get_vertices()[i];
	}
	
	for (int i = 0; i < r; ++i){
		for (int j = i+1; j < r;  ++j){
			if (re[i][j] > 0){
				g.add_edge(i,j,re[i][j]);
			}
		}
	}
	
	auto v = g.neighbours_begin(5);
	v++;
}
