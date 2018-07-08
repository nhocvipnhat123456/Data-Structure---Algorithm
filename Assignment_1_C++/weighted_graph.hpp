#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>

template <typename vertex>
class weighted_graph {

private:

	//A vertex list
	std::vector<vertex> vertices;
	//A matrix to store the weight of edges.
	std::vector<std::vector<int>> adjMatrix;

	//You will need to add some data members here
	//to actually represent the graph internally,
	//and keep track of whatever you need to.

	//The graph_iterator class provides an iterator
	//over the vertices of the graph.
	//This is one of the harder parts, so if you're
	//not too comfortable with C++ leave this for last.
	//If you are, there are many ways of doing this,
	//as long as it passes the tests, it's okay.
	class graph_iterator {

	private:

		const weighted_graph<vertex> * owner;		
		size_t curr_position;						//Used to store current position
		//You may need data members here.

	public:
		graph_iterator(const weighted_graph &);
		graph_iterator(const weighted_graph &, size_t);
		~graph_iterator();
		graph_iterator operator=(const graph_iterator&);
		bool operator==(const graph_iterator&) const;
		bool operator!=(const graph_iterator&) const;
		graph_iterator operator++();
		graph_iterator operator++(int);
		const vertex operator*();
		const vertex* operator->();
	};

	//The neighbour_iterator class provides an iterator
	//over the neighbours of a given vertex. This is
	//probably harder (conceptually) than the graph_iterator.
	//Unless you know how iterators work.
	class neighbour_iterator {

	private:
		
		const weighted_graph<vertex> * owner;
		size_t curr_position;
		vertex curr_vertex;						//Used to store a current vertex that iterator is pointing to.
		//You may need data members here.

	public:
		neighbour_iterator(const neighbour_iterator&);
		neighbour_iterator(const weighted_graph &, const vertex&);
		neighbour_iterator(const weighted_graph &, const vertex&, size_t);
		~neighbour_iterator();
		neighbour_iterator operator=(const neighbour_iterator& it);
		bool operator==(const neighbour_iterator&) const;
		bool operator!=(const neighbour_iterator&) const;
		neighbour_iterator operator++();
		neighbour_iterator operator++(int);
		const std::pair<vertex, int> operator*();
		const std::pair<const vertex, int>* operator->();
	};

public:


	weighted_graph(); //A constructor for weighted_graph. It should start empty.
	weighted_graph(const weighted_graph &);//A copy constructor for weighted_graph
	~weighted_graph(); //A destructor. Depending on how you do things, this may
					   //not be necessary.
	int get_vertex_index(const vertex&) const; //Get position of a vertex.

	bool are_adjacent(const vertex&, const vertex&) const; //Returns true if the two vertices are
														   //adjacent, false otherwise.
	bool has_vertex(const vertex&) const; //Returns true if the passed in vertex is 
										  //a vertex of the graph, false otherwise.

	void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const vertex&, const vertex&, const int&); //Adds an edge between the two vertices
															 //with the given weight (as an int).

	void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.
	void set_edge_weight(const vertex&, const vertex&, const int&); //Changes the edge weight between the two
																	//vertices to the new weight (the int).

	int get_edge_weight(const vertex&, const vertex&) const; //Returns the weight on the edge between the two vertices.
	int degree(const vertex&) const; //Returns the degree of the vertex.
	int weighted_degree(const vertex&); //Returns the sum of the weights on all the edges incident to the vertex.
	int num_vertices() const; //Returns the total number of vertices in the graph.
	int num_edges() const; //Returns the total number of edges in the graph (just the count, not the weight).
	int total_weight(); //Returns the sum of all the edge weights in the graph.

	std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
	std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

	graph_iterator begin(); //Returns a graph_iterator pointing to the start of the vertex set.
	graph_iterator end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

	neighbour_iterator neighbours_begin(const vertex&); //Returns a neighbour_iterator pointing to the start
														//of the neighbour set for the given vertex.
	neighbour_iterator neighbours_end(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end
													  //of the neighbour set for the given vertex.

	std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they
													//are visited in by a depth-first traversal starting at
													//the given vertex.
	std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they
													  //are visisted in by a breadth-first traversal starting
													  //at the given vertex.

	bool check_vertex(std::vector<vertex>, const vertex&) const; //Check if a vertex is in a vector(a bit different to has_vertex function).
	vertex get_vertex_by_index(const int &); //Return a vertex by a given index.
	weighted_graph<vertex> mst(); //Returns a minimum spanning tree of the graph.

};

//Define all your methods down here (or move them up into the header, but be careful you don't double up).
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.


template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g) {
	owner = &g;
	curr_position = 0;
}

template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g, size_t start_pos) {
	owner = &g;
	curr_position = start_pos;
}

template <typename vertex> weighted_graph<vertex>::graph_iterator::~graph_iterator() {}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator=(const graph_iterator& it) {
	curr_position = it.curr_position;	//Assign current position equals to graph iterator position.
	return *this;
}

template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator==(const graph_iterator& it) const {
	return ((owner == it.owner) && (curr_position == it.curr_position));	//Decide if every variable of graph iterator is equal.
}

template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator!=(const graph_iterator& it) const {
	return ((owner != it.owner) || (curr_position != it.curr_position));	//Decide if every variable of graph iterator is not equal.
}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++() {
	curr_position++;
	return *this;
}
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(int n) {
	curr_position += n;
	return *this;
}
template <typename vertex> const vertex weighted_graph<vertex>::graph_iterator::operator*() {
	weighted_graph<vertex> g = weighted_graph(*owner);
	return g.get_vertices()[curr_position];
}
template <typename vertex> const vertex* weighted_graph<vertex>::graph_iterator::operator->() {
	auto v = new vertex();
	v = owner->get_vertices()[curr_position];
	return v;
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const neighbour_iterator &){}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u) {
	owner = &g;
	curr_position = 0;
	curr_vertex = u;
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u, size_t start_pos) {
	owner = &g;
	curr_position = start_pos;
	curr_vertex = u;
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::~neighbour_iterator() {}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator=(const neighbour_iterator& it) {
	curr_position = it.curr_position;
	curr_vertex = it.curr_vertex;
	return *this;
}

template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator==(const neighbour_iterator& it) const {
	return ((owner == it.owner) && (curr_position == it.curr_position) && (curr_vertex == it.curr_vertex));
}

template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator!=(const neighbour_iterator& it) const {
	return ((owner != it.owner) || (curr_position != it.curr_position) || (curr_vertex != it.curr_vertex));
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++() {
	curr_position++;
	return *this;
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++(int n) {
	curr_position += n;
	return *this;
}

template <typename vertex> const std::pair<vertex, int> weighted_graph<vertex>::neighbour_iterator::operator*() {
	//A const pointer cannot call its non const member function,
	//so we have to create a non const weight graph.
	auto g = weighted_graph(*owner);
	vertex temp = g.get_neighbours(curr_vertex)[curr_position];
	auto p = new std::pair<const vertex, int>(temp, curr_position);
	return p;
}
template <typename vertex> const std::pair<const vertex, int>* weighted_graph<vertex>::neighbour_iterator::operator->() {
	//A const pointer cannot call its non const member function,
	//so we have to create a non const weight graph.
	auto g = weighted_graph(*owner);
	vertex temp = g.get_neighbours(curr_vertex)[curr_position];
	auto p = new std::pair<const vertex, int>(temp, curr_position);
	return p;
}

template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::begin() {
	return graph_iterator(*this);
}
template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::end() {
	return graph_iterator(*this, vertices.size());
}


template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_begin(const vertex& u) {
	return neighbour_iterator(*this, u);
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_end(const vertex& u) {
	return neighbour_iterator(*this, u, get_neighbours(u).size());
}

template <typename vertex> weighted_graph<vertex>::weighted_graph() {
}

template <typename vertex> weighted_graph<vertex>::weighted_graph(const weighted_graph &g) {
	vertices = g.vertices;
	adjMatrix = g.adjMatrix;
}

template <typename vertex> weighted_graph<vertex>::~weighted_graph() {

}

template <typename vertex> int weighted_graph<vertex>::get_vertex_index(vertex const &u) const {
	int pos = 0;
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i] == u) {
			pos = i;
			break;
		}
	}
	return pos;
}


template <typename vertex> bool weighted_graph<vertex>::has_vertex(const vertex& u) const {
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i] == u) {
			return true;
		}
	}
	return false;
}

template <typename vertex> bool weighted_graph<vertex>::are_adjacent(const vertex& u, const vertex& v) const {
	if (has_vertex(u) && has_vertex(v)) {
		int u_index = get_vertex_index(u);		//Get position of each vertex then return true or false 
		int v_index = get_vertex_index(v);		//by getting the value in the adjacency matrix.
		return adjMatrix[u_index][v_index];
	}
	return false;
}

template <typename vertex> void weighted_graph<vertex>::add_vertex(const vertex& v) {
	if (!has_vertex(v)) {
		for (int i = 0; i < vertices.size(); i++) {
			adjMatrix[i].push_back(0);	//Insert a collumn full of "0" value in the matrix.
		}
		std::vector<int> last_row_vector(num_vertices() + 1, 0);	//then create a vector of the last row with full of "0" value and add it to adjacency matrix.
		adjMatrix.push_back(last_row_vector);
		vertices.push_back(v);
	}
}

template <typename vertex> void weighted_graph<vertex>::add_edge(const vertex& u, const vertex& v, const int& weight) {
	if (has_vertex(u) && has_vertex(v)) {
		int u_index = get_vertex_index(u);
		int v_index = get_vertex_index(v);
		adjMatrix[u_index][v_index] = weight;	//Obviously, when adjMatrix[u][v] = adjMatrix[v][u] and both are not 0 mean those two vertices are adjacent,
		adjMatrix[v_index][u_index] = weight;	//there is an edge between those vertices, and adding the weight for that edge.
	}

}

template <typename vertex> void weighted_graph<vertex>::remove_vertex(const vertex& u) {
	int index = get_vertex_index(u);						//First step is to get the position of the input vertex,
	for (int i = 0; i < vertices.size(); i++) {
		adjMatrix[i].erase(adjMatrix[i].begin() + index); 	//then use that index to find out where it is in the matrix and erase the collumn of it.
	}
	adjMatrix.erase(adjMatrix.begin() + index);				//After that, it will erase the row of that vertex, so the vertex is gone forever.
	vertices.erase(vertices.begin() + index);				//And it also be deleted in the vertices vector
}


template <typename vertex> void weighted_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	if (has_vertex(u) && has_vertex(v)) {
		int u_index = get_vertex_index(u);
		int v_index = get_vertex_index(v);
		adjMatrix[u_index][v_index] = 0;		//Replace the data of the two vertices in adjacency matrix by 0 will make them not to be adjacent and there will be no edge between them.
		adjMatrix[v_index][u_index] = 0;
	}

}

template <typename vertex> void weighted_graph<vertex>::set_edge_weight(const vertex& u, const vertex& v, const int& weight) {
	if (has_vertex(u) && has_vertex(v)) {
		int u_index = get_vertex_index(u);
		int v_index = get_vertex_index(v);
		adjMatrix[u_index][v_index] = weight;
		adjMatrix[v_index][u_index] = weight;
	}
}

template <typename vertex> int weighted_graph<vertex>::get_edge_weight(const vertex& u, const vertex& v) const {
	int u_index = get_vertex_index(u);
	int v_index = get_vertex_index(v);
	return adjMatrix[u_index][v_index];
}

template <typename vertex> int weighted_graph<vertex>::degree(const vertex& u) const {
	int count = 0;
	int u_index = get_vertex_index(u);
	for (int j = 0; j < vertices.size(); j++) {
		if (adjMatrix[u_index][j] > 0) {
			count++;								//This "count" variable is used to count the numbers of edges.
		}
	}
	return count;
}

template <typename vertex> int weighted_graph<vertex>::weighted_degree(const vertex& u) {
	int sum = 0;
	int u_index = get_vertex_index(u);
	for (int j = 0; j < vertices.size(); j++) {
		sum += adjMatrix[u_index][j];						//Sum of all weights of every edge in adjacency matrix.
	}
	return sum;
}

template <typename vertex> int weighted_graph<vertex>::num_vertices() const {
	return vertices.size();
}

template <typename vertex> int weighted_graph<vertex>::num_edges() const {
	int count = 0;
	for (int i = 0; i < vertices.size(); i++) {
		for (int j = i + 1; j < vertices.size(); j++) {	//j = i+1 is used to scan the half of adjacency matrix, because we do not need to scan all over it.
			if (adjMatrix[i][j] != 0) {
				count++;
			}
		}
	}
	return count;
}

template <typename vertex> int weighted_graph<vertex>::total_weight() {
	int sum = 0;
	for (int i = 0; i < vertices.size(); i++) {
		for (int j = i + 1; j < vertices.size(); j++) {	//j = i+1 is used to scan the half of adjacency matrix, because we do not need to scan all over it.
			sum += adjMatrix[i][j];
		}
	}
	return sum;
}

template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_vertices() {
	return vertices;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::get_neighbours(const vertex& u) {
	std::vector<vertex> neighbours;
	for (int i = 0; i < num_vertices(); i++) {
		if (are_adjacent(u, vertices[i])) {
			neighbours.push_back(vertices[i]);			//After checking if vertex u and vertex in the vertices vector are adjacent, vertex in the vertices vector will be added in neighbours vector.
		}
	}
	return neighbours;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::depth_first(const vertex& start_vertex) {
	bool visited[vertices.size()];
	for (unsigned i = 0; i < vertices.size(); i++) {
		visited[i] = false;							//Put every vertex in the vertices to be not visited.
	}
	//Create stuffs to store data
	std::stack<int> unprocessed;
	unprocessed.push(start_vertex);

	std::vector<int> ordered;

	while (!unprocessed.empty()) {
		int n = unprocessed.top();
		unprocessed.pop();
		//Check if the vertex has been visited or not
		if (!visited[n]) {
			visited[n] = true;
			ordered.push_back(n);
			for (unsigned i = vertices.size(); i != 0; i--) {
				if (adjMatrix[n][i - 1]) {
					unprocessed.push(i - 1);
				}
			}
		}
	}

	return ordered;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::breadth_first(const vertex& start_vertex) {
	bool visited[vertices.size()];

	for (unsigned i = 0; i < vertices.size(); i++) {
		visited[i] = false;							//Put every vertex in the vertices to be not visited.
	}
	//Create stuffs to store data
	std::queue<int> unprocessed;
	unprocessed.push(start_vertex);

	std::vector<int> ordered;

	while (!unprocessed.empty()) {
		int n = unprocessed.front();
		unprocessed.pop();
		//Check if the vertex has been visited or not
		if (!visited[n]) {
			visited[n] = true;
			ordered.push_back(n);
			for (unsigned i = 0; i < vertices.size(); i++) {
				if (adjMatrix[n][i]) {
					unprocessed.push(i);
				}
			}
		}
	}
	return ordered;
}

template <typename vertex> bool weighted_graph<vertex>::check_vertex(std::vector<vertex> a, const vertex &v) const {
	for (int i = 0; i < a.size(); i++) {
		if (a[i] == v) {
			return true;
		}
	}
	return false;
}

template <typename vertex> vertex weighted_graph<vertex>::get_vertex_by_index(const int &n) {
	return vertices[n];	//Return a vertex
}

template <typename vertex>	weighted_graph<vertex> weighted_graph<vertex>::mst() {
	//Initialize a minimum spanning tree graph.
	//A B C D E -> 0 1 2 3 4
	int size = vertices.size();

	weighted_graph<vertex> mst;
	for (int i = 0; i < size; i++) {
		mst.add_vertex(vertices[i]);
	}
	//Creating some vectors to store visited, unvisited and parents of vertices.
	std::vector<int> visited_vertices;
	std::vector<int> parent_vertices;
	std::vector<int> unvisited_vertices;
	for (int i = 0; i < size; i++) {
		unvisited_vertices.push_back(0);
	}
	int count_edge = 0, min;
	//I would like to choose the first vertex in vector vertices rather than choose randomly.
	visited_vertices.push_back(0);
	parent_vertices.push_back(0);
	unvisited_vertices[0] = 1;
	//Create a vector to store weights for mst.
	std::vector<int> weight_mst;
	//Run until getting enough edge for the vertices.
	while (count_edge < size - 1) {
		std::vector<int>::iterator iter = visited_vertices.begin();
		int first_min = 1;	//This is created to make min equals to weight as min is now as a random number and not assigned to anything.
		std::pair<int, int> chosen_candidate;
		while (iter != visited_vertices.end()) {
			for (int i = 0; i < size; i++) {
				int weight = adjMatrix[*iter][i];
				if (!unvisited_vertices[i] && weight) {
					if (first_min) {
						min = weight;
						chosen_candidate.first = *iter;
						chosen_candidate.second = i;
						first_min = 0;
					}
					else {
						if (weight < min) {
							min = weight;
							chosen_candidate.first = *iter;
							chosen_candidate.second = i;
						}
					}
				}
			}
			// switch to next iter
			iter++;
		}
		visited_vertices.push_back(chosen_candidate.second);
		parent_vertices.push_back(chosen_candidate.first);
		unvisited_vertices[chosen_candidate.second] = 1;
		weight_mst.push_back(min);
		count_edge++;
	}
	for (int i = 1; i < size; i++) {
		mst.add_edge(get_vertex_by_index(visited_vertices[i]), get_vertex_by_index(parent_vertices[i]), weight_mst[i - 1]);	//Adding the edges into mst.
	}

	return mst;
}

#endif
