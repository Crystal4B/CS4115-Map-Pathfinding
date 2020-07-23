#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <math.h>
#include <limits>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::string;
using std::vector;
using std::queue;
using std::sort;
using std::sqrt;
using std::pow;
using std::numeric_limits;

// Definitions
#define INF numeric_limits<float>::max()

// Structs
struct Node {int vertex, distance;}; // Used by the BFS
struct NodePath{vector<int> path; int vertex, distance;}; // Used by the BFS to find the path
struct NodeDij{vector<int> path; float distance;}; // Used by the dijkstra
struct Coord {float x, y, z;}; // Used to describe a vertex's xyz coordinates

// Type definitions
typedef vector<vector<int>> Graph;
typedef vector<Coord> XYZ;

// Functions
void readGraph(const string graphname, Graph& graph, int& order, int& size);
void readXYZ(const string xyzname, const int order, XYZ& xyz);
vector<int> getNeighbours(const Graph graph, const int vertex, const int distance);
NodeDij shortestDijkstra(const XYZ xyz, const Graph graph, const int start, const int destination);
NodePath shortestBFS(const Graph graph, const int start, const int destination);

int main()
{
	// get name of our file
	string basename, graphname, xyzname;
	cin >> basename;

	graphname = "./Data/" + basename + ".osm.graph";
	xyzname += "./Data/" + basename + ".osm.xyz";

	// initalize our data
	Graph graph;
	int order, size;
	readGraph(graphname, graph, order, size);

	int query;
	cin >> query;

	// Check which query is being requested
	switch(query)
	{
		// QUERY 1: What is the Order and Size ? (10 marks)
		case 1:
		{
			cout << "n= " << order << "; m= " << size << "." << endl;
			break;
		}
		// QUERY 2: Whats the largest amount of edges that a vertex has ? (5 marks)
		case 2:
		{
			int largest = 0;

			for (unsigned int i = 0; i < graph.size(); i++)
			{
				largest = graph[i].size() > graph[largest].size() ? i : largest;
			}
			cout << "v= " << largest+1 << "; |N(v)|= " << graph[largest].size() << "." << endl;
			break;
		}
		// QUERY 3: Whats the average amount of edges that the vertices have ? (5 marks)
		case 3:
		{
			float average = 0;

			for (vector<int> neighbours : graph)
			{
				average += neighbours.size();
			}
			average /= graph.size();
			printf("avg |N(v)|= %.6f.\n", average);
			break;
		}
		// QUERY 4: List the Neibhours 1 edge away (5 marks)
		case 4:
		{
			int vertex;
			cin >> vertex;

			cout << "N(" << vertex << ")=";
			for (int neighbour : graph[vertex-1])
			{
				cout << " " << neighbour;
			}
			cout << "." << endl;
			break;
		}
		// QUERY 5: List the Neibhours k edges away (15 marks)
		case 5:
		{
			int vertex, distance;
			cin >> vertex >> distance;

			vector<int> neighbours = getNeighbours(graph, vertex, distance);
			cout << "N(" << vertex << "," << distance << ")=";
			for (int neighbour : neighbours)
			{
				cout << " " << neighbour;
			}
			cout << "." << endl;
			break;
		}
		// QUERY 6: Find the shortest path using xyz (15 marks)
		case 6:
		{
			XYZ xyz;
			readXYZ(xyzname, order, xyz);

			int start, destination;
			cin >> start >> destination;

			NodeDij finalChild = shortestDijkstra(xyz, graph, start, destination);
			if (finalChild.distance == -1)
			{
				cout << "No path from " << start << " to " << destination << endl;
			}
			else if (!finalChild.path.empty())
			{
				cout << "Path: " << finalChild.path[0];
				for (unsigned int i = 1; i < finalChild.path.size(); i++)
				{
					cout << " - " << finalChild.path[i];
				}
				cout << " - " << destination << " // optional output" << endl;
			}
			printf("d(%d,%d)= %.6f.\n", start, destination, finalChild.distance);
			break;
		}
		// QUERY 7: Find the shortest path using edges (15 marks)
		case 7:
		{
			int start, destination;
			cin >> start >> destination;

			NodePath finalChild = shortestBFS(graph, start, destination);
			if (finalChild.vertex == -1)
			{
				cout << "No path from " << start << " to " << destination << endl;
			}
			else
			{
				cout << "Path: " << finalChild.path[0];
				for (unsigned int i = 1; i < finalChild.path.size(); i++)
				{
					cout << " - " << finalChild.path[i];
				}
				cout << " - " << destination << " // optional output" << endl;
			}
			cout << "ed(" << start << "," << destination << ")= " << finalChild.distance << "." << endl;
			break;
		}
	}
}

// Reads the graph file (country.osm.graph) and stores all the data in Graph vector, order and size
void readGraph(const string graphname, Graph& graph, int& order, int& size)
{
	// prepare for reading
	ifstream file(graphname);
	string line;

	while (getline(file, line))
	{
		if (line[0] != '%'){break;}
	}

	istringstream lstream(line);
	lstream >> order >> size; // order and size

	Graph temp(order);
	graph = temp; // initialize our graph

	int index = 0;
	// with order and size we can finally read the rest of the graph
	while (getline(file, line)) // get next vertex and neighbours
	{
		istringstream lstream(line);
		int neighbour;
		while (lstream >> neighbour)
		{
			// NOTE: vertex to index conversion takes place here where index = vertex-1
			graph[index].push_back(neighbour);
		}
		index++;
	}
}

// Reads the xyz file (country.osm.xyz) and stores all the data in the XYZ vector
void readXYZ(const string xyzname, const int order, XYZ& xyz)
{
	// prepare for reading
	ifstream file(xyzname);
	string line;

	XYZ temp(order); // initialize our xyz
	xyz = temp;

	int index = 0;
	while (getline(file, line)) // get next vertex and its coordinates
	{
		istringstream lstream(line);
		vector<float> coords;
		float coord;
		while (lstream >> coord)
		{
			coords.push_back(coord);
		}
		// NOTE: vertex to index conversion takes place here where index = vertex-1
		xyz[index] = {coords[0], coords[1], coords[2]};
		index++;
	}
}

// Returns the neighbours of vertex at a specific distance using BFS
vector<int> getNeighbours(const Graph graph, const int vertex, const int distance)
{
	// validate input
	if (distance < 0)
	{
		return {};
	}
	else if (distance == 0)
	{
		return {vertex};
	}

	// list for neighbours who fit our criteria
	vector<int> results;

	// prepare open and closed
	queue<Node> open;
	vector<bool> discovered(graph.size(), false);

	open.push({vertex,0});

	// main loop
	while (!open.empty())
	{
		// get node your looking at
		Node currentNode = open.front();
		open.pop();

		// find its neighbours
		for (int vertex : graph[currentNode.vertex-1])
		{
			Node childNode{vertex, currentNode.distance+1};
			// check if goal condition is met
			if (childNode.distance == distance)
			{
				if (!discovered[vertex-1])
				{
					results.push_back(vertex);
					discovered[vertex-1] = true;
				}
			}
			else
			{
				open.push(childNode);
			}
		}
	}
	// due to specification requirements we must sort the vector here
	sort(results.begin(), results.end());
	return results;
}

// Finds the shortest distance from start to destination using Dijkstra
NodeDij shortestDijkstra(const XYZ xyz, const Graph graph, const int start, const int destination)
{
	// validate input
	if (start == destination)
	{
		return {{start}, 0};
	}

	// prepare necessary datastructures
	vector<NodeDij> distances(graph.size(), {{},INF});
	vector<bool> discovered(graph.size(), false);

	distances[start-1].distance = 0;

	// main loop
	for (unsigned int i = 0; i < graph.size(); i++)
	{
		int index;

		// get nearest index
		float min = INF;
		for (unsigned int i = 0; i < graph.size(); i++)
		{
			if (distances[i].distance <= min && !discovered[i])
			{
				index = i;
				min = distances[i].distance;
			}
		}

		// check if the chosen path is out of reach
		if (distances[index].distance == INF)
		{
			return {{}, -1};
		}

		// check if goal condition has been met
		if (index+1 == destination)
		{
			return distances[index];
		}

		// mark index as complete (nearest to start index)
		discovered[index] = true;

		// update vertexes neighbours dists
		for (int neighbour : graph[index])
		{
			if (!discovered[neighbour-1])
			{
				float distance = sqrt(pow(xyz[index].x - xyz[neighbour-1].x, 2) + pow(xyz[index].y - xyz[neighbour-1].y, 2) + pow(xyz[index].z - xyz[neighbour-1].z, 2));
				if (distances[neighbour-1].distance > distances[index].distance + distance && distances[index].distance != INF)
				{
					vector<int> path = distances[index].path;
					path.push_back(index+1);
					distances[neighbour-1].distance = distance + distances[index].distance;
					distances[neighbour-1].path = path;
				}

			}
		}
	}
	return {{}, -1};
}

// Finds the shortest path from start to destination using BFS
NodePath shortestBFS(const Graph graph, const int start, const int destination)
{
	// validate input
	if (start == destination)
	{
		return {{start}, 0, 0};
	}

	// prepare open and closed lists
	queue<NodePath> open;
	vector<bool> discovered(graph.size(), false);

	open.push({{},start,0});
	discovered[start-1] = true;

	// main loop
	while (!open.empty())
	{
		// get node your looking at
		NodePath currentNode = open.front();
		open.pop();

		// find its neighbours
		for (int vertex : graph[currentNode.vertex-1])
		{
			// prepare path for child
			vector<int> path = currentNode.path;
			path.push_back(currentNode.vertex);

			// prepare child
			NodePath childNode{path, vertex, currentNode.distance+1};
			// check if goal condition is met
			if (childNode.vertex == destination)
			{
				return childNode;
			}
			else if(!discovered[vertex-1])
			{
				open.push(childNode);
				discovered[vertex-1] = true;
			}
		}
	}
	return {{}, -1, -1};
}