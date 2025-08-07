#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <regex>
#include <numeric>
#include <chrono>
#include <random>
#include <climits>
#include <lib/multipleProcessors/multiple_processors.h>

void DAG::count_Nodes_Edges(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open the DAG file for counting");
  }
  std::string line;
  numNodes = 0; 
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string token;
    // Node line
    if (line.find("[weight=") != std::string::npos) {
      numNodes++;
    }

    // Edge line
    if (line.find("->") != std::string::npos) {
      numEdges++;
    }
  }

  nodeWeights.resize(numNodes);
  adjList_inNeighbors.resize(numNodes);
  adjList_outNeighbors.resize(numNodes);

  file.close();
}

void DAG::readGraph(const std::string& filename) {
  std::ifstream file(filename); 
  if (!file.is_open()) {
    throw std::invalid_argument("Could not open the DAG file for reading");
  }

  std::string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string token; 

    if (line.find("[weight=") != std::string::npos) {
      // Node line
      std::string node; 
      ss >> node >> token; 
      int nodeID = std::stoi(node);
      int nodeWeight = -1; 
      std::regex num_regex(R"(\d+)"); 
      std::smatch match; 
      if (std::regex_search(token, match, num_regex)) {
        nodeWeight = std::stoi(match.str()); 
      } else {
        throw::std::runtime_error("No weight specified for the vertex!"); 
      }
      nodeWeights[nodeID] = nodeWeight; 
    } else if (line.find("->") != std::string::npos) {
      std::string from, to; 
      ss >> from >> token >> to >> token;
      int sourceID = std::stoi(from);
      int targetID = std::stoi(to);
      int edgeWeight = -1; 
      std::regex num_regex(R"(\d+)"); 
      std::smatch match; 
      if (std::regex_search(token, match, num_regex)) {
        edgeWeight = std::stoi(match.str()); 
      } else {
        throw::std::runtime_error("No weight specified for the vertex!"); 
      }
        vertex v_out = {targetID, edgeWeight, "comp", std::make_pair(-1, -1)};
        vertex v_in = {sourceID, edgeWeight, "comp", std::make_pair(-1, -1)};
        adjList_outNeighbors[sourceID].push_back(v_out);
        adjList_inNeighbors[targetID].push_back(v_in);
    }
  }
  file.close();
}

void DAG::printGraph() const {
  std::cout << "==================================" << std::endl;
  std::cout << "            PRINT DAG             " << std::endl; 
  std::cout << "==================================" << std::endl;
  std::cout << "Number of nodes: " << numNodes << ". Number of edges: " << numEdges << std::endl;
  std::cout << "Node weights: " << std::endl;
  for (size_t i = 0; i < nodeWeights.size(); i++) {
    std::cout << "Node " << i << ": " << nodeWeights[i] << std::endl;
  }

  std::cout << "Adjacency list (out): " << std::endl;
  for (size_t i = 0; i < numNodes; i++) {
    std::cout << i << ": ";
    for (size_t j = 0; j < adjList_outNeighbors[i].size(); j++) {
      std::cout << adjList_outNeighbors[i][j].id << " "; 
    }
    std::cout << std::endl; 
  }

  std::cout << "Adjacency list (in): " << std::endl;
  for (size_t i = 0; i < numNodes; i++) {
    for (size_t j = 0; j < adjList_inNeighbors[i].size(); j++) {
      std::cout << adjList_inNeighbors[i][j].id << " "; 
    }
    std::cout << std::endl; 
  }
}

////////////////////////////////////////////////////////////////
// SCHEDULER //
///////////////////////////////////////////////////////////////

MultiMachineSchedulingWithoutFreezing::MultiMachineSchedulingWithoutFreezing(const std::string& DAG_file, const std::string& mapping_file, const std::string& cluster_setup_file, int communication_mode) {    
  // read graph
  G.count_Nodes_Edges(DAG_file);
  G.readGraph(DAG_file);

  // read mapping
  readMapping(mapping_file);

  // read setup
  readSetup(cluster_setup_file);

  // Track preprocessing time
  auto start_setup = std::chrono::high_resolution_clock::now();

  intervals.resize(numSubintervals);
  for (int i = 0; i < numSubintervals; i++) {
    interval I; 
    I.interval_id = i; 
    I.start = splitPoints[i]; 
    I.end = splitPoints[i+1]; 
    I.interval_length = subintervalLengths[i]; 
    I.interval_budget = subintervalBudgets[i]; 
    intervals[i] = I;
  }

  addMappingPaths();
  if (communication_mode == 0) {
    startingTimes.resize(G.getNumNodes(), -1);
    already_scheduled.resize(G.getNumNodes(), false);

    compute_EST_LST();

    // for cost computations later on
    init_basePower(); 
  } else {
    // construct communication enhanced DAG, also initializes the base power
    setup_communication_setting();
  }

  auto end_setup = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_setup - start_setup);
  setupTimeCommunicationSetting = duration.count();
}


MultiMachineSchedulingWithoutFreezing::MultiMachineSchedulingWithoutFreezing(const std::string& DAG_file, const std::string& mapping_file, int communication_mode) {  
  // read graph
  G.count_Nodes_Edges(DAG_file);
  G.readGraph(DAG_file);

  // read mapping
  readMapping(mapping_file);

  // setup the graph to get the makespan only
  addMappingPaths();
  
  if (communication_mode == 0) {
    startingTimes.resize(G.getNumNodes(), -1);
    already_scheduled.resize(G.getNumNodes(), false);
    compute_EST_only(); 
  } else {
    // EST is computed here anyway
    setup_communication_setting_for_EST_only(); 
  }
}


MultiMachineSchedulingWithoutFreezing::MultiMachineSchedulingWithoutFreezing(const std::string& DAG_file) {
  G.count_Nodes_Edges(DAG_file);
  G.readGraph(DAG_file);
  G.printGraph();
}

void MultiMachineSchedulingWithoutFreezing::setup_communication_setting_for_EST_only() {
  int number_of_vertices = G.getNumNodes();

  // count additional communication vertices
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    for (size_t j = 0; j < G.adjList_outNeighbors[i].size(); j++) {
      if (processors_of_tasks[i] != processors_of_tasks[G.adjList_outNeighbors[i][j].id]) { 
        number_of_vertices++;
      }
    }
  }

  G_comm.numNodes = number_of_vertices;
  number_of_computation_tasks = G.getNumNodes();

  G_comm.nodeWeights.resize(number_of_vertices, 0);
  G_comm.adjList_inNeighbors.resize(number_of_vertices);
  G_comm.adjList_outNeighbors.resize(number_of_vertices);
    
  // set back for construction
  number_of_vertices = number_of_computation_tasks;
    
  // add all communication tasks and the corresponding edges + the paths on the computation processors
  for (int i = 0; i < G.getNumNodes(); i++) {
    G_comm.nodeWeights[i] = G.nodeWeights[i];
    vertex ve = {i, 0, "comp", std::make_pair(-1, -1)};
    vertices.push_back(ve); 
    for (int j = 0; j < G.adjList_outNeighbors[i].size(); j++) {
      if (processors_of_tasks[i] != processors_of_tasks[G.adjList_outNeighbors[i][j].id]) { 
        G_comm.nodeWeights[number_of_vertices] = G.adjList_outNeighbors[i][j].weight;
        vertex comm_vertex = {number_of_vertices, 0, "comm", std::make_pair(i,G.adjList_outNeighbors[i][j].id)};
        vertices.push_back(comm_vertex);
        vertex v = {i, 0, "comp", std::make_pair(-1, -1)};

        G_comm.adjList_outNeighbors[i].push_back(comm_vertex);
        G_comm.adjList_inNeighbors[comm_vertex.id].push_back(v);
        G_comm.adjList_outNeighbors[comm_vertex.id].push_back(G.adjList_outNeighbors[i][j]);
        G_comm.adjList_inNeighbors[G.adjList_outNeighbors[i][j].id].push_back(comm_vertex);
                    
        number_of_vertices++;
      } else if (processors_of_tasks[i] == processors_of_tasks[G.adjList_outNeighbors[i][j].id]) { 
        vertex v_in = {i, 0, "comp", std::make_pair(-1, -1)};
        G_comm.adjList_outNeighbors[i].push_back(G.adjList_outNeighbors[i][j]);
        G_comm.adjList_inNeighbors[G.adjList_outNeighbors[i][j].id].push_back(v_in); 
      }
    }
  }
    
  compute_EST_only();
  G = std::move(G_comm); 

  /*
  * Assign tasks to processors for computation and communication
  * 0 = PE0, 1 = PE0->PE1 ,... 
  */
  std::vector<int> processors_of_tasks_comm(G.getNumNodes()); 
  processors_of_tasks_comm.resize(G.getNumNodes());
  for (vertex v : vertices) {
    if (v.type == "comp") {
      processors_of_tasks_comm[v.id] = processors_of_tasks[v.id];
    } else {
      if (processors_of_tasks[v.comm_id.first] == 0) {
        processors_of_tasks_comm[v.id] = numProcessors + processors_of_tasks[v.comm_id.first]*(numProcessors-1) + processors_of_tasks[v.comm_id.second]-1;
      } else {
        if (processors_of_tasks[v.comm_id.second] < processors_of_tasks[v.comm_id.first]) {
          processors_of_tasks_comm[v.id] = numProcessors + processors_of_tasks[v.comm_id.first]*(numProcessors-1) + processors_of_tasks[v.comm_id.second];
        } else if (processors_of_tasks[v.comm_id.first] < processors_of_tasks[v.comm_id.second]) {
          processors_of_tasks_comm[v.id] = numProcessors + processors_of_tasks[v.comm_id.first]*(numProcessors-1) + processors_of_tasks[v.comm_id.second]-1;
        }
      }
    }
  }
    
  processors_of_tasks.resize(G.getNumNodes());
  processors_of_tasks = std::move(processors_of_tasks_comm);
    
  // construct mapping data structure for communication tasks
  unsigned int numProcessors_comm = numProcessors + numProcessors*(numProcessors-1);
  std::vector<std::vector<vertex>> mapping_comm(numProcessors_comm);
  std::vector<std::vector<int>> mapping_comm_id(numProcessors_comm);
  for (vertex v : vertices) {
    if (v.type == "comm") {
      mapping_comm[processors_of_tasks[v.id]].push_back(v);
      mapping_comm_id[processors_of_tasks[v.id]].push_back(v.id);
    }
  }

  // sort tasks on communication processors
  for (auto& vec : mapping_comm) {
    if (vec.size() > 1) {
      std::sort(vec.begin(), vec.end(), [this](const vertex& a, const vertex& b) {
        if (EST_LST[a.comm_id.first].first < EST_LST[b.comm_id.first].first) {
          return true;
        } else if (EST_LST[a.comm_id.first].first == EST_LST[b.comm_id.first].first) {
          return (EST_LST[a.comm_id.second].first < EST_LST[b.comm_id.second].first);
        }
        // without this: error because undefined behaviour is possible
        return false;
      });
    }
  }
    
  EST_LST.resize(0);

  // add paths on communication processors
  for (size_t i = 0; i < mapping_comm.size(); i++) {
    if (mapping_comm[i].size() > 1) { 
      for (size_t j = 0; j <= mapping_comm[i].size()-2; j++) {
        vertex source = mapping_comm[i][j];
        vertex target = mapping_comm[i][j+1];
        G.adjList_outNeighbors[source.id].push_back(target);
        G.adjList_inNeighbors[target.id].push_back(source);
      }
    }
  }
    
  numProcessors = numProcessors_comm;
  mapping.resize(0); 
  mapping.resize(numProcessors);
    
  for (vertex v : vertices) {
    mapping[processors_of_tasks[v.id]].push_back(v.id);
  }
    
  compute_EST_only();
  startingTimes.resize(G.getNumNodes(), -1);    
  already_scheduled.resize(G.getNumNodes(), false);  
}

void MultiMachineSchedulingWithoutFreezing::setup_communication_setting() {
  int number_of_vertices = G.getNumNodes();

  // count total number of vertices
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    for (size_t j = 0; j < G.adjList_outNeighbors[i].size(); j++) {
      if (processors_of_tasks[i] != processors_of_tasks[G.adjList_outNeighbors[i][j].id]) { 
        number_of_vertices++;
      }
    }
  }

  G_comm.numNodes = number_of_vertices;
  number_of_computation_tasks = G.getNumNodes();

  G_comm.nodeWeights.resize(number_of_vertices, 0);
  G_comm.adjList_inNeighbors.resize(number_of_vertices);
  G_comm.adjList_outNeighbors.resize(number_of_vertices);
    
  // set back for construction
  number_of_vertices = G.getNumNodes();

  // add all communication tasks and the corresponding edges + the paths on the computation processors
  for (int i = 0; i < G.getNumNodes(); i++) {
    G_comm.nodeWeights[i] = G.nodeWeights[i];
    vertex ve = {i, 0, "comp", std::make_pair(-1, -1)};
    vertices.push_back(ve); 
    for (int j = 0; j < G.adjList_outNeighbors[i].size(); j++) {
      if (processors_of_tasks[i] != processors_of_tasks[G.adjList_outNeighbors[i][j].id]) { 
        G_comm.nodeWeights[number_of_vertices] = G.adjList_outNeighbors[i][j].weight;
        vertex comm_vertex = {number_of_vertices, 0, "comm", std::make_pair(i,G.adjList_outNeighbors[i][j].id)};
        vertices.push_back(comm_vertex);
        vertex v = {i, 0, "comp", std::make_pair(-1, -1)};
        
        G_comm.adjList_outNeighbors[i].push_back(comm_vertex);
        G_comm.adjList_inNeighbors[comm_vertex.id].push_back(v);
        G_comm.adjList_outNeighbors[comm_vertex.id].push_back(G.adjList_outNeighbors[i][j]);
        G_comm.adjList_inNeighbors[G.adjList_outNeighbors[i][j].id].push_back(comm_vertex);
                    
        number_of_vertices++;
      } else if (processors_of_tasks[i] == processors_of_tasks[G.adjList_outNeighbors[i][j].id]) { 
        vertex v_in = {i, 0, "comp", std::make_pair(-1, -1)};
        G_comm.adjList_outNeighbors[i].push_back(G.adjList_outNeighbors[i][j]);
        G_comm.adjList_inNeighbors[G.adjList_outNeighbors[i][j].id].push_back(v_in); 
      }
    }
  }
    
  // for communication order we need EST_LST 
  compute_EST_LST();
  G = std::move(G_comm); 

  /*
  * Assign tasks to processors for computation and communication
  * 0 = PE0, 1 = PE0->PE1 ,... 
  */
  processors_of_tasks_comm.resize(G.getNumNodes());
  for (vertex v : vertices) {
    if (v.type == "comp") {
      processors_of_tasks_comm[v.id] = processors_of_tasks[v.id];
    } else {
      if (processors_of_tasks[v.comm_id.first] == 0) {
        processors_of_tasks_comm[v.id] = numProcessors + processors_of_tasks[v.comm_id.first]*(numProcessors-1) + processors_of_tasks[v.comm_id.second]-1;
      } else {
        if (processors_of_tasks[v.comm_id.second] < processors_of_tasks[v.comm_id.first]) {
          processors_of_tasks_comm[v.id] = numProcessors + processors_of_tasks[v.comm_id.first]*(numProcessors-1) + processors_of_tasks[v.comm_id.second];
        } else if (processors_of_tasks[v.comm_id.first] < processors_of_tasks[v.comm_id.second]) {
          processors_of_tasks_comm[v.id] = numProcessors + processors_of_tasks[v.comm_id.first]*(numProcessors-1) + processors_of_tasks[v.comm_id.second]-1;
        }
      }
    }
  }

  processors_of_tasks.resize(G.getNumNodes());
  processors_of_tasks = std::move(processors_of_tasks_comm);
    
  // construct mapping data structure for communication tasks
  unsigned int numProcessors_comm = numProcessors + numProcessors*(numProcessors-1);
  std::vector<std::vector<vertex>> mapping_comm(numProcessors_comm);
  std::vector<std::vector<int>> mapping_comm_id(numProcessors_comm);
  for (vertex v : vertices) {
    if (v.type == "comm") {
      mapping_comm[processors_of_tasks[v.id]].push_back(v);
      mapping_comm_id[processors_of_tasks[v.id]].push_back(v.id);
    }
  }

  // sort tasks on communication processors
  for (auto& vec : mapping_comm) {
    if (vec.size() > 1) {
      std::sort(vec.begin(), vec.end(), [this](const vertex& a, const vertex& b) {
        if (EST_LST[a.comm_id.first].first < EST_LST[b.comm_id.first].first) {
          return true;
        } else if (EST_LST[a.comm_id.first].first == EST_LST[b.comm_id.first].first) {
          return (EST_LST[a.comm_id.second].first < EST_LST[b.comm_id.second].first);
        }
        // without this: error because undefined behaviour is possible
        return false;
      });
    }
  }

  EST_LST.resize(0);

  // add paths on communication processors.
  for (size_t i = 0; i < mapping_comm.size(); i++) {
    if (mapping_comm[i].size() > 1) { 
      for (size_t j = 0; j <= mapping_comm[i].size()-2; j++) {
        vertex source = mapping_comm[i][j];
        vertex target = mapping_comm[i][j+1];
        G.adjList_outNeighbors[source.id].push_back(target);
        G.adjList_inNeighbors[target.id].push_back(source);        
      }
    }
  }
    
  numProcessors = numProcessors_comm;
  mapping.resize(0); 
  mapping.resize(numProcessors);

  for (vertex v : vertices) {
    mapping[processors_of_tasks[v.id]].push_back(v.id);
  }
    
  compute_EST_LST();
  startingTimes.resize(G.getNumNodes(), -1);
  already_scheduled.resize(G.getNumNodes(), false);    
  init_basePower(); 
}

void MultiMachineSchedulingWithoutFreezing::readMapping(const std::string& filename) {
  std::ifstream file(filename); 
  if (!file.is_open()) {
    throw std::runtime_error("Could not open the mapping file");
  }

  std::string line; 
  std::getline(file, line);
  numProcessors = std::stoi(line);
  mapping.resize(numProcessors);

  int processor_counter = 0; 
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    int task;
    while(ss >> task) {
      mapping[processor_counter].push_back(task);
    }
    processor_counter++;
  }

  processors_of_tasks.resize(G.getNumNodes());
  for (size_t i = 0; i < numProcessors; i++) {
    for (size_t j = 0; j < mapping[i].size(); j++) {
      processors_of_tasks[mapping[i][j]] = i;
    }
  }
  file.close();
}

void MultiMachineSchedulingWithoutFreezing::readSetup(const std::string& filename) {
  std::ifstream file(filename); 
  if (!file.is_open()) {
    throw std::runtime_error("Could not open the setup file");
  }

  std::string line; 
  // read deadline
  std::getline(file, line);
  deadline = std::stoi(line);

  // read number of subintervals
  std::getline(file, line);
  numSubintervals = std::stoi(line);

  // read PStatic
  std::getline(file, line);
  std::stringstream ss_PStatic(line);
  int static_value;
  while (ss_PStatic >> static_value) {
    PStatic.push_back(static_value);
  }

  // read PDynamic
  std::getline(file, line);
  std::stringstream ss_PDynamic(line);
  int dynamic_value;
  while (ss_PDynamic >> dynamic_value) {
    PDynamic.push_back(dynamic_value);
  }

  // read Bcc
  std::getline(file, line);
  BCC = std::stoi(line);

  // read subinterval budgets
  std::getline(file, line);
  std::stringstream ss_budgets(line);
  int budget;
  while (ss_budgets >> budget) {
    subintervalBudgets.push_back(budget);
  }

  // read split points of time horizon
  std::getline(file, line);
  std::stringstream ss_split(line);
  int point;
  while (ss_split >> point) {
    splitPoints.push_back(point);
  }

  file.close();

  // compute subinterval lengths   
  for (size_t i = 0; i < splitPoints.size()-1; i++) {
    subintervalLengths.push_back(splitPoints[i+1] - splitPoints[i]);
  }
}

void MultiMachineSchedulingWithoutFreezing::addMappingPaths() {
  for (size_t i = 0; i < mapping.size(); i++) {
    if (mapping[i].size() > 1) { 
      for (size_t j = 0; j <= mapping[i].size()-2; j++) {
        int source = mapping[i][j];
        int target = mapping[i][j+1];
        vertex v_out = {target, 0, "comp", std::make_pair(-1, -1)}; 
        vertex v_in = {source, 0, "comp", std::make_pair(-1, -1)};

        bool edge_exists = false;
        for (size_t k = 0; k < G.adjList_outNeighbors[source].size(); k++) {
          if (G.adjList_outNeighbors[source][k].id == target) {
            edge_exists = true;
            break;
          }
        }
            
        if (!edge_exists) {
          G.adjList_outNeighbors[source].push_back(v_out);
          G.adjList_inNeighbors[target].push_back(v_in);
        }
      }
    }
  }
}

void MultiMachineSchedulingWithoutFreezing::compute_EST_LST() {
  EST_LST.resize(G.getNumNodes(), std::make_pair(0, deadline+1));

  // compute earliest start time by using an analogy to Kahn's algorithm
  std::queue<int> Q; 
  std::vector<int> in_degree(G.getNumNodes(), 0);
    
  // compute in-degree
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    in_degree[i] = G.adjList_inNeighbors[i].size();
  }
    
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    if (G.adjList_inNeighbors[i].size() == 0) {
      Q.push(i);
    }
  }

  while (!Q.empty()) {
    int node = Q.front(); 
    Q.pop();

    int dependency_constraint = 0; 
    for (size_t k = 0; k < G.adjList_inNeighbors[node].size(); k++) {
      int neighbor = G.adjList_inNeighbors[node][k].id;
      dependency_constraint = std::max(dependency_constraint, EST_LST[neighbor].first + G.nodeWeights[neighbor]);
    }
        
    EST_LST[node].first = dependency_constraint;

    // process children
    for (size_t k = 0; k < G.adjList_outNeighbors[node].size(); k++) {
      int neighbor = G.adjList_outNeighbors[node][k].id;
      in_degree[neighbor]--;
      if (in_degree[neighbor] == 0) {
        Q.push(neighbor);
      }
    }
  }

  // now do Kahn's Algorithm kind of reversed for latest start time
  std::queue<int> Q_LST;

  for (size_t i = 0; i < G.getNumNodes(); i++) {
    if (G.adjList_outNeighbors[i].size() == 0) {
      EST_LST[i].second = deadline - G.nodeWeights[i];
      Q_LST.push(i);
    }
  }

  std::vector<int> out_degree(G.getNumNodes(), 0);
    
  // compute out-degree
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    out_degree[i] = G.adjList_outNeighbors[i].size();
  }
  while (!Q_LST.empty()) {
    int node = Q_LST.front(); 
    Q_LST.pop();

    for (size_t k = 0; k < G.adjList_inNeighbors[node].size(); k++) {
      int neighbor = G.adjList_inNeighbors[node][k].id;
      if (EST_LST[neighbor].second == deadline+1) { // not yet visited
        EST_LST[neighbor].second = EST_LST[node].second - G.nodeWeights[neighbor];
      } else {
        EST_LST[neighbor].second = std::min(EST_LST[neighbor].second, EST_LST[node].second - G.nodeWeights[neighbor]);
      }
            
      out_degree[neighbor]--;
      if (out_degree[neighbor] == 0) {
        Q_LST.push(neighbor);
      }
    }
  }
}

void MultiMachineSchedulingWithoutFreezing::compute_EST_only() {
  EST_LST.resize(G.getNumNodes(), std::make_pair(0, deadline+1));

  // compute earliest start time by using an analogy to Kahn's algorithm
  std::queue<int> Q; 
  std::vector<int> in_degree(G.getNumNodes(), 0);
    
  // compute in-degree
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    in_degree[i] = G.adjList_inNeighbors[i].size();
  }
    
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    if (G.adjList_inNeighbors[i].size() == 0) {
      Q.push(i);
    }
  }

  while (!Q.empty()) {
    int node = Q.front(); 
    Q.pop();

    int dependency_constraint = 0; 
    for (size_t k = 0; k < G.adjList_inNeighbors[node].size(); k++) {
      int neighbor = G.adjList_inNeighbors[node][k].id;
      dependency_constraint = std::max(dependency_constraint, EST_LST[neighbor].first + G.nodeWeights[neighbor]);
    }
        
    EST_LST[node].first = dependency_constraint;

    // process children
    for (size_t k = 0; k < G.adjList_outNeighbors[node].size(); k++) {
      int neighbor = G.adjList_outNeighbors[node][k].id;
      in_degree[neighbor]--;
      if (in_degree[neighbor] == 0) {
        Q.push(neighbor);
      }
    }
  }
}

void MultiMachineSchedulingWithoutFreezing::init_basePower() {
  basePower = 0; 
  for (size_t l = 0; l < numProcessors; l++) {
    basePower += PStatic[l];
  }
}

bool MultiMachineSchedulingWithoutFreezing::is_schedule_valid(const bool experiments) {
  // for every task check deadline and precedence constraints
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    if (startingTimes[i] < 0) {
      if (!experiments) {
        std::cout << "Task " << i << " is not scheduled" << std::endl;
      }
      return false; 
    }

    if (G.nodeWeights[i] <= 0) {
      if (!experiments) {
        std::cout << "Task " << i << " has (less than) zero weight" << std::endl;    
      }
      return false;
    }

    if (startingTimes[i] + G.nodeWeights[i]-1 >= deadline) {
      if (!experiments) {
        std::cout << "Task " << i << " exceeds deadline" << std::endl;
        std::cout << "starts at " << startingTimes[i] << "; ends at " << startingTimes[i]+G.nodeWeights[i]-1 << std::endl;   
      }
      return false; 
    }
        
    for (auto neighbor : G.adjList_outNeighbors[i]) {
      if (startingTimes[i] + G.nodeWeights[i]-1 >= startingTimes[neighbor.id]) {
        if (!experiments) {
          std::cout << "Task " << i << " violates precedence constraint with successor " << neighbor.id << std::endl;
        }
      return false; 
      }
    }
  }

  return true; 
}

int MultiMachineSchedulingWithoutFreezing::get_tight_deadline() {
  int EFT = 0; 
  compute_EST_schedule(); 
  for (int task = 0; task < G.getNumNodes(); task++) {
    EFT = std::max(EFT, startingTimes[task] + G.nodeWeights[task]-1);
  }

  return EFT+1;   
}


int MultiMachineSchedulingWithoutFreezing::get_cost_of_schedule(const bool experiments) {
  unsigned int total_brown_power = 0;

  for (size_t l = 0; l < numSubintervals; l++) {
    std::vector<int> power_per_time_unit(subintervalLengths[l], basePower); 
    for (size_t k = 0; k < numProcessors; k++) {
      for (size_t i = 0; i < mapping[k].size(); i++) {
        int task = mapping[k][i]; 
        if (startingTimes[task] >= splitPoints[l+1]) {continue;} // we can not brake since the order is not necessarily correct in the communication case.
        if (!(startingTimes[task] + G.nodeWeights[task] - 1 < splitPoints[l])){
          int start = std::max(splitPoints[l], startingTimes[task]);
          int end = std::min(splitPoints[l+1], startingTimes[task] + G.nodeWeights[task]) - 1;
          for (size_t c = start; c <= end; c++) {
            power_per_time_unit[c-splitPoints[l]] += PDynamic[k];
          }
        }
      }
    }
    // print power per time unit
    for (size_t t = 0; t < power_per_time_unit.size(); t++) {
      int excess_power = power_per_time_unit[t] - subintervalBudgets[l];
      if (excess_power > 0) {
        total_brown_power += excess_power;
      }
    }
  }

  if (!experiments) {
    std::cout << "=====================================================================" << std::endl;
    std::cout << "            ADDITIONAL AMOUNT OF BROWN POWER = " << total_brown_power * BCC << "             " << std::endl; 
    std::cout << "=====================================================================" << std::endl;
    return total_brown_power * BCC;
  } else {
    return total_brown_power * BCC;
  }
}    

unsigned int MultiMachineSchedulingWithoutFreezing::get_cost_of_schedule_polynomial(const bool experiments) {

  unsigned int total_brown_power = 0;
    
  for (auto I : intervals) {
    std::vector<unsigned int> interval_tasks; 
    std::vector<int> refinedPoints; 

    refinedPoints.push_back(I.start);
    refinedPoints.push_back(I.end);

    for (size_t task = 0; task < G.getNumNodes(); ++task) {
      bool task_intersects = false;
      if (startingTimes[task] >= I.start && startingTimes[task] <= I.end) {
        refinedPoints.push_back(startingTimes[task]);
        task_intersects = true;
      } 
      if (startingTimes[task] + G.nodeWeights[task] >= I.start && startingTimes[task] + G.nodeWeights[task] <= I.end) {
        refinedPoints.push_back(startingTimes[task] + G.nodeWeights[task]);
        task_intersects = true;
      }
      // also add the task if it runs during the interval, but no extra points
      if (startingTimes[task] < I.start && startingTimes[task] + G.nodeWeights[task] >= I.end) {
        task_intersects = true;
      }
      if (task_intersects) {
        interval_tasks.push_back(task);
      }
    }

    // remove duplicates
    std::sort(refinedPoints.begin(), refinedPoints.end());
    refinedPoints.erase(std::unique(refinedPoints.begin(), refinedPoints.end()), refinedPoints.end());

    std::vector<std::pair<int, int>> inlcusive_intervals;
    for (size_t i = 0; i < refinedPoints.size()-1; ++i) {
      std::pair<int, int> inclInterval = std::make_pair(refinedPoints[i], refinedPoints[i+1]-1);
      inlcusive_intervals.push_back(inclInterval);
    }

    // +1 because the endpoints are now inclusive 
    std::vector<int> lengths(inlcusive_intervals.size());
    for (size_t i = 0; i < inlcusive_intervals.size(); ++i) {
      lengths[i] = inlcusive_intervals[i].second - inlcusive_intervals[i].first +1;
    }

    std::vector<unsigned int> power_per_interval(inlcusive_intervals.size(), basePower); 
    for (size_t i = 0; i < inlcusive_intervals.size(); ++i) {
      for (size_t k = 0; k < numProcessors; ++k) {
        for (size_t j = 0; j < interval_tasks.size(); ++j) {
          if (processors_of_tasks[interval_tasks[j]] == k) {
            int int_start = std::max(inlcusive_intervals[i].first, startingTimes[interval_tasks[j]]);
            int int_end = std::min(inlcusive_intervals[i].second, startingTimes[interval_tasks[j]] + G.nodeWeights[interval_tasks[j]]-1);
            if (int_start <= int_end) {
              power_per_interval[i] += PDynamic[k];
              break; // if there are more tasks we do not want to add it again we just want to observe whether to processor computes during this time because if that is the case it will for the whole interval
            }
          }
        }
      }
      int excess_power = power_per_interval[i] - I.interval_budget;
      if (excess_power > 0) {
         total_brown_power += excess_power * lengths[i];
      }
    }
  }
  if (!experiments) {
    std::cout << "=====================================================================" << std::endl;
    std::cout << "            ADDITIONAL AMOUNT OF BROWN POWER = " << total_brown_power * BCC << "             " << std::endl; 
    std::cout << "=====================================================================" << std::endl;
    return total_brown_power * BCC;
  } else {
    return total_brown_power * BCC;
  }
}


int MultiMachineSchedulingWithoutFreezing::get_cost_of_schedule_moreMem() {
  brownPower = 0;
  powerPerTimeUnit.resize(numSubintervals); 
  for (size_t l = 0; l < numSubintervals; l++) {
    powerPerTimeUnit[l].resize(subintervalLengths[l], basePower); 
    for (size_t k = 0; k < numProcessors; k++) {
      for (size_t i = 0; i < mapping[k].size(); i++) {
        int task = mapping[k][i]; 
        if (startingTimes[task] >= splitPoints[l+1]) {continue;} // we can not brake since the order is not necessarily correct in the communication case.
        if (!(startingTimes[task] + G.nodeWeights[task] - 1 < splitPoints[l])){
          int start = std::max(splitPoints[l], startingTimes[task]);
          int end = std::min(splitPoints[l+1], startingTimes[task] + G.nodeWeights[task]) - 1;
          for (size_t c = start; c <= end; c++) {
            powerPerTimeUnit[l][c-splitPoints[l]] += PDynamic[k];
          }
        }
      }
    }

    // print power per time unit
    for (size_t t = 0; t < powerPerTimeUnit[l].size(); t++) {
      int excess_power = powerPerTimeUnit[l][t] - subintervalBudgets[l];
      if (excess_power > 0) {
        brownPower += excess_power;
      }
    }
  }
  return brownPower * BCC;
}    

void MultiMachineSchedulingWithoutFreezing::compute_EST_schedule() {
  for (size_t task = 0; task < G.getNumNodes(); task++) {
    startingTimes[task] = EST_LST[task].first;
  }
}

void MultiMachineSchedulingWithoutFreezing::compute_LST_schedule() {
  for (size_t task = 0; task < G.getNumNodes(); task++) {
    startingTimes[task] = EST_LST[task].second;
  }
}

void MultiMachineSchedulingWithoutFreezing::print_schedule() {
  for (size_t task = 0; task < G.getNumNodes(); task++) {
    std::cout << "Task " << task << " runs: " << startingTimes[task] << " - " << startingTimes[task] + G.nodeWeights[task] -1 << " on processor " << processors_of_tasks[task] << std::endl;
  }
}

int MultiMachineSchedulingWithoutFreezing::getNumProcessors() {
  return numProcessors;
}

int MultiMachineSchedulingWithoutFreezing::getDeadline() {
  return deadline;
}

int MultiMachineSchedulingWithoutFreezing::getNumberOfSubIntervals() {
  return numSubintervals;
}

double MultiMachineSchedulingWithoutFreezing::getSetupTime() {
  return setupTimeCommunicationSetting;
}

int MultiMachineSchedulingWithoutFreezing::getBCC() {
  return BCC;
}

bool MultiMachineSchedulingWithoutFreezing::assign_block_beginning(block B, int point) {
  if (B.block_tasks.size() > 0) {
    if (point + B.running_time(*this) < deadline) {
      return true; 
    }
  }
  return false; 
}

bool MultiMachineSchedulingWithoutFreezing::assign_block_end(block B, int point) {
  if (B.block_tasks.size() > 0) {
    if (point - B.running_time(*this) >= 0) {
      return true; 
    }
  }
  return false; 
}

void MultiMachineSchedulingWithoutFreezing::compute_mixed_scores(bool weighted) {
  score_ordered_tasks.resize(G.getNumNodes());
  std::iota(score_ordered_tasks.begin(), score_ordered_tasks.end(), 0); 
    
  slack.resize(G.getNumNodes(), 0.0f);
  pressure.resize(G.getNumNodes(), 0.0f);
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    slack[i] = static_cast<float>(EST_LST[i].second) - static_cast<float>(EST_LST[i].first);
    pressure[i] = static_cast<float>(G.nodeWeights[i]) / (static_cast<float>(EST_LST[i].second + G.nodeWeights[i] - EST_LST[i].first));
  }

  unsigned int max_power = 0; 
  if (weighted) {
    for (int p = 0; p < numProcessors; p++) {
      max_power = std::max(max_power, PStatic[p] + PDynamic[p]);
    }
  }

  std::vector<float> mixed_score(G.getNumNodes(), 0.0f);
 
  for (size_t task = 0; task < G.getNumNodes(); task++) {
    if (weighted) {
      int pe = processors_of_tasks[task];
      if (slack[task] != 0) {
        mixed_score[task] = (0.3*pressure[task] + 0.7*(1/slack[task])) * ((static_cast<float>(PStatic[pe] + PDynamic[pe])) / static_cast<float>(max_power));
      } else {
        mixed_score[task] = pressure[task] * ((static_cast<float>(PStatic[pe] + PDynamic[pe])) / static_cast<float>(max_power));
      }
    } else {
      if (slack[task] != 0) {
        mixed_score[task] = (0.3*pressure[task] + 0.7*(1/slack[task]));
      } else {
        mixed_score[task] = pressure[task];
      }
    }
  }
  std::sort(score_ordered_tasks.begin(), score_ordered_tasks.end(), [&](int a, int b) {
    return mixed_score[a] > mixed_score[b];
  });
}


void MultiMachineSchedulingWithoutFreezing::compute_slack_machine_focused(bool weighted) {
  score_ordered_tasks.resize(G.getNumNodes());
  std::iota(score_ordered_tasks.begin(), score_ordered_tasks.end(), 0);

  // compute slack
  slack.resize(G.getNumNodes(), 0.0f);

  for (size_t i = 0; i < G.getNumNodes(); i++) {
    slack[i] = static_cast<float>(EST_LST[i].second) - EST_LST[i].first;
  }

  if (weighted) {
    unsigned int max_power = 0; 
    for (int p = 0; p < numProcessors; p++) {
      max_power = std::max(max_power, PStatic[p] + PDynamic[p]);
    }
    for (size_t i = 0; i < G.getNumNodes(); i++) {
      int pe = processors_of_tasks[i];
      slack[i] = slack[i] * (static_cast<float>(max_power) / static_cast<float>((PStatic[pe] + PDynamic[pe]))); 
    }    
  } 
  std::sort(score_ordered_tasks.begin(), score_ordered_tasks.end(), [&](int a, int b) {
    return slack[a] < slack[b];
  });
}

void MultiMachineSchedulingWithoutFreezing::compute_pressure_machine_focused(bool weighted) {
  score_ordered_tasks.resize(G.getNumNodes()); 
  std::iota(score_ordered_tasks.begin(), score_ordered_tasks.end(), 0); 

  // compute pressure
  pressure.resize(G.getNumNodes());

  for (size_t i = 0; i < G.getNumNodes(); i++) {
    pressure[i] = static_cast<float>(G.nodeWeights[i]) / (static_cast<float>(EST_LST[i].second + G.nodeWeights[i] - EST_LST[i].first));
  }
  
  if (weighted) {
    unsigned int max_power = 0; 
    for (int p = 0; p < numProcessors; p++) {
      max_power = std::max(max_power, PStatic[p] + PDynamic[p]);
    }  
    for (size_t i = 0; i < G.getNumNodes(); i++) {
      int pe = processors_of_tasks[i];
      pressure[i] = pressure[i] * ((static_cast<float>(PStatic[pe] + PDynamic[pe])) / static_cast<float>(max_power));
    }
  }
  std::sort(score_ordered_tasks.begin(), score_ordered_tasks.end(), [&](int a, int b) {
    return pressure[a] > pressure[b];
  });
}

void MultiMachineSchedulingWithoutFreezing::greedy_scheduling(bool refined) {
  // compute a topological order for updates using Kahn's algorithm
  std::queue<int> Q;
  std::vector<int> in_degree(G.getNumNodes(), 0);
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    in_degree[i] = G.adjList_inNeighbors[i].size();
    if (in_degree[i] == 0) {
      Q.push(i);
    }
  }
  while (!Q.empty()) {
    int node = Q.front();
    Q.pop();
    topological_order.push_back(node);
    for (auto neighbor : G.adjList_outNeighbors[node]) {
      in_degree[neighbor.id]--;
      if (in_degree[neighbor.id] == 0) {
        Q.push(neighbor.id);
      }
    }
  }
  
  std::vector<interval> refined_intervals;
  if (refined) {
    std::vector<std::vector<block>> blocks_per_processor(numProcessors); 
    std::vector<uint64_t> possible_start_times_for_each_task;
    for (int k = 0; k < numProcessors; k++) {
      for (int t = 0; t < mapping[k].size(); t++) {
        for (int b = 0; (t + b) < mapping[k].size() && b < 3; b++) {
          block B; 
          for (int j = t; j <= t+b; j++) {
            B.block_tasks.push_back(mapping[k][j]);
          }
          blocks_per_processor[k].push_back(B);
        }
      }
    }

    for (int k = 0; k < numProcessors; k++) {
      for (auto B : blocks_per_processor[k]) {
        for (auto point : splitPoints) {
          if (assign_block_beginning(B, point)) {
            int task_start_time_b = point;
            for (auto task : B.block_tasks) {
              possible_start_times_for_each_task.push_back(task_start_time_b);
              task_start_time_b += G.nodeWeights[task]; 
            }
          }
          if (assign_block_end(B, point)) {
            int task_start_time_e = point - G.nodeWeights[B.block_tasks[B.block_tasks.size()-1]];
            for (auto task : B.block_tasks) {
              possible_start_times_for_each_task.push_back(task_start_time_e);
              task_start_time_e -= G.nodeWeights[task]; 
            }
          }
        }
      }
    }

    // add splitting points for edge cases
    for (size_t j = 0; j < splitPoints.size(); ++j) {
      possible_start_times_for_each_task.push_back(splitPoints[j]); 
    }
           
    // apparently three tasks are enough to get a really large amount of possible start times
    std::sort(possible_start_times_for_each_task.begin(), possible_start_times_for_each_task.end());
    possible_start_times_for_each_task.erase(std::unique(possible_start_times_for_each_task.begin(), possible_start_times_for_each_task.end()), possible_start_times_for_each_task.end());
            
    refined_intervals.resize(possible_start_times_for_each_task.size()-1);
    for (int i = 0; i < possible_start_times_for_each_task.size()-1; i++) {
      interval I; 
      I.interval_id = i; 
      I.start = possible_start_times_for_each_task[i]; 
      I.end = possible_start_times_for_each_task[i+1]; 
      I.interval_length = possible_start_times_for_each_task[i+1] - possible_start_times_for_each_task[i]; 
      for (auto inter : intervals) {
        // other cases do not occur since original splitting points are contained by construction for workflows that can be scheduled
        if (I.start >= inter.start && I.end <= inter.end) {
          I.interval_budget = inter.interval_budget;
          break; 
        }
      }     
      refined_intervals[i] = I;
    }
  } else {
    // use input intervals
    refined_intervals = intervals;
  }

  // FIFO is correct here
  std::queue<int> taskQ; 
  for (auto task : score_ordered_tasks) {
    taskQ.push(task);
  }

  while (!taskQ.empty()) {
    int task = taskQ.front();
    taskQ.pop();

    // note that refined intervals are sorted by start time
    std::vector<interval> possible_intervals;
    for (const auto I : refined_intervals) {
      const auto point = I.start; 
      if (point <= EST_LST[task].second) {
        if (point >= EST_LST[task].first) {
          possible_intervals.push_back(I);
        }
      } else {break;}
    }

    std::sort(possible_intervals.begin(), possible_intervals.end(), [](const interval &a, const interval &b) {
      return a.interval_budget > b.interval_budget; 
    });

    std::vector<interval> best_intervals;
    if (possible_intervals.size() > 0) {
      interval max_int = possible_intervals[0];
      best_intervals.push_back(max_int);
      for (size_t i = 1; i < possible_intervals.size(); i++) {
        if (possible_intervals[i].interval_budget == max_int.interval_budget) {
          best_intervals.push_back(possible_intervals[i]);
        } else {
          break;
        }
      }
    }

    std::sort(best_intervals.begin(), best_intervals.end(), [](const interval &a, const interval &b) {
      return a.start < b.start; 
    });
        
    if (best_intervals.size() > 0) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<int> dis(0, best_intervals.size()-1);
      interval max_int = best_intervals[0];
            
      startingTimes[task] = max_int.start; 
      EST_LST[task].first = max_int.start;
      EST_LST[task].second = max_int.start;
      already_scheduled[task] = true;
    } else {         
      // use EST for now.
      startingTimes[task] = EST_LST[task].first;
      EST_LST[task].first = startingTimes[task]; 
      EST_LST[task].second = startingTimes[task];
      already_scheduled[task] = true;
    }

    // update EST
    bool skip_EST = true; 
    for (int i = 0; i < topological_order.size(); i++) {
      int node = topological_order[i];
      if (node == task && skip_EST == true) {skip_EST = false;}
      if (skip_EST || already_scheduled[node]) {continue;}
      int max = 0; 
      for (auto neighbor : G.adjList_inNeighbors[node]) {
        max = std::max(max, EST_LST[neighbor.id].first + G.nodeWeights[neighbor.id]);
      }
      EST_LST[node].first = max;
    }

    // update LST
    bool skip_LST = true;
    for (int i = topological_order.size()-1; i >= 0; i--) {
      int node = topological_order[i];
      if (node == task && skip_LST == true) {skip_LST = false;}
      if (skip_LST || already_scheduled[node]) {continue;}
      int min = deadline; 
      for (auto neighbor : G.adjList_outNeighbors[node]) {
        min = std::min(min, EST_LST[neighbor.id].second);
      }
      EST_LST[node].second = min - G.nodeWeights[node];
    }
  }
}

void MultiMachineSchedulingWithoutFreezing::greedy_scheduling_with_budget_maintenance(bool refined) {
  // compute a topological order for updates
  std::queue<int> Q;
  std::vector<int> in_degree(G.getNumNodes(), 0);
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    in_degree[i] = G.adjList_inNeighbors[i].size();
    if (in_degree[i] == 0) {
      Q.push(i);
    }
  }
    
  while (!Q.empty()) {
    int node = Q.front();
    Q.pop();
    topological_order.push_back(node);
    for (auto neighbor : G.adjList_outNeighbors[node]) {
      in_degree[neighbor.id]--;
      if (in_degree[neighbor.id] == 0) {
        Q.push(neighbor.id);
      }
    }
  }

  std::vector<interval> refined_intervals;
  std::vector<int> budget_per_interval; 
  if (refined) {
    std::vector<std::vector<block>> blocks_per_processor(numProcessors); 
    std::vector<uint64_t> possible_start_times_for_each_task;
    for (int k = 0; k < numProcessors; k++) {
      for (int t = 0; t < mapping[k].size(); t++) {
        for (int b = 0; (t + b) < mapping[k].size() && b < 3; b++) {
          block B; 
          for (int j = t; j <= t+b; j++) {
            B.block_tasks.push_back(mapping[k][j]);
          }
          blocks_per_processor[k].push_back(B);
        }
      }
    }

    for (int k = 0; k < numProcessors; k++) {
      for (auto B : blocks_per_processor[k]) {
        for (auto point : splitPoints) {
          if (assign_block_beginning(B, point)) {
            int task_start_time_b = point;
            for (auto task : B.block_tasks) {
              possible_start_times_for_each_task.push_back(task_start_time_b);
              task_start_time_b += G.nodeWeights[task]; 
            }
          }

          if (assign_block_end(B, point)) {
            int task_start_time_e = point - G.nodeWeights[B.block_tasks[B.block_tasks.size()-1]];
            for (auto task : B.block_tasks) {
              possible_start_times_for_each_task.push_back(task_start_time_e);
              task_start_time_e -= G.nodeWeights[task]; 
            }
          }
        }
      }
    }

    // add splitting points for edge cases
    for (size_t j = 0; j < splitPoints.size(); ++j) {
      possible_start_times_for_each_task.push_back(splitPoints[j]); 
    }
  
    // apparently three tasks are enough to get a really large amount of possible start times
    std::sort(possible_start_times_for_each_task.begin(), possible_start_times_for_each_task.end());
    possible_start_times_for_each_task.erase(std::unique(possible_start_times_for_each_task.begin(), possible_start_times_for_each_task.end()), possible_start_times_for_each_task.end());

    refined_intervals.resize(possible_start_times_for_each_task.size()-1);
    for (int i = 0; i < possible_start_times_for_each_task.size()-1; i++) {
      interval I; 
      I.interval_id = i; 
      I.start = possible_start_times_for_each_task[i]; 
      I.end = possible_start_times_for_each_task[i+1]; 
      I.interval_length = possible_start_times_for_each_task[i+1] - possible_start_times_for_each_task[i]; 
      for (auto inter : intervals) {
        if (I.start >= inter.start && I.end <= inter.end) {
          I.interval_budget = inter.interval_budget;
          break; 
        }
      }     
      refined_intervals[i] = I;
    }
  } else {
    refined_intervals = intervals;
  }

  // FIFO is correct here
  std::queue<int> taskQ; 
  for (auto task : score_ordered_tasks) {
    taskQ.push(task);
  }

  while (!taskQ.empty()) {
    int task = taskQ.front();
    taskQ.pop();
    std::vector<interval> possible_intervals;
    for (int r = 0; r < refined_intervals.size(); r++) {
      const auto I = refined_intervals[r];
      const auto point = I.start; 
      if (point <= EST_LST[task].second) {
        if (point >= EST_LST[task].first) {
          possible_intervals.push_back(I);  
        }
      } else {break;}
    }

    // note that the order in possible_intervals
    std::sort(possible_intervals.begin(), possible_intervals.end(), [](const interval &a, const interval &b) {
      return a.interval_budget > b.interval_budget; 
    });

    std::vector<interval> best_intervals;
    if (possible_intervals.size() > 0) {
      interval max_int = possible_intervals[0];
      best_intervals.push_back(max_int);
      for (size_t i = 1; i < possible_intervals.size(); i++) {
        if (possible_intervals[i].interval_budget == possible_intervals[0].interval_budget) {
          best_intervals.push_back(possible_intervals[i]);
        } else {
          break;
        }
      }
    }   

    std::sort(best_intervals.begin(), best_intervals.end(), [](const interval &a, const interval &b) {
      return a.start < b.start; 
    });
    
    if (best_intervals.size() > 0) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dis(0, best_intervals.size()-1);
      interval max_int = best_intervals[0];
      startingTimes[task] = max_int.start; 
      EST_LST[task].first = max_int.start;
      EST_LST[task].second = max_int.start;
      already_scheduled[task] = true;
    } else {
      // use EST for now.
      startingTimes[task] = EST_LST[task].first;
      EST_LST[task].first = startingTimes[task]; 
      EST_LST[task].second = startingTimes[task];
      already_scheduled[task] = true;
    }

    // update budget per interval --> can be more precise
    for (size_t i = 0; i < refined_intervals.size(); i++) {
      if (startingTimes[task] >= refined_intervals[i].start && startingTimes[task] < refined_intervals[i].end) {
        refined_intervals[i].interval_budget -= (PStatic[processors_of_tasks[task]] + PDynamic[processors_of_tasks[task]]);
      }
    }

    // update EST
    bool skip_EST = true; 
    for (int i = 0; i < topological_order.size(); i++) {
      int node = topological_order[i];
      if (node == task && skip_EST == true) {skip_EST = false;}
      if (skip_EST || already_scheduled[node]) {continue;}
      int max = 0; 
      for (auto neighbor : G.adjList_inNeighbors[node]) {
        max = std::max(max, EST_LST[neighbor.id].first + G.nodeWeights[neighbor.id]);
      }
      EST_LST[node].first = max;
    }

    bool skip_LST = true;
    for (int i = topological_order.size()-1; i >= 0; i--) {
      int node = topological_order[i];
      if (node == task && skip_LST == true) {skip_LST = false;}
      if (skip_LST || already_scheduled[node]) {continue;}
      int min = deadline; 
      for (auto neighbor : G.adjList_outNeighbors[node]) {
        min = std::min(min, EST_LST[neighbor.id].second);
      }
      EST_LST[node].second = min - G.nodeWeights[node];
    }
  }
}

void MultiMachineSchedulingWithoutFreezing::greedy_scheduling_with_budget_maintenance_enhanced(bool refined) {
  // compute a topological order for updates
  std::queue<int> Q;
  std::vector<int> in_degree(G.getNumNodes(), 0);
  for (size_t i = 0; i < G.getNumNodes(); i++) {
    in_degree[i] = G.adjList_inNeighbors[i].size();
    if (in_degree[i] == 0) {
      Q.push(i); 
    }
  }

  while (!Q.empty()) {
    int node = Q.front();
    Q.pop();
    topological_order.push_back(node);
    for (auto neighbor : G.adjList_outNeighbors[node]) {
      in_degree[neighbor.id]--;
      if (in_degree[neighbor.id] == 0) {
        Q.push(neighbor.id);
      }
    }
  }

  dynamicIntervalVector refined_intervals;
  std::vector<int> budget_per_interval; 
  if (refined) {
    std::vector<std::vector<block>> blocks_per_processor(numProcessors); 
    std::vector<uint64_t> possible_start_times_for_each_task;
    for (int k = 0; k < numProcessors; k++) {
      for (int t = 0; t < mapping[k].size(); t++) {
        for (int b = 0; (t + b) < mapping[k].size() && b < 3; b++) {
          block B; 
          for (int j = t; j <= t+b; j++) {
            B.block_tasks.push_back(mapping[k][j]);
          }
          blocks_per_processor[k].push_back(B);
        }
      }
    }

    for (int k = 0; k < numProcessors; k++) {
      for (auto B : blocks_per_processor[k]) {
        for (auto point : splitPoints) {
          if (assign_block_beginning(B, point)) {
            int task_start_time_b = point;
            for (auto task : B.block_tasks) {
              possible_start_times_for_each_task.push_back(task_start_time_b);
              task_start_time_b += G.nodeWeights[task]; 
            }
          }
          
          if (assign_block_end(B, point)) {
            int task_start_time_e = point - G.nodeWeights[B.block_tasks[B.block_tasks.size()-1]];
            for (auto task : B.block_tasks) {
              possible_start_times_for_each_task.push_back(task_start_time_e);
              task_start_time_e -= G.nodeWeights[task]; 
            }
          }
        }
      }
    }
  
    // add splitting points for edge cases
    for (size_t j = 0; j < splitPoints.size(); ++j) {
      possible_start_times_for_each_task.push_back(splitPoints[j]); 
    }

    // apparently three tasks are enough to get a really large amount of possible start times
    std::sort(possible_start_times_for_each_task.begin(), possible_start_times_for_each_task.end());
    possible_start_times_for_each_task.erase(std::unique(possible_start_times_for_each_task.begin(), possible_start_times_for_each_task.end()), possible_start_times_for_each_task.end());

    for (int i = 0; i < possible_start_times_for_each_task.size()-1; i++) {
      interval I; 
      I.interval_id = i; 
      I.start = possible_start_times_for_each_task[i]; 
      I.end = possible_start_times_for_each_task[i+1]; 
      I.interval_length = possible_start_times_for_each_task[i+1] - possible_start_times_for_each_task[i]; 
      for (auto inter : intervals) {
        if (I.start >= inter.start && I.end <= inter.end) {
          I.interval_budget = inter.interval_budget;
          break; 
        }
      }     
      refined_intervals.push_back(I);
    }
  } else {
    for (auto interv : intervals) {
      interval I;
      I.interval_id = interv.interval_id;
      I.start = interv.start;
      I.end = interv.end;
      I.interval_length = interv.interval_length;
      I.interval_budget = interv.interval_budget;
      refined_intervals.push_back(I);
    }
  }
 
  // FIFO is correct here
  std::queue<int> taskQ; 
  for (auto task : score_ordered_tasks) {
    taskQ.push(task);
  }

  while (!taskQ.empty()) {
    int task = taskQ.front();
    taskQ.pop();
    std::vector<interval> possible_intervals;
    // Note that refined intervals is also sorted automatically
    for (auto it = refined_intervals.begin(); it != refined_intervals.end(); ++it) {
      const interval I = *it;
      const auto point = I.start; 
      if (point <= EST_LST[task].second) {
        if (point >= EST_LST[task].first) {
          possible_intervals.push_back(I);  // now decision on the fly or store the interval for all tasks first?
        }
      } else {break;}
    }

    std::sort(possible_intervals.begin(), possible_intervals.end(), [](const interval &a, const interval &b) {
      return a.interval_budget > b.interval_budget; 
    });

    std::vector<interval> best_intervals;
    if (possible_intervals.size() > 0) {
      interval max_int = possible_intervals[0];
      best_intervals.push_back(max_int);
      for (size_t i = 1; i < possible_intervals.size(); i++) {
        if (possible_intervals[i].interval_budget == possible_intervals[0].interval_budget) {
          best_intervals.push_back(possible_intervals[i]);
        } else {
          break;
        }
      }
    }   
    
    std::sort(best_intervals.begin(), best_intervals.end(), [](const interval &a, const interval &b) {
      return a.start < b.start; 
    });
     
    if (best_intervals.size() > 0) {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dis(0, best_intervals.size()-1);
      interval max_int = best_intervals[0];
      
      startingTimes[task] = max_int.start; 
      EST_LST[task].first = max_int.start;
      EST_LST[task].second = max_int.start;
      already_scheduled[task] = true;
    } else {     
      // use EST for now.
      startingTimes[task] = EST_LST[task].first;
      EST_LST[task].first = startingTimes[task]; 
      EST_LST[task].second = startingTimes[task];
      already_scheduled[task] = true;
    }

    int splitPoint1 = startingTimes[task]; // not sure whether this is needed
    int splitPoint2 = startingTimes[task] + G.nodeWeights[task];

    refined_intervals.split_interval(splitPoint1);
    refined_intervals.split_interval(splitPoint2);
  
    for (auto it = refined_intervals.begin(); it != refined_intervals.end(); ++it) {
      if (startingTimes[task] == it->start || it->end == startingTimes[task]+G.nodeWeights[task]) {
        interval* modif = const_cast<interval*>(&(*it)); 
        modif->interval_budget -= (PStatic[processors_of_tasks[task]] + PDynamic[processors_of_tasks[task]]);
      }
    }

    // update EST
    bool skip_EST = true; 
    for (int i = 0; i < topological_order.size(); i++) {
      int node = topological_order[i];
      if (node == task && skip_EST == true) {skip_EST = false;}
      if (skip_EST || already_scheduled[node]) {continue;}
      int max = 0; 
      for (auto neighbor : G.adjList_inNeighbors[node]) {
        max = std::max(max, EST_LST[neighbor.id].first + G.nodeWeights[neighbor.id]);
      }
      EST_LST[node].first = max;
    }

    // update LST   
    bool skip_LST = true;
    for (int i = topological_order.size()-1; i >= 0; i--) {
      int node = topological_order[i];
      if (node == task && skip_LST == true) {skip_LST = false;}
      if (skip_LST || already_scheduled[node]) {continue;}
      int min = deadline; 
      for (auto neighbor : G.adjList_outNeighbors[node]) {
        min = std::min(min, EST_LST[neighbor.id].second);
      }
      EST_LST[node].second = min - G.nodeWeights[node];
    }
  }
}

bool MultiMachineSchedulingWithoutFreezing::is_move_possible(int task, int pot_time) {  
  if (pot_time < 0 || pot_time + G.nodeWeights[task]-1 >= deadline) {return false;}
  bool possible = true;
  for (int i = 0; i < G.adjList_inNeighbors[task].size(); i++) {
    int neighbor = G.adjList_inNeighbors[task][i].id;
    if (startingTimes[neighbor] + G.nodeWeights[neighbor] > pot_time) {
      possible = false;
      break;
    }
  }
  if (!possible) {return false;}
  for (int i = 0; i < G.adjList_outNeighbors[task].size(); i++) {
    int neighbor = G.adjList_outNeighbors[task][i].id;
    if (pot_time + G.nodeWeights[task] > startingTimes[neighbor]) {
      possible = false;
      break;
    }
  }
  if (!possible) {return false;}

  return true; 
}

void MultiMachineSchedulingWithoutFreezing::LS_machine_focused_fastUpdates(int scope) {
  // initialize data structures
  get_cost_of_schedule_moreMem(); 

  std::vector<int> sorted_processors(numProcessors); 
  std::iota(sorted_processors.begin(), sorted_processors.end(), 0); 

  std::sort(sorted_processors.begin(), sorted_processors.end(), [&](int a, int b) {
    return PDynamic[a] + PStatic[a] > PDynamic[b] + PStatic[b];
  });
  
  // for random order of tasks
  std::random_device rd; 
  std::mt19937 gen(rd()); 

  bool found_gain = true; 
  while (found_gain) {
    found_gain = false;
    for (auto pe : sorted_processors) {
      std::vector<int> random_order_of_tasks;
      for (int i = 0; i < mapping[pe].size(); i++) {
        random_order_of_tasks.push_back(mapping[pe][i]);
      }

      // Optional: Use a random order
      // std::shuffle(random_order_of_tasks.begin(), random_order_of_tasks.end(), gen);

      for (auto task : random_order_of_tasks) {
        int time = startingTimes[task];
        int best_move = time;
        int best_gain_so_far = 0;
        int bestCost = 0; 
        for (int t = time - scope; t <= time + scope; t++) {
          int new_time = t;
          
          if (!is_move_possible(task, new_time)) {continue;}

          // now we know the move is possible and we check the gain 
          std::pair<int, int> gain_newCost = is_gain_fast(task, new_time, time, processors_of_tasks[task]);
          
          if (gain_newCost.first > 0) {
            if (gain_newCost.first > best_gain_so_far) {
              best_gain_so_far = gain_newCost.first;
              best_move = new_time;
              bestCost = gain_newCost.second; 
            }
            break; 
          } else {
            continue;
          }
        }
        if (best_gain_so_far > 0) { 
          startingTimes[task] = best_move;
          brownPower = bestCost; 
          found_gain = true;
          apply_move(task, best_move, time, pe); 
        }
      }
    }
  }
}

std::pair<int,int> MultiMachineSchedulingWithoutFreezing::is_gain_fast(int task, int new_time, int old_time, int pe) {
  int gain = 0;
  int oldBrownPower = static_cast<int>(brownPower); // int und unsigned int grrrrr... 

  int newBrownPower = 0; 

  // further optimization: identify only the intervals that are used. 

  for (size_t l = 0; l < numSubintervals; l++) {
    if (old_time >= splitPoints[l+1]) {continue;} 
    if (!(old_time + G.nodeWeights[task] - 1 < splitPoints[l])){
      int start = std::max(splitPoints[l], old_time);
      int end = std::min(splitPoints[l+1], old_time + G.nodeWeights[task]) - 1;
      for (size_t c = start; c <= end; c++) {
        powerPerTimeUnit[l][c-splitPoints[l]] -= PDynamic[pe];
      }
    }
  }

  for (size_t l = 0; l < numSubintervals; l++) {
    if (new_time >= splitPoints[l+1]) {continue;} 
    if (!(new_time + G.nodeWeights[task] - 1 < splitPoints[l])){
      int start = std::max(splitPoints[l], new_time);
      int end = std::min(splitPoints[l+1], new_time + G.nodeWeights[task]) - 1;
      for (size_t c = start; c <= end; c++) {
        powerPerTimeUnit[l][c-splitPoints[l]] += PDynamic[pe];
      }
    }
  }
    
  for (size_t l = 0; l < numSubintervals; l++) {
    for (size_t t = 0; t < powerPerTimeUnit[l].size(); t++) {
      int excess_power = powerPerTimeUnit[l][t] - subintervalBudgets[l];
      if (excess_power > 0) {
        newBrownPower += excess_power;
      }
    }
  }

  // discard changes
  for (size_t l = 0; l < numSubintervals; l++) {
    if (new_time >= splitPoints[l+1]) {continue;} 
    if (!(new_time + G.nodeWeights[task] - 1 < splitPoints[l])){
      int start = std::max(splitPoints[l], new_time);
      int end = std::min(splitPoints[l+1], new_time + G.nodeWeights[task]) - 1;
      for (size_t c = start; c <= end; c++) {
        powerPerTimeUnit[l][c-splitPoints[l]] -= PDynamic[pe];
      }
    }
  }

  for (size_t l = 0; l < numSubintervals; l++) {
    if (old_time >= splitPoints[l+1]) {continue;} 
    if (!(old_time + G.nodeWeights[task] - 1 < splitPoints[l])){
      int start = std::max(splitPoints[l], old_time);
      int end = std::min(splitPoints[l+1], old_time + G.nodeWeights[task]) - 1;
      for (size_t c = start; c <= end; c++) {
        powerPerTimeUnit[l][c-splitPoints[l]] += PDynamic[pe];
      }
    }
  }
  
  gain = oldBrownPower - newBrownPower;

  return std::make_pair(gain, newBrownPower); 
}

unsigned int MultiMachineSchedulingWithoutFreezing::getCost() {
  return brownPower * BCC; 
}

void MultiMachineSchedulingWithoutFreezing::apply_move(int task, int new_time, int old_time, int pe) {
  for (size_t l = 0; l < numSubintervals; l++) {
    if (old_time >= splitPoints[l+1]) {continue;} 
    if (!(old_time + G.nodeWeights[task] - 1 < splitPoints[l])){
      int start = std::max(splitPoints[l], old_time);
      int end = std::min(splitPoints[l+1], old_time + G.nodeWeights[task]) - 1;
      for (size_t c = start; c <= end; c++) {
        powerPerTimeUnit[l][c-splitPoints[l]] -= PDynamic[pe];
      }
    }
  }

  for (size_t l = 0; l < numSubintervals; l++) {
    if (new_time >= splitPoints[l+1]) {continue;}
    if (!(new_time + G.nodeWeights[task] - 1 < splitPoints[l])){
      int start = std::max(splitPoints[l], new_time);
      int end = std::min(splitPoints[l+1], new_time + G.nodeWeights[task]) - 1;
      for (size_t c = start; c <= end; c++) {
        powerPerTimeUnit[l][c-splitPoints[l]] += PDynamic[pe];
      }
    }
  }
}
