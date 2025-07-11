#include <unordered_map>
#include <queue>
#include <vector>
#include <functional>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <set>
#include <iterator>
#include <optional>


/*!
* Interval data structure for the subdivision of the time horizon. 
*/
struct interval {
  int interval_id;
  int interval_budget; 
  int interval_length; 
  int start; 
  int end;  

  bool operator==(const interval& other) const {
    return start == other.start && end == other.end;
  } 

  /*!
  * Intervals are sorted according to their start time
  */
  bool operator<(const interval& other) const {
    return start < other.start;
  }
};

/*!
* Stores a vertex with weight and type 
* In case it is a communication vertex, we store the original edge
*/
struct vertex {
  int id; 
  int weight;
  std::string type; // "comp" or "comm"
  std::pair<int, int> comm_id; // = -1,-1 if type = comp
}; 

struct costInterval {
  std::pair<int, int> Int; // start and end of the interval
  std::vector<int> active_processors; 
};

/*!
* Data structure to support the dynamic interval splitting in our CaWoSched framework.
* The data structure makes sure the intervals are sorted according to their start time at every point in time.
*/
class dynamicIntervalVector {
  public: 
    // Iterator
    using iterator = std::set<interval>::iterator;
    iterator begin() { return intervalSet.begin(); }
    iterator end() { return intervalSet.end(); }

    // Const iterator
    using const_iterator = std::set<interval>::const_iterator;
    const_iterator begin() const { return intervalSet.begin(); }
    const_iterator end() const { return intervalSet.end(); }

    void push_back(const interval& new_interval) {
      intervalSet.insert(new_interval); 
    }

    /*!
    * This function splits the interval at splitPoint, where splitPoint becomes the end of the first and start of the second interval
    * @param splitPoint Determines where the interval is split
    */
    void split_interval(int splitPoint) {
      auto it = intervalSet.begin();
      auto it_next = intervalSet.begin(); 

      std::advance(it_next, 1); 

      for (; it_next != intervalSet.end() && it != intervalSet.end(); ++it, ++it_next) {
        interval current = *it; 
        interval next = *it_next;

        // if split point is start or end -> no action
        if (current.start == splitPoint || current.end == splitPoint) {
          break; 
        } 

        // The budget is adapted by the algo
        if (current.start < splitPoint && next.start > splitPoint) {
          interval refined_first = {current.interval_id, current.interval_budget, splitPoint - current.start, current.start, splitPoint}; 
          interval refined_second = {current.interval_id, current.interval_budget, current.end - splitPoint, splitPoint, current.end};
          it = intervalSet.erase(it); 
          intervalSet.insert(refined_first); 
          intervalSet.insert(refined_second); 
        }
      }
    }

  private: 
    std::set<interval> intervalSet;  
};

/*!
* Data structure used to store the input workflow 
*/
class DAG {
  public: 
    void readGraph(const std::string& filename); 

    void count_Nodes_Edges(const std::string& filename); 

    void printGraph() const;

    unsigned int getNumNodes() const {
      return numNodes;
    }

    friend struct vertex;
    friend class MultiMachineSchedulingWithoutFreezing;

  private:
    std::vector<std::vector<vertex>> adjList_inNeighbors;
    std::vector<std::vector<vertex>> adjList_outNeighbors; 
    std::vector<int> nodeWeights;

    unsigned int numNodes; 
    unsigned int numEdges; 
};

////////////////////////////////////////////////////////////////
// SCHEDULER //
///////////////////////////////////////////////////////////////

/*!
* This is the scheduler class where the actual carbon-aware scheduling happens. This implements the CaWoSched framework.
*/
class MultiMachineSchedulingWithoutFreezing {
  public: 
    /*!
    * This is the data structure for the task blocks that are used to refine the initial interval set
    */
    struct block {
      std::vector<int> block_tasks;

      void print() const {
        for (const auto task : block_tasks) {
          std::cout << task << " ";
        }
        std::cout << std::endl;
      }

      // return the total running time of the block
      int running_time(MultiMachineSchedulingWithoutFreezing& obj) {
        unsigned int rt = 0;
        for (auto task : block_tasks) {
          rt += obj.G.nodeWeights[task];
        } 
        return rt;
      }
    };

    friend struct block;

    /*
    *
    * INITIAL SETUP
    *
    */

    /*!
    * This is the constructor for the carbon aware scheduler
    */
    MultiMachineSchedulingWithoutFreezing(const std::string& DAG_file, const std::string& mapping_file, const std::string& cluster_setup_file, int communication_mode);
        
    /*!
    * This is the constructor for the HEFT baseline, i.e. to extract the tight deadline
    */
    MultiMachineSchedulingWithoutFreezing(const std::string& DAG_file, const std::string& mapping_file, int communication_mode);
        
    /*!
    * This is the constructor to use when looking at the workflow is the only goal
    */
    MultiMachineSchedulingWithoutFreezing(const std::string& DAG_file);

    /*!
    * The input file should contain the number of processors and the mapping 
    */
    void readMapping(const std::string& filename);

    /*!
    * The input file should contain:
    * - deadline
    * - number of subintervals
    * - values for P_idle
    * - values for P_work
    * - the brown carbon cost BCC
    * - subinterval budgets
    * - split points of the time horizon
    */
    void readSetup(const std::string& filename);

    /*!
    * This function adds paths on the communication processors to enforce the ordering constraint on the processors
    */
    void addMappingPaths();

    /*!
    * This function computes for every task the earliest start time and the earliest finish time
    */
    void compute_EST_LST();

    void compute_EST_only();

    /*!
    * This function constructs the communication enhanced DAG
    */
    void setup_communication_setting(); 

    /*!
    * Use this function only if you want the thight deadline using the communication enhanced DAG
    */
    void setup_communication_setting_for_EST_only();

    /*!
    * Computes the power which the cluster consumes at any point in time 
    */
    void init_basePower(); 

    /*
    *
    * COST, CORRECTNESS AND DEADLINE
    * 
    */

    /*!
    * This function checks the validity of a schedule 
    * @param experiments If the parameter is true, the function just returns a boolean without printing 
    */
    bool is_schedule_valid(const bool experiments = false); 

    /*!
    * This function computes the makespan of the schedule for the given mapping
    * To be used with the corresponding constructor
    */
    int get_tight_deadline(); 

    /*!
    * This function computes the carbon cost of the computed schedule
    * @param experiments If the parameter is true, the function just returns the cost without printing
    */
    int get_cost_of_schedule(const bool experiments = false);

    /*!
    * This function computes the carbon cost of the computed schedule in polynomial time using the interval formulation for the cost function
    * @param experiments If the parameter is true, the function just returns the cost without printing
    */
    unsigned int get_cost_of_schedule_polynomial(const bool experiments = false); 

    /*!
    * This function computes the carbon cost before the local search starts and stores more information
    * for faster updates
    */
    int get_cost_of_schedule_moreMem(); 

    /*!
    * This function prints the computed schedule along with the mapping
    */
    void print_schedule(); 

    /*
    *
    * GREEDY SCHEDULING
    * 
    */

    ////////////////////////////
    ////////// SCORES //////////
    //////////////////////////// 
        
    /*!
    * This function computes the slack for each task and sorts the tasks accordingly 
    * @param weighted If true, the scores is multiplied by a weight the accounts for the power heterogenity of the cluster
    */
    void compute_pressure_machine_focused(bool weighted);

    /*!
    * This function computes the pressure for each task and sorts the tasks accordingly
    * @param weighted If true, the scores is multiplied by a weight the accounts for the power heterogenity of the cluster
    */
    void compute_slack_machine_focused(bool weighted);        
        
    /*! 
    * This function computes a mix of slack and pressure and sorts the tasks accordingly
    * This was only used in preliminary experiments 
    * @param weighted If true, the score is multiplied by a weight factor that accounts for the power heterogenity of the cluster
    */
    void compute_mixed_scores(bool weighted);
    
    ////////////////////////////////
    ////////// SCHEDULING //////////
    ////////////////////////////////

    /*!
    * This function computes the schedule corresponding to the fixed mapping given by HEFT
    * Hence, this is the schedule given by the baseline HEFT
    */
    void compute_EST_schedule();

    /*!
    * This function computes a schedule by scheduling each task as late as possible
    */
    void compute_LST_schedule();

    /*!
    * This function tests whether block B can be scheduled starting from time unit 'point', i.e., does not exceed the deadline
    * The function is needed for the interval refinement
    * @param B Describes the block which should be tentatively scheduled
    * @param point This is the time unit at which the block should start
    */
    bool assign_block_beginning(block B, int point);

    /*!
    * This function tests whether block B can be scheduled such that the last occupied time unit is point-1, i.e., does not have to start before time unit 0
    * The function is needed for the interval refinement
    * @param B describes the block which should bhe tentatively scheduled
    * @param point point-1 is the last occupied time unit by the block
    */
    bool assign_block_end(block B, int point);

    /*!
    * This is the greedy scheduling heuristic without budget maintenance 
    * @param refined If true, the algorithm first refines the input intervals
    */
    void greedy_scheduling(bool refined = false);

    /*!
    * This is the greedy scheduling heuristic with a coarse budget maintenance, i.e. the interval is not further split
    * @param refined If true, the algorithm first refines the input intervals
    */
    void greedy_scheduling_with_budget_maintenance(bool refined = false);

    /*!
    * This is the most sophisticated form of the greedy scheduling heuristic with budget maintenance where we furhter split the interval if appropriate
    * @param refined If true, the algorithm first refines the input intervals
    */
    void greedy_scheduling_with_budget_maintenance_enhanced(bool refined = false);  

    /*
    *
    * LOCAL SEARCH
    *
    */

    /*!
    * This is the local search for the scheduling heuristic. It is designed as a hill climber.
    * @param scope This parameter determines the size of the neighborhood that the hill climber looks at.
    */
    void LS_machine_focused_fastUpdates(int scope); 

    /*!
    * This function checks whether moving task to new_time on processor pe brings a gain
    * The first value is the gain, the second value is the new carbon cost
    * @param task The task the should be moved
    * @param new_time The tentative new time for the task
    * @param old_time The time the task is currently scheduled at
    * @param pe This is the processor to which the task is mapped
    */
    std::pair<int,int> is_gain_fast(int task, int new_time, int old_time, int pe);

    /*!
    * This function checks whether task can be scheduled to pot_time
    * @param task This is the task that should be scheduled
    * @param pot_time This is the time unit to which task should be moved
    */
    bool is_move_possible(int task, int pot_time); 

    /*!
    * This function schedules task to new_time
    * @param task This is the task that is moved
    * @param new_time This is the new time unit where task starts
    * @param old_time This is the time unit to which task was scheduled before
    * @param pe This is the processor to which the task is mapped
    */
    void apply_move(int task, int new_time, int old_time, int pe); 

    /*!
    * This function returns the final carbon cost after applying the local search
    */
    unsigned int getCost(); 

    /*
    *
    * ILP COMPARISON
    * 
    */

    /*!
    * This function computes the cost of the schedule that is computed by our ILP formulation of the scheduling problem
    * @param experiments If true, the carbon cost of the schedule is printed, otherwise the carbon cost is only returned 
    */
    int print_cost_of_ILP_schedule_zeroComm(bool experiments) {
      
      unsigned int total_brown_power = 0;
      int base_power = 0; 
      for (size_t l = 0; l < numProcessors; l++) {
        base_power += PStatic[l];
      }
        
      for (size_t l = 0; l < numSubintervals; l++) {
        std::vector<int> power_per_time_unit(subintervalLengths[l], base_power); 
        std::vector<unsigned int> covered_time_units(numProcessors, 0);
        for (size_t k = 0; k < numProcessors; k++) {
          for (size_t i = 0; i < mapping_ILP[k].size(); i++) {
            int task = mapping_ILP[k][i]; 
            if (startingTimes_ILP[task] >= splitPoints[l+1]) {continue;}
            if (!(startingTimes_ILP[task] + nodeWeights_ILP[task] - 1 < splitPoints[l])){
              int start = std::max(splitPoints[l], startingTimes_ILP[task]);
              int end = std::min(splitPoints[l+1], startingTimes_ILP[task] + nodeWeights_ILP[task]) - 1;
              for (size_t c = start; c <= end; c++) {
                power_per_time_unit[c-splitPoints[l]] += PDynamic[k];
              }
            }
          }
        }

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

    /*!
    * This function reads the schedule computed by the ILP which is given csv format
    * @param filename This is the csv file where the schedule is stored
    */
    void read_ILP_schedule(const std::string& filename) {
      startingTimes_ILP.resize(G.getNumNodes(), -1);
      nodeWeights_ILP.resize(G.getNumNodes(), -1);
      mapping_ILP.resize(numProcessors);
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
      }
    
      std::string line;
      // Skip the header line
      std::getline(file, line);
    
      while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;
    
        int task, processor, start, end, duration;
    
        // Read task
        std::getline(ss, token, ',');
        task = std::stoi(token);
    
        // Read processor
        std::getline(ss, token, ',');
        processor = std::stoi(token);
    
        // Read start
        std::getline(ss, token, ',');
        start = std::stoi(token);
    
        // Read end
        std::getline(ss, token, ',');
        end = std::stoi(token);
    
        // Read duration
        std::getline(ss, token, ',');
        duration = std::stoi(token);
    
        // Store the values in the appropriate data structures
        mapping_ILP[processor].push_back(task);
        startingTimes_ILP[task] = start;
        nodeWeights_ILP[task] = duration;
      }
    
      file.close();
    }

    /*
    *
    * MISCELLANEOUS
    *
    */
 
    int getNumProcessors();
    
    int getDeadline();
    
    int getNumberOfSubIntervals();
    
    /*!
    * This function returns the time in milliseconds for the initial communication setup before the scheduling begins
    */
    double getSetupTime();
    
    int getBCC();

    /*!
    * This function prints the Power values for each processor
    */
    void print_power() {
      std::cout << "Static Power: ";
      for (int i = 0; i < PStatic.size(); i++) {
        std::cout << PStatic[i] << " ";
      }
      std::cout << std::endl;
            
      std::cout << "Dynamic Power: ";
      for (int i = 0; i < PDynamic.size(); i++) {
        std::cout << PDynamic[i] << " ";
      }
      std::cout << std::endl;
    }

  private: 

    /********** Input data  **********/
    unsigned int deadline; 
    unsigned int numSubintervals; 
    unsigned int numProcessors;
    std::vector<unsigned int> PStatic; 
    std::vector<unsigned int> PDynamic; 
    unsigned int BCC; 
    std::vector<std::vector<int>> mapping; 
    DAG G; // input DAG 
    std::vector<int> subintervalBudgets; 
    std::vector<int> splitPoints; 
    std::vector<int> subintervalLengths; 
    std::vector<interval> intervals; 
    std::vector<int> processors_of_tasks; // processor_of_tasks[i] = j means task i is assigned to processor j
    unsigned int brownPower; 
    
    /********** Data for communication setup  **********/
    DAG G_comm; // for communication mode
    unsigned int number_of_computation_tasks;
    std::vector<vertex> vertices; 
    std::vector<int> processors_of_tasks_comm;
    double setupTimeCommunicationSetting;

    /********** Algorithm data  **********/
    std::vector<std::pair<int, int>> EST_LST; // earliest start time and latest start time for every task
    std::vector<int> topological_order; 
    std::vector<float> slack; 
    std::vector<float> pressure;
    std::vector<int> score_ordered_tasks;
    std::vector<int> startingTimes; 
    std::vector<bool> already_scheduled; 
    unsigned int basePower; 
    std::vector<std::vector<unsigned int>> powerPerTimeUnit; 

    /********** Data for ILP comparison **********/ 
    std::vector<std::vector<int>> mapping_ILP;
    std::vector<int> startingTimes_ILP;
    std::vector<int> nodeWeights_ILP;
};
