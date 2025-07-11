#!/usr/bin/env python3

import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from functools import cmp_to_key
import argparse
import timeit
import re
import csv
import matplotlib.pyplot as plt

class SchedulingProblem: 
  def __init__(self, graph_file, mapping_file, setup_file): 

    # READ THE SETUP ################################################
    content = []
    with open(setup_file, 'r') as file: 
      for line in file: 
        line = line.strip()
        numbers = list(map(int, line.split()))
        content.append(numbers)

    match = re.search(r'_(.*?)_', setup_file)
    integer_match = re.findall(r'\d+', setup_file)
    self.seed = int(integer_match[-1])
    self.workflow = match.group(1)
    self.deadline = 0
    self.numSubintervals = 0
    self.PStatic = []
    self.PDynamic = []
    self.bcc = 0
    self.subinterval_budgets = []
    self.splitting_points = []

    counter = 0
    for numbers in content: 
      if (counter == 0):  
        self.deadline = int(numbers[0])
      elif (counter == 1):
        self.numSubintervals = int(numbers[0])
      elif (counter == 2):
        self.PStatic = numbers
      elif (counter == 3):
        self.PDynamic = numbers
      elif (counter == 4):
        self.bcc = numbers[0]
      elif (counter == 5):
        self.subinterval_budgets = numbers
      elif (counter == 6):
        self.splitting_points = numbers
      counter += 1
    ################################################################

    # READ THE MAPPING ################################################
    mapping_info = []
    with open(mapping_file, 'r') as file: 
      for line in file: 
        line = line.strip()
        numbers = list(map(int, line.split()))
        mapping_info.append(numbers)

    mapping_counter = 0

    self.numProcessors = 0
    self.processor_tasks = []
    for con in mapping_info: 
      if mapping_counter == 0: 
        self.numProcessors = int(con[0])
      else: 
        self.processor_tasks.append(con)
      mapping_counter += 1
    ################################################################

    # READ THE GRAPH ################################################
    # Read the DOT file
    self.graph = nx.drawing.nx_pydot.read_dot(graph_file)

    # Extract vertices and edges
    self.vertices = list(self.graph.nodes(data=True))
    self.edges = list(self.graph.edges(data=True))

    self.vertex_ids = []
    self.weights = {}
    for vertex in self.vertices:
      self.vertex_ids.append(int(vertex[0]))
      self.weights[int(vertex[0])] = int(vertex[1]['weight'])

    self.edge_list = []
    self.edge_weights = {}

    for edge in self.edges: 
      source = int(edge[0])
      target = int(edge[1])
      self.edge_list.append((source, target))
      self.edge_weights[(source, target)] = int(edge[2]['size'])
    ################################################################

  def solve_problem(self): 
    try:
      # CREATE THE MODEL 
      self.model = gp.Model("Single-Machine-No-Freezing")
      self.model.setParam("OutputFlag", 0)

      # setup of the communication graph
      # USE CONSTANTS FROM ABOVE ######################################################
      self.number_of_machines = self.numProcessors
      eps = 0.01
      M = 1000000000
      T = self.deadline

      # compute interval lengths from switching points
      self.interval_lengths = []
      for i in range(len(self.splitting_points)-1):
        self.interval_lengths.append(self.splitting_points[i+1] - self.splitting_points[i])

      number_of_intervals = self.numSubintervals
      self.ge = []
      budget_counter = 0
      for p in range(len(self.splitting_points)-1):
        start_of_interval = self.splitting_points[p]
        end_of_interval_excl = self.splitting_points[p+1]
        for i in range(start_of_interval, end_of_interval_excl):
          self.ge.append(self.subinterval_budgets[budget_counter])
        budget_counter += 1

      self.number_of_tasks = len(self.vertex_ids)

      P_static = self.PStatic
      P_dynamic = self.PDynamic

      # Define mapping and graph ####################################

      # construct mapping where we store pairs : (task, processor)
      self.mapping = {}
      proc_counter = 0
      for tasks in self.processor_tasks:
        #print(type(self.processor_tasks))
        for task in tasks:
          self.mapping[task] = proc_counter
        proc_counter += 1

      # this is just the mapping with the correct order
      mapping_order_counter = 0
      self.mapping_order = {}
      for tasks in self.processor_tasks:
        self.mapping_order[(mapping_order_counter)] = self.processor_tasks[mapping_order_counter]
        mapping_order_counter += 1

      new_task_id = self.number_of_tasks
      self.communication_tasks = {}
      for edge in self.edge_list: 
        source = edge[0]
        target = edge[1]
        if self.mapping[source] != self.mapping[target]:
          self.communication_tasks[(new_task_id)] = edge
          self.weights[new_task_id] = self.edge_weights[(source, target)]
          print(f"ID: {new_task_id}, Edge: {source} -> {target}")
          new_task_id += 1

      # new number of tasks
      self.number_of_tasks = self.number_of_tasks + len(self.communication_tasks)

      # introduce new processors for communication tasks
      # for each pair of processors (i,j) for i != j we introduce a new processor
      for task in self.communication_tasks:
        source = self.communication_tasks[task][0]
        target = self.communication_tasks[task][1]
        if self.mapping[source] == 0 : 
          self.mapping[task] = self.numProcessors + self.mapping[source] * (self.numProcessors - 1) + self.mapping[target]-1
        else: 
          if self.mapping[target] < self.mapping[source]: 
            self.mapping[task] = self.numProcessors + self.mapping[source] * (self.numProcessors - 1) + self.mapping[target]
          elif self.mapping[target] > self.mapping[source]: 
            self.mapping[task] = self.numProcessors + self.mapping[source] * (self.numProcessors - 1) + self.mapping[target]-1

      self.numProcessors = self.numProcessors + self.numProcessors * (self.numProcessors - 1)  

      # now construct the rest of the mapping order.
      for i in range(self.number_of_tasks): 
        if i < self.number_of_tasks - len(self.communication_tasks): 
          continue
        processor = self.mapping[i]
        if processor not in self.mapping_order:
          self.mapping_order[processor] = []
        self.mapping_order[processor].append(i)

        def compare_tasks(task1, task2):
          source_a, target_a = self.communication_tasks[task1]
          source_b, target_b = self.communication_tasks[task2]

          processor_tasks_a = self.processor_tasks[self.mapping[source_a]]
          processor_tasks_b = self.processor_tasks[self.mapping[target_a]] 

          if processor_tasks_a.index(source_a) < processor_tasks_a.index(source_b): 
            return -1
          elif processor_tasks_a.index(source_a) > processor_tasks_a.index(source_b): 
            return 1
          else: 
            if processor_tasks_b.index(target_a) < processor_tasks_b.index(target_b): 
              return -1
            else:
              return 1

          return 0

        # now I have to sort the mapping order for the communication tasks
        for proce, tasks in self.mapping_order.items(): 
          if proce >= self.number_of_machines:
            tasks.sort(key=cmp_to_key(compare_tasks))

        # CREATE THE VARIABLES ##########################################

        # green energy usage
        self.gu = self.model.addVars(self.deadline, vtype=GRB.INTEGER, name="green energy usage")

        # brown energy usage
        self.bu = self.model.addVars(self.deadline, vtype=GRB.INTEGER, name="brown energy usage")

        # beginning of tasks
        self.begin_comp = self.model.addVars(
          [(task, processor, t) for task, processor in self.mapping.items() for t in range(T)], vtype=GRB.BINARY, name="beginning of tasks")

        # end of tasks
        self.end_comp = self.model.addVars(
          [(task, processor, t) for task, processor in self.mapping.items() for t in range(T)], vtype=GRB.BINARY, name="end of tasks")

        # running tasks
        self.run_comp = self.model.addVars(
          [(task, processor, t) for task, processor in self.mapping.items() for t in range(T)], vtype=GRB.BINARY, name="running tasks")

        # total power usage
        self.gamma = self.model.addVars(self.deadline, vtype=GRB.INTEGER, name="total energy usage combined per interval")

        # auxiliary variables to check gamma[l] > ge[l]
        self.a = self.model.addVars(self.deadline, vtype=GRB.BINARY, name="auxiliary variables indicating whether we need brown energy")

        ###################################################################

        # SET THE OBJECTIVE ################################################
        self.model.setObjective(
          sum(self.bu[l] for l in range(self.deadline)),
          GRB.MINIMIZE
        )

        ####################################################################
        # CREATE THE CONSTRAINTS ##########################################

        # 'scheduling' constraints

        # constraints for computation
        for task_id, pe in self.mapping.items(): 
          self.model.addConstr(sum(self.begin_comp[task_id, pe, i] for i in range(T-self.weights[task_id]+1)) == 1)

        for task_id, pe in self.mapping.items(): 
          self.model.addConstr(sum(self.begin_comp[task_id, pe, i] for i in range(T-self.weights[task_id]+1, T)) == 0)

        for task_id, pe in self.mapping.items(): 
          self.model.addConstr(sum(self.end_comp[task_id, pe, i] for i in range(self.weights[task_id]-1)) == 0)

        for task_id, pe in self.mapping.items(): 
          self.model.addConstr(sum(self.end_comp[task_id, pe, i]  for i in range(self.weights[task_id]-1, T)) == 1)

        for task_id, pe in self.mapping.items(): 
          for t in range(T-self.weights[task_id] + 1):
            self.model.addConstr(self.begin_comp[task_id, pe, t] == self.end_comp[task_id, pe, t + self.weights[task_id] - 1])

        # THIS HAS TO BE CHECKED FOR CORRECTNESS
        # keep the internal order on each of the processors
        for pe, tasks in self.mapping_order.items(): 
          for i in range(len(tasks)-1): 
            for t in range(T):
              self.model.addConstr(self.begin_comp[tasks[i+1], pe, t] <= sum(self.end_comp[tasks[i], pe, l] for l in range(t)))

        for task_id, pe in self.mapping.items():
          self.model.addConstr(sum(self.run_comp[task_id, pe, l] for l in range(T)) == self.weights[task_id])

        for task_id, pe in self.mapping.items(): 
          for t in range(T-self.weights[task_id]+1):
            for k in range(t, t + self.weights[task_id]):
              self.model.addConstr(self.run_comp[task_id, pe, k] >= self.begin_comp[task_id, pe, t])

        for task_id, edge in self.communication_tasks.items():
          source = edge[0]
          target = edge[1]
          source_processor = self.mapping[source]
          target_processor = self.mapping[target]
          processor = self.mapping[task_id]
          for t in range(T): 
            self.model.addConstr(self.begin_comp[task_id, processor, t] <= sum(self.end_comp[source, source_processor, l] for l in range(t)))
            self.model.addConstr(self.begin_comp[target, target_processor, t] <= sum(self.end_comp[task_id, processor, l] for l in range(t)))


        # energy constraints'
        for l in range(self.deadline): #
          self.model.addConstr(self.gu[l] + self.bu[l] == self.gamma[l], f"green and brown energy usage should sum up to total energy usage during interval {l}")

        for l in range(self.deadline): 
          self.model.addConstr(self.gu[l] <= self.ge[l], f"Cannot use more green energy than available in interval {l}")

        for l in range(self.deadline): 
          self.model.addConstr(self.gu[l] >= 0)

        for l in range(self.deadline): 
          self.model.addConstr(self.gu[l] <= self.gamma[l])

        for l in range(self.deadline): 
          self.model.addConstr(self.gu[l] >= (self.gamma[l] - self.bu[l]))

        for l in range(self.deadline): 
          self.model.addConstr(self.bu[l] >= 0)

        for l in range(self.deadline): 
          self.model.addConstr(self.bu[l] >= self.gamma[l] - self.ge[l])

        for l in range(self.deadline): 
          self.model.addConstr(self.bu[l] <= self.gamma[l] - self.ge[l] + M*(1-self.a[l]))

        for l in range(self.deadline): 
          self.model.addConstr(self.bu[l] <= M*self.a[l])

        for l in range(self.deadline): 
          self.model.addConstr((self.gamma[l]-self.ge[l]) <= M*self.a[l])

        for l in range(self.deadline): 
          self.model.addConstr((self.gamma[l]-self.ge[l]) >= (eps - M*(1 - self.a[l])))

        for l in range(self.deadline): 
          self.model.addConstr(self.gamma[l] >= 0)

        for l in range(self.deadline): 
          self.model.addConstr(self.gamma[l] == sum(P_static[p] for p in range(self.numProcessors)) + sum((self.run_comp[task_id, pe, l]*P_dynamic[pe]) for task_id, pe in self.mapping.items()))
        ########################################################################################

        # OPTIMIZE MODEL ##########################################################
        self.model.setParam("Threads", 1)
        self.model.setParam("TimeLimit", 3600)
        self.model.optimize()
            
        carbon_intensity = 1  # e.g. 500 gCO₂ per unit of brown energy
        carbon_cost = sum(self.bu[l].X for l in range(self.deadline)) * carbon_intensity
        print(f"Total Carbon Cost: {carbon_cost:.2f} gCO₂")

        # Collect values
        gu_vals = [self.gu[l].X for l in range(self.deadline)]
        bu_vals = [self.bu[l].X for l in range(self.deadline)]
        ge_vals = [self.ge[l] for l in range(self.deadline)]  # assuming ge is a fixed list/array

        total_energy_used = [g + b for g, b in zip(gu_vals, bu_vals)]

        # Plot
        plt.figure(figsize=(12, 6))
        plt.plot(total_energy_used, label='Total Energy Used (green + brown)', linewidth=2)
        plt.plot(ge_vals, label='Green Energy Available (ge)', linestyle='--')
        #plt.plot(bu_vals, label='Brown Energy Used (bu)', linestyle=':', linewidth=2)

        plt.xlabel('Time Interval')
        plt.ylabel('Energy')
        plt.title('Energy Consumption vs Available Green Energy')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("test.png")
        print(bu_vals)

    except gp.GurobiError as e: 
      print(f"Error code {e.errno}: {e}")

    except AttributeError:
      print("Encountered an attribute error.")

  def check_schedule_and_solver_status(self): 
    # check the solver status
    print('')
    print("SOLVER STATUS #############################################################################")
    if self.model.Status == GRB.OPTIMAL:
      print("Optimal solution found.")
    elif self.model.Status == GRB.INFEASIBLE:
      print("The model is infeasible.")
    elif self.model.Status == GRB.UNBOUNDED:
      print("The model is unbounded.")
    elif self.model.Status == GRB.TIME_LIMIT:
      print("Time limit reached, solution found may not be optimal.")
    else:
      print(f"Solver ended with status: {self.model.Status}")
    print("############################################################################################")
    print('')

    #########################################################################
    # compute the schedule
    self.begin_time = {}
    self.end_time = {}
    for task_id, pe in self.mapping.items(): 
      for t in range(self.deadline): 
        if self.begin_comp[task_id, pe, t].X > 0.5: 
          self.begin_time[task_id] = t
        if self.end_comp[task_id, pe, t].X > 0.5: 
          self.end_time[task_id] = t
    #########################################################################

    #########################################################################
    # check validity of the schedule.
    # first check edges
    self.invalid = False 
    for vertex in self.vertices:
      neighbors = list(map(int, self.graph.adj[vertex[0]]))

      for neigh in neighbors: 
        if self.end_time[int(vertex[0])] >= self.begin_time[neigh]:
          self.invalid = True
    # then check communication tasks
    for task_id, edge in self.communication_tasks.items(): 
      source = edge[0]
      target = edge[1]
      if self.end_time[source] >= self.begin_time[task_id]:
        self.invalid = True
      if self.end_time[task_id] >= self.begin_time[target]:
        self.invalid = True

    # Now go through mappings and check the order 
    for pe, tasks in self.mapping_order.items(): 
      for i in range(len(tasks)-1): 
        if self.end_time[tasks[i]] >= self.begin_time[tasks[i+1]]:
          self.invalid = True

    if self.invalid:
      print("The schedule is invalid")
    else:  
      print("The schedule is valid")
      print(f"Carbon cost: {self.model.ObjVal:g}")
    #########################################################################

  def print_schedule(self):
    for task_id, pe in self.mapping.items(): 
      print(f"Task {task_id} is scheduled on processor {pe} from {self.begin_time[task_id]} to {self.end_time[task_id]}")

  def print_schedule_csv(self):
    data = [
      ["Task", "Processor", "Start", "End", "Duration"],
    ]

    for task_id, pe in self.mapping.items():
      data.append([task_id, pe, self.begin_time[task_id], self.end_time[task_id], self.end_time[task_id] - self.begin_time[task_id] + 1])

    output_file = f"ILP_{self.workflow}_{self.deadline}_{self.numProcessors}_{self.bcc}_{self.seed}_schedule.csv"

    with open(output_file, 'w', newline="") as file:
      writer = csv.writer(file)
      writer.writerows(data)


  def print_csv(self, time_in_ms): 
    import csv

    valid = not self.invalid

    data = [
      ["Variant","Workflow","PEs","Deadline","Subintervals","BCC","Seed","Cost","Time[ms]","Validity"],
      ["ILP",self.workflow,self.numProcessors,self.deadline,self.numSubintervals,self.bcc,self.seed,self.model.ObjVal,time_in_ms,valid]
    ]

    # write 
    output_file = f"ILP_{self.workflow}_{self.deadline}_{self.numProcessors}_{self.bcc}_{self.seed}.csv"
    with open(output_file, 'w', newline="") as file:
      writer = csv.writer(file)
      writer.writerows(data)
    
def getCost(self): 
  return self.model.ObjVal

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="ILP for scheduling problem with multiple machines and no freezing")
  parser.add_argument("--graph", type=str, help="Path to the graph file in .dot format", required=True)
  parser.add_argument("--mapping", type=str, help="Path to the mapping file", required=True)
  parser.add_argument("--setup", type=str, help="Path to the setup specs file", required=True)

  args = parser.parse_args()

  schedule = SchedulingProblem(args.graph, args.mapping, args.setup)
  execution_time = timeit.timeit(schedule.solve_problem, number=1) * 1000
  schedule.check_schedule_and_solver_status()
  finalCost = schedule.getCost()
  schedule.print_schedule()
  schedule.print_schedule_csv()
  print(f"{execution_time:.2f},{finalCost}")