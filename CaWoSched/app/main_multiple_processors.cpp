#include <iostream>
#include <lib/multipleProcessors/multiple_processors.h>
#include <chrono>
#include <string_view>
#include <filesystem>


struct result {
    std::string variant;
    std::string workflow;
    int PEs;
    int deadline;
    int subIntervals;
    int seed;
    int BCC;
    int cost;
    double time;
    bool validity;
};

bool has_flag(int argc, char** argv, const std::string_view flag) {
  for (int i = 0; i < argc; ++i) {
    if (argv[i] == flag) return true;
  }
  return false;
}

void slackR(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "slackR";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;
  
  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_slack_machine_focused(false);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(true);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
      Result.variant = "slackR-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

void slack(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();
  
  Result.variant = "slack";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_slack_machine_focused(false);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(false);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "slack-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

void slackWR(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "slackWR";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_slack_machine_focused(true);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(true);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "slackWR-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

void slackW(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "slackW";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_slack_machine_focused(true);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(false);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "slackW-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}



void pressR(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "pressR";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_pressure_machine_focused(false);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(true);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "pressR-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

void press(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "press";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_pressure_machine_focused(false);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(false);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "press-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

void pressWR(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "pressWR";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_pressure_machine_focused(true);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(true);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "pressWR-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

void pressW(bool LS, result& Result, const std::filesystem::path& wf_file, const std::string& mapping_file, const std::string& setup_file) {
  MultiMachineSchedulingWithoutFreezing scheduler(wf_file.c_str(), mapping_file, setup_file, 1);
  int PEs = scheduler.getNumProcessors();
  int deadl = scheduler.getDeadline();
  int subInt = scheduler.getNumberOfSubIntervals();
  int seed = 1;
  int BCC = scheduler.getBCC();

  Result.variant = "pressW";
  Result.workflow = wf_file.filename();
  Result.PEs = PEs;
  Result.deadline = deadl;
  Result.subIntervals = subInt;
  Result.seed = seed;
  Result.BCC = BCC;

  auto start_time = std::chrono::high_resolution_clock::now();
  scheduler.compute_pressure_machine_focused(true);
  scheduler.greedy_scheduling_with_budget_maintenance_enhanced(false);
  if (LS) {
    scheduler.LS_machine_focused_fastUpdates(10);
    Result.variant = "pressW-LS";
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
  double algo_time = duration.count() + scheduler.getSetupTime();
  int carbon_cost = -1;
  bool valid = scheduler.is_schedule_valid(false);
  if (valid) {
    if (LS) {
      carbon_cost = scheduler.getCost();
    } else {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }
  }
  Result.cost = carbon_cost;
  Result.time = algo_time;
  Result.validity = valid;
}

int main (int argc, char** argv) {    
  if (argc < 4 || argc > 6) {
    std::cout << "Usage: " << argv[0] << " <DAG file> <Mapping file> <Setup file> [--baseline-only] [--no_LS]" << std::endl;
    return 1;
  }

  std::filesystem::path const graphPath(std::move(argv[1])); 
  std::string const wf = graphPath.filename(); 

  std::string const mapping_file = std::move(argv[2]);
  std::string const setup_file = std::move(argv[3]);

  bool baseline_only = has_flag(argc, argv, "--baseline_only");
  bool no_LS = has_flag(argc, argv, "--no_LS");

  result result_slack;
  result result_slackR;
  result result_slackW;
  result result_slackWR;

  result result_press;
  result result_pressR;
  result result_pressW;
  result result_pressWR;

  if (baseline_only) {
    result result_baseline;
    MultiMachineSchedulingWithoutFreezing scheduler(graphPath.c_str(), mapping_file, setup_file, 1);
    
    int PEs = scheduler.getNumProcessors();
    int deadl = scheduler.getDeadline();
    int subInt = scheduler.getNumberOfSubIntervals();
    int seed = 1;
    int BCC = scheduler.getBCC();
    
    result_baseline.variant = "baseline";
    result_baseline.workflow = wf;
    result_baseline.PEs = PEs;
    result_baseline.deadline = deadl;
    result_baseline.subIntervals = subInt;
    result_baseline.seed = seed;
    result_baseline.BCC = BCC;

    auto start_time = std::chrono::high_resolution_clock::now();
    scheduler.compute_EST_schedule();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    double algo_time = duration.count() + scheduler.getSetupTime();
    int carbon_cost = -1;
    bool valid = scheduler.is_schedule_valid(false);
    if (valid) {
      carbon_cost = scheduler.get_cost_of_schedule(true);
    }

    result_baseline.cost = carbon_cost;
    result_baseline.time = algo_time;
    result_baseline.validity = valid;

    // write header
    std::cout << "Variant,Workflow,PEs,Deadline,Subintervals,BCC,Seed,Cost,Time[ms],Validity\n";
    // write results
    std::cout << result_baseline.variant << "," << result_baseline.workflow << "," << result_baseline.PEs << "," << result_baseline.deadline << "," << result_baseline.subIntervals << "," << result_baseline.BCC << "," << result_baseline.seed << "," << result_baseline.cost << "," << result_baseline.time << "," << result_baseline.validity << "\n";
    
    return 0;
  }

  if (no_LS) {
    slack(false, result_slack, graphPath, mapping_file, setup_file);
    slackR(false, result_slackR, graphPath, mapping_file, setup_file);
    slackW(false, result_slackW, graphPath, mapping_file, setup_file);
    slackWR(false, result_slackWR, graphPath, mapping_file, setup_file);
    press(false, result_press, graphPath, mapping_file, setup_file);
    pressR(false, result_pressR, graphPath, mapping_file, setup_file);
    pressW(false, result_pressW, graphPath, mapping_file, setup_file);
    pressWR(false, result_pressWR, graphPath, mapping_file, setup_file);
  } else {
    slack(true, result_slack, graphPath, mapping_file, setup_file);
    slackR(true, result_slackR, graphPath, mapping_file, setup_file);
    slackW(true, result_slackW, graphPath, mapping_file, setup_file);
    slackWR(true, result_slackWR, graphPath, mapping_file, setup_file);
    press(true, result_press, graphPath, mapping_file, setup_file);
    pressR(true, result_pressR, graphPath, mapping_file, setup_file);
    pressW(true, result_pressW, graphPath, mapping_file, setup_file);
    pressWR(true, result_pressWR, graphPath, mapping_file, setup_file);
  }

  /****************** WRITE OUTPUT *******************************/
  // Structure: Variant, Workflow, PEs, Deadline, Subintervals, BCC, Seed, Cost Variant, Time Variant[ms], Validity Variant
  // write header
  std::cout << "Variant,Workflow,PEs,Deadline,Subintervals,BCC,Seed,Cost,Time[ms],Validity\n";

  // write results
  std::cout << result_slackR.variant << "," << result_slackR.workflow << "," << result_slackR.PEs << "," << result_slackR.deadline << "," << result_slackR.subIntervals << "," << result_slackR.BCC << "," << result_slackR.seed << "," << result_slackR.cost << "," << result_slackR.time << "," << result_slackR.validity << "\n";
  std::cout << result_slack.variant << "," << result_slack.workflow << "," << result_slack.PEs << "," << result_slack.deadline << "," << result_slack.subIntervals << "," << result_slack.BCC << "," << result_slack.seed << "," << result_slack.cost << "," << result_slack.time << "," << result_slack.validity << "\n";
  std::cout << result_slackWR.variant << "," << result_slackWR.workflow << "," << result_slackWR.PEs << "," << result_slackWR.deadline << "," << result_slackWR.subIntervals << "," << result_slackWR.BCC << "," << result_slackWR.seed << "," << result_slackWR.cost << "," << result_slackWR.time << "," << result_slackWR.validity << "\n";
  std::cout << result_slackW.variant << "," << result_slackW.workflow << "," << result_slackW.PEs << "," << result_slackW.deadline << "," << result_slackW.subIntervals << "," << result_slackW.BCC << "," << result_slackW.seed << "," << result_slackW.cost << "," << result_slackW.time << "," << result_slackW.validity << "\n";

  std::cout << result_pressR.variant << "," << result_pressR.workflow << "," << result_pressR.PEs << "," << result_pressR.deadline << "," << result_pressR.subIntervals << "," << result_pressR.BCC << "," << result_pressR.seed << "," << result_pressR.cost << "," << result_pressR.time << "," << result_pressR.validity << "\n";
  std::cout << result_press.variant << "," << result_press.workflow << "," << result_press.PEs << "," << result_press.deadline << "," << result_press.subIntervals << "," << result_press.BCC << "," << result_press.seed << "," << result_press.cost << "," << result_press.time << "," << result_press.validity << "\n";
  std::cout << result_pressWR.variant << "," << result_pressWR.workflow << "," << result_pressWR.PEs << "," << result_pressWR.deadline << "," << result_pressWR.subIntervals << "," << result_pressWR.BCC << "," << result_pressWR.seed << "," << result_pressWR.cost << "," << result_pressWR.time << "," << result_pressWR.validity << "\n";
  std::cout << result_pressW.variant << "," << result_pressW.workflow << "," << result_pressW.PEs << "," << result_pressW.deadline << "," << result_pressW.subIntervals << "," << result_pressW.BCC << "," << result_pressW.seed << "," << result_pressW.cost << "," << result_pressW.time << "," << result_pressW.validity << "\n";

  return 0;
}
