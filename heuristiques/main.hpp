#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <filesystem>
#include <string>
#include <chrono>
#include <functional>

namespace fs = std::filesystem;

using namespace std;

// --- Data Structures ---

struct Clause {
    int weight;
    vector<int> literals;
};

// Function type aliases
using ReductionFunc = std::function<vector<Clause>(const vector<Clause>&)>;
using SolverFunc = std::function<map<int, bool>(const vector<Clause>&)>;


// --- CPLEX Solvers ---
map<int, bool> solve_max2sat_subgradient(const vector<Clause>& wcnf_array, int max_iters = 1, double initial_step = 1.0);

// --- Generators and I/O ---
vector<Clause> generate_max_2sat(int num_vars, int num_clauses, int max_weight = 10, bool allow_unit_clauses = false);
void write_wcnf(const vector<Clause>& clauses, const string& file_path);

// --- Utilities ---
int get_max_var(const vector<Clause>& clauses);
pair<int, int> evaluate_max2sat_solution(const vector<Clause>& formula, const map<int, bool>& solution);

// --- Heuristics and Reduction ---
vector<Clause> reduce_max2sat_to_implications_v1_1(const vector<Clause>& clauses);
vector<Clause> reduce_max2sat_to_implications_v1_2(const vector<Clause>& clauses);

// --- CPLEX Solvers ---
map<int, bool> solve_max2sat_implications_lp(const vector<Clause>& wcnf_array);

vector<Clause> read_wcnf(const string& file_path, int& out_total_weight);
void process_folder_to_csv(const string& folder_path, const string& csv_path, double timeout_seconds, ReductionFunc reduction_func, SolverFunc solver_func);

#endif // MAIN_HPP