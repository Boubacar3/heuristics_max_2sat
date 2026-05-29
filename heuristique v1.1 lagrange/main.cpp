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

namespace fs = std::filesystem;

using namespace std;

// --- Data Structures ---

struct Clause {
    int weight;
    vector<int> literals;
};

// --- Generators and I/O ---

vector<Clause> generate_max_2sat(int num_vars, int num_clauses, int max_weight = 10, bool allow_unit_clauses = false) {
    vector<Clause> clauses;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> weight_dist(1.0, max_weight);
    uniform_int_distribution<> bool_dist(0, 1);
    uniform_int_distribution<> var_dist(1, num_vars);

    for (int i = 0; i < num_clauses; ++i) {
        int weight = std::round(weight_dist(gen));
        bool is_unit = allow_unit_clauses && bool_dist(gen);

        if (is_unit) {
            int var1 = var_dist(gen);
            int sign1 = bool_dist(gen) ? 1 : -1;
            clauses.push_back({weight, {sign1 * var1}});
        } else {
            int var1 = var_dist(gen);
            int var2 = var_dist(gen);
            while (var1 == var2) {
                var2 = var_dist(gen); // Ensure distinct variables
            }
            int sign1 = bool_dist(gen) ? 1 : -1;
            int sign2 = bool_dist(gen) ? 1 : -1;
            clauses.push_back({weight, {sign1 * var1, sign2 * var2}});
        }
    }
    return clauses;
}

void write_wcnf(const vector<Clause>& clauses, const string& file_path) {
    ofstream f(file_path);
    if (!f.is_open()) {
        cerr << "Error opening file for writing." << endl;
        return;
    }
    for (const auto& clause : clauses) {
        f << clause.weight << " ";
        for (int lit : clause.literals) {
            f << lit << " ";
        }
        f << "0\n";
    }
    f.close();
}

// --- Utilities ---

int get_max_var(const vector<Clause>& clauses) {
    int max_var = 0;
    for (const auto& clause : clauses) {
        for (int lit : clause.literals) {
            max_var = max(max_var, abs(lit));
        }
    }
    return max_var;
}
// Returns a pair: {satisfied_weight, violated_weight}
pair<int, int> evaluate_max2sat_solution(const vector<Clause>& formula, const map<int, bool>& solution) {
    int total_satisfied_weight = 0;
    int total_violated_weight = 0;

    for (const auto& item : formula) {
        int weight = item.weight;
        bool is_satisfied = false;

        for (int lit : item.literals) {
            int var = abs(lit);
            bool var_value = false;
            
            auto it = solution.find(var);
            if (it != solution.end()) {
                var_value = it->second;
            }

            if ((lit > 0 && var_value) || (lit < 0 && !var_value)) {
                is_satisfied = true;
                break;
            }
        }
        
        // Directly accumulate the exact weights
        if (is_satisfied) {
            total_satisfied_weight += weight;
        } else {
            total_violated_weight += weight; 
        }
    }
    return {total_satisfied_weight, total_violated_weight};
}
// --- Heuristics and Reduction ---

vector<Clause> reduce_max2sat_to_implications(const vector<Clause>& clauses) {
    vector<Clause> new_clauses;

    for (const auto& clause : clauses) {
        int w = clause.weight;

        if (clause.literals.size() == 1) {
            new_clauses.push_back(clause);
        } else if (clause.literals.size() == 2) {
            int l1 = clause.literals[0];
            int l2 = clause.literals[1];

            if (l1 > 0 && l2 > 0) {
                new_clauses.push_back({w, {l2}});
                new_clauses.push_back({-w, {-l1, l2}});
            } else if (l1 < 0 && l2 < 0) {
                new_clauses.push_back({w, {l1}});
                new_clauses.push_back({-w, {l1, -l2}});
            } else {
                new_clauses.push_back(clause);
            }
        }
    }
    return new_clauses;
}    

// --- CPLEX Solvers ---
/**
 * Solves the Lagrangian relaxation of MAX-2SAT by optimizing
 * the multipliers via a subgradient algorithm.
 * @param wcnf_array The clauses of the problem.
 * @param timeout_seconds Maximum total time for the subgradient loop in seconds.
 * @param initial_step Initial step for the update (learning rate).
 * @return The boolean assignment from the best relaxation.
 */
map<int, bool> solve_max2sat_subgradient(const vector<Clause>& wcnf_array, double timeout_seconds = 1.0, double initial_step = 1.0) {
    IloEnv env;
    map<int, bool> best_assignment;
    
    // Initialize random number generator for random rounding
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> prob_dist(0.0, 1.0);

    try {
        // 1. Deduce the total number of variables and identify the implications
        int num_vars = 0;
        vector<pair<int, int>> D_edges; // Stores the arcs i => j
        
        for (const auto& clause : wcnf_array) {
            for (int lit : clause.literals) {
                num_vars = max(num_vars, abs(lit));
            }
            if (clause.literals.size() == 2) {
                int lit1 = clause.literals[0];
                int lit2 = clause.literals[1];
                int i = (lit1 < 0) ? -lit1 : -lit2;
                int j = (lit1 > 0) ? lit1 : lit2;
                D_edges.push_back({i, j});
            }
        }

        // 2. Initialization of Lagrangian multipliers to zero
        map<pair<int, int>, double> lambda;
        map<pair<int, int>, double> mu;
        for (auto edge : D_edges) {
            lambda[edge] = 0.0;
            mu[edge] = 0.0;
        }

        // 3. Construction of the static architecture of the CPLEX model
        IloModel model(env);
        IloNumVarArray y(env, num_vars + 1, 0.0, 1.0, ILOFLOAT);
        map<pair<int, int>, IloNumVar> s;

        // Creation of s_ij variables and addition of the strict constraint (1)
        for (auto edge : D_edges) {
            if (s.find(edge) == s.end()) {
                s[edge] = IloNumVar(env, 0.0, 1.0, ILOFLOAT);
                // Constraint (1): sij + yi - yj <= 1
                model.add(s[edge] + y[edge.first] - y[edge.second] <= 1.0);
            }
        }

        // Objective (empty for now, will be updated dynamically)
        IloObjective objective = IloMaximize(env);
        model.add(objective);

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); // Silent mode

        double best_dual_bound = IloInfinity; // Objective to minimize the relaxation

        // =================================================================
        // 4. SUBGRADIENT LOOP
        // =================================================================
        auto loop_start = chrono::high_resolution_clock::now();
        for (int iter = 1; ; ++iter) {
            IloExpr obj_expr(env);

            // Construction of the objective function with current lambda and mu
            for (const auto& clause : wcnf_array) {
                if (clause.literals.size() == 2) {
                    int lit1 = clause.literals[0];
                    int lit2 = clause.literals[1];
                    int i = (lit1 < 0) ? -lit1 : -lit2;
                    int j = (lit1 > 0) ? lit1 : lit2;
                    pair<int, int> edge = {i, j};

                    double w_ij = clause.weight;
                    double l_ij = lambda[edge];
                    double m_ij = mu[edge];

                    obj_expr += (w_ij + l_ij + m_ij) * s[edge] - l_ij * y[j] + m_ij * y[i] - m_ij;

                } else if (clause.literals.size() == 1) {
                    int lit = clause.literals[0];
                    if (lit > 0) obj_expr += clause.weight * y[lit];
                    else obj_expr += clause.weight * (1.0 - y[-lit]);
                }
            }

            // Update of the objective in CPLEX
            objective.setExpr(obj_expr);
            obj_expr.end(); // Memory release

            // Resolution of the iteration
            if (!cplex.solve()) {
                cerr << "Error: CPLEX could not solve the iteration " << iter << endl;
                break;
            }

            double current_dual = cplex.getObjValue();

            // Recording of the best solution (the smallest upper bound)
            if (current_dual < best_dual_bound) {
                best_dual_bound = current_dual;
                best_assignment.clear();
                for (int v = 1; v <= num_vars; ++v) {
                    double value = cplex.getValue(y[v]);
                    // Random rounding: set to 1 with probability equal to the fractional value
                    best_assignment[v] = (prob_dist(gen) < value);
                }
            }

            // Calculation of the step (Step size): decreases with the square root of iterations
            double step = initial_step / sqrt(iter);

            // Update of the multipliers via the subgradients (violations)
            for (auto edge : D_edges) {
                int i = edge.first;
                int j = edge.second;

                double val_y_i = cplex.getValue(y[i]);
                double val_y_j = cplex.getValue(y[j]);
                double val_s_ij = cplex.getValue(s[edge]);

                // Calculation of violations
                double violation_2 = val_y_j - val_s_ij;             // For s_ij >= y_j
                double violation_3 = 1.0 - val_y_i - val_s_ij;       // For s_ij >= 1 - y_i

                // Updating of the multipliers (projected on >= 0)
                lambda[edge] = max(0.0, lambda[edge] + step * violation_2);
                mu[edge] = max(0.0, mu[edge] + step * violation_3);
            }

            auto loop_now = chrono::high_resolution_clock::now();
            double elapsed = chrono::duration<double>(loop_now - loop_start).count();
            if (elapsed >= timeout_seconds) {
                break;
            }
            
            // Optional line to track convergence:
            // cout << "Iter " << iter << " | Dual Bound: " << current_dual << endl;
        }

    } catch (IloException& e) {
        cerr << "Concert Exception: " << e << endl;
    } catch (...) {
        cerr << "Unknown error." << endl;
    }

    env.end(); // Properly releases the entire CPLEX environment
    return best_assignment;
}
vector<Clause> read_wcnf(const string& file_path, int& out_total_weight) {
    vector<Clause> clauses;
    out_total_weight = 0; 
    
    ifstream f(file_path);
    if (!f.is_open()) {
        cerr << "Error opening file: " << file_path << endl;
        return clauses;
    }

    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;

        istringstream iss(line);
        string token;
        
        if (!(iss >> token)) continue;

        if (token == "c") continue; 
        
        // Parse the header line: p wcnf <vars> <clauses> [<total_weight>]
        if (token == "p") {
            string format;
            int vars, num_clauses;
            int sum_weight = 0;
            // The total weight is optional in some formats, so we handle it gracefully
            if (iss >> format >> vars >> num_clauses) {
                if (iss >> sum_weight) {
                    out_total_weight = sum_weight;
                }
            }
            continue; 
        }

        // Since there are no hard clauses, we attempt to read the weight directly as a int
        int weight = 0;
        try {
            weight = stod(token);
        } catch (const std::invalid_argument& e) {
            // If it's not a valid number, we skip this malformed line
            continue; 
        }

        vector<int> current_clause;
        while (iss >> token) {
            try {
                int lit = stoi(token);
                if (lit == 0) {
                    clauses.push_back({weight, current_clause});
                    current_clause.clear();
                    break;
                } else {
                    current_clause.push_back(lit);
                }
            } catch (const std::invalid_argument& e) {}
        }
    }
    return clauses;
}void process_folder_to_csv(const string& folder_path, const string& csv_path, double timeout_seconds = 1.0) {
    ofstream csv_file(csv_path);
    
    csv_file << "Filename,Num_Variables,Num_Clauses,Total_Weight,Satisfied_Weight,Violated_Weight,Elapsed_Time\n";

    if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
        cerr << "Error: Directory does not exist -> " << folder_path << endl;
        return;
    }

    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".wcnf") {
            string filepath = entry.path().string();
            string filename = entry.path().filename().string();
            
            cout << "Processing: " << filename << "..." << endl;

            int total_weight_from_header = 0;
            vector<Clause> original_clauses = read_wcnf(filepath, total_weight_from_header);
            if (original_clauses.empty()) continue;

            int num_vars = get_max_var(original_clauses);
            int num_clauses = original_clauses.size();

            vector<Clause> reduced_clauses = reduce_max2sat_to_implications(original_clauses);
            // --- START TIMER ---
            auto start_time = chrono::high_resolution_clock::now();

            // Solve the Subgradient Relaxation
            map<int, bool> assignment = solve_max2sat_subgradient(reduced_clauses, timeout_seconds);

            // --- STOP TIMER ---
            auto end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_seconds = end_time - start_time;
            int satisfied_weight = 0;
            int violated_weight = 0;

            if (!assignment.empty()) {
                // Get the exact weights directly from the assignment
                pair<int, int> weights = evaluate_max2sat_solution(original_clauses, assignment);
                satisfied_weight = weights.first;
                violated_weight = weights.second;
            }

            // We calculate the actual total weight dynamically as a sanity check
            int actual_total_weight = satisfied_weight + violated_weight;

            csv_file << filename << "," 
                     << num_vars << "," 
                     << num_clauses << "," 
                     << actual_total_weight << ","
                     << satisfied_weight << "," 
                     << violated_weight << ","
                     << elapsed_seconds.count() << "\n";
            cout << "  -> Solved! Violated Weight: " << violated_weight << "\n";
        }
    }
    
    cout << "\nAll files processed. Results saved to " << csv_path << endl;
}
int main(int argc, char* argv[]) {
    string input_folder = argv[1];   // Change to your folder path
    string output_csv = argv[2];      // The desired output CSV name
    
    // Create the folder for testing purposes if it doesn't exist
    if (!fs::exists(input_folder)) {
        fs::create_directory(input_folder);
        cout << "Created directory: " << input_folder << ". Please put your .wcnf files there." << endl;
    } else {
        process_folder_to_csv(input_folder, output_csv, atoi(argv[3]));
    }

    return 0;
}