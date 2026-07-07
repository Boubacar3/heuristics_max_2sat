#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <ilcplex/ilocplex.h>
#include <filesystem>
#include <string>
#include <chrono>
long long int total_weights_shift = 0;
long long int total_weight = 0;
namespace fs = std::filesystem;

using namespace std;

// --- Data Structures ---

struct Clause {
    long long int weight;
    vector<long long int> literals;
};

// --- Generators and I/O ---

vector<Clause> generate_max_2sat(long long int num_vars, long long int num_clauses, long long int max_weight = 10, bool allow_unit_clauses = false) {
    vector<Clause> clauses;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> weight_dist(1.0, max_weight);
    uniform_int_distribution<> bool_dist(0, 1);
    uniform_int_distribution<> var_dist(1, num_vars);

    for (long long int i = 0; i < num_clauses; ++i) {
        long long int weight = std::round(weight_dist(gen));
        bool is_unit = allow_unit_clauses && bool_dist(gen);

        if (is_unit) {
            long long int var1 = var_dist(gen);
            long long int sign1 = bool_dist(gen) ? 1 : -1;
            clauses.push_back({weight, {sign1 * var1}});
        } else {
            long long int var1 = var_dist(gen);
            long long int var2 = var_dist(gen);
            while (var1 == var2) {
                var2 = var_dist(gen); // Ensure distinct variables
            }
            long long int sign1 = bool_dist(gen) ? 1 : -1;
            long long int sign2 = bool_dist(gen) ? 1 : -1;
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
        for (long long int lit : clause.literals) {
            f << lit << " ";
        }
        f << "0\n";
    }
    f.close();
}

// --- Utilities ---

long long int get_max_var(const vector<Clause>& clauses) {
    long long int max_var = 0;
    for (const auto& clause : clauses) {
        for (long long int lit : clause.literals) {
            max_var = max(max_var, abs(lit));
        }
    }
    return max_var;
}
// Returns a pair: {satisfied_weight, violated_weight}
pair<long long int, long long int> evaluate_max2sat_solution(const vector<Clause>& formula, const map<long long int, bool>& solution) {
    long long int total_satisfied_weight = 0;
    long long int total_violated_weight = 0;

    for (const auto& item : formula) {
        long long int weight = item.weight;
        bool is_satisfied = false;

        for (long long int lit : item.literals) {
            long long int var = abs(lit);
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
    total_weights_shift = 0;
    total_weight = 0;
    for (const auto& clause : clauses) {
        long long int w = clause.weight;
        total_weight += w;
        if (clause.literals.size() == 1) {
            new_clauses.push_back(clause);
        } else if (clause.literals.size() == 2) {
            long long int l1 = clause.literals[0];
            long long int l2 = clause.literals[1];

            if (l1 > 0 && l2 > 0) {
                new_clauses.push_back({w, {l2}});
                new_clauses.push_back({-w, {-l1, l2}});
                total_weights_shift -= w;
            } else if (l1 < 0 && l2 < 0) {
                new_clauses.push_back({w, {l1}});
                new_clauses.push_back({-w, {l1, -l2}});
                total_weights_shift -= w;
            } else {
                new_clauses.push_back(clause);
            }
        }
    }
    return new_clauses;
}

// --- CPLEX Solvers ---

// 1. Solves the ILP with binary variable values and a timeout parameter
map<long long int, bool> solve_max2sat_implications_ilp(const vector<Clause>& wcnf_array, long long int timeout_seconds) {
    IloEnv env;
    map<long long int, bool> assignments;
    long long int max_var = get_max_var(wcnf_array);

    try {
        IloModel model(env);
        IloBoolVarArray y(env, max_var + 1);
        for (long long int i = 0; i <= max_var; ++i) {
            y[i] = IloBoolVar(env);
        }
        IloExpr objExpr(env);

        for (size_t k = 0; k < wcnf_array.size(); ++k) {
            long long int weight = wcnf_array[k].weight;
            IloNum cplex_weight = static_cast<IloNum>(weight);
            const vector<long long int>& lits = wcnf_array[k].literals;

            if (lits.size() == 1) {
                long long int lit = lits[0];
                long long int var_idx = abs(lit);
                if (lit > 0) {
                    objExpr += cplex_weight * y[var_idx];
                } else {
                    objExpr += cplex_weight * (1 - y[var_idx]);
                }
            } else if (lits.size() == 2) {
                long long int lit1 = lits[0], lit2 = lits[1];
                long long int premise = 0, conclusion = 0;

                if (lit1 < 0 && lit2 > 0) { premise = abs(lit1); conclusion = lit2; }
                else if (lit2 < 0 && lit1 > 0) { premise = abs(lit2); conclusion = lit1; }
                else {
                    cerr << "Error: Clause is not a valid implication." << endl;
                    continue;
                }

                // Continuous slack variable s_k
                IloNumVar s_k(env, 0, 1.0, ILOFLOAT);
                objExpr += cplex_weight * s_k;

                model.add(s_k + y[premise] - y[conclusion] <= 1.0);
                if (weight < 0) {
                    model.add(s_k >= y[conclusion]);
                    model.add(s_k >= 1 - y[premise]);
               }
            }
        }

        model.add(IloMaximize(env, objExpr));
        IloCplex cplex(model);
        // Force CPLEX to use a single core/thread
        cplex.setParam(IloCplex::Threads, 1);
        cplex.setParam(IloCplex::TiLim, timeout_seconds);
        cplex.setParam(IloCplex::EpAGap, 0);
        cplex.setParam(IloCplex::EpGap, 0);
        cplex.setParam(IloCplex::EpInt, 0);
        cplex.setParam(IloCplex::EpRHS, 1e-9);
        cplex.setOut(env.getNullStream()); // Suppress CPLEX logging

        if (cplex.solve()) {
            for (long long int i = 1; i <= max_var; ++i) {
                assignments[i] = (cplex.getValue(y[i]) > 0.5);
            }
        } else {
            cout << "ILP Solver Status: " << cplex.getStatus() << endl;
        }

    } catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    } catch (...) {
        cerr << "Unknown exception caught" << endl;
    }

    env.end();
    return assignments;
}

// 2. Exports the linear relaxation of the MAX-SAT problem to .lp format
long long int export_linear_relaxation_to_lp(const vector<Clause>& wcnf_array, const string& output_file_path) {
    IloEnv env;
    long long int max_var = get_max_var(wcnf_array);
    long long int result = std::numeric_limits<double>::quiet_NaN();
    long long int objective_value = std::numeric_limits<double>::quiet_NaN();

    try {
        IloModel model(env);
        
        // Create continuous variables [0, 1] for the relaxation
        IloNumVarArray y(env, max_var + 1, 0, 1, ILOFLOAT);
        IloExpr objExpr(env);
        
        // Vector to keep track of slack variables
        vector<IloNumVar> slack_vars;
        

        for (size_t k = 0; k < wcnf_array.size(); ++k) {
            long long int weight = wcnf_array[k].weight;
            IloNum cplex_weight = static_cast<IloNum>(weight);
            const vector<long long int>& lits = wcnf_array[k].literals;

            if (lits.size() == 1) {
                // Unit clause
                long long int lit = lits[0];
                long long int var_idx = abs(lit);
                if (lit > 0) {
                    objExpr += cplex_weight * y[var_idx];
                } else {
                    objExpr += cplex_weight * (1 - y[var_idx]);
                }
            } else if (lits.size() == 2) {
                // Implication clause
                long long int lit1 = lits[0], lit2 = lits[1];
                long long int premise = 0, conclusion = 0;

                if (lit1 < 0 && lit2 > 0) { 
                    premise = abs(lit1); 
                    conclusion = lit2; 
                }
                else if (lit2 < 0 && lit1 > 0) { 
                    premise = abs(lit2); 
                    conclusion = lit1; 
                }
                else {
                    cerr << "Error: Clause is not a valid implication." << endl;
                    continue;
                }

                // Continuous slack variable s_k
                IloNumVar s_k(env, 0, 1.0, ILOFLOAT);
                slack_vars.push_back(s_k);
                objExpr += cplex_weight * s_k;

                // Add constraint: s_k + y[premise] - y[conclusion] <= 1
                model.add(s_k + y[premise] - y[conclusion] <= 1.0);
                
                // Additional constraints for negative weights
                if (weight < 0) {
                    model.add(s_k >= y[conclusion]);
                    model.add(s_k >= 1 - y[premise]);
                }
            }
        }
        model.add(IloMaximize(env, objExpr));
        
        // Export the model to LP format and solve
        IloCplex cplex(model);
        // Force CPLEX to use a single core/thread for deterministic behavior
        cplex.setParam(IloCplex::Threads, 1);
        cplex.exportModel(output_file_path.c_str());
        
        // Solve the model
        if (cplex.solve()) {
            objective_value = cplex.getObjValue();
            result = total_weight - objective_value + total_weights_shift;
            cout << "Linear relaxation exported to: " << output_file_path << endl;
            cout << "Total Weight: " << total_weight << endl;
            cout << "Objective Value: " << objective_value << endl;
            cout << "Total Weight - Objective Value: " << total_weight - objective_value << endl;
            cout << "Total Weights Shift: " << total_weights_shift << endl;
            cout << "Computed Difference (Total Weight - Objective Value - Total Weights Shift): "<< result << endl;
        } else {
            cout << "Linear relaxation exported to: " << output_file_path << endl;
            cout << "Warning: Could not solve the model to compute the difference." << endl;
        }

    } catch (IloException& e) {
        cerr << "Concert exception caught during LP export: " << e << endl;
    } catch (...) {
        cerr << "Unknown exception caught during LP export" << endl;
    }

    env.end();
    return result;
}

vector<Clause> read_wcnf(const string& file_path, long long int& out_total_weight) {
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
            long long int vars, num_clauses;
            long long int sum_weight = 0;
            // The total weight is optional in some formats, so we handle it gracefully
            if (iss >> format >> vars >> num_clauses) {
                if (iss >> sum_weight) {
                    out_total_weight = sum_weight;
                }
            }
            continue; 
        }

        // Since there are no hard clauses, we attempt to read the weight directly as a long long int
        long long int weight = 0;
        try {
            weight = stod(token);
        } catch (const std::invalid_argument& e) {
            // If it's not a valid number, we skip this malformed line
            continue; 
        }

        vector<long long int> current_clause;
        while (iss >> token) {
            try {
                long long int lit = stoi(token);
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
}

void process_folder_to_csv(const string& folder_path, const string& csv_path, long long int timeout_seconds, bool export_lp = true) {
    ofstream csv_file(csv_path);
    
    csv_file << "Filename,Num_Variables,Num_Clauses,Total_Weight,Satisfied_Weight,Violated_Weight,Elapsed_Time,Linear_Relaxation_Objective\n";

    if (!fs::exists(folder_path) || !fs::is_directory(folder_path)) {
        cerr << "Error: Directory does not exist -> " << folder_path << endl;
        return;
    }

    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".wcnf") {
            string filepath = entry.path().string();
            string filename = entry.path().filename().string();
            
            cout << "Processing: " << filename << "..." << endl;

            long long int total_weight_from_header = 0;
            vector<Clause> original_clauses = read_wcnf(filepath, total_weight_from_header);
            if (original_clauses.empty()) continue;

            long long int num_vars = get_max_var(original_clauses);
            long long int num_clauses = original_clauses.size();

            vector<Clause> reduced_clauses = reduce_max2sat_to_implications(original_clauses);
            
            long long int linear_relaxation_objective = std::numeric_limits<long long int>::quiet_NaN();
            // Export linear relaxation if requested
            if (export_lp) {
                string lp_path = filepath;
                size_t pos = lp_path.find(".wcnf");
                if (pos != string::npos) {
                    lp_path.replace(pos, 5, ".lp");
                } else {
                    lp_path += ".lp";
                }
                linear_relaxation_objective = export_linear_relaxation_to_lp(reduced_clauses, lp_path);
            }
            
            // --- START TIMER ---
            auto start_time = chrono::high_resolution_clock::now();

            // Solve the ILP formulation with binary values
            map<long long int, bool> assignment = solve_max2sat_implications_ilp(reduced_clauses, timeout_seconds);

            // --- STOP TIMER ---
            auto end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed_seconds = end_time - start_time;
            long long int satisfied_weight = 0;
            long long int violated_weight = 0;

            if (!assignment.empty()) {
                // Get the exact weights directly from the assignment
                pair<long long int, long long int> weights = evaluate_max2sat_solution(original_clauses, assignment);
                satisfied_weight = weights.first;
                violated_weight = weights.second;
            }

            // We calculate the actual total weight dynamically as a sanity check
            long long int actual_total_weight = satisfied_weight + violated_weight;

            csv_file << filename << "," 
                     << num_vars << "," 
                     << num_clauses << "," 
                     << actual_total_weight << ","
                     << satisfied_weight << "," 
                     << violated_weight << ","
                     << elapsed_seconds.count() << ",";
            if (!std::isnan(linear_relaxation_objective)) {
                csv_file << linear_relaxation_objective;
            }
            csv_file << "\n";
            cout << "  -> Solved! Violated Weight: " << violated_weight << "\n";
        }
    }
    
    cout << "\nAll files processed. Results saved to " << csv_path << endl;
}
// --- Main Execution ---
int main(int argc, char* argv[]) {
    string input_folder = "";   // Default input folder
    string output_csv = "";      // Default output CSV name
    long long int timeout_seconds; // Default CPLEX time limit in seconds
    bool export_lp = true; // Whether to export the linear relaxation to .lp format
    if (argc > 1) {
        input_folder = argv[1];
    }
    else {
        cout << "No input folder specified." << endl;
        return 0;
    }
    if (argc > 2) {
        output_csv = argv[2];
    }
    else {
        cout << "No output CSV specified." << endl;
        return 0;
    }
    if (argc > 3) {
        timeout_seconds = stoi(argv[3]);
    }
    else {
        cout << "No timeout specified." << endl;
        return 0;
    }
    // Create the folder for testing purposes if it doesn't exist
    if (!fs::exists(input_folder)) {
        fs::create_directory(input_folder);
        cout << "Created directory: " << input_folder << ". Please put your .wcnf files there." << endl;
    } else {
        process_folder_to_csv(input_folder, output_csv, timeout_seconds, export_lp);
    }

    return 0;
}