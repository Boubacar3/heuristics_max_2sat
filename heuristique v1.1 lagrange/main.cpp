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
 * Résout la relaxation Lagrangienne du MAX-2SAT en optimisant 
 * les multiplicateurs via un algorithme de sous-gradient.
 * * @param wcnf_array Les clauses du problème.
 * @param max_iters Nombre maximum d'itérations du sous-gradient.
 * @param initial_step Pas initial pour la mise à jour (learning rate).
 * @return L'assignement booléen issu de la meilleure relaxation.
 */
map<int, bool> solve_max2sat_subgradient(const vector<Clause>& wcnf_array, int max_iters = 1, double initial_step = 1.0) {
    IloEnv env;
    map<int, bool> best_assignment;

    try {
        // 1. Déduire le nombre total de variables et identifier les implications
        int num_vars = 0;
        vector<pair<int, int>> D_edges; // Stocke les arcs i => j
        
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

        // 2. Initialisation des multiplicateurs Lagrangiens à zéro
        map<pair<int, int>, double> lambda;
        map<pair<int, int>, double> mu;
        for (auto edge : D_edges) {
            lambda[edge] = 0.0;
            mu[edge] = 0.0;
        }

        // 3. Construction de l'architecture statique du modèle CPLEX
        IloModel model(env);
        IloNumVarArray y(env, num_vars + 1, 0.0, 1.0, ILOFLOAT);
        map<pair<int, int>, IloNumVar> s;

        // Création des variables s_ij et ajout de la contrainte (1) stricte
        for (auto edge : D_edges) {
            if (s.find(edge) == s.end()) {
                s[edge] = IloNumVar(env, 0.0, 1.0, ILOFLOAT);
                // Contrainte (1) : sij + yi - yj <= 1
                model.add(s[edge] + y[edge.first] - y[edge.second] <= 1.0);
            }
        }

        // Objectif (vide pour l'instant, sera mis à jour dynamiquement)
        IloObjective objective = IloMaximize(env);
        model.add(objective);

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); // Mode silencieux

        double best_dual_bound = IloInfinity; // Objectif de minimiser la relaxation

        // =================================================================
        // 4. BOUCLE DU SOUS-GRADIENT
        // =================================================================
        for (int iter = 1; iter <= max_iters; ++iter) {
            IloExpr obj_expr(env);

            // Construction de la fonction objectif avec lambda et mu actuels
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

            // Mise à jour de l'objectif dans CPLEX
            objective.setExpr(obj_expr);
            obj_expr.end(); // Libération mémoire

            // Résolution de l'itération
            if (!cplex.solve()) {
                cerr << "Erreur: CPLEX n'a pas pu résoudre l'itération " << iter << endl;
                break;
            }

            double current_dual = cplex.getObjValue();

            // Enregistrement de la meilleure solution (la plus petite borne supérieure)
            if (current_dual < best_dual_bound) {
                best_dual_bound = current_dual;
                best_assignment.clear();
                for (int v = 1; v <= num_vars; ++v) {
                    best_assignment[v] = (cplex.getValue(y[v]) >= 0.5);
                }
            }

            // Calcul du pas (Step size) : décroît avec la racine carrée des itérations
            double step = initial_step / sqrt(iter);

            // Mise à jour des multiplicateurs via les sous-gradients (violations)
            for (auto edge : D_edges) {
                int i = edge.first;
                int j = edge.second;

                double val_y_i = cplex.getValue(y[i]);
                double val_y_j = cplex.getValue(y[j]);
                double val_s_ij = cplex.getValue(s[edge]);

                // Calcul des violations
                double violation_2 = val_y_j - val_s_ij;             // Pour s_ij >= y_j
                double violation_3 = 1.0 - val_y_i - val_s_ij;       // Pour s_ij >= 1 - y_i

                // Actualisation des multiplicateurs (projetés sur >= 0)
                lambda[edge] = max(0.0, lambda[edge] + step * violation_2);
                mu[edge] = max(0.0, mu[edge] + step * violation_3);
            }
            
            // Ligne optionnelle pour tracker la convergence :
            // cout << "Iter " << iter << " | Dual Bound: " << current_dual << endl;
        }

    } catch (IloException& e) {
        cerr << "Exception Concert : " << e << endl;
    } catch (...) {
        cerr << "Erreur inconnue." << endl;
    }

    env.end(); // Libère proprement tout l'environnement CPLEX
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
}void process_folder_to_csv(const string& folder_path, const string& csv_path) {
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
            map<int, bool> assignment = solve_max2sat_subgradient(reduced_clauses);

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
// --- Main Execution ---
int main() {
    string input_folder = "./test_wcnf_files_20000_100000";   // Change to your folder path
    string output_csv = "results_random_v1.1.csv";      // The desired output CSV name
    
    // Create the folder for testing purposes if it doesn't exist
    if (!fs::exists(input_folder)) {
        fs::create_directory(input_folder);
        cout << "Created directory: " << input_folder << ". Please put your .wcnf files there." << endl;
    } else {
        process_folder_to_csv(input_folder, output_csv);
    }

    return 0;
}