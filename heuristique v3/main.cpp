#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <map>

// CPLEX Header
#include <ilcplex/ilocplex.h>

using namespace std;

// =====================================================================
// DATA STRUCTURES
// =====================================================================

struct Clause {
    vector<int> literals;
};

struct OrClause {
    int u, v, w;
};

// =====================================================================
// 1. ORIGINAL SCORE EVALUATOR
// =====================================================================

int get_original_score(const vector<bool>& assignment, const vector<Clause>& clauses, const vector<int>& weights) {
    int score = 0;
    for (size_t idx = 0; idx < clauses.size(); ++idx) {
        bool satisfied = false;
        for (int lit : clauses[idx].literals) {
            int var_idx = abs(lit) - 1;
            bool is_true = (lit > 0) ? assignment[var_idx] : !assignment[var_idx];
            if (is_true) {
                satisfied = true;
                break;
            }
        }
        if (satisfied) {
            score += weights[idx];
        }
    }
    return score;
}

// =====================================================================
// 2. HEURISTIC FUNCTIONS
// =====================================================================

int get_score(const vector<bool>& assignment, const vector<OrClause>& or_clauses, const vector<int>& neg_weights) {
    int score = 0;
    for (const auto& clause : or_clauses) {
        if (assignment[clause.u] || assignment[clause.v]) {
            score += clause.w;
        }
    }
    for (size_t i = 0; i < neg_weights.size(); ++i) {
        if (!assignment[i]) {
            score += neg_weights[i];
        }
    }
    return score;
}
struct Edge {
    int to;
    int weight;
};
pair<vector<bool>, int> solve_heuristic_cpp(int N, const vector<OrClause>& or_clauses, const vector<int>& neg_weights) {
    // 1. Use uint8_t instead of vector<bool> for much faster memory access
    vector<uint8_t> assignment(N, 0);

    // 2. Build an Adjacency List to avoid O(M) loop lookups
    vector<vector<Edge>> adj(N);
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    // 3. Delta Function: Flips a sequence of variables, updates assignment in-place,
    // and returns the exact change in score in O(degree) time.
    auto apply_flips = [&](const vector<int>& vars) -> int {
        int delta = 0;
        for (int v : vars) {
            if (assignment[v]) { // Transitioning True -> False
                delta += neg_weights[v];
                for (const auto& edge : adj[v]) {
                    if (!assignment[edge.to]) delta -= edge.weight;
                }
            } else { // Transitioning False -> True
                delta -= neg_weights[v];
                for (const auto& edge : adj[v]) {
                    if (!assignment[edge.to]) delta += edge.weight;
                }
            }
            assignment[v] ^= 1; // Toggle state in-place
        }
        return delta;
    };

    // Initialize base score (all zeros)
    int current_score = 0;
    for(int w : neg_weights) current_score += w;

    // =========================================================
    // PHASE 1: Constructive Group Lookahead
    // =========================================================
    while (true) {
        int best_gain = 0;
        vector<int> best_group;

        for (int v = 0; v < N; ++v) {
            if (assignment[v]) continue;

            // Try flipping just 'v'
            int gain1 = apply_flips({v});
            if (gain1 > best_gain) {
                best_gain = gain1;
                best_group = {v};
            }
            apply_flips({v}); // Revert immediately

            // Try flipping 'v' and all of its currently false neighbors
            vector<int> group = {v};
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) {
                    group.push_back(edge.to);
                }
            }

            if (group.size() > 1) {
                int gain2 = apply_flips(group);
                if (gain2 > best_gain) {
                    best_gain = gain2;
                    best_group = group;
                }
                // Revert in reverse order to perfectly undo the delta calculation
                vector<int> rev_group = group;
                reverse(rev_group.begin(), rev_group.end());
                apply_flips(rev_group); 
            }
        }

        if (best_gain > 0) {
            current_score += apply_flips(best_group); // Commit the best move
        } else {
            break;
        }
    }

    // =========================================================
    // PHASES 2, 3, & 4: Refinement, Reversals, and Proactive Insertions
    // =========================================================
    bool changed = true;
    for(int i=0;i<10;i++) { // Limit to 10 iterations to avoid infinite loops}
        changed = false;

        // 1. Standard 1-opt Pruning
        for (int v = 0; v < N; ++v) {
            int gain = apply_flips({v});
            if (gain > 0) {
                current_score += gain;
                changed = true;
            } else {
                apply_flips({v}); // Revert
            }
        }
        if (changed) continue;

        // 2. Cascading Swap (Reversal with smart repair)
        for (int v = 0; v < N; ++v) {
            if (assignment[v]) {
                vector<int> to_flip = {v};

                for (const auto& edge : adj[v]) {
                    int other = edge.to;
                    if (!assignment[other]) {
                        if (edge.weight > neg_weights[other]) {
                            to_flip.push_back(other);
                        }
                    }
                }

                int gain = apply_flips(to_flip);
                if (gain > 0) {
                    current_score += gain;
                    changed = true;
                } else {
                    reverse(to_flip.begin(), to_flip.end());
                    apply_flips(to_flip); // Revert
                }
            }
        }
        if (changed) continue;

        // 3. Proactive Insertion
        for (int v = 0; v < N; ++v) {
            if (!assignment[v]) {
                int initial_gain = apply_flips({v}); // Tentatively commit v=1
                
                // Copy state for the sub-search to avoid complicated backtracking logic
                vector<uint8_t> test_assign = assignment;
                int test_score = current_score + initial_gain;

                bool sub_changed = true;
                while (sub_changed) {
                    sub_changed = false;

                    for (int u = 0; u < N; ++u) {
                        if (test_assign[u] && u != v) {
                            vector<int> drop_flips = {u};
                            
                            for (const auto& edge : adj[u]) {
                                int other = edge.to;
                                if (!test_assign[other] && edge.weight > neg_weights[other]) {
                                    drop_flips.push_back(other);
                                }
                            }

                            // Inline fast delta for the local test array
                            int delta = 0;
                            for (int idx : drop_flips) {
                                if (test_assign[idx]) {
                                    delta += neg_weights[idx];
                                    for (const auto& e : adj[idx]) if (!test_assign[e.to]) delta -= e.weight;
                                } else {
                                    delta -= neg_weights[idx];
                                    for (const auto& e : adj[idx]) if (!test_assign[e.to]) delta += e.weight;
                                }
                                test_assign[idx] ^= 1; 
                            }

                            if (delta > 0) {
                                test_score += delta;
                                sub_changed = true;
                            } else {
                                // Revert local test flips
                                for (int i = drop_flips.size() - 1; i >= 0; --i) {
                                    test_assign[drop_flips[i]] ^= 1;
                                }
                            }
                        }
                    }
                }

                if (test_score > current_score) {
                    assignment = test_assign;
                    current_score = test_score;
                    changed = true;
                } else {
                    apply_flips({v}); // Undo initial tentative v=1
                }
            }
        }
    }

    // Convert fast uint8_t array back to vector<bool> for the return signature
    vector<bool> final_assign(N);
    for (int i = 0; i < N; ++i) final_assign[i] = assignment[i];

    return {final_assign, current_score};
}
// =====================================================================
// 3. GENERATOR & L-REDUCTION
// =====================================================================

void generate_standard_max2sat(int num_vars, int num_clauses, int max_weight, vector<Clause>& clauses, vector<int>& weights) {
    mt19937 rng(random_device{}());
    uniform_int_distribution<int> size_dist(1, 2);
    uniform_int_distribution<int> var_dist(1, num_vars);
    uniform_int_distribution<int> sign_dist(0, 1);
    uniform_int_distribution<int> weight_dist(1, max_weight);

    clauses.clear();
    weights.clear();

    for (int i = 0; i < num_clauses; ++i) {
        int size = size_dist(rng);
        Clause c;
        int v1 = var_dist(rng);
        c.literals.push_back(sign_dist(rng) ? v1 : -v1);
        
        if (size == 2) {
            int v2 = var_dist(rng);
            while (v2 == v1) v2 = var_dist(rng); // Ensure distinct variables
            c.literals.push_back(sign_dist(rng) ? v2 : -v2);
        }
        clauses.push_back(c);
        weights.push_back(weight_dist(rng));
    }
}

void reduce_to_weighted_restricted(int num_vars, const vector<Clause>& original_clauses, const vector<int>& weights, 
                                   int& N_new, vector<OrClause>& or_clauses_arr, vector<int>& neg_weights) {
    int V = num_vars;
    N_new = 2 * V;
    neg_weights.assign(N_new, 0);
    map<pair<int, int>, int> or_clauses_dict;

    auto get_var_idx = [V](int literal) {
        int var_idx = abs(literal) - 1;
        return (literal > 0) ? var_idx : var_idx + V;
    };

    // 1. Translate
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        int w = weights[idx];
        if (original_clauses[idx].literals.size() == 2) {
            int u = get_var_idx(original_clauses[idx].literals[0]);
            int v = get_var_idx(original_clauses[idx].literals[1]);
            if (u > v) swap(u, v);
            or_clauses_dict[{u, v}] += w;
        } else if (original_clauses[idx].literals.size() == 1) {
            int literal = original_clauses[idx].literals[0];
            if (literal > 0) neg_weights[abs(literal) - 1 + V] += w;
            else neg_weights[abs(literal) - 1] += w;
        }
    }

    // 2. Weighted variable degrees
    vector<int> degrees(V + 1, 0);
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        int w = weights[idx];
        for (int literal : original_clauses[idx].literals) {
            degrees[abs(literal)] += w;
        }
    }

    // 3. Consistency Gadgets
    for (int i = 1; i <= V; ++i) {
        int d_i = degrees[i];
        if (d_i == 0) continue;
        int u = i - 1;
        int v = i - 1 + V;
        or_clauses_dict[{u, v}] += (2 * d_i);
        neg_weights[u] += d_i;
        neg_weights[v] += d_i;
    }

    or_clauses_arr.clear();
    for (const auto& pair_kv : or_clauses_dict) {
        or_clauses_arr.push_back({pair_kv.first.first, pair_kv.first.second, pair_kv.second});
    }
}

// =====================================================================
// 4. CPLEX EXACT SOLVER
// =====================================================================

pair<vector<bool>, int> solve_classical_max2sat_cplex(int num_vars, const vector<Clause>& clauses, const vector<int>& weights) {
    IloEnv env;
    vector<bool> assignment(num_vars, false);
    int opt_score = 0;

    try {
        IloModel model(env);
        int M = clauses.size();
        
        // Variables
        IloBoolVarArray x(env, num_vars); 
        IloBoolVarArray z(env, M);        

        // Objective: Maximize sum(weight_j * z_j)
        IloExpr objExpr(env);
        for (int j = 0; j < M; ++j) {
            objExpr += weights[j] * z[j];
        }
        model.add(IloMaximize(env, objExpr));
        objExpr.end();

        // Constraints: z_j <= sum(literals)
        for (int j = 0; j < M; ++j) {
            IloExpr clauseExpr(env);
            for (int literal : clauses[j].literals) {
                int var_idx = abs(literal) - 1;
                if (literal > 0) {
                    clauseExpr += x[var_idx];
                } else {
                    clauseExpr += (1 - x[var_idx]);
                }
            }
            model.add(z[j] <= clauseExpr);
            clauseExpr.end();
        }

        // Solve
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); // Suppress CPLEX output
        cplex.setWarning(env.getNullStream());

        if (cplex.solve()) {
            opt_score = round(cplex.getObjValue());
            for (int i = 0; i < num_vars; ++i) {
                assignment[i] = (round(cplex.getValue(x[i])) == 1.0);
            }
        } else {
            cerr << "CPLEX Failed to find a solution!" << endl;
        }

    } catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    } catch (...) {
        cerr << "Unknown exception caught" << endl;
    }

    env.end();
    return {assignment, opt_score};
}

// =====================================================================
// 5. UNIFIED TEST RUNNER
// =====================================================================

void test_unified_pipeline(int num_tests) {
    cout << "\n--- UNIFIED CPLEX EXACT vs C++ HEURISTIC ---\n";
    cout << left << setw(15) << "Graph Size" << " | "
         << setw(10) << "CPLEX Exact" << " | "
         << setw(10) << "Heur Score" << " | "
         << setw(10) << "Time CPLEX" << " | "
         << "Time Heur\n";
    cout << string(65, '-') << "\n";

    for (int t = 0; t < num_tests; ++t) {
        int V = 2000;
        int M = 10000;


        vector<Clause> orig_clauses;
        vector<int> orig_weights;
        generate_standard_max2sat(V, M, 1, orig_clauses, orig_weights);

        // 2. Solve exactly via CPLEX
        auto t0 = chrono::high_resolution_clock::now();
        auto cplex_res = solve_classical_max2sat_cplex(V, orig_clauses, orig_weights);
        auto t1 = chrono::high_resolution_clock::now();

        // 3. L-Reduce and Solve via Heuristic
        int N_red;
        vector<OrClause> or_clauses;
        vector<int> neg_weights;
        reduce_to_weighted_restricted(V, orig_clauses, orig_weights, N_red, or_clauses, neg_weights);

        auto heur_res = solve_heuristic_cpp(N_red, or_clauses, neg_weights);
        auto t2 = chrono::high_resolution_clock::now();

        // 4. Evaluate Heuristic's assignment back on Original Graph
        vector<bool> heur_orig_assign(V);
        for(int i = 0; i < V; ++i) heur_orig_assign[i] = heur_res.first[i];
        
        int heur_true_score = get_original_score(heur_orig_assign, orig_clauses, orig_weights);

        // 5. Format Output
        double cplex_time = chrono::duration<double, milli>(t1 - t0).count();
        double heur_time = chrono::duration<double, milli>(t2 - t1).count();

        string graph_size = "V:" + to_string(V) + " E:" + to_string(M);
        
        cout << left << setw(15) << graph_size << " | "
             << setw(10) << cplex_res.second << " | "
             << setw(10) << heur_true_score << " | "
             << fixed << setprecision(1) << cplex_time << "ms   | "
             << heur_time << "ms\n";
    }
}

int main() {
    test_unified_pipeline(15);
    return 0;
}