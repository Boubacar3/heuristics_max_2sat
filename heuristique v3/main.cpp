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
    vector<long long int> literals;
};

struct OrClause {
    long long int u, v, w;
};

// =====================================================================
// 1. ORIGINAL SCORE EVALUATOR
// =====================================================================

long long int get_original_score(const vector<bool>& assignment, const vector<Clause>& clauses, const vector<long long int>& weights) {
    long long int score = 0;
    for (size_t idx = 0; idx < clauses.size(); ++idx) {
        bool satisfied = false;
        for (long long int lit : clauses[idx].literals) {
            long long int var_idx = abs(lit) - 1;
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

long long int get_score(const vector<bool>& assignment, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights) {
    long long int score = 0;
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
    long long int to;
    long long int weight;
};
pair<vector<bool>, long long int> solve_heuristic_cpp(long long int N, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights) {
    // 1. Use uint8_t instead of vector<bool> for much faster memory access
    vector<uint8_t> assignment(N, 0);

    // 2. Build an Adjacency List to avoid O(M) loop lookups
    vector<vector<Edge>> adj(N);
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    // 3a. Single-flip Delta Function (Zero heap allocations)
    auto apply_flip_single = [&](long long int v) -> long long int {
        long long int delta = 0;
        if (assignment[v]) { // True -> False
            delta += neg_weights[v];
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) delta -= edge.weight;
            }
        } else { // False -> True
            delta -= neg_weights[v];
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) delta += edge.weight;
            }
        }
        assignment[v] ^= 1; 
        return delta;
    };

    // 3b. Multi-flip Delta Function
    auto apply_flips = [&](const vector<long long int>& vars) -> long long int {
        long long int delta = 0;
        for (long long int v : vars) {
            delta += apply_flip_single(v);
        }
        return delta;
    };

    // Initialize base score
    long long int current_score = 0;
    for(long long int w : neg_weights) current_score += w;

    // Pre-allocate vectors outside the tight loops to prevent millions of allocations
    vector<long long int> group;         group.reserve(100);
    vector<long long int> rev_group;     rev_group.reserve(100);
    vector<long long int> to_flip;       to_flip.reserve(100);
    vector<long long int> drop_flips;    drop_flips.reserve(100);

    // =========================================================
    // PHASE 1: Constructive Group Lookahead (Active Queue)
    // =========================================================
    vector<long long int> queue_vars;
    queue_vars.reserve(N * 2); // Avoid reallocation during pushes
    vector<bool> in_queue(N, true);
    for(long long int i = 0; i < N; ++i) queue_vars.push_back(i);

    long long int head = 0;
    while(head < queue_vars.size()) {
        long long int v = queue_vars[head++];
        in_queue[v] = false;

        if (assignment[v]) continue; // Only check False variables initially

        // Try 1-opt
        long long int gain1 = apply_flip_single(v);
        if (gain1 > 0) {
            current_score += gain1;
            for (const auto& edge : adj[v]) {
                if (!in_queue[edge.to]) {
                    in_queue[edge.to] = true;
                    queue_vars.push_back(edge.to);
                }
            }
            continue; 
        }
        apply_flip_single(v); // Revert

        // Try group flip
        group.clear();
        group.push_back(v);
        for (const auto& edge : adj[v]) {
            if (!assignment[edge.to]) {
                group.push_back(edge.to);
            }
        }

        if (group.size() > 1) {
            long long int gain2 = apply_flips(group);
            if (gain2 > 0) {
                current_score += gain2;
                for (long long int u : group) {
                    for (const auto& edge : adj[u]) {
                        if (!in_queue[edge.to]) {
                            in_queue[edge.to] = true;
                            queue_vars.push_back(edge.to);
                        }
                    }
                }
            } else {
                // Revert
                rev_group = group;
                reverse(rev_group.begin(), rev_group.end());
                apply_flips(rev_group); 
            }
        }
    }

    // =========================================================
    // PHASES 2, 3, & 4: Refinement, Reversals, and Proactive Insertions
    // =========================================================
    bool changed = true;
    while (changed) {
        changed = false;

        // 1. Standard 1-opt Pruning
        for (long long int v = 0; v < N; ++v) {
            long long int gain = apply_flip_single(v);
            if (gain > 0) {
                current_score += gain;
                changed = true;
            } else {
                apply_flip_single(v); // Revert
            }
        }
        if (changed) continue;

        // 2. Cascading Swap (Reversal with smart repair)
        for (long long int v = 0; v < N; ++v) {
            if (assignment[v]) {
                to_flip.clear();
                to_flip.push_back(v);

                for (const auto& edge : adj[v]) {
                    long long int other = edge.to;
                    if (!assignment[other] && edge.weight > neg_weights[other]) {
                        to_flip.push_back(other);
                    }
                }

                long long int gain = apply_flips(to_flip);
                if (gain > 0) {
                    current_score += gain;
                    changed = true;
                } else {
                    rev_group = to_flip;
                    reverse(rev_group.begin(), rev_group.end());
                    apply_flips(rev_group); // Revert
                }
            }
        }
        if (changed) continue;

        // 3. Proactive Insertion (Restricted Neighborhood Search)
        for (long long int v = 0; v < N; ++v) {
            if (!assignment[v]) {
                long long int initial_gain = apply_flip_single(v); // Tentatively commit v=1
                long long int best_sub_gain = 0;
                vector<long long int> best_sub_flips;
                
                // Only scan NEIGHBORS, not all N variables
                for (const auto& edge : adj[v]) {
                    long long int u = edge.to;
                    if (assignment[u] && u != v) {
                        drop_flips.clear();
                        drop_flips.push_back(u);
                        
                        for (const auto& edge_u : adj[u]) {
                            long long int other = edge_u.to;
                            if (!assignment[other] && edge_u.weight > neg_weights[other]) {
                                drop_flips.push_back(other);
                            }
                        }

                        long long int sub_gain = apply_flips(drop_flips);
                        if (sub_gain > best_sub_gain) {
                            best_sub_gain = sub_gain;
                            best_sub_flips = drop_flips;
                        }
                        
                        // Undo immediately to fairly test the next neighbor
                        rev_group = drop_flips;
                        reverse(rev_group.begin(), rev_group.end());
                        apply_flips(rev_group);
                    }
                }

                if (initial_gain + best_sub_gain > 0) {
                    current_score += (initial_gain + best_sub_gain);
                    if (!best_sub_flips.empty()) {
                        apply_flips(best_sub_flips); // Commit the chosen repair
                    }
                    changed = true;
                } else {
                    apply_flip_single(v); // Undo initial tentative v=1
                }
            }
        }
    }

    // Convert fast uint8_t array back to vector<bool> for the return signature
    vector<bool> final_assign(N);
    for (long long int i = 0; i < N; ++i) final_assign[i] = assignment[i];

    return {final_assign, current_score};
}
// =====================================================================
// 3. GENERATOR & L-REDUCTION
// =====================================================================

void generate_standard_max2sat(long long int num_vars, long long int num_clauses, long long int max_weight, vector<Clause>& clauses, vector<long long int>& weights) {
    mt19937 rng(random_device{}());
    uniform_int_distribution<long long int> size_dist(1, 2);
    uniform_int_distribution<long long int> var_dist(1, num_vars);
    uniform_int_distribution<long long int> sign_dist(0, 1);
    uniform_int_distribution<long long int> weight_dist(1, max_weight);

    clauses.clear();
    weights.clear();

    for (long long int i = 0; i < num_clauses; ++i) {
        long long int size = size_dist(rng);
        Clause c;
        long long int v1 = var_dist(rng);
        c.literals.push_back(sign_dist(rng) ? v1 : -v1);
        
        if (size == 2) {
            long long int v2 = var_dist(rng);
            while (v2 == v1) v2 = var_dist(rng); // Ensure distinct variables
            c.literals.push_back(sign_dist(rng) ? v2 : -v2);
        }
        clauses.push_back(c);
        weights.push_back(weight_dist(rng));
    }
}

void reduce_to_weighted_restricted(long long int num_vars, const vector<Clause>& original_clauses, const vector<long long int>& weights, 
                                   long long int& N_new, vector<OrClause>& or_clauses_arr, vector<long long int>& neg_weights) {
    long long int V = num_vars;
    N_new = 2 * V;
    neg_weights.assign(N_new, 0);
    map<pair<long long int, long long int>, long long int> or_clauses_dict;

    auto get_var_idx = [V](long long int literal) {
        long long int var_idx = abs(literal) - 1;
        return (literal > 0) ? var_idx : var_idx + V;
    };

    // 1. Translate
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        long long int w = weights[idx];
        if (original_clauses[idx].literals.size() == 2) {
            long long int u = get_var_idx(original_clauses[idx].literals[0]);
            long long int v = get_var_idx(original_clauses[idx].literals[1]);
            if (u > v) swap(u, v);
            or_clauses_dict[{u, v}] += w;
        } else if (original_clauses[idx].literals.size() == 1) {
            long long int literal = original_clauses[idx].literals[0];
            if (literal > 0) neg_weights[abs(literal) - 1 + V] += w;
            else neg_weights[abs(literal) - 1] += w;
        }
    }

    // 2. Weighted variable degrees
    vector<long long int> degrees(V + 1, 0);
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        long long int w = weights[idx];
        for (long long int literal : original_clauses[idx].literals) {
            degrees[abs(literal)] += w;
        }
    }

    // 3. Consistency Gadgets
    for (long long int i = 1; i <= V; ++i) {
        long long int d_i = degrees[i];
        if (d_i == 0) continue;
        long long int u = i - 1;
        long long int v = i - 1 + V;
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

pair<vector<bool>, long long int> solve_classical_max2sat_cplex(long long int num_vars, const vector<Clause>& clauses, const vector<long long int>& weights) {
    IloEnv env;
    vector<bool> assignment(num_vars, false);
    long long int opt_score = 0;

    try {
        IloModel model(env);
        long long int M = clauses.size();
        
        // Variables
        IloBoolVarArray x(env, num_vars); 
        IloBoolVarArray z(env, M);        

        // Objective: Maximize sum(weight_j * z_j)
        IloExpr objExpr(env);
        for (long long int j = 0; j < M; ++j) {
            objExpr += IloNum(weights[j]) * z[j];
        }
        model.add(IloMaximize(env, objExpr));
        objExpr.end();

        // Constraints: z_j <= sum(literals)
        for (long long int j = 0; j < M; ++j) {
            IloExpr clauseExpr(env);
            for (long long int literal : clauses[j].literals) {
                long long int var_idx = abs(literal) - 1;
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
            for (long long int i = 0; i < num_vars; ++i) {
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

void test_unified_pipeline(long long int num_tests) {
    cout << "\n--- UNIFIED CPLEX EXACT vs C++ HEURISTIC ---\n";
    cout << left << setw(15) << "Graph Size" << " | "
         << setw(10) << "CPLEX Exact" << " | "
         << setw(10) << "Heur Score" << " | "
         << setw(10) << "Time CPLEX" << " | "
         << "Time Heur\n";
    cout << string(65, '-') << "\n";

    for (long long int t = 0; t < num_tests; ++t) {
        long long int V = 2000;
        long long int M = 10000;
        long long int w = 100;

        vector<Clause> orig_clauses;
        vector<long long int> orig_weights;
        generate_standard_max2sat(V, M, w, orig_clauses, orig_weights);

        // 2. Solve exactly via CPLEX
        auto t0 = chrono::high_resolution_clock::now();
        auto cplex_res = solve_classical_max2sat_cplex(V, orig_clauses, orig_weights);
        auto t1 = chrono::high_resolution_clock::now();

        // 3. L-Reduce and Solve via Heuristic
        long long int N_red;
        vector<OrClause> or_clauses;
        vector<long long int> neg_weights;
        reduce_to_weighted_restricted(V, orig_clauses, orig_weights, N_red, or_clauses, neg_weights);

        auto heur_res = solve_heuristic_cpp(N_red, or_clauses, neg_weights);
        auto t2 = chrono::high_resolution_clock::now();

        // 4. Evaluate Heuristic's assignment back on Original Graph
        vector<bool> heur_orig_assign(V);
        for(long long int i = 0; i < V; ++i) heur_orig_assign[i] = heur_res.first[i];
        
        long long int heur_true_score = get_original_score(heur_orig_assign, orig_clauses, orig_weights);

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