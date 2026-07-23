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

struct EdgeCUT {
    long long int u, v;
};

// =====================================================================
// 1. ORIGINAL SCORE EVALUATORS
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
        if (satisfied) score += weights[idx];
    }
    return score;
}

long long int get_maxcut_score(const vector<bool>& assignment, const vector<EdgeCUT>& edges) {
    long long int score = 0;
    for (const auto& edge : edges) {
        // Cut if they are in different sets (true vs false)
        if (assignment[edge.u - 1] != assignment[edge.v - 1]) {
            score++;
        }
    }
    return score;
}

// =====================================================================
// 2. YOUR HEURISTIC
// =====================================================================
struct Edge{
    long long int to;
    long long int weight;
};
pair<vector<bool>, long long int> solve_heuristic_cpp(int N, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights) {
    // 1. Use uint8_t instead of vector<bool> for much faster memory access
    vector<uint8_t> assignment(N, 0);

    // 2. Build an Adjacency List to avoid O(M) loop lookups
    vector<vector<Edge>> adj(N);
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    // 3a. Single-flip Delta Function (Zero heap allocations, using long long int to prevent overflow)
    auto apply_flip_single = [&](int v) -> long long int {
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
    auto apply_flips = [&](const vector<int>& vars) -> long long int {
        long long int delta = 0;
        for (int v : vars) {
            delta += apply_flip_single(v);
        }
        return delta;
    };

    // Initialize base score
    long long int current_score = 0;
    for(long long int w : neg_weights) current_score += w;

    // Pre-allocate vectors outside the tight loops to prevent millions of allocations
    vector<int> group;         group.reserve(100);
    vector<int> rev_group;     rev_group.reserve(100);
    vector<int> to_flip;       to_flip.reserve(100);
    vector<int> drop_flips;    drop_flips.reserve(100);

    // =========================================================
    // PHASE 1: Constructive Group Lookahead (Active Queue)
    // =========================================================
    vector<int> queue_vars;
    queue_vars.reserve(N * 2); // Avoid reallocation during pushes
    vector<bool> in_queue(N, true);
    for(int i = 0; i < N; ++i) queue_vars.push_back(i);

    int head = 0;
    while(head < queue_vars.size()) {
        int v = queue_vars[head++];
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
                for (int u : group) {
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
    int max_passes = 15; // Prevent long-tail stalling on massive graphs
    int pass = 0;
    
    while (changed && pass < max_passes) { 
        changed = false;
        pass++;

        // 1. Standard 1-opt Pruning
        for (int v = 0; v < N; ++v) {
            if (assignment[v]) {
                long long int gain = apply_flip_single(v);
                if (gain > 0) {
                    current_score += gain;
                    changed = true;
                } else {
                    apply_flip_single(v); // Revert
                }
            }
        }
        
        // 2. Cascading Swap (Reversal with smart repair)
        for (int v = 0; v < N; ++v) {
            if (assignment[v]) {
                to_flip.clear();
                to_flip.push_back(v);

                for (const auto& edge : adj[v]) {
                    int other = edge.to;
                    if (!assignment[other] && edge.weight > neg_weights[other]) {
                        to_flip.push_back(other);
                    }
                }

                long long int gain = apply_flips(to_flip);
                if (gain > 0) {
                    current_score += gain;
                    changed = true;
                } else {
                    // Fast revert without std::reverse
                    for (int i = to_flip.size() - 1; i >= 0; --i) {
                        apply_flip_single(to_flip[i]); 
                    }
                }
            }
        }
        
        // 3. Proactive Insertion (Neighborhood-Restricted Cascade)
        for (int v = 0; v < N; ++v) {
            if (!assignment[v]) {
                long long int initial_gain = apply_flip_single(v); 
                long long int accumulated_sub_gain = 0;
                
                // Track all committed flips in this cascade for fast rollback
                vector<int> committed_sub_flips; 
                committed_sub_flips.reserve(20);
                
                bool sub_changed = true;
                while (sub_changed) {
                    sub_changed = false;
                    
                    // FAST SCAN: Only check neighbors of v, not the whole graph (O(Degree) instead of O(N))
                    for (const auto& neighbor_edge : adj[v]) {
                        int u = neighbor_edge.to;
                        if (assignment[u]) {
                            drop_flips.clear();
                            drop_flips.push_back(u);
                            
                            // Check for necessary repairs
                            for (const auto& edge_u : adj[u]) {
                                int other = edge_u.to;
                                if (!assignment[other] && edge_u.weight > neg_weights[other]) {
                                    drop_flips.push_back(other);
                                }
                            }

                            long long int sub_gain = apply_flips(drop_flips);
                            
                            if (sub_gain > 0) { 
                                accumulated_sub_gain += sub_gain;
                                for(int flip_var : drop_flips) committed_sub_flips.push_back(flip_var);
                                sub_changed = true;
                            } else {
                                // Fast revert
                                for (int i = drop_flips.size() - 1; i >= 0; --i) {
                                    apply_flip_single(drop_flips[i]); 
                                }
                            }
                        }
                    }
                }

                // Final Evaluation
                if (initial_gain + accumulated_sub_gain > 0) {
                    current_score += (initial_gain + accumulated_sub_gain);
                    changed = true;
                } else {
                    // Rollback the entire cascade in reverse order
                    for (int i = committed_sub_flips.size() - 1; i >= 0; --i) {
                        apply_flip_single(committed_sub_flips[i]);
                    }
                    apply_flip_single(v); // Undo the initial insertion
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
// 3. GENERATORS & L-REDUCTIONS
// =====================================================================

void generate_random_maxcut(long long int num_vars, long long int num_edges, vector<EdgeCUT>& edges) {
    mt19937 rng(random_device{}());
    uniform_int_distribution<long long int> var_dist(1, num_vars);

    edges.clear();
    for (long long int i = 0; i < num_edges; ++i) {
        EdgeCUT edge;
        edge.u = var_dist(rng);
        edge.v = var_dist(rng);
        while (edge.u == edge.v) edge.v = var_dist(rng);
        edges.push_back(edge);
    }
}

void convert_maxcut_to_max2sat(const vector<EdgeCUT>& edges, vector<Clause>& clauses, vector<long long int>& weights) {
    clauses.clear();
    weights.clear();

    for (const auto& edge : edges) {
        // (u V v)
        clauses.push_back({{edge.u, edge.v}});
        weights.push_back(1);
        
        // (!u V !v)
        clauses.push_back({{-edge.u, -edge.v}});
        weights.push_back(1);
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
// 4. CPLEX EXACT SOLVERS
// =====================================================================

pair<vector<bool>, long long int> solve_maxcut_cplex(long long int num_vars, const vector<EdgeCUT>& edges) {
    IloEnv env;
    vector<bool> assignment(num_vars, false);
    long long int opt_score = 0;

    try {
        IloModel model(env);
        long long int M = edges.size();
        
        IloBoolVarArray x(env, num_vars); 
        IloBoolVarArray y(env, M); // y[j] = 1 if edge j is cut

        IloExpr objExpr(env);
        for (long long int j = 0; j < M; ++j) objExpr += y[j];
        model.add(IloMaximize(env, objExpr));
        objExpr.end();

        // Edge is cut if x[u] != x[v]
        for (long long int j = 0; j < M; ++j) {
            long long int u = edges[j].u - 1;
            long long int v = edges[j].v - 1;
            
            // Linearization of XOR
            model.add(y[j] <= x[u] + x[v]);
            model.add(y[j] <= 2 - x[u] - x[v]);
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); 
        cplex.setWarning(env.getNullStream());

        if (cplex.solve()) {
            opt_score = round(cplex.getObjValue());
            for (long long int i = 0; i < num_vars; ++i) {
                assignment[i] = (round(cplex.getValue(x[i])) == 1.0);
            }
        }
    } catch (IloException& e) { cerr << "Concert exception: " << e << endl; }
    
    env.end();
    return {assignment, opt_score};
}

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

void test_maxcut_pipeline(long long int num_tests) {
    cout << "\n--- MAX CUT EXACT vs MAX 2-SAT GADGET HEURISTIC ---\n";
    cout << left << setw(13) << "Vars / Edges" << " | "
         << setw(10) << "Opt CUT" << " | "
         << setw(10) << "Heur CUT" << " | "
         << setw(10) << "Opt 2SAT" << " | "
         << setw(10) << "Heur 2SAT" << " | "
         << "Approximability Ratios (CUT vs 2SAT)\n";
    cout << string(95, '-') << "\n";

    for (long long int t = 0; t < num_tests; ++t) {
        long long int V = 200;
        long long int E = 400; // Dense enough to be challenging

        vector<EdgeCUT> edges;
        generate_random_maxcut(V, E, edges);

        // 1. Exact CPLEX solve for MAX CUT
        auto cplex_cut = solve_maxcut_cplex(V, edges);
        long long int opt_cut = cplex_cut.second;

        // 2. Convert to MAX 2-SAT 
        vector<Clause> sat_clauses;
        vector<long long int> sat_weights;
        convert_maxcut_to_max2sat(edges, sat_clauses, sat_weights);
        long long int total_sat_vars = V; // No auxiliary variables needed!

        // 3. Exact CPLEX solve for MAX 2-SAT
        auto cplex_2sat = solve_classical_max2sat_cplex(total_sat_vars, sat_clauses, sat_weights);
        long long int opt_2sat = cplex_2sat.second;

        // 4. Heuristic solve via your custom mapping
        long long int N_red;
        vector<OrClause> or_clauses;
        vector<long long int> neg_weights;
        reduce_to_weighted_restricted(total_sat_vars, sat_clauses, sat_weights, N_red, or_clauses, neg_weights);

        auto heur_res = solve_heuristic_cpp(N_red, or_clauses, neg_weights);
        
        vector<bool> heur_sat_assign(total_sat_vars);
        for(long long int i = 0; i < total_sat_vars; ++i) heur_sat_assign[i] = heur_res.first[i];
        long long int heur_2sat = get_original_score(heur_sat_assign, sat_clauses, sat_weights);

        // 5. Evaluate the Heuristic's assignment back on MAX CUT
        long long int heur_cut = get_maxcut_score(heur_sat_assign, edges);

        // Format Output
        string graph_size = "V:" + to_string(V) + " E:" + to_string(E);
        double ratio_cut = (opt_cut > 0) ? (double)heur_cut / opt_cut : 1.0;
        double ratio_2sat = (opt_2sat > 0) ? (double)heur_2sat / opt_2sat : 1.0;

        cout << left << setw(13) << graph_size << " | "
             << setw(10) << opt_cut << " | "
             << setw(10) << heur_cut << " | "
             << setw(10) << opt_2sat << " | "
             << setw(10) << heur_2sat << " | "
             << fixed << setprecision(3) << "CUT: " << ratio_cut << "  ->  2SAT: " << ratio_2sat << "\n";
    }
}

int main() {
    test_maxcut_pipeline(10);
    return 0;
}