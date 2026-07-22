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

struct Lin2Equation {
    long long int v1, v2, v3;
    int rhs; // 0 or 1
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
        if (satisfied) {
            score += weights[idx];
        }
    }
    return score;
}

long long int get_e3lin2_score(const vector<bool>& assignment, const vector<Lin2Equation>& equations) {
    long long int score = 0;
    for (const auto& eq : equations) {
        int v1_val = assignment[eq.v1 - 1] ? 1 : 0;
        int v2_val = assignment[eq.v2 - 1] ? 1 : 0;
        int v3_val = assignment[eq.v3 - 1] ? 1 : 0;
        if ((v1_val ^ v2_val ^ v3_val) == eq.rhs) {
            score++;
        }
    }
    return score;
}

// =====================================================================
// 2. HEURISTIC FUNCTIONS
// =====================================================================

struct Edge {
    long long int to;
    long long int weight;
};

pair<vector<bool>, long long int> solve_heuristic_cpp(long long int N, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights) {
    vector<uint8_t> assignment(N, 0);
    vector<vector<Edge>> adj(N);
    
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    auto apply_flip_single = [&](long long int v) -> long long int {
        long long int delta = 0;
        if (assignment[v]) { 
            delta += neg_weights[v];
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) delta -= edge.weight;
            }
        } else { 
            delta -= neg_weights[v];
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) delta += edge.weight;
            }
        }
        assignment[v] ^= 1; 
        return delta;
    };

    auto apply_flips = [&](const vector<long long int>& vars) -> long long int {
        long long int delta = 0;
        for (long long int v : vars) delta += apply_flip_single(v);
        return delta;
    };

    long long int current_score = 0;
    for(long long int w : neg_weights) current_score += w;

    vector<long long int> group, rev_group, to_flip, drop_flips;
    group.reserve(100); rev_group.reserve(100); to_flip.reserve(100); drop_flips.reserve(100);

    vector<long long int> queue_vars;
    queue_vars.reserve(N * 2);
    vector<bool> in_queue(N, true);
    for(long long int i = 0; i < N; ++i) queue_vars.push_back(i);

    long long int head = 0;
    while(head < queue_vars.size()) {
        long long int v = queue_vars[head++];
        in_queue[v] = false;

        if (assignment[v]) continue; 

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
        apply_flip_single(v); 

        group.clear();
        group.push_back(v);
        for (const auto& edge : adj[v]) {
            if (!assignment[edge.to]) group.push_back(edge.to);
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
                rev_group = group;
                reverse(rev_group.begin(), rev_group.end());
                apply_flips(rev_group); 
            }
        }
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (long long int v = 0; v < N; ++v) {
            long long int gain = apply_flip_single(v);
            if (gain > 0) {
                current_score += gain;
                changed = true;
            } else {
                apply_flip_single(v);
            }
        }
        if (changed) continue;

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
                    apply_flips(rev_group); 
                }
            }
        }
    }

    vector<bool> final_assign(N);
    for (long long int i = 0; i < N; ++i) final_assign[i] = assignment[i];
    return {final_assign, current_score};
}

// =====================================================================
// 3. GENERATORS & L-REDUCTIONS
// =====================================================================

void generate_random_max_e3lin2(long long int num_vars, long long int num_eqs, vector<Lin2Equation>& equations) {
    mt19937 rng(random_device{}());
    uniform_int_distribution<long long int> var_dist(1, num_vars);
    uniform_int_distribution<int> rhs_dist(0, 1);

    equations.clear();
    for (long long int i = 0; i < num_eqs; ++i) {
        Lin2Equation eq;
        eq.v1 = var_dist(rng);
        eq.v2 = var_dist(rng);
        while (eq.v2 == eq.v1) eq.v2 = var_dist(rng);
        eq.v3 = var_dist(rng);
        while (eq.v3 == eq.v1 || eq.v3 == eq.v2) eq.v3 = var_dist(rng);
        eq.rhs = rhs_dist(rng);
        equations.push_back(eq);
    }
}

void convert_e3lin2_to_max2sat(long long int num_vars, const vector<Lin2Equation>& equations, 
                               vector<Clause>& clauses, vector<long long int>& weights) {
    clauses.clear();
    weights.clear();

    for (size_t i = 0; i < equations.size(); ++i) {
        const auto& eq = equations[i];
        long long int l1 = eq.v1;
        long long int l2 = eq.v2;
        // If rhs is 1, x1 ^ x2 ^ x3 = 1 is equivalent to x1 ^ x2 ^ !x3 = 0
        long long int l3 = (eq.rhs == 0) ? eq.v3 : -eq.v3;
        
        // Allocate a fresh auxiliary variable w for the 12-clause gadget
        long long int w = num_vars + i + 1; 

        vector<vector<long long int>> raw_clauses = {
            {l1, w}, {l2, w}, {l3, w},
            {-l1, -w}, {-l2, -w}, {-l3, -w},
            {l1, l2}, {l2, l3}, {l3, l1},
            {-l1, -l2}, {-l2, -l3}, {-l3, -l1}
        };

        for (const auto& lits : raw_clauses) {
            Clause c;
            c.literals = lits;
            clauses.push_back(c);
            weights.push_back(1); 
        }
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

    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        long long int w = weights[idx];
        if (original_clauses[idx].literals.size() == 2) {
            long long int u = get_var_idx(original_clauses[idx].literals[0]);
            long long int v = get_var_idx(original_clauses[idx].literals[1]);
            if (u > v) swap(u, v);
            or_clauses_dict[{u, v}] += w;
        }
    }

    vector<long long int> degrees(V + 1, 0);
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        long long int w = weights[idx];
        for (long long int literal : original_clauses[idx].literals) degrees[abs(literal)] += w;
    }

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

pair<vector<bool>, long long int> solve_max_e3lin2_cplex(long long int num_vars, const vector<Lin2Equation>& equations) {
    IloEnv env;
    vector<bool> assignment(num_vars, false);
    long long int opt_score = 0;

    try {
        IloModel model(env);
        long long int M = equations.size();
        
        IloBoolVarArray x(env, num_vars); 
        IloBoolVarArray s(env, M); // s[j] = 1 if equation j is satisfied

        IloExpr objExpr(env);
        for (long long int j = 0; j < M; ++j) objExpr += s[j];
        model.add(IloMaximize(env, objExpr));
        objExpr.end();

        // Exact linear formulations for the 3-XOR parity constraints
        for (long long int j = 0; j < M; ++j) {
            long long int v1 = equations[j].v1 - 1;
            long long int v2 = equations[j].v2 - 1;
            long long int v3 = equations[j].v3 - 1;
            
            if (equations[j].rhs == 0) {
                model.add(s[j] <= 1 - x[v1] + x[v2] + x[v3]);
                model.add(s[j] <= 1 + x[v1] - x[v2] + x[v3]);
                model.add(s[j] <= 1 + x[v1] + x[v2] - x[v3]);
                model.add(s[j] <= 3 - x[v1] - x[v2] - x[v3]);
            } else {
                model.add(s[j] <= x[v1] + x[v2] + x[v3]);
                model.add(s[j] <= 2 - x[v1] - x[v2] + x[v3]);
                model.add(s[j] <= 2 - x[v1] + x[v2] - x[v3]);
                model.add(s[j] <= 2 + x[v1] - x[v2] - x[v3]);
            }
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); 
        cplex.setWarning(env.getNullStream());
        //cplex.setParam(IloCplex::Param::TimeLimit, 120.0); // E3-LIN-2 is very hard

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
        IloBoolVarArray x(env, num_vars); 
        IloBoolVarArray z(env, M);        

        IloExpr objExpr(env);
        for (long long int j = 0; j < M; ++j) objExpr += IloNum(weights[j]) * z[j];
        model.add(IloMaximize(env, objExpr));
        objExpr.end();

        for (long long int j = 0; j < M; ++j) {
            IloExpr clauseExpr(env);
            for (long long int literal : clauses[j].literals) {
                long long int var_idx = abs(literal) - 1;
                if (literal > 0) clauseExpr += x[var_idx];
                else clauseExpr += (1 - x[var_idx]);
            }
            model.add(z[j] <= clauseExpr);
            clauseExpr.end();
        }

        IloCplex cplex(model);
        cplex.setOut(env.getNullStream()); 
        cplex.setWarning(env.getNullStream());
        //cplex.setParam(IloCplex::Param::TimeLimit, 120.0);

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

// =====================================================================
// 5. UNIFIED TEST RUNNER
// =====================================================================

void test_unified_pipeline(long long int num_tests) {
    cout << "\n--- E3-LIN-2 EXACT vs MAX 2-SAT GADGET HEURISTIC ---\n";
    cout << left << setw(13) << "Vars / Eqs" << " | "
         << setw(10) << "Opt LIN" << " | "
         << setw(10) << "Heur LIN" << " | "
         << setw(10) << "Opt 2SAT" << " | "
         << setw(10) << "Heur 2SAT" << " | "
         << "Approximability Ratios (LIN vs 2SAT)\n";
    cout << string(95, '-') << "\n";

    for (long long int t = 0; t < num_tests; ++t) {
        // Kept moderately sized because 3-XOR causes heavy CPLEX branching
        long long int V = 40;
        long long int M = 100; 

        vector<Lin2Equation> equations;
        generate_random_max_e3lin2(V, M, equations);

        // 1. Exact CPLEX solve for E3-LIN-2
        auto cplex_lin2 = solve_max_e3lin2_cplex(V, equations);
        long long int opt_e3lin2 = cplex_lin2.second;

        // 2. Convert to MAX 2-SAT using 12-clause gadget
        vector<Clause> sat_clauses;
        vector<long long int> sat_weights;
        convert_e3lin2_to_max2sat(V, equations, sat_clauses, sat_weights);
        long long int total_sat_vars = V + M; // V base + M auxiliary vars

        // 3. Exact CPLEX solve for MAX 2-SAT representation
        auto cplex_2sat = solve_classical_max2sat_cplex(total_sat_vars, sat_clauses, sat_weights);
        long long int opt_2sat = cplex_2sat.second;

        // 4. Heuristic solve via L-reduction
        long long int N_red;
        vector<OrClause> or_clauses;
        vector<long long int> neg_weights;
        reduce_to_weighted_restricted(total_sat_vars, sat_clauses, sat_weights, N_red, or_clauses, neg_weights);

        auto heur_res = solve_heuristic_cpp(N_red, or_clauses, neg_weights);
        
        vector<bool> heur_sat_assign(total_sat_vars);
        for(long long int i = 0; i < total_sat_vars; ++i) heur_sat_assign[i] = heur_res.first[i];
        long long int heur_2sat = get_original_score(heur_sat_assign, sat_clauses, sat_weights);

        // 5. Evaluate the Heuristic's assignment back on the E3-LIN-2 equations
        vector<bool> heur_lin2_assign(V);
        for(long long int i = 0; i < V; ++i) heur_lin2_assign[i] = heur_res.first[i];
        long long int heur_e3lin2 = get_e3lin2_score(heur_lin2_assign, equations);

        // Format Output
        string graph_size = "V:" + to_string(V) + " M:" + to_string(M);
        double ratio_lin2 = (opt_e3lin2 > 0) ? (double)heur_e3lin2 / opt_e3lin2 : 1.0;
        double ratio_2sat = (opt_2sat > 0) ? (double)heur_2sat / opt_2sat : 1.0;

        cout << left << setw(13) << graph_size << " | "
             << setw(10) << opt_e3lin2 << " | "
             << setw(10) << heur_e3lin2 << " | "
             << setw(10) << opt_2sat << " | "
             << setw(10) << heur_2sat << " | "
             << fixed << setprecision(3) << "LIN: " << ratio_lin2 << "  ->  2SAT: " << ratio_2sat << "\n";
    }
}

int main() {
    test_unified_pipeline(10);
    return 0;
}