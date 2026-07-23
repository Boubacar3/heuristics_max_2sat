#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>

using namespace std;

// Required structs based on the heuristic's signature and internal usage
struct OrClause {
    int u;
    int v;
    long long int w;
};

struct Edge {
    int to;
    long long int weight;
};

// ============================================================================
// The original heuristic provided in the prompt
// ============================================================================
pair<vector<bool>, long long int> solve_heuristic_cpp(int N, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights) {
    // 1. Use uint8_t instead of vector<bool> for much faster memory access
    vector<uint8_t> assignment(N, 0);

    // 2. Build an Adjacency List to avoid O(M) loop lookups
    vector<vector<Edge>> adj(N);
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    // 3a. Single-flip Delta Function (Zero heap allocations)
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
    auto apply_flips = [&](const vector<int>& vars) -> int {
        int delta = 0;
        for (int v : vars) {
            delta += apply_flip_single(v);
        }
        return delta;
    };

    // Initialize base score
    int current_score = 0;
    for(int w : neg_weights) current_score += w;

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
        int gain1 = apply_flip_single(v);
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
            int gain2 = apply_flips(group);
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
    for(int i = 0; i < 5; i++) { 
        changed = false;

        // 1. Standard 1-opt Pruning
        for (int v = 0; v < N; ++v) {
            int gain = apply_flip_single(v);
            if (gain > 0) {
                current_score += gain;
                changed = true;
            } else {
                apply_flip_single(v); // Revert
            }
        }
        if (changed) continue;

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

                int gain = apply_flips(to_flip);
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
        for (int v = 0; v < N; ++v) {
            if (!assignment[v]) {
                int initial_gain = apply_flip_single(v); // Tentatively commit v=1
                int best_sub_gain = 0;
                vector<int> best_sub_flips;
                
                // Only scan NEIGHBORS, not all N variables
                for (const auto& edge : adj[v]) {
                    int u = edge.to;
                    if (assignment[u] && u != v) {
                        drop_flips.clear();
                        drop_flips.push_back(u);
                        
                        for (const auto& edge_u : adj[u]) {
                            int other = edge_u.to;
                            if (!assignment[other] && edge_u.weight > neg_weights[other]) {
                                drop_flips.push_back(other);
                            }
                        }

                        int sub_gain = apply_flips(drop_flips);
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
    for (int i = 0; i < N; ++i) final_assign[i] = assignment[i];

    return {final_assign, current_score};
}

// ============================================================================
// Test harness
// ============================================================================
int main() {
    int N = 5;
    
    // A straight line chain: A - B - C - D - E
    vector<OrClause> or_clauses = {
        {0, 1, 100},
        {1, 2, 100},
        {2, 3, 100},
        {3, 4, 100}
    };
    
    // Penalties for variables being assigned '0' (False)
    vector<long long int> neg_weights = {10, 20, 10, 20, 10};

    // Run the heuristic
    auto result = solve_heuristic_cpp(N, or_clauses, neg_weights);
    vector<bool> assignment = result.first;
    long long int score = result.second;

    // Output the results
    cout << "--- HEURISTIC TEST ---" << endl;
    cout << "Optimal Target Assignment: [1, 0, 1, 0, 1]" << endl;
    cout << "Optimal Target Score:      440" << endl;
    cout << "----------------------" << endl;
    cout << "Heuristic Assignment:      [";
    for (int i = 0; i < N; ++i) {
        cout << (assignment[i] ? "1" : "0") << (i < N - 1 ? ", " : "");
    }
    cout << "]" << endl;
    cout << "Heuristic Score:           " << score << endl;
    cout << "----------------------" << endl;

    if (score == 430) {
        cout << "RESULT: Test Passed! The heuristic successfully failed and got trapped at the suboptimal plateau (Score 430)." << endl;
    } else if (score == 440) {
        cout << "RESULT: Test Failed. The heuristic unexpectedly found the optimal solution." << endl;
    } else {
        cout << "RESULT: Unexpected outcome." << endl;
    }

    return 0;
}