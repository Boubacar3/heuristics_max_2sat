#!/bin/bash
./maxsat_solver test_wcnf_files_20000_100000 results_rand_20000_v1.1_1200.csv 1200
./maxsat_solver test_wcnf_files results_rand_200000_v1.1_1200.csv 1200
./maxsat_solver wcnf_files results_grid_v1.1_1200.csv 1200
./maxsat_solver max_cut results_max_cut_v1.1_1200.csv 1200