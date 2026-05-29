#!/bin/bash
./script_olllp.sh test_wcnf_files_20000_100000 results_rand_20000_olllp_1200.csv 1200
./script_olllp.sh test_wcnf_files results_rand_200000_olllp_1200.csv 1200
./script_olllp.sh wcnf_files results_grid_olllp_1200.csv 1200
./script_olllp.sh max_cut results_max_cut_olllp_1200.csv 1200