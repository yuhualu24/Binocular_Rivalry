#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <random>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

// This code simulates rivalry between two networks with adaption.

// To run the code, need to create the following files:
// matrix_1.txt matrix_2.txt matrix_3.txt matrix_4.txt 
// f_1.txt f_2.txt image_1.txt image_2.txt start_1.txt start_2.txt 
// metric.txt time_step.txt fire_num_1.txt fire_num_2.txt indices.txt 
// re.txt iter_step.txt mean_thresh_1.txt mean_thresh_2.txt

// To run the code, need to have the following files:
// stripes_matrix.txt vertical.txt

int main(){
    int num_1 = 2000;
    int num_e_1 = 1000;
    int num_i_1 = 1000;
    int bin_size = 300;

    int num_2 = 2000;
    int num_e_2 = 1000;
    int num_i_2 = 1000;

    int num = num_1 + num_2;
    int number = 0;

    int i, j, l, a, b, c;
    double time_step = 0;
    double alter = 0.0;
    int iter = 0;
    double trace = 0.0;
    int check_1 = 0;
    int check_2 = 0;
    int random = 0;
    int checkpoint = 0;
    int mean = 0;
    int height = 2000;
    int width = 10000;
    int temp_3 = 0;
    int temp_4 = 0;
    int temp_5 = 0;
    int temp_6 = 0;
    double trace_step = 0.0;
    double vary = 0.04;

    double thresh_sum_1 = 0.0;
    double thresh_sum_2 = 0.0;
    double thresh_avg_1 = 0.0;
    double thresh_avg_2 = 0.0;
    double re1 = 0.0;
    double diff = 0.0;
    double t = 0.00;
    double run_time = 2100.0;
    double lambda = 0.00625;
    double phi = 0.005;
    double value;
    double add = 0;
    double average_1 = 0;
    double average_2 = 0;
    double temp = 0;
    double temp_1 = 0;
    double temp_2 = 0;
    double metric = 0;
    string line;
    ofstream file;
    ofstream file_1;
    ofstream file_2;
    ofstream file_3;
    ofstream file_8;
    ofstream file_9;
    // file 10 and 11 are for storing dominance durations
    // ofstream file_10;
    // ofstream file_11;
    bool flag = false;
    bool alert = false;

    string fname = "firingRates";
    string fname_2 = "forceAvg";
    string fname_3 = "threshAvg";
    string result;
    string result2;
    string result3;
    ostringstream convert;
    ostringstream convert2;
    ostringstream convert3;

    double r_ee = 1.0;
    double r_ie = 1.0;
    double r_ii = -1.8;
    double r_ei = -2.0;
    double K = 40;
    double f_e = 1.25;
    double f_i = 1.0;
    double m_0 = 1 / 750;
    double threshold_e = 1;
    double threshold_i = 0.8;

    double *firing_num_1 = new double[num_1];
    double *input_1 = new double[num_1];
    double *avg_1 = new double[num_1];
    double *firing_num_2 = new double[num_2];
    double *input_2 = new double[num_2];
    double *avg_2 = new double[num_2];
    double * avg_thresh_1 = new double[num_1];
    double * avg_thresh_2 = new double[num_2];

    for (i = 0; i < num_1; i++){
        avg_thresh_1[i] = 0;
        avg_thresh_2[i] = 0;
    }

    vector<int> neighbors_1;
    vector<int> neighbors_2;
    vector<int> check_switch_1;
    vector<int> check_switch_2;

    // initialize fire state to zero
    int *fire_1 = new int[num_1];
    for (i = 0; i < num_1; i++){
        fire_1[i] = 0;
        input_1[i] = 0;
        input_2[i] = 0;
        avg_1[i] = 0;
        avg_2[i] = 0;
        firing_num_1[i] = 0;
        firing_num_2[i] = 0;
    }

    int *index = new int[1000];
    for (i = 0; i < 20; i++) {
        index[i] = 0;
    }
    int dom_num = 0;

    int *fire_2 = new int[num_2];
    for (i = 0; i < num_2; i++){
        fire_2[i] = 0;
    }

    double *threshold_1 = new double[num_1];
    double *threshold_2 = new double[num_2];
    for (i = 0; i < num_e_1; i++){
        threshold_1[i] = threshold_e;
        threshold_2[i] = threshold_e;
    }
    for (i = num_e_1; i < num_1; i++){
        threshold_1[i] = threshold_i;
        threshold_2[i] = threshold_i;
    }

    // initialize inner connectivity of network 1
    double** matrix_1 = new double*[num_1];
    for (i = 0; i < num_1; i++){
        matrix_1[i] = new double[num_1];
    }

    double KE = K / num_e_1;
    // divide K by the population of neurons that are firing (sending the signal)
    double KI = K / num_i_1;
    // KI = K / num_i_1;
    srand(time(0));
    for (i = 0; i < num_1; i++){
        for (j = 0; j < num_1; j++){
            value = (double)rand() / RAND_MAX;
            if (i < num_e_1 && j < num_e_1){
                if (value < KE){
                    matrix_1[i][j] = r_ee / sqrt(K);
                }
                else {
                    matrix_1[i][j] = 0;
                }
            }
            else if (i >= num_e_1 && j < num_e_1){
                if (value < KE ){
                    matrix_1[i][j] = r_ie / sqrt(K);
                }
                else {
                    matrix_1[i][j] = 0;
                }
            }
            else if (i < num_e_1 && j >= num_e_1){
                if (value < KI ){
                    matrix_1[i][j] = r_ei / sqrt(K);
                }
                else {
                    matrix_1[i][j] = 0;
                }
            }
            else if (i >= num_e_1 && j >= num_e_1){
                if (value < KI){
                    matrix_1[i][j] = r_ii / sqrt(K);
                }
                else {
                    matrix_1[i][j] = 0;
                }
            }
        }
    }
    
    // diagonal matrix entries are zero
    for (i = 0; i < num_1; i++){
        matrix_1[i][i] = 0;
    }

    file.open("matrix_1.txt");
    for (i = 0; i < num_1; i++){
        for (j = 0; j < num_1; j++){
            if (j == num_1 - 1){
                file << matrix_1[i][j] << endl;
            }
            else{
                file << matrix_1[i][j] << "\t";
            }
        }
    }
    file.close();

    // initialize inner connectivity of network 2
    double** matrix_2 = new double*[num_2];
    for (i = 0; i < num_2; i++){
        matrix_2[i] = new double[num_2];
    }

    for (i = 0; i < num_2; i++){
        for (j = 0; j < num_2; j++){
            value = (double)rand() / RAND_MAX;
            if (i < num_e_2 && j < num_e_2){
                if (value < KE){
                    matrix_2[i][j] = r_ee / sqrt(K);
                }
                else {
                    matrix_2[i][j] = 0;
                }
            }
            else if (i >= num_e_2 && j < num_e_2){
                if (value < KE ){
                    matrix_2[i][j] = r_ie / sqrt(K);
                }
                else {
                    matrix_2[i][j] = 0;
                }
            }
            else if (i < num_e_2 && j >= num_e_2){
                if (value < KI ){
                    matrix_2[i][j] = r_ei / sqrt(K);
                }
                else {
                    matrix_2[i][j] = 0;
                }
            }
            else if (i >= num_e_2 && j >= num_e_2){
                if (value < KI){
                    matrix_2[i][j] = r_ii / sqrt(K);
                }
                else {
                    matrix_2[i][j] = 0;
                }
            }
        }
    }
    
    // diagonal matrix entries are zero
    for (i = 0; i < num_2; i++){
        matrix_2[i][i] = 0;
    }

    file.open("matrix_2.txt");
    for (i = 0; i < num_2; i++){
        for (j = 0; j < num_2; j++){
            if (j == num_2 - 1){
                file << matrix_2[i][j] << endl;
            }
            else{
                file << matrix_2[i][j] << "\t";
            }
        }
    }
    file.close();

    // initialize connectivity between two matrix
    // network 1 E to network 2 I
    double** matrix_3 = new double*[num_2];
    for (i = 0; i < num_2; i++){
        matrix_3[i] = new double[num_1];
    }

    for (i = 0; i < num_2; i++){
        for (j = 0; j < num_1; j++){
            value = (double)rand() / RAND_MAX;
            if (i >= num_e_2 && j < num_e_1){
                if (value < KE){
                    matrix_3[i][j] = r_ie / sqrt(K);
                }
                else {
                    matrix_3[i][j] = 0;
                }
            }
            else {
                matrix_3[i][j] = 0;
            }
        }
    }

    file.open("matrix_3.txt");
    for (i = 0; i < num_2; i++){
        for (j = 0; j < num_1; j++){
            if (j == num_1 - 1){
                file << matrix_3[i][j] << endl;
            }
            else{
                file << matrix_3[i][j] << "\t";
            }
        }
    }
    file.close();

    // network 2 E to network 1 I
    double** matrix_4 = new double*[num_1];
    for (i = 0; i < num_1; i++){
        matrix_4[i] = new double[num_2];
    }

    for (i = 0; i < num_1; i++){
        for (j = 0; j < num_2; j++){
            value = (double)rand() / RAND_MAX;
            if (i >= num_e_1 && j < num_e_2){
                if (value < KE){
                    matrix_4[i][j] = r_ie / sqrt(K);
                }
                else {
                    matrix_4[i][j] = 0;
                }
            }
            else {
                matrix_4[i][j] = 0;
            }
        }
    }

    file.open("matrix_4.txt");
    for (i = 0; i < num_1; i++){
        for (j = 0; j < num_2; j++){
            if (j == num_2 - 1){
                file << matrix_4[i][j] << endl;
            }
            else{
                file << matrix_4[i][j] << "\t";
            }
        }
    }
    file.close();

    double *p1 = new double[width];
    ifstream infile("stripes_matrix.txt");
    for (i = 0; i < width; i++){
        getline(infile, line);
        p1[i] = stod(line);
    }
    infile.close();

    double *p2 = new double[width];
    ifstream infile_2("vertical.txt");
    for (i = 0; i < width; i++){
        getline(infile_2, line);
        p2[i] = stod(line);
    }
    infile_2.close();

    for (i = 0; i < width; i++){
        add = add + p1[i];
    }
    
    average_1 = add / width;
    add = 0;

    for (i = 0; i < width; i++){
        add = add + p2[i];
    }
    average_2 = add / width;
    add = 0;

    // initialize the F matrix
    int** matrix_f_1 = new int*[height];
    for (i = 0; i < height; i++){
        matrix_f_1[i] = new int[width];
    }

    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            matrix_f_1[i][j] = 0;
        }
    }

    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            value = double(rand()) / double(RAND_MAX);
            if (value < 0.001){
                matrix_f_1[i][j] = 1;
            }
            else{
                matrix_f_1[i][j] = 0;
            }
        }
    }

    file.open("f_1.txt");
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            if (j == width - 1){
                file << matrix_f_1[i][j] << endl;
            }
            else{
                file << matrix_f_1[i][j] << "\t";
            }
        }
    }
    file.close();

    double *image_1 = new double[num_1];
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            if (i < num_e_1){
                add = add + matrix_f_1[i][j] * p1[j] / average_1 * f_e * vary;
            }
            else {
                add = add + matrix_f_1[i][j] * p1[j] / average_1 * f_i * vary;
            }
        }
        image_1[i] = add;
        add = 0;
    }

    file.open("image_1.txt");
    for (i = 0; i < 1000; i++){
        file << image_1[i] << endl;
    }
    file.close();

    int** matrix_f_2 = new int*[height];
    for (i = 0; i < height; i++){
        matrix_f_2[i] = new int[width];
    }

    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            matrix_f_2[i][j] = 0;
        }
    }

    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            value = double(rand()) / double(RAND_MAX);
            if (value < 0.001){
                matrix_f_2[i][j] = 1;
            }
            else{
                matrix_f_2[i][j] = 0;
            }
        }
    }

    file.open("f_2.txt");
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            if (j == width - 1){
                file << matrix_f_2[i][j] << endl;
            }
            else{
                file << matrix_f_2[i][j] << "\t";
            }
        }
    }
    file.close();

    double *image_2 = new double[num_2];
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            if (i < num_e_2){
                add = add + matrix_f_2[i][j] * p2[j] / average_2 * f_e * vary;
            }
            else {
                add = add + matrix_f_2[i][j] * p2[j] / average_2 * f_i * vary;
            }
        }
        image_2[i] = add;
        add = 0;
    }

    file.open("image_2.txt");
    for (i = 0; i < height; i++){
        file << image_2[i] << endl;
    }
    file.close();

    double *start_1 = new double[num_1];
    double *end_1 = new double[num_1];
    double *start_2 = new double[num_2];
    double *end_2 = new double[num_2];

    file.open("start_1.txt");
    for (i = 0; i < num_1; i++){
        start_1[i] = rand() % 10000 / 10000.0;
        file << start_1[i] << endl;
    }
    file.close();

    file.open("start_2.txt");
    for (i = 0; i < num_2; i++){
        start_2[i] = rand() % 10000 / 10000.0;
        file << start_2[i] << endl;
    }
    file.close();

    double *k_1 = new double[num_1];
    double *k_2 = new double[num_1];
    double *k = new double[num_1];

    for (i = 0; i < num_1; i++){
        k_1[i] = 0;
        k_2[i] = 0;
        k[i] = 0;
    }

    ofstream file_12;
    ofstream file_13;
    ofstream file_14;

    file.open("metric.txt");
    file_1.open("time_step.txt");
    file_2.open("fire_num_1.txt");
    file_3.open("fire_num_2.txt");
    file_8.open("indices.txt");
    file_9.open("re.txt");

    // file_10.open("duration.txt", fstream::app);
    // file_11.open("duration_2.txt", fstream::app);

    file_12.open("mean_thresh_1.txt");
    file_13.open("mean_thresh_2.txt");
    file_14.open("iter_step.txt");

    ofstream file_4;
    ofstream file_5;
    ofstream file_6;
    ofstream file_7;

    while (t < run_time){

        for (i = 0; i < num_1; i++){
            avg_1[i] = avg_1[i] + start_1[i];
            k_1[i] = (- start_1[i] + image_1[i]) * 0.01;
            k_2[i] = (- (start_1[i] + k_1[i]) + image_1[i]) * 0.01;
            k[i] = (k_1[i] + k_2[i]) / 2;
            end_1[i] = start_1[i] + k[i];

            avg_2[i] = avg_2[i] + start_2[i];
            k_1[i] = (- start_2[i] + image_2[i]) * 0.01;
            k_2[i] = (- (start_2[i] + k_1[i]) + image_2[i]) * 0.01;
            k[i] = (k_1[i] + k_2[i]) / 2;
            end_2[i] = start_2[i] + k[i];
        }

        // for each neuron, we run the following steps
        for (i = 0; i < num_1; i++){
            if (end_1[i] >= threshold_1[i] && fire_1[i] == 0) {
                end_1[i] = 0;
                fire_1[i] = 1;
                firing_num_1[i] ++;
                input_1[i] = input_1[i] + 1;
                threshold_1[i] = threshold_1[i] + phi;
                for (j = 0; j < num_1; j++){
                    if (matrix_1[j][i] != 0 && fire_1[j] == 0){
                        end_1[j] = end_1[j] + matrix_1[j][i];
                        neighbors_1.push_back(j);
                    }
                }
                for (l = num_e_2; l < num_2; l++){
                    if (matrix_3[l][i] != 0 && fire_2[l] == 0){
                        end_2[l] = end_2[l] + matrix_3[l][i];
                        neighbors_2.push_back(l);
                    }
                }
            }

            if (end_2[i] >= threshold_2[i] && fire_2[i] == 0) {
                end_2[i] = 0;
                fire_2[i] = 1;
                firing_num_2[i] ++;
                input_2[i] = input_2[i] + 1;
                threshold_2[i] = threshold_2[i] + phi;
                for (j = 0; j < num_2; j++){
                    if (matrix_2[j][i] != 0 && fire_2[j] == 0){
                        end_2[j] = end_2[j] + matrix_2[j][i];
                        neighbors_2.push_back(j);
                    }
                }
                for (l = num_e_1; l < num_1; l++){
                    if (matrix_4[l][i] != 0 && fire_1[l] == 0){
                        end_1[l] = end_1[l] + matrix_4[l][i];
                        neighbors_1.push_back(l);
                    }
                }
            }

            // check if neighbors spike
            while (flag == false) {
                for (l = 0; l < neighbors_1.size(); l++){
                    if (end_1[neighbors_1[l]] >= threshold_1[neighbors_1[l]] && fire_1[neighbors_1[l]] == 0){
                        b = neighbors_1[l];
                        end_1[b] = 0;
                        fire_1[b] = 1;
                        firing_num_1[b] ++;
                        input_1[b] = input_1[b] + 1;
                        threshold_1[b] = threshold_1[b] + phi;
                        for (j = 0; j < num_1; j++){
                            if (matrix_1[j][b] != 0 && fire_1[j] == 0){
                                end_1[j] = end_1[j] + matrix_1[j][b];
                                neighbors_1.push_back(j);
                            }
                        }
                        for (a = num_e_2; a < num_2; a++){
                            if (matrix_3[a][b] != 0 && fire_2[a] == 0){
                                end_2[a] = end_2[a] + matrix_3[a][b];
                                neighbors_2.push_back(a);
                            }
                        }
                    }
                }

                for (l = 0; l < neighbors_2.size(); l++){
                    if (end_2[neighbors_2[l]] >= threshold_2[neighbors_2[l]] && fire_2[neighbors_2[l]] == 0){
                        b = neighbors_2[l];
                        end_2[b] = 0;
                        fire_2[b] = 1;
                        firing_num_2[b] ++;
                        input_2[b] = input_2[b] + 1;
                        threshold_2[b] = threshold_2[b] + phi;
                        for (j = 0; j < num_2; j++){
                            if (matrix_2[j][b] != 0 && fire_2[j] == 0){
                                end_2[j] = end_2[j] + matrix_2[j][b];
                                neighbors_2.push_back(j);
                            }
                        }
                        for (a = num_e_1; a < num_1; a++){
                            if (matrix_4[a][b] != 0 && fire_1[a] == 0){
                                end_1[a] = end_1[a] + matrix_4[a][b];
                                neighbors_1.push_back(a);
                                temp ++;
                            }
                        }
                    }
                }

                if (temp != 0){
                    flag = false;
                    temp = 0;
                }

                else {
                    flag = true;
                }

            }

            // clear neighbors vector
            neighbors_1.clear();
            neighbors_2.clear();

        }
        // end of each iteration of neurons

        // adaption
        for (l = 0; l < num_e_1; l++) {
            if (fire_1[l] == 0){
                temp = - 0.01 * lambda * (threshold_1[l] - threshold_e) + threshold_1[l];
                threshold_1[l] = threshold_1[l] - 0.5 * 0.01 * (lambda * (threshold_1[l] - threshold_e) + lambda * (temp - threshold_e));
            }

            if (fire_2[l] == 0){
                temp = - 0.01 * lambda * (threshold_2[l] - threshold_e) + threshold_2[l];
                threshold_2[l] = threshold_2[l] - 0.5 * 0.01 * (lambda * (threshold_2[l] - threshold_e) + lambda * (temp - threshold_e));
            }
        }

        for (l = num_e_1; l < num_1; l++) {
            if (fire_1[l] == 0){
                temp = - 0.01 * lambda * (threshold_1[l] - threshold_i) + threshold_1[l];
                threshold_1[l] = threshold_1[l] - 0.5 * 0.01 * (lambda * (threshold_1[l] - threshold_i) + lambda * (temp - threshold_i));
            }

            if (fire_2[l] == 0){
                temp = - 0.01 * lambda * (threshold_2[l] - threshold_i) + threshold_2[l];
                threshold_2[l] = threshold_2[l] - 0.5 * 0.01 * (lambda * (threshold_2[l] - threshold_i) + lambda * (temp - threshold_i));
            }
        }

        // calculate mean threshold
        for (l = 0; l < num_1; l++) {
            thresh_sum_1 = thresh_sum_1 + threshold_1[l];
            thresh_sum_2 = thresh_sum_2 + threshold_2[l];
        }

        thresh_avg_1 = thresh_sum_1 / num_1;
        thresh_avg_2 = thresh_sum_2 / num_1;
        file_12 << thresh_avg_1 << endl;
        file_13 << thresh_avg_2 << endl;
        thresh_sum_1 = 0;
        thresh_sum_2 = 0;

        // after each time step, fire state resets to zero for all neurons
        // starting voltage is set to the ending voltage of last round
        for (l = 0; l < num_1; l++){
            fire_1[l] = 0;
            fire_2[l] = 0;
            start_1[l] = end_1[l];
            start_2[l] = end_2[l];
            avg_thresh_1[l] = avg_thresh_1[l] + threshold_1[l];
            avg_thresh_2[l] = avg_thresh_2[l] + threshold_2[l];
        }

        t = t + 0.01;
        iter ++;
        trace = trace + 1;
        trace_step = trace / 100.0;

        file_14 << trace_step << endl;

        // check metric every bin
        if (iter == bin_size){

            for (l = 0; l < num_e_1; l++){
                temp_1 = temp_1 + firing_num_1[l];
            }

            file_2 << temp_1 << endl;

            for (l = 0; l < num_1; l++){
                firing_num_1[l] = 0;
            }

            for (l = 0; l < num_e_2; l++){
                temp_2 = temp_2 + firing_num_2[l];
            }

            file_3 << temp_2 << endl;

            for (l = 0; l < num_1; l++){
                firing_num_2[l] = 0;
            }

            metric = (temp_1 - temp_2) / (temp_1 + temp_2);
            file << metric << endl;

            iter = 0;
            temp_1 = 0;
            temp_2 = 0;
            time_step = time_step + 0.01 * bin_size;
            file_1 << time_step << endl;

            // check when first reconstruction starts
            for (l = random; l < 5 + random; l++){
                if (time_step == l * 0.01 * bin_size){
                    if (metric > 0.3) {
                        check_1 ++;
                    }
                    else if (metric < -0.3) {
                        check_2 ++;
                    }
                }
            }

            if (time_step == (4 + random) * 0.01 * bin_size) {
                if (check_1 >= 4){
                    checkpoint = time_step - 0.01 * bin_size * 4;
                    file_8 << checkpoint << endl;
                    index[0] = checkpoint;
                    re1 = 0.3;
                    file_9 << re1 << endl;
                    check_1 = 0;
                }
                else if (check_2 >= 4){
                    checkpoint = time_step - 0.01 * bin_size * 4;
                    file_8 << checkpoint << endl;
                    index[0] = checkpoint;
                    re1 = - 0.3;
                    file_9 << re1 << endl;
                    check_2 = 0;
                }
                else {
                    // if first reconstruction has not started
                    // set all values to zero
                    random = random + 5;
                    for (l = 0; l < num_1; l++){
                        input_1[i] = 0;
                        input_2[i] = 0;
                        avg_1[i] = 0;
                        avg_2[i] = 0;
                        avg_thresh_1[i] = 0;
                        avg_thresh_2[i] = 0;
                    }
                }
            }

            // check when switching happens
            if (re1 > 0.0 && metric < 0.3){
                check_switch_1.push_back(time_step);
            }
            else if (re1 < 0.0 && metric > -0.3){
                check_switch_2.push_back(time_step);
            }

            if (check_switch_1.size() == 5 && check_switch_1[4] - check_switch_1[0] == 0.01 * bin_size * 4){
                checkpoint = check_switch_1[0] - 0.01 * bin_size;
                dom_num ++;
                alter = alter + 1.0;
                cout << dom_num << endl;
                temp_4 = dom_num * 2 - 1;
                index[temp_4] = checkpoint;
                file_8 << checkpoint << endl;
                // start index of the next dominance
                checkpoint = check_switch_1[0] + 0.01 * bin_size;
                index[temp_4 + 1] = checkpoint;
                file_8 << checkpoint << endl;
                re1 = - re1;
                file_9 << re1 << endl;
                check_switch_1.clear();
                // save the data for the previous dominance
                // calculate the runtime of the previous dominance
                temp_5 = index[temp_4] - index[temp_4 - 1];
                // cout << temp_4 << endl;
                cout << temp_5 << endl;
                temp_5 = double(temp_5);
                // file_10 << temp_5 << endl;
                for (l = 0; l < num_1; l++){
                    input_1[l] = input_1[l] / temp_5;
                    input_2[l] = input_2[l] / temp_5;
                    avg_1[l] = avg_1[l] / temp_5 / 100.0;
                    avg_2[l] = avg_2[l] / temp_5 / 100.0;
                    avg_thresh_1[l] = avg_thresh_1[l] / temp_5 / 100.0;
                    avg_thresh_2[l] = avg_thresh_2[l] / temp_5 / 100.0;
                }

                // steps to create files and save information for reconstruction

                // convert << dom_num;
                // result = convert.str();
                // fname.append(result);
                // fname.append("_1.txt");
                // file_4.open(fname);
                // for (l = 0; l < num_1; l++){
                //     file_4 << input_1[l] << endl;
                // }
                // file_4.close();

                // for (l = 0; l < 5; l++) {
                //     fname.pop_back();
                // }
                // fname.append("2.txt");

                // file_4.open(fname);
                // for (l = 0; l < num_1; l++){
                //     file_4 << input_2[l] << endl;
                // }
                // file_4.close();

                // for (l = 0; l < 7; l++) {
                //     fname.pop_back();
                // }

                // fname_2.append(result);
                // fname_2.append("_1.txt");
                // file_4.open(fname_2);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_1[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 5; l++) {
                //     fname_2.pop_back();
                // }
                // fname_2.append("2.txt");
                // file_4.open(fname_2);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname_2.pop_back();
                // }
                
                // fname_3.append(result);
                // fname_3.append("_1.txt");
                // file_4.open(fname_3);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_thresh_1[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 5; l++) {
                //     fname_3.pop_back();
                // }
                // fname_3.append("2.txt");
                // file_4.open(fname_3);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_thresh_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname_3.pop_back();
                // }
                // convert.str("");
                // convert.clear();

                // reset values to zero for the next reconstruction
                for (l = 0; l < num_1; l++){
                    input_1[i] = 0;
                    input_2[i] = 0;
                    avg_1[i] = 0;
                    avg_2[i] = 0;
                    avg_thresh_1[i] = 0;
                    avg_thresh_2[i] = 0;
                }
            }

            else if (check_switch_1.size() == 3 && check_switch_1[2] - check_switch_1[0] > 0.01 * bin_size * 2){
                check_switch_1.erase(check_switch_1.begin());
            }

            else if (check_switch_2.size() == 5 && check_switch_2[4] - check_switch_2[0] == 0.01 * bin_size * 4){
                checkpoint = check_switch_2[0] - 0.01 * bin_size;
                dom_num ++;
                alter = alter + 1.0;
                cout << dom_num << endl;
                temp_4 = dom_num * 2 - 1;
                index[temp_4] = checkpoint;
                file_8 << checkpoint << endl;
                // start index of the next dominance
                checkpoint = check_switch_2[0] + 0.01 * bin_size;
                index[temp_4 + 1] = checkpoint;
                file_8 << checkpoint << endl;
                re1 = - re1;
                file_9 << re1 << endl;
                check_switch_2.clear();
                // save the data for the previous dominance
                // calculate the runtime of the previous dominance
                temp_5 = index[temp_4] - index[temp_4 - 1];
                // cout << temp_4 << endl;
                cout << temp_5 << endl;
                temp_5 = double(temp_5);
                // file_11 << temp_5 << endl;
                for (l = 0; l < num_1; l++){
                    input_1[l] = input_1[l] / temp_5;
                    input_2[l] = input_2[l] / temp_5;
                    avg_1[l] = avg_1[l] / temp_5 / 100.0;
                    avg_2[l] = avg_2[l] / temp_5 / 100.0;
                    avg_thresh_1[l] = avg_thresh_1[l] / temp_5 / 100.0;
                    avg_thresh_2[l] = avg_thresh_2[l] / temp_5 / 100.0;
                }

                // steps to create files and save information for reconstruction

                // convert2 << dom_num;
                // result = convert2.str();
                // fname.append(result);
                // fname.append("_1.txt");
                // file_4.open(fname);
                // for (l = 0; l < num_1; l++){
                //     file_4 << input_1[l] << endl;
                // }
                // file_4.close();

                // for (l = 0; l < 5; l++) {
                //     fname.pop_back();
                // }
                // fname.append("2.txt");

                // file_4.open(fname);
                // for (l = 0; l < num_1; l++){
                //     file_4 << input_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname.pop_back();
                // }

                // fname_2.append(result);
                // fname_2.append("_1.txt");
                // file_4.open(fname_2);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_1[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 5; l++) {
                //     fname_2.pop_back();
                // }
                // fname_2.append("2.txt");
                // file_4.open(fname_2);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname_2.pop_back();
                // }
                
                // fname_3.append(result);
                // fname_3.append("_1.txt");
                // file_4.open(fname_3);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_thresh_1[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 5; l++) {
                //     fname_3.pop_back();
                // }
                // fname_3.append("2.txt");
                // file_4.open(fname_3);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_thresh_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname_3.pop_back();
                // }
                // convert2.str("");
                // convert2.clear();

                // reset values to zero for the next reconstruction
                for (l = 0; l < num_1; l++){
                    input_1[i] = 0;
                    input_2[i] = 0;
                    avg_1[i] = 0;
                    avg_2[i] = 0;
                    avg_thresh_1[i] = 0;
                    avg_thresh_2[i] = 0;
                }
            }

            else if (check_switch_2.size() == 3 && check_switch_2[2] - check_switch_2[0] > 0.01 * bin_size * 2){
                check_switch_2.erase(check_switch_2.begin());
            }

            // check when last reconstruction ends
            // look at the last 10 data points
            // assume that at least one data point within the last 5 datapoints is within the previous dominance
            for (l = 9; l > - 1; l--){
                if (time_step == run_time - l * 0.01 * run_time){
                    if (re1 > 0.0){
                        if (metric >= 0.3 && time_step > temp_3){
                            temp_3 = time_step;
                        }
                    }
                    else if (re1 < 0.0) {
                        if (metric <= - 0.3 && time_step > temp_3){
                            temp_3 = time_step;
                        }
                    }
                }
            }

            if (time_step == run_time){
                file_8 << temp_3 << endl;
                dom_num ++;
                cout << dom_num << endl;
                temp_4 = dom_num * 2 - 1;
                index[temp_4] = temp_3;

                // calculate the runtime of the previous dominance
                temp_5 = index[temp_4] - index[temp_4 - 1];
                // cout << temp_4 << endl;
                cout << temp_5 << endl;
                temp_5 = double(temp_5);
                // if (re1 > 0){
                //     file_10 << temp_5 << endl;
                // }
                // else {
                //     file_11 << temp_5 << endl;
                // }
                for (l = 0; l < num_1; l++){
                    input_1[l] = input_1[l] / temp_5;
                    input_2[l] = input_2[l] / temp_5;
                    avg_1[l] = avg_1[l] / temp_5 / 100.0;
                    avg_2[l] = avg_2[l] / temp_5 / 100.0;
                    avg_thresh_1[l] = avg_thresh_1[l] / temp_5 / 100.0;
                    avg_thresh_2[l] = avg_thresh_2[l] / temp_5 / 100.0;
                }

                // save data for the last dominance

                // convert3 << dom_num;
                // result = convert3.str();
                // fname.append(result);
                // fname.append("_1.txt");
                // file_4.open(fname);
                // for (l = 0; l < num_1; l++){
                //     file_4 << input_1[l] << endl;
                // }
                // file_4.close();

                // for (l = 0; l < 5; l++) {
                //     fname.pop_back();
                // }
                // fname.append("2.txt");

                // file_4.open(fname);
                // for (l = 0; l < num_1; l++){
                //     file_4 << input_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname.pop_back();
                // }

                // fname_2.append(result);
                // fname_2.append("_1.txt");
                // file_4.open(fname_2);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_1[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 5; l++) {
                //     fname_2.pop_back();
                // }
                // fname_2.append("2.txt");
                // file_4.open(fname_2);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname_2.pop_back();
                // }
                
                // fname_3.append(result);
                // fname_3.append("_1.txt");
                // file_4.open(fname_3);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_thresh_1[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 5; l++) {
                //     fname_3.pop_back();
                // }
                // fname_3.append("2.txt");
                // file_4.open(fname_3);
                // for (l = 0; l < num_1; l++){
                //     file_4 << avg_thresh_2[l] << endl;
                // }
                // file_4.close();
                // for (l = 0; l < 7; l++) {
                //     fname_3.pop_back();
                // }
                // convert3.str("");
                // convert3.clear();

            }

        }

    }

    file.close();
    file_1.close();
    file_2.close();
    file_3.close();
    file_8.close();
    file_9.close();
    // file_10.close();
    // file_11.close();
    file_12.close();
    file_13.close();
    file_14.close();

    // delete arrays
    delete [] firing_num_1;
    delete [] firing_num_2;
    delete [] input_1;
    delete [] input_2;
    delete [] avg_1;
    delete [] avg_2;
    delete [] fire_1;
    delete [] fire_2;
    delete [] p1;
    delete [] p2;
    delete [] image_1;
    delete [] image_2;
    delete [] start_1;
    delete [] start_2;
    delete [] end_1;
    delete [] end_2;
    delete [] k;
    delete [] k_1;
    delete [] k_2;
    delete [] threshold_1;
    delete [] threshold_2;
    delete [] avg_thresh_1;
    delete [] avg_thresh_2;
    delete [] index;

    for (i = 0; i < num_1; i++){
        delete [] matrix_1[i];
        delete [] matrix_2[i];
    }

    for (i = 0; i < num_1; i++){
        delete [] matrix_3[i];
        delete [] matrix_4[i];
    }

    for (i = 0; i < height; i++){
        delete [] matrix_f_1[i];
        delete [] matrix_f_2[i];
    }

    return 0;

}

// notes about binocular rivalry

// binocular rivalry, two networks
// metric = (mu_1 - mu_2) / (mu_1 + mu_2)
// 1 time unit is 20 millisecond; 50 time units is 1 second
// mean activity: number of neurons fired divided by the population

// expects little volatility with no adaptation
// reconstruct the image where firing rate saturates
// weak adap: takes time for switching to happen
// max matric diff(phi, lambda)
// test values of phi and lambda so that switching is clear and periodic

// Levelt's law 3
// fixed pool 1 strength, increase pool 2 strength, plot pool 2 strength v.s. alternation rate
// alternation rate is number of dominance durations divided by time

// vary the strenth of pool 1 and 2 by varying m_0
// Levelt's law 2
// increase strength of pool 1 (m_0), plot pool 1's average dominance duration, plot pool 2's average dominance duration (separately)

// law 4
// have the same strength for pool 1 and 2, vary that strength, plot strength v.s. average dominance duration (or alternation rate)
// law 1
// increase the strength of pool 1, plot strength v.s. percent of time pool 1 is dominant (predominance)
