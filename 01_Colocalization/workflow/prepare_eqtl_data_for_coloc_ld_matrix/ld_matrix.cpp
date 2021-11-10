#include <iostream>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>


int main(int argc, char** argv) {

    // Ensure correct arguments
    if (argc != 2) {

        std::cerr << "The program must be run with the following syntax: ld_matrix <snps_file>" << std::endl;
        std::cerr << std::endl;
        std::cerr << "\tsnps_file - The list of snps being tested for LD. Should be a tab-delimited file containing chromosome and position." << std::endl;
        std::cerr << std::endl;

        return EXIT_FAILURE;
    }

    // Open SNPs file and create a mapping from SNP position to array index
    std::ifstream snps_file (argv[1]);

    if (!snps_file.is_open()) {

        std::cerr << "Unable to open the SNPs file" << std::endl;

        return EXIT_FAILURE;
    }

    std::unordered_map<int, size_t> position_map;
    std::unordered_map<size_t, int> index_map;
    std::string snp_line;
    size_t i = 0;
    while (getline(snps_file, snp_line)) {

        size_t tab_pos = snp_line.find('\t');
        std::string position = snp_line.substr(tab_pos + 1);

        if (position != "") {
            int position_int = std::stoi(position);
            position_map[position_int] = i;
            index_map[i] = position_int;
            i += 1;
        }
    }

    snps_file.close();

    // Allocate memory for SNP LD Matrix
    double** ld_matrix = new double*[position_map.size()];
    if (ld_matrix == NULL) {
        std::cerr << "Unable to allocate memory for LD matrix" << std::endl;
    }
    for (size_t i = 0; i < position_map.size(); i++) {
        ld_matrix[i] = new double[position_map.size()];
        if (ld_matrix[i] == NULL) {
            std::cerr << "Unable to allocate memory for LD matrix" << std::endl;
        }
    }

    // Fill LD Matrix with -1
    for (size_t i = 0; i < position_map.size(); i++) {
        for (size_t j = 0; j < position_map.size(); j++) {
            ld_matrix[i][j] = -1.0;
        }
    }

    // Fill diagonal with 1 (LD of SNP with self is trivially 1)
    for (size_t i = 0; i < position_map.size(); i++) {
        ld_matrix[i][i] = 1.0;
    }

    // Ignore first line
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Parse all the other lines
    std::string input_line;
    while (getline(std::cin, input_line)) {
    
        size_t tab_1_pos = input_line.find('\t');
        size_t tab_2_pos = input_line.find('\t', tab_1_pos + 1);
        size_t tab_3_pos = input_line.find('\t', tab_2_pos + 1);
        size_t tab_4_pos = input_line.find('\t', tab_3_pos + 1);
        std::string snp_1_pos = input_line.substr(tab_1_pos + 1, tab_2_pos - tab_1_pos - 1);
        std::string snp_2_pos = input_line.substr(tab_2_pos + 1, tab_3_pos - tab_2_pos - 1);
        std::string snps_r2 = input_line.substr(tab_4_pos + 1);

        int snp_1_pos_int = std::stoi(snp_1_pos);
        int snp_2_pos_int = std::stoi(snp_2_pos);
        if (position_map.count(snp_1_pos_int) > 0 && position_map.count(snp_2_pos_int) > 0) {
            size_t i = position_map[snp_1_pos_int];
            size_t j = position_map[snp_2_pos_int];
            double r2 = std::stod(snps_r2);
            ld_matrix[i][j] = r2;
            ld_matrix[j][i] = r2;
        }
    }

    // Print LD Matrix Header
    for (size_t i = 0; i < index_map.size(); i++) {
        std::cout << '\t' << index_map[i];
    }
    std::cout << '\n';

    // Print each line of LD matrix
    for (size_t i = 0; i < index_map.size(); i++) {
        std::cout << index_map[i];
        for (size_t j = 0; j < index_map.size(); j++) {
            std::cout << '\t' << ld_matrix[i][j];
        }
        std::cout << '\n';
    }

    // Deallocate memory
    for (size_t i = 0; i < position_map.size(); i++) {
        delete ld_matrix[i];
    }
    delete ld_matrix;

    return EXIT_SUCCESS;
}