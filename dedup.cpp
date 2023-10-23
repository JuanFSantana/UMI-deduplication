#include <iostream>
#include <filesystem>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cstdlib>
#include <thread>

void process_zip_file(std::map<std::string, std::string>& r2_paths, std::filesystem::path& directory_dedup, std::map<std::string, std::string>& read_umi_map);
void dedup(std::map<std::string, std::string>& bed_paths, std::filesystem::path& directory_dedup, std::map<std::string, std::string>& read_umi_map);

int main(int argc, char* argv[]) {
    // make sure that 6 arguments are passed
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <path bed file> <path to samplekey.csv>" << std::endl << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    // Convert argv[] to a std::filesystem::path object
    std::filesystem::path path_bed(argv[1]);
    std::filesystem::path samplesheet(argv[2]);
    std::filesystem::path directory_dedup = path_bed.parent_path();

    // Check if the path exists
    if (!std::filesystem::exists(path_bed) || path_bed.extension() != ".bed" || !std::filesystem::is_regular_file(path_bed)) {
        std::cerr << "The path to the bed file is invalid." << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    if (!std::filesystem::exists(samplesheet) || samplesheet.extension() != ".csv" || samplesheet.filename() != "sampleskey.csv" || !std::filesystem::is_regular_file(samplesheet)) {
        std::cerr << "The path to the samplesheet file is invalid." << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
    if (!std::filesystem::exists(directory_dedup)){
        std::filesystem::create_directories(directory_dedup);
    }

    // create a map to store bed file name and path
    std::map<std::string, std::string> bed_paths;
    std::string bed_name = path_bed.stem().string();
    bed_paths[bed_name] = path_bed.string(); 

    // read in samplekey.csv
    std::ifstream samplesheet_file(samplesheet);
    if (!samplesheet_file.is_open()) {
        std::cerr << "Failed to open the input file." << std::endl;
        std::exit(1);
    }
    // read in the header: name is the same as the name of the bed files (bed files have extension.bed)
    std::map<std::string, std::string> r2_paths;
    std::string header;
    if (std::getline(samplesheet_file, header)) {
        // Header read and discarded
    } else {
        std::cerr << "Failed to read header line." << std::endl;
        return 1;
    }
    // Variables to store the data
    std::string name, fastq_1, fastq_2, fastq_3;
    // Loop through the data lines
    while (std::getline(samplesheet_file, name, ',') &&
           std::getline(samplesheet_file, fastq_1, ',') &&
           std::getline(samplesheet_file, fastq_2, ',') &&
           std::getline(samplesheet_file, fastq_3)) {
            // paths for R2 files
            if (name == bed_name){
                std::string r2_path = std::filesystem::current_path().string() + "/" + fastq_2;
                // add to map
                r2_paths[name] = r2_path;
            }
        } 
    // close samplesheet file
    samplesheet_file.close();

    // map to store read name and umi barcode
    std::map<std::string, std::string> read_umi_map;
    // process zip files
    process_zip_file(r2_paths, directory_dedup, read_umi_map);
    // dedup
    dedup(bed_paths, directory_dedup, read_umi_map);

    return 0;
}


void process_zip_file(std::map<std::string, std::string>& r2_paths, std::filesystem::path& directory_dedup, std::map<std::string, std::string>& read_umi_map) {
    std::string name = r2_paths.begin()->first;
    std::string r2_zipped_path = r2_paths.begin()->second;
   
    std::string output_unzipped_r2 = directory_dedup.string() + "/" + name + ".txt";

    // Construct the unzip command using the gzip utility
    int ncores = std::thread::hardware_concurrency();
    std::string unzip_command = "unpigz -p " + std::to_string(ncores) + " -c " + r2_zipped_path + " > " + output_unzipped_r2;

    // Execute the unzip command using the system function
    int result = std::system(unzip_command.c_str());

    if (result != 0) {
        std::cerr << "Failed to unzip the file: " << r2_zipped_path << std::endl;
        std::exit(1);
    }

    // Open the unzipped file for reading
    std::ifstream unzipped_file(output_unzipped_r2);

    if (!unzipped_file.is_open()) {
        std::cerr << "Failed to open the unzipped file: " << output_unzipped_r2 << std::endl;
        std::exit(1);
    }

    std::vector<std::string> read_info; // Store read information
    std::string line; 
    while (std::getline(unzipped_file, line)) {
        read_info.push_back(line); // Store the line

        // If we have collected 4 lines, process the read information
        if (read_info.size() == 4) {
            // Split the first line at the space
            std::istringstream iss(read_info[0]);
            std::string read_name;
            iss >> read_name;
            // Remove the "@" symbol from read_name
            if (!read_name.empty() && read_name[0] == '@') {
                read_name = read_name.substr(1); 
            }

            // Extract other information from the vector
            std::string umi_barcode = read_info[1];
            std::string strand = read_info[2];
            std::string other_info = read_info[3];

            // Store the read name and UMI barcode in the map
            read_umi_map[read_name] = umi_barcode;

            // Clear the vector for the next read
            read_info.clear();
        }

    }

    unzipped_file.close();
    // Remove the zipped files
    std::filesystem::remove(output_unzipped_r2);

}

void dedup (std::map<std::string, std::string>& bed_paths, std::filesystem::path& directory_dedup, std::map<std::string, std::string>& read_umi_map) {
    // sort bed file based on chromosome,start,end,id
    std::string bed_name = bed_paths.begin()->first;
    std::string bed_path = bed_paths.begin()->second;
    int ncores = std::thread::hardware_concurrency();

    std::filesystem::path sorted_bed_path = bed_path;
    sorted_bed_path.replace_extension("-sorted.bed");
    std::string sort_command = "sort --parallel=" + std::to_string(ncores) + " -k1,1 -k2,2n -k3,3n -k4,4 " + bed_path + " > " + sorted_bed_path.string();
    int result = std::system(sort_command.c_str());

    if (result != 0) {
        std::cerr << "Failed to sort the file: " << bed_path << std::endl;
        std::exit(1);
    }

    std::ifstream sorted_file(sorted_bed_path);
    if (!sorted_file.is_open()) {
        std::cerr << "Failed to open the sorted file: " << sorted_bed_path << std::endl;
        std::exit(1);
    }

    // loop through the sorted bed file
    std::string below_chromosome, below_start, below_end, below_id, below_quality, below_strand;
    size_t total_reads = 0;
    size_t total_duplicates = 0;
    std::string chromosome, start, end, id, quality, strand;
    std::map<std::string, int> umi_dulplicate_map;
    while (sorted_file >> chromosome >> start >> end >> id >> quality >> strand) {
        std::string above_chromosome = chromosome;
        std::string above_start = start;
        std::string above_end = end;
        std::string above_id = id; 
        std::string above_quality = quality; 
        std::string above_strand = strand;
       
        if(total_reads > 0){
            if((above_chromosome == below_chromosome) && (above_start == below_start) && (above_end == below_end)){
                // copmaring umi barcodes between reads
                std::string above_umi_barcodes = read_umi_map[above_id];
                std::string below_umi_barcodes = read_umi_map[below_id];

                if(above_umi_barcodes.compare(below_umi_barcodes) == 0){ // if umi barcodes are the same
                    umi_dulplicate_map[below_umi_barcodes]++;
                    total_duplicates++;
                }
                else{
                    std::cout << above_chromosome << '\t' << above_start << '\t' << above_end << '\t' << above_id << '\t' << above_strand << '\n';
                }
            }  
            else{
                std::cout << above_chromosome << '\t' << above_start << '\t' << above_end << '\t' << above_id << '\t' << above_strand << '\n';
            }
        }
        else{
            std::cout << above_chromosome << '\t' << above_start << '\t' << above_end << '\t' << above_id << '\t' << above_strand << '\n';
        }

        below_chromosome = above_chromosome;
        below_start = above_start;
        below_end = above_end;
        below_id = above_id;
        below_quality = above_quality;
        below_strand = above_strand;

        total_reads++;
    }
    // calculate percentage of reads removed
    double percentage;
    if(total_duplicates != 0){
        percentage = ((double)total_duplicates / (double)total_reads) * 100;
    }
    else{
        percentage = 0;
    }

    // count times of duplication 
    std::map<int, int> duplicate_map;
    for (auto const& x : umi_dulplicate_map){
        duplicate_map[x.second]++;
    }
    std::ofstream out_file_log(directory_dedup.string() + "/" + bed_name + "-dedup_log.txt");
    out_file_log << "Total reads" << '\t' << "Total duplicates" << '\t' << "Percentage of reads removed" << '\n';
    out_file_log << total_reads << '\t' << total_duplicates << '\t' << percentage << '\n' << '\n';
    out_file_log << "Number of duplications for a given UMI" << '\t' << "Number of UMIs" << '\t' << "Total" << '\n';
    for (auto const& x : duplicate_map){
        out_file_log << x.first << '\t' << x.second << '\t' << x.first*x.second << '\n';
    }

    out_file_log.close();
    sorted_file.close();
    // remove sorted bed file
    std::filesystem::remove(sorted_bed_path);
}


