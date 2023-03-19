#include <chrono>
#include <dirent.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_do.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_scheduler_init.h>
#include <tuple>
#include <vector>

using namespace std;

#define THREAD_COUNT_OPT 'T'
#define TAXONOMY_PATH_OPT 'A'
#define KMER_LENGTH 32

namespace TaxonomicInfo {
uint16_t num_levels = 7;
enum level { KINGDOM = 1, PHYLUM = 2, CLASS = 3, ORDER = 4, FAMILY = 5, GENUS = 6, SPECIES = 7 };

enum Kingdoms { BACTERIA = 2, ARCHAEA = 2157 };

static const Kingdoms AllKingdoms[] = {BACTERIA, ARCHAEA};
} // namespace TaxonomicInfo

struct kmer_match {
  uint16_t dist;
  float vote;
  uint64_t taxID;
};

struct read_info {
  bool isRC;
  string readID;
  tuple<uint64_t, float, float> pred_taxID_info;
  vector<kmer_match> match_vector;
};

vector<string> list_dir(const char *path) {
  vector<string> userString;
  struct dirent *entry;
  DIR *dir = opendir(path);

  while ((entry = readdir(dir)) != NULL) {
    if ((strcmp(entry->d_name, "..") != 0) && (strcmp(entry->d_name, ".") != 0)) {
      userString.push_back(string(path) + "/" + entry->d_name);
    }
  }

  closedir(dir);
  return (userString);
}

void read_taxonomy_lookup(string filepath, tbb::concurrent_unordered_map<uint64_t, vector<uint64_t>> &ancestor_map) {
  ifstream ftable;
  ftable.open(filepath);

  if (!ftable) {
    cout << "Cannot open file for LCA lookup table." << endl;
    exit(1);
  }

  for (string line; getline(ftable, line);) {
    istringstream iss(line);
    string taxIDstr;
    getline(iss, taxIDstr, ' ');
    uint64_t taxID = stoi(taxIDstr);

    vector<uint64_t> ancestors;
    string ancestorID;

    while (getline(iss, ancestorID, ',')) {
      ancestors.push_back(stoi(ancestorID));
    }
    ancestor_map.insert({taxID, ancestors});
  }
}

void read_matches(string filepath, vector<read_info> &all_read_info) {
  ifstream infile(filepath);
  string line;
  string readID;
  uint64_t line_counter = 0;

  uint16_t k = KMER_LENGTH;

  while (getline(infile, line)) {
    istringstream iss(line);

    if ((line_counter % 3) == 0)
      iss >> readID;
    else {
      bool isRC;
      string tmp;
      if ((line_counter % 3) == 1) {
        isRC = false;
        iss >> tmp;
      } else {
        isRC = true;
        iss >> tmp;
      }

      string match_str;
      vector<kmer_match> match_vector;

      while (iss >> match_str) {
        stringstream ss(match_str);
        string taxID_str;
        string dist_str;
        getline(ss, taxID_str, ':');
        getline(ss, dist_str, ':');

        kmer_match curr_match;
        curr_match.dist = stoi(dist_str);
        curr_match.taxID = stoi(taxID_str);
        /* curr_match.vote = pow((1.0 - curr_match.dist / (float)k), k); */
        match_vector.push_back(curr_match);
      }

      read_info curr_read;
      curr_read.readID = readID;
      curr_read.isRC = isRC;
      curr_read.match_vector = match_vector;
      all_read_info.push_back(curr_read);
    }
    line_counter++;
  }
}

tbb::concurrent_unordered_map<uint64_t, float>
get_level_votes(kmer_match amatch, tbb::concurrent_unordered_map<uint64_t, vector<uint64_t>> &ancestor_map) {
  uint16_t k = KMER_LENGTH;
  tbb::concurrent_unordered_map<uint64_t, float> match_votes;
  for (uint64_t a_taxID : ancestor_map[amatch.taxID]) {
    amatch.vote = pow((1.0 - (float)amatch.dist / (float)k), k);
    match_votes[a_taxID] = amatch.vote;
  }
  return match_votes;
}

void aggregate_votes(tbb::concurrent_unordered_map<uint64_t, vector<uint64_t>> ancestor_map,
                     vector<read_info> &all_read_info) {
  auto sum_func = [](const tbb::concurrent_vector<float> &vec, size_t begin, size_t end) {
    float sum = 0;
    for (size_t i = begin; i < end; ++i) {
      sum += vec[i];
    }
    return sum;
  };

  /* for (auto &curr_read : all_read_info) { */
  tbb::parallel_for(0, (int)all_read_info.size(), [&](int i) {
    auto &curr_read = all_read_info[i];
    pair<uint64_t, float> identity(0, 0.0);
    pair<uint64_t, float> maxID(0, 0.0);
    float max_vote = 0.0;
    if (curr_read.match_vector.size() > 0) {
      tbb::concurrent_unordered_map<uint64_t, tbb::concurrent_vector<float>> vote_collector;
      tbb::concurrent_unordered_map<uint64_t, float> final_votes;
      tbb::concurrent_unordered_set<uint64_t> all_taxIDs;
      tbb::concurrent_unordered_map<uint16_t, tbb::concurrent_unordered_set<uint64_t>> taxIDs_by_level;

      tbb::parallel_do(curr_read.match_vector.begin(), curr_read.match_vector.end(), [&](kmer_match &curr_match) {
        tbb::concurrent_unordered_map<uint64_t, float> match_votes = get_level_votes(curr_match, ancestor_map);
        for (auto &vote : match_votes) {
          vote_collector[vote.first].push_back(vote.second);
          all_taxIDs.insert(vote.first);
          taxIDs_by_level[(int)ancestor_map[vote.first].size()].insert(vote.first);
        }
      });

      tbb::concurrent_vector<uint64_t> taxIDs_vec(all_taxIDs.begin(), all_taxIDs.end());

      tbb::parallel_for(0, (int)taxIDs_vec.size(), [&](int i) {
        final_votes[taxIDs_vec[i]] = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, vote_collector[taxIDs_vec[i]].size()), 0.0,
            [&](const tbb::blocked_range<int> &range, float partial_sum) -> float {
              return partial_sum + sum_func(vote_collector[taxIDs_vec[i]], range.begin(), range.end());
            },
            [](float x, float y) -> float { return x + y; });
      });

      uint64_t rootID;
      for (const auto &taxon : TaxonomicInfo::AllKingdoms) {
        if (final_votes[static_cast<int>(taxon)] > max_vote) {
          max_vote = final_votes[static_cast<int>(taxon)];
          rootID = static_cast<int>(taxon);
        }
      }

      float th_vote = 0.5 * max_vote;

      for (uint16_t lvl = TaxonomicInfo::num_levels; lvl >= 1; --lvl) {
        tbb::concurrent_vector<uint64_t> taxIDs_vec_lvl(taxIDs_by_level[lvl].begin(), taxIDs_by_level[lvl].end());
        maxID = tbb::parallel_reduce(
            tbb::blocked_range<int>(0, taxIDs_vec_lvl.size()), identity,
            [&](const tbb::blocked_range<int> &r, pair<uint64_t, float> init) -> pair<uint64_t, float> {
              for (int i = r.begin(); i != r.end(); ++i) {
                float curr_vote = final_votes[taxIDs_vec_lvl[i]];
                if (init.second < curr_vote) {
                  init.first = taxIDs_vec_lvl[i];
                  init.second = curr_vote;
                }
              }
              return init;
            },
            [](pair<uint64_t, float> x, pair<uint64_t, float> y) -> pair<uint64_t, float> {
              if (x.second > y.second) {
                return x;
              } else {
                return y;
              }
            });
        if (maxID.second > th_vote) {
          break;
        }
      }
    }
    get<0>(curr_read.pred_taxID_info) = maxID.first;
    get<1>(curr_read.pred_taxID_info) = maxID.second;
    get<2>(curr_read.pred_taxID_info) = max_vote;
  });
  /* } */
}

void write_predictions_to_file(string filepath, vector<read_info> all_read_info) {
  ofstream outfile(filepath);
  string read_form;
  for (auto &curr_read : all_read_info) {
    if (curr_read.isRC) {
      read_form = "rc";
    } else {
      read_form = "--";
    }
    outfile << curr_read.readID << "\t" << read_form << "\t" << to_string(get<0>(curr_read.pred_taxID_info)) << "\t"
            << "\t" << to_string(get<1>(curr_read.pred_taxID_info)) << "\t"
            << to_string(get<2>(curr_read.pred_taxID_info)) << endl;
  }
  outfile.close();
}

int main(int argc, char *argv[]) {
  uint64_t thread_count = 1;
  char *input_matches_dir = NULL;
  string taxonomy_lookup_path;
  string output_predictions_dir = ".";

  int cf_tmp;
  opterr = 0;

  while (1) {
    static struct option long_options[] = {
        {"input-matches-dir", 1, 0, 'i'},
        {"output-predictions-dir", 1, 0, 'o'},
        {"taxonomy-path", 1, 0, TAXONOMY_PATH_OPT},
        {"thread-count", 1, 0, THREAD_COUNT_OPT},
        {0, 0, 0, 0},
    };

    int option_index = 0;
    cf_tmp = getopt_long(argc, argv, "i:o:", long_options, &option_index);

    if ((optarg != NULL) && (*optarg == '-')) {
      cf_tmp = ':';
    }

    if (cf_tmp == -1)
      break;
    else if (cf_tmp == TAXONOMY_PATH_OPT)
      taxonomy_lookup_path = optarg;
    else if (cf_tmp == THREAD_COUNT_OPT)
      thread_count = atoi(optarg); // Default is 1.
    else {
      switch (cf_tmp) {
      case 'i':
        input_matches_dir = optarg;
        break;
      case 'o':
        output_predictions_dir = optarg;
        break;
      case '?':
        if (optopt == 'i')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (optopt == 'q')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option '-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option '%s'.\n", argv[optind - 1]);
        return 1;
      default:
        abort();
      }
    }
  }
  tbb::task_scheduler_init init(thread_count);
  cout << "Number of threads is " << thread_count << "." << endl;

  struct stat s_input_path;
  vector<string> match_info_path_list;
  if (stat(input_matches_dir, &s_input_path) == 0) {
    if (s_input_path.st_mode & S_IFDIR) {
      match_info_path_list = list_dir(input_matches_dir);
    } else {
      cout << "Filetype in the given query path is not recognized." << endl;
      exit(1);
    }
  } else {
    cout << "Given query path is not valid!" << endl;
    exit(1);
  }

  int total_readTime = 0;
  int total_classificationTime = 0;
  int total_numberOfReads = 0;

  chrono::steady_clock::time_point t1;
  chrono::steady_clock::time_point t2;

  for (string input_path : match_info_path_list) {
    string query_name = input_path.substr(input_path.find_last_of("/") + 1);
    query_name = query_name.substr(0, query_name.find_last_of("."));
    query_name = query_name.substr(query_name.find_first_of("_") + 1);
    string output_path = output_predictions_dir + "/" + "predictions_" + query_name;

    tbb::concurrent_unordered_map<uint64_t, vector<uint64_t>> ancestor_map;
    read_taxonomy_lookup(taxonomy_lookup_path, ancestor_map);

    vector<read_info> all_read_info;
    t1 = chrono::steady_clock::now();
    read_matches(input_path, all_read_info);
    t2 = chrono::steady_clock::now();
    total_readTime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();

    total_numberOfReads += all_read_info.size();

    t1 = chrono::steady_clock::now();
    aggregate_votes(ancestor_map, all_read_info);
    t2 = chrono::steady_clock::now();
    total_classificationTime += chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();

    write_predictions_to_file(output_path, all_read_info);
  }
  cout << "Time past (read_matches) = " << total_readTime << "[ms]" << endl;
  cout << "Time past (aggregate_votes) = " << total_classificationTime << "[ms]" << endl;
  cout << "Total number of reads: " << total_numberOfReads << endl;
}
