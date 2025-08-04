#ifndef IT_FAST_DEFINITIONS_HPP
#define IT_FAST_DEFINITIONS_HPP

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include "definitions.hpp"


typedef std::list<Vertex> path;
typedef std::pair<const path, float> path_info; // pair of the path and its weight

inline std::string path_to_string(const path& p) {
    std::ostringstream oss;
    for (const auto& v : p) {
        oss << v << " ";
    }
    return oss.str();
}

// check two paths are the same
inline bool are_paths_equal(const path& p1, const path& p2) {
    if (p1.size() != p2.size()) return false;
    auto it1 = p1.begin();
    auto it2 = p2.begin();
    while (it1 != p1.end() && it2 != p2.end()) {
        if (*it1 != *it2) return false;
        ++it1;
        ++it2;
    }
    return true;
}

inline path_info* deep_copy_path_info(const path_info* original) {
    // Create a new path_info object and copy the contents
    path_info* new_path_info = new path_info(original->first, original->second);
    return new_path_info;
}

inline std::vector< std::vector< std::vector< path_info* > > > deep_copy_RT(const std::vector< std::vector< std::vector< path_info* > > >& original) {
    std::vector< std::vector< std::vector< path_info* > > > copy;
    // set the size of the copy to match the original
    copy.resize(original.size());

    for (size_t i = 0; i < original.size(); ++i) {
        copy[i].resize(original[i].size());
        for (size_t j = 0; j < original[i].size(); ++j) {
            for (const auto& path_ptr : original[i][j]) {
                copy[i][j].push_back(deep_copy_path_info(path_ptr));
            }
        }
    }

    return copy;
}

inline void set_RT_as_RT(std::vector<std::vector<std::vector<path_info*>>>& target_RT, const std::vector<std::vector<std::vector<path_info*>>>& initial_RT) {
    // set the values in path_info pointers in target_RT to be the same as in initial_RT
    // the two RT should have the same size and contain the same paths, only the float values in path_info pointers should be different
    for (size_t i = 0; i < target_RT.size(); ++i) {
        for (size_t j = 0; j < target_RT[i].size(); ++j) {
            for (size_t k = 0; k < target_RT[i][j].size(); ++k) {
                target_RT[i][j][k]->second = initial_RT[i][j][k]->second; // copy the float value
                assert(are_paths_equal(target_RT[i][j][k]->first, initial_RT[i][j][k]->first) && 
                       "Paths in target_RT and initial_RT should be the same.");
            }
        }
    }
}

using namespace boost::multi_index;
struct by_key {};
struct by_value {};

struct Entry {
    // The key is a pointer to path_info, which contains the path and its weight
    path_info* key;
    // The value is a float representing the flow on this path
    float value;
};

typedef multi_index_container<
    Entry,
    indexed_by<
        hashed_unique<
            tag<by_key>,
            member<Entry, path_info*, &Entry::key>
        >,
        ordered_non_unique<
            tag<by_value>,
            member<Entry, float, &Entry::value>,
            std::greater<float>  // Add this for descending order
        >
    >
> descend_loadcontri_pathinfos;

#endif