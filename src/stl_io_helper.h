#pragma once
#include <array>
#include <map>
#include <ostream>
#include <set>
#include <vector>
#include <utility>

//! Print a set of objects
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T> &set_objs);

//! Print a vector of objects
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &vector_objs);

//! Print a pair of objects
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2> &pair_objs);

//! Print an array containing N objects
template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N> &arr_objs);

//! Print a map
template <typename Tkey, typename Tval>
std::ostream& operator<<(std::ostream& os, const std::map<Tkey, Tval> &map_objs);

//! Print a the keys of a nested map
template <typename Tkey1, typename Tkey2, typename Tval>
void print_keys(std::ostream& os, const std::map<Tkey1, std::map<Tkey2, Tval>> &nested_map);

//! Get total number of elements
template <typename Tkey1, typename Tkey2, typename Tkey3, typename Tval>
void print_keys(std::ostream& os, const std::map<Tkey1, std::map<Tkey2, std::map<Tkey3, Tval>>> &nested_map);

//! Get total number of elements
template <typename Tkey1, typename Tkey2, typename Tval>
int get_num_keys(const std::map<Tkey1, std::map<Tkey2, Tval>> &nested_map);

//! Get total number of elements
template <typename Tkey1, typename Tkey2, typename Tkey3, typename Tval>
int get_num_keys(const std::map<Tkey1, std::map<Tkey2, std::map<Tkey3, Tval>>> &nested_map);

#include "stl_io_helper.hpp"  // IWYU pragma: export
