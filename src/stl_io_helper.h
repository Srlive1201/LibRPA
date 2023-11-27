#pragma once
#include <array>
#include <map>
#include <ostream>
#include <set>
#include <vector>

//! Print a set of objects
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T> &set_objs)
{
    os << "[";
    if (!set_objs.empty())
    {
        auto pobj = set_objs.cbegin();
        for (; pobj != std::prev(set_objs.cend()); pobj++)
            os << *pobj << ",";
        os << *pobj;
    }
    os << "](s" << set_objs.size() << ")";
    return os;
}

//! Print a vector of objects
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &vector_objs)
{
    os << "[";
    if (!vector_objs.empty())
    {
        auto pobj = vector_objs.cbegin();
        for (; pobj != std::prev(vector_objs.cend()); pobj++)
            os << *pobj << ",";
        os << *pobj;
    }
    os << "](v" << vector_objs.size() << ")";
    return os;
}

//! Print a pair of objects
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2> &pair_objs)
{
    os << "{" << pair_objs.first << "," << pair_objs.second << "}";
    return os;
}

//! Print an array containing N objects
template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N> &arr_objs)
{
    os << "[";
    for (int i = 0; i < N - 1; i++)
        os << arr_objs[i] << ",";
    os << arr_objs[N-1] << "](a" << N << ")";
    return os;
}

//! Print a map
template <typename Tkey, typename Tval>
std::ostream& operator<<(std::ostream& os, const std::map<Tkey, Tval> &map_objs)
{
    os << "{";
    for (const auto& kv: map_objs)
    {
        os << "\"" << kv.first << "\":\"" << kv.second << "\"\n";
    }
    os << "}";
    return os;
}

//! Print a the keys of a nested map
template <typename Tkey1, typename Tkey2, typename Tval>
void print_keys(std::ostream& os, const std::map<Tkey1, std::map<Tkey2, Tval>> &nested_map)
{
    for (const auto& k1k2v: nested_map)
        for (const auto& k2v: k1k2v.second)
            os << k1k2v.first << " " << k2v.first << "\n";
}

//! Get total number of elements
template <typename Tkey1, typename Tkey2, typename Tkey3, typename Tval>
void print_keys(std::ostream& os, const std::map<Tkey1, std::map<Tkey2, std::map<Tkey3, Tval>>> &nested_map)
{
    for (const auto& k1k2k3v: nested_map)
        for (const auto& k2k3v: k1k2k3v.second)
            for (const auto& k3v: k2k3v.second)
            os << k1k2k3v.first << " " << k2k3v.first << " " << k3v.first << "\n";
}

//! Get total number of elements
template <typename Tkey1, typename Tkey2, typename Tval>
int get_num_keys(const std::map<Tkey1, std::map<Tkey2, Tval>> &nested_map)
{
    int nkeys = 0;
    for (const auto& k1k2v: nested_map)
        for (const auto& k2v: k1k2v.second)
            nkeys++;
    return nkeys;
}

//! Get total number of elements
template <typename Tkey1, typename Tkey2, typename Tkey3, typename Tval>
int get_num_keys(const std::map<Tkey1, std::map<Tkey2, std::map<Tkey3, Tval>>> &nested_map)
{
    int nkeys = 0;
    for (const auto& k1k2k3v: nested_map)
        for (const auto& k2k3v: k1k2k3v.second)
            for (const auto& k3v: k2k3v.second)
                nkeys++;
    return nkeys;
}
