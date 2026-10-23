#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace librpa_int
{
//! Streaming, uncompressed ZIP64/NumPy archive; no Python runtime or object arrays.
class NpzWriter
{
public:
    explicit NpzWriter(const std::string& filename);
    ~NpzWriter();
    void array(const std::string& name, const std::string& dtype,
               const std::vector<std::size_t>& shape, const void* data, std::size_t item_size);
    void text(const std::string& name, const std::string& value);
    void finish();

private:
    struct Entry
    {
        std::string name;
        std::uint64_t size, offset;
        std::uint32_t crc;
    };
    std::string filename_, temporary_;
    std::ofstream stream_;
    std::vector<Entry> entries_;
    bool finished_ = false;
    void integer(std::uint64_t value, int bytes);
};
}  // namespace librpa_int
