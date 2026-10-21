#include "npz_writer.h"

#include <array>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace librpa_int
{
namespace
{
std::uint32_t update_crc(std::uint32_t crc, const void* data, std::size_t size)
{
    static const auto table = []
    {
        std::array<std::uint32_t, 256> values{};
        for (std::uint32_t i = 0; i < values.size(); ++i)
        {
            auto value = i;
            for (int bit = 0; bit < 8; ++bit)
                value = (value >> 1) ^ ((value & 1) ? 0xedb88320u : 0u);
            values[i] = value;
        }
        return values;
    }();
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) crc = table[(crc ^ bytes[i]) & 255] ^ (crc >> 8);
    return crc;
}
}  // namespace

NpzWriter::NpzWriter(const std::string& filename)
    : filename_(filename), temporary_(filename + ".tmp"), stream_(temporary_, std::ios::binary)
{
    const std::uint32_t endian = 1;
    if (*reinterpret_cast<const unsigned char*>(&endian) != 1)
        throw std::runtime_error("NumPy export requires little-endian floating-point storage");
    if (!stream_) throw std::runtime_error("cannot create model archive " + temporary_);
}

NpzWriter::~NpzWriter()
{
    if (!finished_)
    {
        stream_.close();
        std::error_code error;
        std::filesystem::remove(temporary_, error);
    }
}

void NpzWriter::integer(std::uint64_t value, int bytes)
{
    for (int i = 0; i < bytes; ++i) stream_.put(static_cast<char>((value >> (8 * i)) & 255));
}

void NpzWriter::array(const std::string& name, const std::string& dtype,
                      const std::vector<std::size_t>& shape, const void* data,
                      std::size_t item_size)
{
    if (finished_) throw std::runtime_error("model archive is already closed");
    std::size_t count = 1;
    for (auto extent : shape)
    {
        if (extent && count > std::numeric_limits<std::size_t>::max() / extent)
            throw std::runtime_error("model array size overflow");
        count *= extent;
    }
    if (item_size == 0 || count > std::numeric_limits<std::size_t>::max() / item_size ||
        (count && !data))
        throw std::runtime_error("invalid model array");
    std::ostringstream dictionary;
    dictionary << "{'descr': '" << dtype << "', 'fortran_order': False, 'shape': (";
    for (auto extent : shape) dictionary << extent << ", ";
    dictionary << "), }";
    auto header = dictionary.str();
    header.append((64 - ((10 + header.size() + 1) % 64)) % 64, ' ');
    header += '\n';
    if (header.size() > 65535) throw std::runtime_error("NumPy array header is too long");
    std::string prefix("\x93NUMPY", 6);
    prefix += '\1';
    prefix += '\0';
    prefix += static_cast<char>(header.size() & 255);
    prefix += static_cast<char>((header.size() >> 8) & 255);
    prefix += header;
    const auto bytes = count * item_size;
    const auto crc =
        update_crc(update_crc(0xffffffffu, prefix.data(), prefix.size()), data, bytes) ^
        0xffffffffu;
    Entry entry{name + ".npy", prefix.size() + bytes, static_cast<std::uint64_t>(stream_.tellp()),
                crc};
    if (entry.name.size() > 65535) throw std::runtime_error("NumPy array name is too long");
    // ZIP64 local header, followed by an ordinary NPY v1 payload.
    integer(0x04034b50, 4);
    integer(45, 2);
    integer(0, 2);
    integer(0, 2);
    integer(0, 2);
    integer(0, 2);
    integer(crc, 4);
    integer(0xffffffffu, 4);
    integer(0xffffffffu, 4);
    integer(entry.name.size(), 2);
    integer(20, 2);
    stream_.write(entry.name.data(), entry.name.size());
    integer(1, 2);
    integer(16, 2);
    integer(entry.size, 8);
    integer(entry.size, 8);
    stream_.write(prefix.data(), prefix.size());
    if (bytes) stream_.write(static_cast<const char*>(data), bytes);
    if (!stream_) throw std::runtime_error("failed writing model array " + name);
    entries_.push_back(std::move(entry));
}

void NpzWriter::text(const std::string& name, const std::string& value)
{
    array(name, "|u1", {value.size()}, value.data(), 1);
}

void NpzWriter::finish()
{
    const auto start = static_cast<std::uint64_t>(stream_.tellp());
    for (const auto& entry : entries_)
    {
        integer(0x02014b50, 4);
        integer(45, 2);
        integer(45, 2);
        integer(0, 2);
        integer(0, 2);
        integer(0, 2);
        integer(0, 2);
        integer(entry.crc, 4);
        integer(0xffffffffu, 4);
        integer(0xffffffffu, 4);
        integer(entry.name.size(), 2);
        integer(28, 2);
        integer(0, 2);
        integer(0, 2);
        integer(0, 2);
        integer(0, 4);
        integer(0xffffffffu, 4);
        stream_.write(entry.name.data(), entry.name.size());
        integer(1, 2);
        integer(24, 2);
        integer(entry.size, 8);
        integer(entry.size, 8);
        integer(entry.offset, 8);
    }
    const auto end = static_cast<std::uint64_t>(stream_.tellp());
    integer(0x06064b50, 4);
    integer(44, 8);
    integer(45, 2);
    integer(45, 2);
    integer(0, 4);
    integer(0, 4);
    integer(entries_.size(), 8);
    integer(entries_.size(), 8);
    integer(end - start, 8);
    integer(start, 8);
    integer(0x07064b50, 4);
    integer(0, 4);
    integer(end, 8);
    integer(1, 4);
    integer(0x06054b50, 4);
    integer(0, 2);
    integer(0, 2);
    integer(65535, 2);
    integer(65535, 2);
    integer(0xffffffffu, 4);
    integer(0xffffffffu, 4);
    integer(0, 2);
    stream_.close();
    if (!stream_) throw std::runtime_error("failed closing model archive " + temporary_);
    std::filesystem::rename(temporary_, filename_);
    finished_ = true;
}
}  // namespace librpa_int
