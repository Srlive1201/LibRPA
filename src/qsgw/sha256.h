#pragma once

#include <string>
#include <string_view>

namespace librpa_int
{
namespace qsgw
{

std::string sha256_string(std::string_view input);
std::string sha256_file(const std::string& path);
bool is_sha256_hex(std::string_view digest) noexcept;

} // namespace qsgw
} // namespace librpa_int
