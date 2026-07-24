#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/sha256.h"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using librpa_int::qsgw::is_sha256_hex;
using librpa_int::qsgw::sha256_file;
using librpa_int::qsgw::sha256_string;

namespace
{

template <typename Function>
void assert_throws(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    assert(threw);
}

void test_known_vectors()
{
    assert(sha256_string("") ==
           "e3b0c44298fc1c149afbf4c8996fb924"
           "27ae41e4649b934ca495991b7852b855");
    assert(sha256_string("abc") ==
           "ba7816bf8f01cfea414140de5dae2223"
           "b00361a396177a9cb410ff61f20015ad");
}

void test_file_hash_and_validation()
{
    const std::string path = "test_qsgw_sha256_input.tmp";
    {
        std::ofstream output(path, std::ios::binary);
        output << "abc";
    }
    assert(sha256_file(path) == sha256_string("abc"));
    std::remove(path.c_str());

    assert(is_sha256_hex(
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad"));
    assert(!is_sha256_hex("BA7816BF8F01CFEA414140DE5DAE2223"
                          "B00361A396177A9CB410FF61F20015AD"));
    assert(!is_sha256_hex("abc"));
    assert_throws([&] { (void)sha256_file("missing-qsgw-input.dat"); });
}

} // namespace

int main()
{
    test_known_vectors();
    test_file_hash_and_validation();
    std::cout << "test_qsgw_sha256: all tests passed\n";
    return 0;
}
