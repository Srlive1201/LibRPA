#include <cassert>

#include "../src/io/fs.h"
#include "../src/io/stl_io_helper.h"

int main (int argc, char *argv[])
{
    using namespace librpa_int;

    int myid = 0;
    create_directories("librpa.d", myid);

    std::map<int, std::map<int, std::map<int, double>>> nested_map
    {
        {0, {
                {0, {
                        {1, 1.0},
                        {2, 2.0},
                    }
                },
                {1, {
                        {2, 3.0},
                    }
                },
            }
        },
        {1, {
                {4, {
                        {0, -1.0},
                        {1, 1.0},
                        {2, 2.0},
                    }
                },
            }
        },
    };
    assert(get_num_keys(nested_map) == 6);

    return 0;
}
