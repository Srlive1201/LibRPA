#include "../src/io/fs.h"

int main (int argc, char *argv[])
{
    using namespace librpa_int;

    int myid = 0;
    create_directories("librpa.d", myid);

    return 0;
}
