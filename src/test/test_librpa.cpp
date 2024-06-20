#include "../librpa.h"
#include <cstring>

void test_set_librpa_params()
{
    LibRPAParams params;
    get_default_librpa_params(&params);
    strcpy(params.task, "LibRPA_output.txt"); // OKAY
    strcpy(params.output_file, "LibRPA_output.txt"); // FAIL
    strcpy(params.output_dir, "LibRPA_output.txt"); // OKAY
    params.nfreq = 16;
    set_librpa_params(&params);

}

int main(int argc, char *argv[])
{
    test_set_librpa_params();
}
