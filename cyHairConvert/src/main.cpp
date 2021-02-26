#include <iostream>
#include <cstring>
#include <cstdio>
#include "Convert.h"
#include "cyHairFile.h"

int main(int argc, char **argv)
{
    //The first parameter is input file path,
    // the second parameter is output filepath
    if ((argc < 3) || (argc > 3))
    {
        std::cerr << "Wrong Number of Arguments. \n";
        return 0;
    }

    HairConvert convert(argv[1], argv[2]);

    return 1;
}
