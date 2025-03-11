// https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2023-2/cmake-config-for-onemkl.html

#include <stdio.h>
#include "mkl.h"

int main(void)
{
    MKLVersion mkl_version;
    mkl_get_version(&mkl_version);

    printf("oneMKL %d.%d\n", mkl_version.MajorVersion, mkl_version.UpdateVersion);

    return 0;
}