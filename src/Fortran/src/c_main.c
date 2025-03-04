#include <stdio.h>

/* Declare the Fortran subroutine with appropriate name mangling for Intel Fortran */
extern void run_fortran_main(void);

int main(int argc, char **argv) {
    printf("Starting C main function...\n");
    
    /* Call the Fortran subroutine */
    run_fortran_main();
    
    printf("Fortran subroutine completed successfully.\n");
    
    return 0;
}

