#include <iostream>
#include <omp.h>

int main()
{
printf("wat\n");
#pragma omp parallel                   
{
    printf("Hello World... from thread = %d\n", 
           omp_get_thread_num());
}
}

