#include <stdio.h>
#include <rt/rt_api.h>

int __rt_fpga_fc_frequency = 20000000;// e.g. 20000000 for 20MHz;
int __rt_fpga_periph_frequency = 10000000; // e.g. 10000000 for 10MHz;

int main()
{
    printf("Hello from FPGA\n");
    return 0;
}