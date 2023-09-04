#include<iostream>
#include<omp.h>
void main()
{
    int sum = 0;
    int thread_sum[4] = {0, 0, 0, 0};
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        int ID = omp_get_thread_num();
        for(int i = ID * 25 + 1; i <= ID * 25 + 25; i++){
            thread_sum[ID] += i;
            #pragma omp critical
            {
                std::cout << "thread " << ID << " added " << i << std::endl;
            }
            
        }
    }
    for(int i = 0; i < 4; i++){
        sum += thread_sum[i];
    }
    std::cout << "sum = " << sum << std::endl;
    return;
}