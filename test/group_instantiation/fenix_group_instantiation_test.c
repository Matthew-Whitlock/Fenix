#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <fenix.h>
#include <fenix_data_policy.h>

int main(int argc, char **argv)
{
   int i;
   int isends = 1000;
   int its = 50;
   int it;
   MPI_Comm newcomm;

   MPI_Init(&argc, &argv);

   int fenix_status;
   int error;
   int spare_ranks = 0;
   Fenix_Init(&fenix_status, MPI_COMM_WORLD, &newcomm, &argc, &argv, spare_ranks, 0, MPI_INFO_NULL, &error);

   int rank, size;
   MPI_Comm_rank(newcomm, &rank);
   MPI_Comm_size(newcomm, &size);
   printf("Hello world from rank %d/%d\n", rank, size);

   int my_group = 6;
   int policy_vals[2] = {0,1};
   Fenix_Data_group_create(my_group, newcomm, 0, 0, FENIX_DATA_POLICY_IN_MEMORY_RAID, policy_vals, &error);

   int policy_name;
   int policy_vals_recvd[2];
   int flag;
   Fenix_Data_group_get_redundancy_policy(my_group, &policy_name, &policy_vals_recvd, &flag);

   printf("Policy_name: %d, policy vals: %d %d, flag: %d\n", policy_name, policy_vals_recvd[0], policy_vals_recvd[1], flag);
   
   int member1 = 0;
   int member2 = 1; 
   Fenix_Data_member_create(my_group, member1, &policy_vals, 2, MPI_INT);
   Fenix_Data_member_create(my_group, member2, &policy_vals, 2, MPI_INT);

   printf("Created members\n");

   Fenix_Data_member_delete(my_group, member1);
   
   printf("Deleted member 1\n");

   Fenix_Finalize();
   fprintf(stderr, "Finalized Fenix\n");

   MPI_Finalize();
   return 0;
   }
