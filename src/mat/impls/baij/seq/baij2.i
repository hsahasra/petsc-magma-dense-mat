# 1 "baij2.c"
# 1 "/home/dpnkarthik/petsc-rnet/src/mat/impls/baij/seq//"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "baij2.c"

# 1 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/baij.h" 1



# 1 "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h" 1




# 1 "/home/dpnkarthik/petsc-rnet/include/petscmat.h" 1





# 1 "/home/dpnkarthik/petsc-rnet/include/petscvec.h" 1
# 9 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscis.h" 1






# 1 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 1
# 13 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/petscconf.h" 1
# 14 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/petscfix.h" 1
# 15 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 59 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscversion.h" 1
# 60 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 105 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h" 1
# 23 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int MPI_Datatype;
# 124 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int MPI_Comm;




typedef int MPI_Group;



typedef int MPI_Win;







typedef struct ADIOI_FileD *MPI_File;



typedef int MPI_Op;
# 210 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
extern int MPICH_ATTR_FAILED_PROCESSES;


typedef enum MPIR_Topo_type { MPI_GRAPH=1, MPI_CART=2, MPI_DIST_GRAPH=3 } MPIR_Topo_type;
# 227 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef void (MPI_Handler_function) ( MPI_Comm *, int *, ... );
typedef int (MPI_Comm_copy_attr_function)(MPI_Comm, int, void *, void *,
       void *, int *);
typedef int (MPI_Comm_delete_attr_function)(MPI_Comm, int, void *, void *);
typedef int (MPI_Type_copy_attr_function)(MPI_Datatype, int, void *, void *,
       void *, int *);
typedef int (MPI_Type_delete_attr_function)(MPI_Datatype, int, void *, void *);
typedef int (MPI_Win_copy_attr_function)(MPI_Win, int, void *, void *, void *,
      int *);
typedef int (MPI_Win_delete_attr_function)(MPI_Win, int, void *, void *);

typedef void (MPI_Comm_errhandler_function)(MPI_Comm *, int *, ...);
typedef void (MPI_File_errhandler_function)(MPI_File *, int *, ...);
typedef void (MPI_Win_errhandler_function)(MPI_Win *, int *, ...);

typedef MPI_Comm_errhandler_function MPI_Comm_errhandler_fn;
typedef MPI_File_errhandler_function MPI_File_errhandler_fn;
typedef MPI_Win_errhandler_function MPI_Win_errhandler_fn;
# 255 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int MPI_Errhandler;
# 276 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int MPI_Request;


typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * );


typedef int (MPI_Copy_function) ( MPI_Comm, int, void *, void *, void *, int * );
typedef int (MPI_Delete_function) ( MPI_Comm, int, void *, void * );
# 327 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
enum MPIR_Combiner_enum {
    MPI_COMBINER_NAMED = 1,
    MPI_COMBINER_DUP = 2,
    MPI_COMBINER_CONTIGUOUS = 3,
    MPI_COMBINER_VECTOR = 4,
    MPI_COMBINER_HVECTOR_INTEGER = 5,
    MPI_COMBINER_HVECTOR = 6,
    MPI_COMBINER_INDEXED = 7,
    MPI_COMBINER_HINDEXED_INTEGER = 8,
    MPI_COMBINER_HINDEXED = 9,
    MPI_COMBINER_INDEXED_BLOCK = 10,
    MPI_COMBINER_STRUCT_INTEGER = 11,
    MPI_COMBINER_STRUCT = 12,
    MPI_COMBINER_SUBARRAY = 13,
    MPI_COMBINER_DARRAY = 14,
    MPI_COMBINER_F90_REAL = 15,
    MPI_COMBINER_F90_COMPLEX = 16,
    MPI_COMBINER_F90_INTEGER = 17,
    MPI_COMBINER_RESIZED = 18
};


typedef int MPI_Info;
# 372 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef long MPI_Aint;
typedef int MPI_Fint;
# 385 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef long long MPI_Offset;


typedef struct MPI_Status {
    int count;
    int cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;

} MPI_Status;
# 441 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
extern MPI_Fint * MPI_F_STATUS_IGNORE;
extern MPI_Fint * MPI_F_STATUSES_IGNORE;
# 460 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int (MPI_Grequest_cancel_function)(void *, int);
typedef int (MPI_Grequest_free_function)(void *);
typedef int (MPI_Grequest_query_function)(void *, MPI_Status *);
# 546 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int (MPI_Datarep_conversion_function)(void *, MPI_Datatype, int,
             void *, MPI_Offset, void *);
typedef int (MPI_Datarep_extent_function)(MPI_Datatype datatype, MPI_Aint *,
       void *);
# 569 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Get_count(MPI_Status *, MPI_Datatype, int *);
int MPI_Bsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Ssend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Rsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int MPI_Buffer_attach( void*, int);
int MPI_Buffer_detach( void*, int *);
int MPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Ibsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Issend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Irsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Wait(MPI_Request *, MPI_Status *);
int MPI_Test(MPI_Request *, int *, MPI_Status *);
int MPI_Request_free(MPI_Request *);
int MPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Waitall(int, MPI_Request *, MPI_Status *);
int MPI_Testall(int, MPI_Request *, int *, MPI_Status *);
int MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
int MPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *);
int MPI_Probe(int, int, MPI_Comm, MPI_Status *);
int MPI_Cancel(MPI_Request *);
int MPI_Test_cancelled(MPI_Status *, int *);
int MPI_Send_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int MPI_Bsend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Ssend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Rsend_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Recv_init(void*, int, MPI_Datatype, int,int, MPI_Comm, MPI_Request *);
int MPI_Start(MPI_Request *);
int MPI_Startall(int, MPI_Request *);
int MPI_Sendrecv(void *, int, MPI_Datatype,int, int, void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int MPI_Sendrecv_replace(void*, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status *);
int MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int MPI_Type_indexed(int, int *, int *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *);
int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
int MPI_Address(void*, MPI_Aint *);

int MPI_Type_extent(MPI_Datatype, MPI_Aint *);



int MPI_Type_size(MPI_Datatype, int *);

int MPI_Type_lb(MPI_Datatype, MPI_Aint *);
int MPI_Type_ub(MPI_Datatype, MPI_Aint *);
int MPI_Type_commit(MPI_Datatype *);
int MPI_Type_free(MPI_Datatype *);
int MPI_Get_elements(MPI_Status *, MPI_Datatype, int *);
int MPI_Pack(void*, int, MPI_Datatype, void *, int, int *, MPI_Comm);
int MPI_Unpack(void*, int, int *, void *, int, MPI_Datatype, MPI_Comm);
int MPI_Pack_size(int, MPI_Datatype, MPI_Comm, int *);
int MPI_Barrier(MPI_Comm );
int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
int MPI_Gather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int MPI_Gatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, int, MPI_Comm);
int MPI_Scatter(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int MPI_Scatterv(void* , int *, int *, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int MPI_Allgather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
int MPI_Allgatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
int MPI_Alltoall(void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
int MPI_Alltoallv(void* , int *, int *, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
int MPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
int MPI_Op_free( MPI_Op *);
int MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
int MPI_Reduce_scatter(void* , void*, int *, MPI_Datatype, MPI_Op, MPI_Comm);
int MPI_Scan(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm );
int MPI_Group_size(MPI_Group, int *);
int MPI_Group_rank(MPI_Group, int *);
int MPI_Group_translate_ranks (MPI_Group, int, int *, MPI_Group, int *);
int MPI_Group_compare(MPI_Group, MPI_Group, int *);
int MPI_Comm_group(MPI_Comm, MPI_Group *);
int MPI_Group_union(MPI_Group, MPI_Group, MPI_Group *);
int MPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *);
int MPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *);
int MPI_Group_incl(MPI_Group, int, int *, MPI_Group *);
int MPI_Group_excl(MPI_Group, int, int *, MPI_Group *);
int MPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *);
int MPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *);
int MPI_Group_free(MPI_Group *);
int MPI_Comm_size(MPI_Comm, int *);
int MPI_Comm_rank(MPI_Comm, int *);
int MPI_Comm_compare(MPI_Comm, MPI_Comm, int *);
int MPI_Comm_dup(MPI_Comm, MPI_Comm *);
int MPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *);
int MPI_Comm_split(MPI_Comm, int, int, MPI_Comm *);
int MPI_Comm_free(MPI_Comm *);
int MPI_Comm_test_inter(MPI_Comm, int *);
int MPI_Comm_remote_size(MPI_Comm, int *);
int MPI_Comm_remote_group(MPI_Comm, MPI_Group *);
int MPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm * );
int MPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *);
int MPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void*);
int MPI_Keyval_free(int *);
int MPI_Attr_put(MPI_Comm, int, void*);
int MPI_Attr_get(MPI_Comm, int, void *, int *);
int MPI_Attr_delete(MPI_Comm, int);
int MPI_Topo_test(MPI_Comm, int *);
int MPI_Cart_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
int MPI_Dims_create(int, int, int *);
int MPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
int MPI_Graphdims_get(MPI_Comm, int *, int *);
int MPI_Graph_get(MPI_Comm, int, int, int *, int *);
int MPI_Cartdim_get(MPI_Comm, int *);
int MPI_Cart_get(MPI_Comm, int, int *, int *, int *);
int MPI_Cart_rank(MPI_Comm, int *, int *);
int MPI_Cart_coords(MPI_Comm, int, int, int *);
int MPI_Graph_neighbors_count(MPI_Comm, int, int *);
int MPI_Graph_neighbors(MPI_Comm, int, int, int *);
int MPI_Cart_shift(MPI_Comm, int, int, int *, int *);
int MPI_Cart_sub(MPI_Comm, int *, MPI_Comm *);
int MPI_Cart_map(MPI_Comm, int, int *, int *, int *);
int MPI_Graph_map(MPI_Comm, int, int *, int *, int *);
int MPI_Get_processor_name(char *, int *);
int MPI_Get_version(int *, int *);
int MPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *);
int MPI_Errhandler_set(MPI_Comm, MPI_Errhandler);
int MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *);
int MPI_Errhandler_free(MPI_Errhandler *);
int MPI_Error_string(int, char *, int *);
int MPI_Error_class(int, int *);
double MPI_Wtime(void);
double MPI_Wtick(void);

double PMPI_Wtime(void);
double PMPI_Wtick(void);

int MPI_Init(int *, char ***);
int MPI_Finalize(void);
int MPI_Initialized(int *);
int MPI_Abort(MPI_Comm, int);




int MPI_Pcontrol(const int, ...);

int MPIR_Dup_fn ( MPI_Comm, int, void *, void *, void *, int * );





int MPI_Close_port(char *);
int MPI_Comm_accept(char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int MPI_Comm_connect(char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int MPI_Comm_disconnect(MPI_Comm *);
int MPI_Comm_get_parent(MPI_Comm *);
int MPI_Comm_join(int, MPI_Comm *);
int MPI_Comm_spawn(char *, char *[], int, MPI_Info, int, MPI_Comm, MPI_Comm *,
                   int []);
int MPI_Comm_spawn_multiple(int, char *[], char **[], int [], MPI_Info [], int,
       MPI_Comm, MPI_Comm *, int []);
int MPI_Lookup_name(char *, MPI_Info, char *);
int MPI_Open_port(MPI_Info, char *);
int MPI_Publish_name(char *, MPI_Info, char *);
int MPI_Unpublish_name(char *, MPI_Info, char *);


int MPI_Accumulate(void *, int, MPI_Datatype, int, MPI_Aint, int,
     MPI_Datatype, MPI_Op, MPI_Win);
int MPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
     MPI_Win);
int MPI_Put(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
     MPI_Win);
int MPI_Win_complete(MPI_Win);
int MPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *);
int MPI_Win_fence(int, MPI_Win);
int MPI_Win_free(MPI_Win *);
int MPI_Win_get_group(MPI_Win, MPI_Group *);
int MPI_Win_lock(int, int, int, MPI_Win);
int MPI_Win_post(MPI_Group, int, MPI_Win);
int MPI_Win_start(MPI_Group, int, MPI_Win);
int MPI_Win_test(MPI_Win, int *);
int MPI_Win_unlock(int, MPI_Win);
int MPI_Win_wait(MPI_Win);


int MPI_Alltoallw(void *, int [], int [], MPI_Datatype [], void *, int [],
    int [], MPI_Datatype [], MPI_Comm);
int MPI_Exscan(void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm) ;


int MPI_Add_error_class(int *);
int MPI_Add_error_code(int, int *);
int MPI_Add_error_string(int, char *);
int MPI_Comm_call_errhandler(MPI_Comm, int);
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *,
                           MPI_Comm_delete_attr_function *, int *, void *);
int MPI_Comm_delete_attr(MPI_Comm, int);
int MPI_Comm_free_keyval(int *);
int MPI_Comm_get_attr(MPI_Comm, int, void *, int *);
int MPI_Comm_get_name(MPI_Comm, char *, int *);
int MPI_Comm_set_attr(MPI_Comm, int, void *);
int MPI_Comm_set_name(MPI_Comm, char *);
int MPI_File_call_errhandler(MPI_File, int);
int MPI_Grequest_complete(MPI_Request);
int MPI_Grequest_start(MPI_Grequest_query_function *,
                       MPI_Grequest_free_function *,
                       MPI_Grequest_cancel_function *, void *, MPI_Request *);
int MPI_Init_thread(int *, char ***, int, int *);
int MPI_Is_thread_main(int *);
int MPI_Query_thread(int *);
int MPI_Status_set_cancelled(MPI_Status *, int);
int MPI_Status_set_elements(MPI_Status *, MPI_Datatype, int);
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *,
                           MPI_Type_delete_attr_function *, int *, void *);
int MPI_Type_delete_attr(MPI_Datatype, int);
int MPI_Type_dup(MPI_Datatype, MPI_Datatype *);
int MPI_Type_free_keyval(int *);
int MPI_Type_get_attr(MPI_Datatype, int, void *, int *);
int MPI_Type_get_contents(MPI_Datatype, int, int, int, int [], MPI_Aint [],
                          MPI_Datatype []);
int MPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *);
int MPI_Type_get_name(MPI_Datatype, char *, int *);
int MPI_Type_set_attr(MPI_Datatype, int, void *);
int MPI_Type_set_name(MPI_Datatype, char *);
int MPI_Type_match_size( int, int, MPI_Datatype *);
int MPI_Win_call_errhandler(MPI_Win, int);
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *,
                         MPI_Win_delete_attr_function *, int *, void *);
int MPI_Win_delete_attr(MPI_Win, int);
int MPI_Win_free_keyval(int *);
int MPI_Win_get_attr(MPI_Win, int, void *, int *);
int MPI_Win_get_name(MPI_Win, char *, int *);
int MPI_Win_set_attr(MPI_Win, int, void *);
int MPI_Win_set_name(MPI_Win, char *);
# 823 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
int MPI_Alloc_mem(MPI_Aint, MPI_Info info, void *baseptr);
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *, MPI_Errhandler *);
int MPI_Comm_get_errhandler(MPI_Comm, MPI_Errhandler *);
int MPI_Comm_set_errhandler(MPI_Comm, MPI_Errhandler);
int MPI_File_create_errhandler(MPI_File_errhandler_function *, MPI_Errhandler *);
int MPI_File_get_errhandler(MPI_File, MPI_Errhandler *);
int MPI_File_set_errhandler(MPI_File, MPI_Errhandler);
int MPI_Finalized(int *);
int MPI_Free_mem(void *);
int MPI_Get_address(void *, MPI_Aint *);
int MPI_Info_create(MPI_Info *);
int MPI_Info_delete(MPI_Info, char *);
int MPI_Info_dup(MPI_Info, MPI_Info *);
int MPI_Info_free(MPI_Info *info);
int MPI_Info_get(MPI_Info, char *, int, char *, int *);
int MPI_Info_get_nkeys(MPI_Info, int *);
int MPI_Info_get_nthkey(MPI_Info, int, char *);
int MPI_Info_get_valuelen(MPI_Info, char *, int *, int *);
int MPI_Info_set(MPI_Info, char *, char *);
int MPI_Pack_external(char *, void *, int, MPI_Datatype, void *, MPI_Aint,
                      MPI_Aint *);
int MPI_Pack_external_size(char *, int, MPI_Datatype, MPI_Aint *);
int MPI_Request_get_status(MPI_Request, int *, MPI_Status *);
int MPI_Status_c2f(MPI_Status *, MPI_Fint *);
int MPI_Status_f2c(MPI_Fint *, MPI_Status *);
int MPI_Type_create_darray(int, int, int, int [], int [], int [], int [], int,
                           MPI_Datatype, MPI_Datatype *);
int MPI_Type_create_hindexed(int, int [], MPI_Aint [], MPI_Datatype,
                             MPI_Datatype *);
int MPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int MPI_Type_create_indexed_block(int, int, int [], MPI_Datatype,
                                  MPI_Datatype *);
int MPI_Type_create_resized(MPI_Datatype, MPI_Aint, MPI_Aint, MPI_Datatype *);
int MPI_Type_create_struct(int, int [], MPI_Aint [], MPI_Datatype [],
                           MPI_Datatype *);
int MPI_Type_create_subarray(int, int [], int [], int [], int, MPI_Datatype,
                             MPI_Datatype *);
int MPI_Type_get_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int MPI_Type_get_true_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int MPI_Unpack_external(char *, void *, MPI_Aint, MPI_Aint *, void *, int,
                        MPI_Datatype);
int MPI_Win_create_errhandler(MPI_Win_errhandler_function *, MPI_Errhandler *);
int MPI_Win_get_errhandler(MPI_Win, MPI_Errhandler *);
int MPI_Win_set_errhandler(MPI_Win, MPI_Errhandler);




int MPI_Type_create_f90_integer( int, MPI_Datatype * );
int MPI_Type_create_f90_real( int, int, MPI_Datatype * );
int MPI_Type_create_f90_complex( int, int, MPI_Datatype * );


int MPI_Reduce_local(void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op);
int MPI_Op_commutative(MPI_Op op, int *commute);
int MPI_Reduce_scatter_block(void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, int [], int [], int outdegree, int [], int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, int [], int [], int [], int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int [], int [], int maxoutdegree, int [], int []);
# 891 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
int PMPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
int PMPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int PMPI_Get_count(MPI_Status *, MPI_Datatype, int *);
int PMPI_Bsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int PMPI_Ssend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int PMPI_Rsend(void*, int, MPI_Datatype, int, int, MPI_Comm);
int PMPI_Buffer_attach( void* buffer, int);
int PMPI_Buffer_detach( void* buffer, int *);
int PMPI_Isend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Ibsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Issend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Irsend(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Irecv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Wait(MPI_Request *, MPI_Status *);
int PMPI_Test(MPI_Request *, int *, MPI_Status *);
int PMPI_Request_free(MPI_Request *);
int PMPI_Waitany(int, MPI_Request *, int *, MPI_Status *);
int PMPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Waitall(int, MPI_Request *, MPI_Status *);
int PMPI_Testall(int, MPI_Request *, int *, MPI_Status *);
int PMPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *);
int PMPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *);
int PMPI_Probe(int, int, MPI_Comm, MPI_Status *);
int PMPI_Cancel(MPI_Request *);
int PMPI_Test_cancelled(MPI_Status *, int *);
int PMPI_Send_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Bsend_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Ssend_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Rsend_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Recv_init(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *);
int PMPI_Start(MPI_Request *);
int PMPI_Startall(int, MPI_Request *);
int PMPI_Sendrecv(void *, int, MPI_Datatype, int, int, void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
int PMPI_Sendrecv_replace(void*, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status *);
int PMPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_indexed(int, int *, int *, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
int PMPI_Address(void*, MPI_Aint *);
int PMPI_Type_extent(MPI_Datatype, MPI_Aint *);
int PMPI_Type_size(MPI_Datatype, int *);
int PMPI_Type_lb(MPI_Datatype, MPI_Aint *);
int PMPI_Type_ub(MPI_Datatype, MPI_Aint *);
int PMPI_Type_commit(MPI_Datatype *);
int PMPI_Type_free(MPI_Datatype *);
int PMPI_Get_elements(MPI_Status *, MPI_Datatype, int *);
int PMPI_Pack(void*, int, MPI_Datatype, void *, int, int *, MPI_Comm);
int PMPI_Unpack(void*, int, int *, void *, int, MPI_Datatype, MPI_Comm);
int PMPI_Pack_size(int, MPI_Datatype, MPI_Comm, int *);
int PMPI_Barrier(MPI_Comm );
int PMPI_Bcast(void* buffer, int, MPI_Datatype, int, MPI_Comm );
int PMPI_Gather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int PMPI_Gatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, int, MPI_Comm);
int PMPI_Scatter(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int PMPI_Scatterv(void* , int *, int *displs, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm);
int PMPI_Allgather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
int PMPI_Allgatherv(void* , int, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
int PMPI_Alltoall(void* , int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm);
int PMPI_Alltoallv(void* , int *, int *, MPI_Datatype, void*, int *, int *, MPI_Datatype, MPI_Comm);
int PMPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int PMPI_Op_create(MPI_User_function *, int, MPI_Op *);
int PMPI_Op_free( MPI_Op *);
int PMPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
int PMPI_Reduce_scatter(void* , void*, int *, MPI_Datatype, MPI_Op, MPI_Comm);
int PMPI_Scan(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm );
int PMPI_Group_size(MPI_Group, int *);
int PMPI_Group_rank(MPI_Group, int *);
int PMPI_Group_translate_ranks (MPI_Group, int, int *, MPI_Group, int *);
int PMPI_Group_compare(MPI_Group, MPI_Group, int *);
int PMPI_Comm_group(MPI_Comm, MPI_Group *);
int PMPI_Group_union(MPI_Group, MPI_Group, MPI_Group *);
int PMPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *);
int PMPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *);
int PMPI_Group_incl(MPI_Group, int, int *, MPI_Group *);
int PMPI_Group_excl(MPI_Group, int, int *, MPI_Group *);
int PMPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *);
int PMPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *);
int PMPI_Group_free(MPI_Group *);
int PMPI_Comm_size(MPI_Comm, int *);
int PMPI_Comm_rank(MPI_Comm, int *);
int PMPI_Comm_compare(MPI_Comm, MPI_Comm, int *);
int PMPI_Comm_dup(MPI_Comm, MPI_Comm *);
int PMPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *);
int PMPI_Comm_split(MPI_Comm, int, int, MPI_Comm *);
int PMPI_Comm_free(MPI_Comm *);
int PMPI_Comm_test_inter(MPI_Comm, int *);
int PMPI_Comm_remote_size(MPI_Comm, int *);
int PMPI_Comm_remote_group(MPI_Comm, MPI_Group *);
int PMPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm *);
int PMPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *);
int PMPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void*);
int PMPI_Keyval_free(int *);
int PMPI_Attr_put(MPI_Comm, int, void*);
int PMPI_Attr_get(MPI_Comm, int, void *, int *);
int PMPI_Attr_delete(MPI_Comm, int);
int PMPI_Topo_test(MPI_Comm, int *);
int PMPI_Cart_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
int PMPI_Dims_create(int, int, int *);
int PMPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *);
int PMPI_Graphdims_get(MPI_Comm, int *, int *);
int PMPI_Graph_get(MPI_Comm, int, int, int *, int *);
int PMPI_Cartdim_get(MPI_Comm, int *);
int PMPI_Cart_get(MPI_Comm, int, int *, int *, int *);
int PMPI_Cart_rank(MPI_Comm, int *, int *);
int PMPI_Cart_coords(MPI_Comm, int, int, int *);
int PMPI_Graph_neighbors_count(MPI_Comm, int, int *);
int PMPI_Graph_neighbors(MPI_Comm, int, int, int *);
int PMPI_Cart_shift(MPI_Comm, int, int, int *, int *);
int PMPI_Cart_sub(MPI_Comm, int *, MPI_Comm *);
int PMPI_Cart_map(MPI_Comm, int, int *, int *, int *);
int PMPI_Graph_map(MPI_Comm, int, int *, int *, int *);
int PMPI_Get_processor_name(char *, int *);
int PMPI_Get_version(int *, int *);
int PMPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *);
int PMPI_Errhandler_set(MPI_Comm, MPI_Errhandler);
int PMPI_Errhandler_get(MPI_Comm, MPI_Errhandler *);
int PMPI_Errhandler_free(MPI_Errhandler *);
int PMPI_Error_string(int, char *, int *);
int PMPI_Error_class(int, int *);


int PMPI_Init(int *, char ***);
int PMPI_Finalize(void);
int PMPI_Initialized(int *);
int PMPI_Abort(MPI_Comm, int);

int PMPI_Pcontrol(const int, ...);




int PMPI_Close_port(char *);
int PMPI_Comm_accept(char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int PMPI_Comm_connect(char *, MPI_Info, int, MPI_Comm, MPI_Comm *);
int PMPI_Comm_disconnect(MPI_Comm *);
int PMPI_Comm_get_parent(MPI_Comm *);
int PMPI_Comm_join(int, MPI_Comm *);
int PMPI_Comm_spawn(char *, char *[], int, MPI_Info, int, MPI_Comm, MPI_Comm *,
                   int []);
int PMPI_Comm_spawn_multiple(int, char *[], char **[], int [], MPI_Info [], int,
       MPI_Comm, MPI_Comm *, int []);
int PMPI_Lookup_name(char *, MPI_Info, char *);
int PMPI_Open_port(MPI_Info, char *);
int PMPI_Publish_name(char *, MPI_Info, char *);
int PMPI_Unpublish_name(char *, MPI_Info, char *);


int PMPI_Accumulate(void *, int, MPI_Datatype, int, MPI_Aint, int,
     MPI_Datatype, MPI_Op, MPI_Win);
int PMPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
     MPI_Win);
int PMPI_Put(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
     MPI_Win);
int PMPI_Win_complete(MPI_Win);
int PMPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *);
int PMPI_Win_fence(int, MPI_Win);
int PMPI_Win_free(MPI_Win *);
int PMPI_Win_get_group(MPI_Win, MPI_Group *);
int PMPI_Win_lock(int, int, int, MPI_Win);
int PMPI_Win_post(MPI_Group, int, MPI_Win);
int PMPI_Win_start(MPI_Group, int, MPI_Win);
int PMPI_Win_test(MPI_Win, int *);
int PMPI_Win_unlock(int, MPI_Win);
int PMPI_Win_wait(MPI_Win);


int PMPI_Alltoallw(void *, int [], int [], MPI_Datatype [], void *, int [],
    int [], MPI_Datatype [], MPI_Comm);
int PMPI_Exscan(void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm) ;


int PMPI_Add_error_class(int *);
int PMPI_Add_error_code(int, int *);
int PMPI_Add_error_string(int, char *);
int PMPI_Comm_call_errhandler(MPI_Comm, int);
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *,
                           MPI_Comm_delete_attr_function *, int *, void *);
int PMPI_Comm_delete_attr(MPI_Comm, int);
int PMPI_Comm_free_keyval(int *);
int PMPI_Comm_get_attr(MPI_Comm, int, void *, int *);
int PMPI_Comm_get_name(MPI_Comm, char *, int *);
int PMPI_Comm_set_attr(MPI_Comm, int, void *);
int PMPI_Comm_set_name(MPI_Comm, char *);
int PMPI_File_call_errhandler(MPI_File, int);
int PMPI_Grequest_complete(MPI_Request);
int PMPI_Grequest_start(MPI_Grequest_query_function *,
                       MPI_Grequest_free_function *,
                       MPI_Grequest_cancel_function *, void *, MPI_Request *);
int PMPI_Init_thread(int *, char ***, int, int *);
int PMPI_Is_thread_main(int *);
int PMPI_Query_thread(int *);
int PMPI_Status_set_cancelled(MPI_Status *, int);
int PMPI_Status_set_elements(MPI_Status *, MPI_Datatype, int);
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *,
                           MPI_Type_delete_attr_function *, int *, void *);
int PMPI_Type_delete_attr(MPI_Datatype, int);
int PMPI_Type_dup(MPI_Datatype, MPI_Datatype *);
int PMPI_Type_free_keyval(int *);
int PMPI_Type_get_attr(MPI_Datatype, int, void *, int *);
int PMPI_Type_get_contents(MPI_Datatype, int, int, int, int [], MPI_Aint [],
                          MPI_Datatype []);
int PMPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *);
int PMPI_Type_get_name(MPI_Datatype, char *, int *);
int PMPI_Type_set_attr(MPI_Datatype, int, void *);
int PMPI_Type_set_name(MPI_Datatype, char *);
int PMPI_Type_match_size( int, int, MPI_Datatype *);
int PMPI_Win_call_errhandler(MPI_Win, int);
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function *,
                         MPI_Win_delete_attr_function *, int *, void *);
int PMPI_Win_delete_attr(MPI_Win, int);
int PMPI_Win_free_keyval(int *);
int PMPI_Win_get_attr(MPI_Win, int, void *, int *);
int PMPI_Win_get_name(MPI_Win, char *, int *);
int PMPI_Win_set_attr(MPI_Win, int, void *);
int PMPI_Win_set_name(MPI_Win, char *);




int PMPI_Type_create_f90_integer( int, MPI_Datatype * );
int PMPI_Type_create_f90_real( int, int, MPI_Datatype * );
int PMPI_Type_create_f90_complex( int, int, MPI_Datatype * );


int PMPI_Alloc_mem(MPI_Aint, MPI_Info info, void *baseptr);
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function *, MPI_Errhandler *);
int PMPI_Comm_get_errhandler(MPI_Comm, MPI_Errhandler *);
int PMPI_Comm_set_errhandler(MPI_Comm, MPI_Errhandler);
int PMPI_File_create_errhandler(MPI_File_errhandler_function *, MPI_Errhandler *);
int PMPI_File_get_errhandler(MPI_File, MPI_Errhandler *);
int PMPI_File_set_errhandler(MPI_File, MPI_Errhandler);
int PMPI_Finalized(int *);
int PMPI_Free_mem(void *);
int PMPI_Get_address(void *, MPI_Aint *);
int PMPI_Info_create(MPI_Info *);
int PMPI_Info_delete(MPI_Info, char *);
int PMPI_Info_dup(MPI_Info, MPI_Info *);
int PMPI_Info_free(MPI_Info *info);
int PMPI_Info_get(MPI_Info, char *, int, char *, int *);
int PMPI_Info_get_nkeys(MPI_Info, int *);
int PMPI_Info_get_nthkey(MPI_Info, int, char *);
int PMPI_Info_get_valuelen(MPI_Info, char *, int *, int *);
int PMPI_Info_set(MPI_Info, char *, char *);
int PMPI_Pack_external(char *, void *, int, MPI_Datatype, void *, MPI_Aint,
                      MPI_Aint *);
int PMPI_Pack_external_size(char *, int, MPI_Datatype, MPI_Aint *);
int PMPI_Request_get_status(MPI_Request, int *, MPI_Status *);
int PMPI_Status_c2f(MPI_Status *, MPI_Fint *);
int PMPI_Status_f2c(MPI_Fint *, MPI_Status *);
int PMPI_Type_create_darray(int, int, int, int [], int [], int [], int [], int,
                           MPI_Datatype, MPI_Datatype *);
int PMPI_Type_create_hindexed(int, int [], MPI_Aint [], MPI_Datatype,
                             MPI_Datatype *);
int PMPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *);
int PMPI_Type_create_indexed_block(int, int, int [], MPI_Datatype,
                                  MPI_Datatype *);
int PMPI_Type_create_resized(MPI_Datatype, MPI_Aint, MPI_Aint, MPI_Datatype *);
int PMPI_Type_create_struct(int, int [], MPI_Aint [], MPI_Datatype [],
                           MPI_Datatype *);
int PMPI_Type_create_subarray(int, int [], int [], int [], int, MPI_Datatype,
                             MPI_Datatype *);
int PMPI_Type_get_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int PMPI_Type_get_true_extent(MPI_Datatype, MPI_Aint *, MPI_Aint *);
int PMPI_Unpack_external(char *, void *, MPI_Aint, MPI_Aint *, void *, int,
                        MPI_Datatype);
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function *, MPI_Errhandler *);
int PMPI_Win_get_errhandler(MPI_Win, MPI_Errhandler *);
int PMPI_Win_set_errhandler(MPI_Win, MPI_Errhandler);
int PMPI_Reduce_local(void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype, MPI_Op op);
int PMPI_Op_commutative(MPI_Op op, int *commute);
int PMPI_Reduce_scatter_block(void *sendbuf, void *recvbuf, int recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, int [], int [], int outdegree, int [], int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, int [], int [], int [], int [], MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int [], int [], int maxoutdegree, int [], int []);
# 1184 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h" 1
# 13 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h" 1
# 14 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h" 2
# 119 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h"
int MPI_File_open(MPI_Comm, char *, int, MPI_Info, MPI_File *);
int MPI_File_close(MPI_File *);
int MPI_File_delete(char *, MPI_Info);
int MPI_File_set_size(MPI_File, MPI_Offset);
int MPI_File_preallocate(MPI_File, MPI_Offset);
int MPI_File_get_size(MPI_File, MPI_Offset *);
int MPI_File_get_group(MPI_File, MPI_Group *);
int MPI_File_get_amode(MPI_File, int *);
int MPI_File_set_info(MPI_File, MPI_Info);
int MPI_File_get_info(MPI_File, MPI_Info *);


int MPI_File_set_view(MPI_File, MPI_Offset, MPI_Datatype,
          MPI_Datatype, char *, MPI_Info);
int MPI_File_get_view(MPI_File, MPI_Offset *,
                 MPI_Datatype *, MPI_Datatype *, char *);


int MPI_File_read_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);
int MPI_File_read_at_all(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);
int MPI_File_write_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);
int MPI_File_write_at_all(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);





int MPI_File_iread_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Request *);
int MPI_File_iwrite_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Request *);


int MPI_File_read(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int MPI_File_read_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int MPI_File_write(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int MPI_File_write_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *);





int MPI_File_iread(MPI_File, void *, int, MPI_Datatype, MPI_Request *);
int MPI_File_iwrite(MPI_File, void *, int, MPI_Datatype, MPI_Request *);

int MPI_File_seek(MPI_File, MPI_Offset, int);
int MPI_File_get_position(MPI_File, MPI_Offset *);
int MPI_File_get_byte_offset(MPI_File, MPI_Offset, MPI_Offset *);


int MPI_File_read_shared(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int MPI_File_write_shared(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int MPI_File_iread_shared(MPI_File, void *, int, MPI_Datatype, MPI_Request *);
int MPI_File_iwrite_shared(MPI_File, void *, int,
      MPI_Datatype, MPI_Request *);
int MPI_File_read_ordered(MPI_File, void *, int,
                          MPI_Datatype, MPI_Status *);
int MPI_File_write_ordered(MPI_File, void *, int,
                           MPI_Datatype, MPI_Status *);
int MPI_File_seek_shared(MPI_File, MPI_Offset, int);
int MPI_File_get_position_shared(MPI_File, MPI_Offset *);


int MPI_File_read_at_all_begin(MPI_File, MPI_Offset, void *,
                               int, MPI_Datatype);
int MPI_File_read_at_all_end(MPI_File, void *, MPI_Status *);
int MPI_File_write_at_all_begin(MPI_File, MPI_Offset, void *,
                                int, MPI_Datatype);
int MPI_File_write_at_all_end(MPI_File, void *, MPI_Status *);
int MPI_File_read_all_begin(MPI_File, void *, int, MPI_Datatype);
int MPI_File_read_all_end(MPI_File, void *, MPI_Status *);
int MPI_File_write_all_begin(MPI_File, void *, int, MPI_Datatype);
int MPI_File_write_all_end(MPI_File, void *, MPI_Status *);
int MPI_File_read_ordered_begin(MPI_File, void *, int, MPI_Datatype);
int MPI_File_read_ordered_end(MPI_File, void *, MPI_Status *);
int MPI_File_write_ordered_begin(MPI_File, void *, int, MPI_Datatype);
int MPI_File_write_ordered_end(MPI_File, void *, MPI_Status *);


int MPI_File_get_type_extent(MPI_File, MPI_Datatype, MPI_Aint *);


int MPI_Register_datarep(char *,
    MPI_Datarep_conversion_function *,
    MPI_Datarep_conversion_function *,
    MPI_Datarep_extent_function *,
    void *);


int MPI_File_set_atomicity(MPI_File, int);
int MPI_File_get_atomicity(MPI_File, int *);
int MPI_File_sync(MPI_File);
# 248 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h"
MPI_File MPI_File_f2c(MPI_Fint);
MPI_Fint MPI_File_c2f(MPI_File);
# 305 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h"
int PMPI_File_open(MPI_Comm, char *, int, MPI_Info, MPI_File *);
int PMPI_File_close(MPI_File *);
int PMPI_File_delete(char *, MPI_Info);
int PMPI_File_set_size(MPI_File, MPI_Offset);
int PMPI_File_preallocate(MPI_File, MPI_Offset);
int PMPI_File_get_size(MPI_File, MPI_Offset *);
int PMPI_File_get_group(MPI_File, MPI_Group *);
int PMPI_File_get_amode(MPI_File, int *);
int PMPI_File_set_info(MPI_File, MPI_Info);
int PMPI_File_get_info(MPI_File, MPI_Info *);


int PMPI_File_set_view(MPI_File, MPI_Offset,
    MPI_Datatype, MPI_Datatype, char *, MPI_Info);
int PMPI_File_get_view(MPI_File, MPI_Offset *,
      MPI_Datatype *, MPI_Datatype *, char *);


int PMPI_File_read_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);
int PMPI_File_read_at_all(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);
int PMPI_File_write_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);
int PMPI_File_write_at_all(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *);





int PMPI_File_iread_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Request *);
int PMPI_File_iwrite_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Request *);


int PMPI_File_read(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_read_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_write(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_write_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *);





int PMPI_File_iread(MPI_File, void *, int, MPI_Datatype, MPI_Request *);
int PMPI_File_iwrite(MPI_File, void *, int, MPI_Datatype, MPI_Request *);

int PMPI_File_seek(MPI_File, MPI_Offset, int);
int PMPI_File_get_position(MPI_File, MPI_Offset *);
int PMPI_File_get_byte_offset(MPI_File, MPI_Offset, MPI_Offset *);


int PMPI_File_read_shared(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_write_shared(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_iread_shared(MPI_File, void *, int,
      MPI_Datatype, MPI_Request *);
int PMPI_File_iwrite_shared(MPI_File, void *, int,
       MPI_Datatype, MPI_Request *);
int PMPI_File_read_ordered(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_write_ordered(MPI_File, void *, int, MPI_Datatype, MPI_Status *);
int PMPI_File_seek_shared(MPI_File, MPI_Offset, int);
int PMPI_File_get_position_shared(MPI_File, MPI_Offset *);


int PMPI_File_read_at_all_begin(MPI_File, MPI_Offset, void *,
                               int, MPI_Datatype);
int PMPI_File_read_at_all_end(MPI_File, void *, MPI_Status *);
int PMPI_File_write_at_all_begin(MPI_File, MPI_Offset, void *,
                                int, MPI_Datatype);
int PMPI_File_write_at_all_end(MPI_File, void *, MPI_Status *);
int PMPI_File_read_all_begin(MPI_File, void *, int, MPI_Datatype);
int PMPI_File_read_all_end(MPI_File, void *, MPI_Status *);
int PMPI_File_write_all_begin(MPI_File, void *, int, MPI_Datatype);
int PMPI_File_write_all_end(MPI_File, void *, MPI_Status *);
int PMPI_File_read_ordered_begin(MPI_File, void *, int, MPI_Datatype);
int PMPI_File_read_ordered_end(MPI_File, void *, MPI_Status *);
int PMPI_File_write_ordered_begin(MPI_File, void *, int, MPI_Datatype);
int PMPI_File_write_ordered_end(MPI_File, void *, MPI_Status *);


int PMPI_File_get_type_extent(MPI_File, MPI_Datatype, MPI_Aint *);


int PMPI_Register_datarep(char *,
    MPI_Datarep_conversion_function *,
    MPI_Datarep_conversion_function *,
    MPI_Datarep_extent_function *,
    void *);


int PMPI_File_set_atomicity(MPI_File, int);
int PMPI_File_get_atomicity(MPI_File, int *);
int PMPI_File_sync(MPI_File);
# 419 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpio.h"
MPI_File PMPI_File_f2c(MPI_Fint);
MPI_Fint PMPI_File_c2f(MPI_File);
# 1185 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h" 2
# 1208 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/mpi.h"
typedef int (MPIX_Grequest_poll_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_wait_function)(int, void **, double, MPI_Status *);

typedef int MPIX_Grequest_class;
int MPIX_Grequest_class_create(MPI_Grequest_query_function *,
                       MPI_Grequest_free_function *,
                       MPI_Grequest_cancel_function *,
         MPIX_Grequest_poll_function *,
         MPIX_Grequest_wait_function *,
         MPIX_Grequest_class *);

int MPIX_Grequest_class_allocate(MPIX_Grequest_class,
         void *,
         MPI_Request *);

int MPIX_Grequest_start(MPI_Grequest_query_function *,
                       MPI_Grequest_free_function *,
                       MPI_Grequest_cancel_function *,
         MPIX_Grequest_poll_function *,
         MPIX_Grequest_wait_function *, void *, MPI_Request *);

int PMPIX_Grequest_class_create(MPI_Grequest_query_function *,
                       MPI_Grequest_free_function *,
                       MPI_Grequest_cancel_function *,
         MPIX_Grequest_poll_function *,
         MPIX_Grequest_wait_function *,
         MPIX_Grequest_class *);

int PMPIX_Grequest_class_allocate(MPIX_Grequest_class,
         void *,
         MPI_Request *);
int PMPIX_Grequest_start(MPI_Grequest_query_function *,
                       MPI_Grequest_free_function *,
                       MPI_Grequest_cancel_function *,
         MPIX_Grequest_poll_function *,
         MPIX_Grequest_wait_function *, void *, MPI_Request *);
# 106 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2







# 1 "/usr/include/stdio.h" 1 3 4
# 28 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/features.h" 1 3 4
# 361 "/usr/include/features.h" 3 4
# 1 "/usr/include/sys/cdefs.h" 1 3 4
# 365 "/usr/include/sys/cdefs.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 366 "/usr/include/sys/cdefs.h" 2 3 4
# 362 "/usr/include/features.h" 2 3 4
# 385 "/usr/include/features.h" 3 4
# 1 "/usr/include/gnu/stubs.h" 1 3 4



# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 5 "/usr/include/gnu/stubs.h" 2 3 4




# 1 "/usr/include/gnu/stubs-64.h" 1 3 4
# 10 "/usr/include/gnu/stubs.h" 2 3 4
# 386 "/usr/include/features.h" 2 3 4
# 29 "/usr/include/stdio.h" 2 3 4





# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 211 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 3 4
typedef long unsigned int size_t;
# 35 "/usr/include/stdio.h" 2 3 4

# 1 "/usr/include/bits/types.h" 1 3 4
# 28 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 29 "/usr/include/bits/types.h" 2 3 4


typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;


typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;







typedef long int __quad_t;
typedef unsigned long int __u_quad_t;
# 131 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/typesizes.h" 1 3 4
# 132 "/usr/include/bits/types.h" 2 3 4


typedef unsigned long int __dev_t;
typedef unsigned int __uid_t;
typedef unsigned int __gid_t;
typedef unsigned long int __ino_t;
typedef unsigned long int __ino64_t;
typedef unsigned int __mode_t;
typedef unsigned long int __nlink_t;
typedef long int __off_t;
typedef long int __off64_t;
typedef int __pid_t;
typedef struct { int __val[2]; } __fsid_t;
typedef long int __clock_t;
typedef unsigned long int __rlim_t;
typedef unsigned long int __rlim64_t;
typedef unsigned int __id_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;

typedef int __daddr_t;
typedef long int __swblk_t;
typedef int __key_t;


typedef int __clockid_t;


typedef void * __timer_t;


typedef long int __blksize_t;




typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;


typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;


typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;

typedef long int __ssize_t;



typedef __off64_t __loff_t;
typedef __quad_t *__qaddr_t;
typedef char *__caddr_t;


typedef long int __intptr_t;


typedef unsigned int __socklen_t;
# 37 "/usr/include/stdio.h" 2 3 4
# 45 "/usr/include/stdio.h" 3 4
struct _IO_FILE;



typedef struct _IO_FILE FILE;





# 65 "/usr/include/stdio.h" 3 4
typedef struct _IO_FILE __FILE;
# 75 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/libio.h" 1 3 4
# 32 "/usr/include/libio.h" 3 4
# 1 "/usr/include/_G_config.h" 1 3 4
# 15 "/usr/include/_G_config.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 16 "/usr/include/_G_config.h" 2 3 4




# 1 "/usr/include/wchar.h" 1 3 4
# 83 "/usr/include/wchar.h" 3 4
typedef struct
{
  int __count;
  union
  {

    unsigned int __wch;



    char __wchb[4];
  } __value;
} __mbstate_t;
# 21 "/usr/include/_G_config.h" 2 3 4

typedef struct
{
  __off_t __pos;
  __mbstate_t __state;
} _G_fpos_t;
typedef struct
{
  __off64_t __pos;
  __mbstate_t __state;
} _G_fpos64_t;
# 53 "/usr/include/_G_config.h" 3 4
typedef int _G_int16_t __attribute__ ((__mode__ (__HI__)));
typedef int _G_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int _G_uint16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int _G_uint32_t __attribute__ ((__mode__ (__SI__)));
# 33 "/usr/include/libio.h" 2 3 4
# 53 "/usr/include/libio.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stdarg.h" 1 3 4
# 40 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stdarg.h" 3 4
typedef __builtin_va_list __gnuc_va_list;
# 54 "/usr/include/libio.h" 2 3 4
# 170 "/usr/include/libio.h" 3 4
struct _IO_jump_t; struct _IO_FILE;
# 180 "/usr/include/libio.h" 3 4
typedef void _IO_lock_t;





struct _IO_marker {
  struct _IO_marker *_next;
  struct _IO_FILE *_sbuf;



  int _pos;
# 203 "/usr/include/libio.h" 3 4
};


enum __codecvt_result
{
  __codecvt_ok,
  __codecvt_partial,
  __codecvt_error,
  __codecvt_noconv
};
# 271 "/usr/include/libio.h" 3 4
struct _IO_FILE {
  int _flags;




  char* _IO_read_ptr;
  char* _IO_read_end;
  char* _IO_read_base;
  char* _IO_write_base;
  char* _IO_write_ptr;
  char* _IO_write_end;
  char* _IO_buf_base;
  char* _IO_buf_end;

  char *_IO_save_base;
  char *_IO_backup_base;
  char *_IO_save_end;

  struct _IO_marker *_markers;

  struct _IO_FILE *_chain;

  int _fileno;



  int _flags2;

  __off_t _old_offset;



  unsigned short _cur_column;
  signed char _vtable_offset;
  char _shortbuf[1];



  _IO_lock_t *_lock;
# 319 "/usr/include/libio.h" 3 4
  __off64_t _offset;
# 328 "/usr/include/libio.h" 3 4
  void *__pad1;
  void *__pad2;
  void *__pad3;
  void *__pad4;
  size_t __pad5;

  int _mode;

  char _unused2[15 * sizeof (int) - 4 * sizeof (void *) - sizeof (size_t)];

};


typedef struct _IO_FILE _IO_FILE;


struct _IO_FILE_plus;

extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;
# 364 "/usr/include/libio.h" 3 4
typedef __ssize_t __io_read_fn (void *__cookie, char *__buf, size_t __nbytes);







typedef __ssize_t __io_write_fn (void *__cookie, __const char *__buf,
     size_t __n);







typedef int __io_seek_fn (void *__cookie, __off64_t *__pos, int __w);


typedef int __io_close_fn (void *__cookie);
# 416 "/usr/include/libio.h" 3 4
extern int __underflow (_IO_FILE *);
extern int __uflow (_IO_FILE *);
extern int __overflow (_IO_FILE *, int);
# 460 "/usr/include/libio.h" 3 4
extern int _IO_getc (_IO_FILE *__fp);
extern int _IO_putc (int __c, _IO_FILE *__fp);
extern int _IO_feof (_IO_FILE *__fp) __attribute__ ((__nothrow__));
extern int _IO_ferror (_IO_FILE *__fp) __attribute__ ((__nothrow__));

extern int _IO_peekc_locked (_IO_FILE *__fp);





extern void _IO_flockfile (_IO_FILE *) __attribute__ ((__nothrow__));
extern void _IO_funlockfile (_IO_FILE *) __attribute__ ((__nothrow__));
extern int _IO_ftrylockfile (_IO_FILE *) __attribute__ ((__nothrow__));
# 490 "/usr/include/libio.h" 3 4
extern int _IO_vfscanf (_IO_FILE * __restrict, const char * __restrict,
   __gnuc_va_list, int *__restrict);
extern int _IO_vfprintf (_IO_FILE *__restrict, const char *__restrict,
    __gnuc_va_list);
extern __ssize_t _IO_padn (_IO_FILE *, int, __ssize_t);
extern size_t _IO_sgetn (_IO_FILE *, void *, size_t);

extern __off64_t _IO_seekoff (_IO_FILE *, __off64_t, int, int);
extern __off64_t _IO_seekpos (_IO_FILE *, __off64_t, int);

extern void _IO_free_backup_area (_IO_FILE *) __attribute__ ((__nothrow__));
# 76 "/usr/include/stdio.h" 2 3 4




typedef __gnuc_va_list va_list;
# 91 "/usr/include/stdio.h" 3 4
typedef __off_t off_t;
# 103 "/usr/include/stdio.h" 3 4
typedef __ssize_t ssize_t;







typedef _G_fpos_t fpos_t;




# 161 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/bits/stdio_lim.h" 1 3 4
# 162 "/usr/include/stdio.h" 2 3 4



extern struct _IO_FILE *stdin;
extern struct _IO_FILE *stdout;
extern struct _IO_FILE *stderr;









extern int remove (__const char *__filename) __attribute__ ((__nothrow__));

extern int rename (__const char *__old, __const char *__new) __attribute__ ((__nothrow__));




extern int renameat (int __oldfd, __const char *__old, int __newfd,
       __const char *__new) __attribute__ ((__nothrow__));








extern FILE *tmpfile (void) ;
# 208 "/usr/include/stdio.h" 3 4
extern char *tmpnam (char *__s) __attribute__ ((__nothrow__)) ;





extern char *tmpnam_r (char *__s) __attribute__ ((__nothrow__)) ;
# 226 "/usr/include/stdio.h" 3 4
extern char *tempnam (__const char *__dir, __const char *__pfx)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;








extern int fclose (FILE *__stream);




extern int fflush (FILE *__stream);

# 251 "/usr/include/stdio.h" 3 4
extern int fflush_unlocked (FILE *__stream);
# 265 "/usr/include/stdio.h" 3 4






extern FILE *fopen (__const char *__restrict __filename,
      __const char *__restrict __modes) ;




extern FILE *freopen (__const char *__restrict __filename,
        __const char *__restrict __modes,
        FILE *__restrict __stream) ;
# 294 "/usr/include/stdio.h" 3 4

# 305 "/usr/include/stdio.h" 3 4
extern FILE *fdopen (int __fd, __const char *__modes) __attribute__ ((__nothrow__)) ;
# 318 "/usr/include/stdio.h" 3 4
extern FILE *fmemopen (void *__s, size_t __len, __const char *__modes)
  __attribute__ ((__nothrow__)) ;




extern FILE *open_memstream (char **__bufloc, size_t *__sizeloc) __attribute__ ((__nothrow__)) ;






extern void setbuf (FILE *__restrict __stream, char *__restrict __buf) __attribute__ ((__nothrow__));



extern int setvbuf (FILE *__restrict __stream, char *__restrict __buf,
      int __modes, size_t __n) __attribute__ ((__nothrow__));





extern void setbuffer (FILE *__restrict __stream, char *__restrict __buf,
         size_t __size) __attribute__ ((__nothrow__));


extern void setlinebuf (FILE *__stream) __attribute__ ((__nothrow__));








extern int fprintf (FILE *__restrict __stream,
      __const char *__restrict __format, ...);




extern int printf (__const char *__restrict __format, ...);

extern int sprintf (char *__restrict __s,
      __const char *__restrict __format, ...) __attribute__ ((__nothrow__));





extern int vfprintf (FILE *__restrict __s, __const char *__restrict __format,
       __gnuc_va_list __arg);




extern int vprintf (__const char *__restrict __format, __gnuc_va_list __arg);

extern int vsprintf (char *__restrict __s, __const char *__restrict __format,
       __gnuc_va_list __arg) __attribute__ ((__nothrow__));





extern int snprintf (char *__restrict __s, size_t __maxlen,
       __const char *__restrict __format, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
        __const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 0)));

# 416 "/usr/include/stdio.h" 3 4
extern int vdprintf (int __fd, __const char *__restrict __fmt,
       __gnuc_va_list __arg)
     __attribute__ ((__format__ (__printf__, 2, 0)));
extern int dprintf (int __fd, __const char *__restrict __fmt, ...)
     __attribute__ ((__format__ (__printf__, 2, 3)));








extern int fscanf (FILE *__restrict __stream,
     __const char *__restrict __format, ...) ;




extern int scanf (__const char *__restrict __format, ...) ;

extern int sscanf (__const char *__restrict __s,
     __const char *__restrict __format, ...) __attribute__ ((__nothrow__));
# 447 "/usr/include/stdio.h" 3 4
extern int fscanf (FILE *__restrict __stream, __const char *__restrict __format, ...) __asm__ ("" "__isoc99_fscanf")

                               ;
extern int scanf (__const char *__restrict __format, ...) __asm__ ("" "__isoc99_scanf")
                              ;
extern int sscanf (__const char *__restrict __s, __const char *__restrict __format, ...) __asm__ ("" "__isoc99_sscanf")

                          __attribute__ ((__nothrow__));
# 467 "/usr/include/stdio.h" 3 4








extern int vfscanf (FILE *__restrict __s, __const char *__restrict __format,
      __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;





extern int vscanf (__const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;


extern int vsscanf (__const char *__restrict __s,
      __const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__scanf__, 2, 0)));
# 498 "/usr/include/stdio.h" 3 4
extern int vfscanf (FILE *__restrict __s, __const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vfscanf")



     __attribute__ ((__format__ (__scanf__, 2, 0))) ;
extern int vscanf (__const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vscanf")

     __attribute__ ((__format__ (__scanf__, 1, 0))) ;
extern int vsscanf (__const char *__restrict __s, __const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vsscanf")



     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__scanf__, 2, 0)));
# 526 "/usr/include/stdio.h" 3 4









extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);





extern int getchar (void);

# 554 "/usr/include/stdio.h" 3 4
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);
# 565 "/usr/include/stdio.h" 3 4
extern int fgetc_unlocked (FILE *__stream);











extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);





extern int putchar (int __c);

# 598 "/usr/include/stdio.h" 3 4
extern int fputc_unlocked (int __c, FILE *__stream);







extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);






extern int getw (FILE *__stream);


extern int putw (int __w, FILE *__stream);








extern char *fgets (char *__restrict __s, int __n, FILE *__restrict __stream)
     ;






extern char *gets (char *__s) ;

# 660 "/usr/include/stdio.h" 3 4
extern __ssize_t __getdelim (char **__restrict __lineptr,
          size_t *__restrict __n, int __delimiter,
          FILE *__restrict __stream) ;
extern __ssize_t getdelim (char **__restrict __lineptr,
        size_t *__restrict __n, int __delimiter,
        FILE *__restrict __stream) ;







extern __ssize_t getline (char **__restrict __lineptr,
       size_t *__restrict __n,
       FILE *__restrict __stream) ;








extern int fputs (__const char *__restrict __s, FILE *__restrict __stream);





extern int puts (__const char *__s);






extern int ungetc (int __c, FILE *__stream);






extern size_t fread (void *__restrict __ptr, size_t __size,
       size_t __n, FILE *__restrict __stream) ;




extern size_t fwrite (__const void *__restrict __ptr, size_t __size,
        size_t __n, FILE *__restrict __s) ;

# 732 "/usr/include/stdio.h" 3 4
extern size_t fread_unlocked (void *__restrict __ptr, size_t __size,
         size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite_unlocked (__const void *__restrict __ptr, size_t __size,
          size_t __n, FILE *__restrict __stream) ;








extern int fseek (FILE *__stream, long int __off, int __whence);




extern long int ftell (FILE *__stream) ;




extern void rewind (FILE *__stream);

# 768 "/usr/include/stdio.h" 3 4
extern int fseeko (FILE *__stream, __off_t __off, int __whence);




extern __off_t ftello (FILE *__stream) ;
# 787 "/usr/include/stdio.h" 3 4






extern int fgetpos (FILE *__restrict __stream, fpos_t *__restrict __pos);




extern int fsetpos (FILE *__stream, __const fpos_t *__pos);
# 810 "/usr/include/stdio.h" 3 4

# 819 "/usr/include/stdio.h" 3 4


extern void clearerr (FILE *__stream) __attribute__ ((__nothrow__));

extern int feof (FILE *__stream) __attribute__ ((__nothrow__)) ;

extern int ferror (FILE *__stream) __attribute__ ((__nothrow__)) ;




extern void clearerr_unlocked (FILE *__stream) __attribute__ ((__nothrow__));
extern int feof_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;
extern int ferror_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;








extern void perror (__const char *__s);






# 1 "/usr/include/bits/sys_errlist.h" 1 3 4
# 27 "/usr/include/bits/sys_errlist.h" 3 4
extern int sys_nerr;
extern __const char *__const sys_errlist[];
# 849 "/usr/include/stdio.h" 2 3 4




extern int fileno (FILE *__stream) __attribute__ ((__nothrow__)) ;




extern int fileno_unlocked (FILE *__stream) __attribute__ ((__nothrow__)) ;
# 868 "/usr/include/stdio.h" 3 4
extern FILE *popen (__const char *__command, __const char *__modes) ;





extern int pclose (FILE *__stream);





extern char *ctermid (char *__s) __attribute__ ((__nothrow__));
# 908 "/usr/include/stdio.h" 3 4
extern void flockfile (FILE *__stream) __attribute__ ((__nothrow__));



extern int ftrylockfile (FILE *__stream) __attribute__ ((__nothrow__)) ;


extern void funlockfile (FILE *__stream) __attribute__ ((__nothrow__));
# 938 "/usr/include/stdio.h" 3 4

# 114 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 128 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef int PetscErrorCode;
# 141 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef int PetscClassId;
# 171 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef int PetscBLASInt;
# 190 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef int PetscMPIInt;
# 205 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum { ENUM_DUMMY } PetscEnum;
# 220 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef int PetscInt;


typedef long long Petsc64bitInt;
# 233 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum { PETSC_PRECISION_SINGLE=4,PETSC_PRECISION_DOUBLE=8 } PetscPrecision;
extern const char *PetscPrecisions[];
# 256 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern FILE* PETSC_STDOUT;





extern FILE* PETSC_STDERR;





extern FILE* PETSC_ZOPEFD;
# 386 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscmath.h" 1
# 13 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
# 1 "/usr/include/math.h" 1 3 4
# 30 "/usr/include/math.h" 3 4




# 1 "/usr/include/bits/huge_val.h" 1 3 4
# 35 "/usr/include/math.h" 2 3 4

# 1 "/usr/include/bits/huge_valf.h" 1 3 4
# 37 "/usr/include/math.h" 2 3 4
# 1 "/usr/include/bits/huge_vall.h" 1 3 4
# 38 "/usr/include/math.h" 2 3 4


# 1 "/usr/include/bits/inf.h" 1 3 4
# 41 "/usr/include/math.h" 2 3 4


# 1 "/usr/include/bits/nan.h" 1 3 4
# 44 "/usr/include/math.h" 2 3 4



# 1 "/usr/include/bits/mathdef.h" 1 3 4
# 26 "/usr/include/bits/mathdef.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 27 "/usr/include/bits/mathdef.h" 2 3 4




typedef float float_t;
typedef double double_t;
# 48 "/usr/include/math.h" 2 3 4
# 71 "/usr/include/math.h" 3 4
# 1 "/usr/include/bits/mathcalls.h" 1 3 4
# 53 "/usr/include/bits/mathcalls.h" 3 4


extern double acos (double __x) __attribute__ ((__nothrow__)); extern double __acos (double __x) __attribute__ ((__nothrow__));

extern double asin (double __x) __attribute__ ((__nothrow__)); extern double __asin (double __x) __attribute__ ((__nothrow__));

extern double atan (double __x) __attribute__ ((__nothrow__)); extern double __atan (double __x) __attribute__ ((__nothrow__));

extern double atan2 (double __y, double __x) __attribute__ ((__nothrow__)); extern double __atan2 (double __y, double __x) __attribute__ ((__nothrow__));


extern double cos (double __x) __attribute__ ((__nothrow__)); extern double __cos (double __x) __attribute__ ((__nothrow__));

extern double sin (double __x) __attribute__ ((__nothrow__)); extern double __sin (double __x) __attribute__ ((__nothrow__));

extern double tan (double __x) __attribute__ ((__nothrow__)); extern double __tan (double __x) __attribute__ ((__nothrow__));




extern double cosh (double __x) __attribute__ ((__nothrow__)); extern double __cosh (double __x) __attribute__ ((__nothrow__));

extern double sinh (double __x) __attribute__ ((__nothrow__)); extern double __sinh (double __x) __attribute__ ((__nothrow__));

extern double tanh (double __x) __attribute__ ((__nothrow__)); extern double __tanh (double __x) __attribute__ ((__nothrow__));

# 87 "/usr/include/bits/mathcalls.h" 3 4


extern double acosh (double __x) __attribute__ ((__nothrow__)); extern double __acosh (double __x) __attribute__ ((__nothrow__));

extern double asinh (double __x) __attribute__ ((__nothrow__)); extern double __asinh (double __x) __attribute__ ((__nothrow__));

extern double atanh (double __x) __attribute__ ((__nothrow__)); extern double __atanh (double __x) __attribute__ ((__nothrow__));







extern double exp (double __x) __attribute__ ((__nothrow__)); extern double __exp (double __x) __attribute__ ((__nothrow__));


extern double frexp (double __x, int *__exponent) __attribute__ ((__nothrow__)); extern double __frexp (double __x, int *__exponent) __attribute__ ((__nothrow__));


extern double ldexp (double __x, int __exponent) __attribute__ ((__nothrow__)); extern double __ldexp (double __x, int __exponent) __attribute__ ((__nothrow__));


extern double log (double __x) __attribute__ ((__nothrow__)); extern double __log (double __x) __attribute__ ((__nothrow__));


extern double log10 (double __x) __attribute__ ((__nothrow__)); extern double __log10 (double __x) __attribute__ ((__nothrow__));


extern double modf (double __x, double *__iptr) __attribute__ ((__nothrow__)); extern double __modf (double __x, double *__iptr) __attribute__ ((__nothrow__));

# 127 "/usr/include/bits/mathcalls.h" 3 4


extern double expm1 (double __x) __attribute__ ((__nothrow__)); extern double __expm1 (double __x) __attribute__ ((__nothrow__));


extern double log1p (double __x) __attribute__ ((__nothrow__)); extern double __log1p (double __x) __attribute__ ((__nothrow__));


extern double logb (double __x) __attribute__ ((__nothrow__)); extern double __logb (double __x) __attribute__ ((__nothrow__));






extern double exp2 (double __x) __attribute__ ((__nothrow__)); extern double __exp2 (double __x) __attribute__ ((__nothrow__));


extern double log2 (double __x) __attribute__ ((__nothrow__)); extern double __log2 (double __x) __attribute__ ((__nothrow__));








extern double pow (double __x, double __y) __attribute__ ((__nothrow__)); extern double __pow (double __x, double __y) __attribute__ ((__nothrow__));


extern double sqrt (double __x) __attribute__ ((__nothrow__)); extern double __sqrt (double __x) __attribute__ ((__nothrow__));





extern double hypot (double __x, double __y) __attribute__ ((__nothrow__)); extern double __hypot (double __x, double __y) __attribute__ ((__nothrow__));






extern double cbrt (double __x) __attribute__ ((__nothrow__)); extern double __cbrt (double __x) __attribute__ ((__nothrow__));








extern double ceil (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __ceil (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern double fabs (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __fabs (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern double floor (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __floor (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern double fmod (double __x, double __y) __attribute__ ((__nothrow__)); extern double __fmod (double __x, double __y) __attribute__ ((__nothrow__));




extern int __isinf (double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern int __finite (double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));





extern int isinf (double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern int finite (double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern double drem (double __x, double __y) __attribute__ ((__nothrow__)); extern double __drem (double __x, double __y) __attribute__ ((__nothrow__));



extern double significand (double __x) __attribute__ ((__nothrow__)); extern double __significand (double __x) __attribute__ ((__nothrow__));





extern double copysign (double __x, double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __copysign (double __x, double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));






extern double nan (__const char *__tagb) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __nan (__const char *__tagb) __attribute__ ((__nothrow__)) __attribute__ ((__const__));





extern int __isnan (double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern int isnan (double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern double j0 (double) __attribute__ ((__nothrow__)); extern double __j0 (double) __attribute__ ((__nothrow__));
extern double j1 (double) __attribute__ ((__nothrow__)); extern double __j1 (double) __attribute__ ((__nothrow__));
extern double jn (int, double) __attribute__ ((__nothrow__)); extern double __jn (int, double) __attribute__ ((__nothrow__));
extern double y0 (double) __attribute__ ((__nothrow__)); extern double __y0 (double) __attribute__ ((__nothrow__));
extern double y1 (double) __attribute__ ((__nothrow__)); extern double __y1 (double) __attribute__ ((__nothrow__));
extern double yn (int, double) __attribute__ ((__nothrow__)); extern double __yn (int, double) __attribute__ ((__nothrow__));






extern double erf (double) __attribute__ ((__nothrow__)); extern double __erf (double) __attribute__ ((__nothrow__));
extern double erfc (double) __attribute__ ((__nothrow__)); extern double __erfc (double) __attribute__ ((__nothrow__));
extern double lgamma (double) __attribute__ ((__nothrow__)); extern double __lgamma (double) __attribute__ ((__nothrow__));






extern double tgamma (double) __attribute__ ((__nothrow__)); extern double __tgamma (double) __attribute__ ((__nothrow__));





extern double gamma (double) __attribute__ ((__nothrow__)); extern double __gamma (double) __attribute__ ((__nothrow__));






extern double lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__)); extern double __lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__));







extern double rint (double __x) __attribute__ ((__nothrow__)); extern double __rint (double __x) __attribute__ ((__nothrow__));


extern double nextafter (double __x, double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __nextafter (double __x, double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));

extern double nexttoward (double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __nexttoward (double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern double remainder (double __x, double __y) __attribute__ ((__nothrow__)); extern double __remainder (double __x, double __y) __attribute__ ((__nothrow__));



extern double scalbn (double __x, int __n) __attribute__ ((__nothrow__)); extern double __scalbn (double __x, int __n) __attribute__ ((__nothrow__));



extern int ilogb (double __x) __attribute__ ((__nothrow__)); extern int __ilogb (double __x) __attribute__ ((__nothrow__));




extern double scalbln (double __x, long int __n) __attribute__ ((__nothrow__)); extern double __scalbln (double __x, long int __n) __attribute__ ((__nothrow__));



extern double nearbyint (double __x) __attribute__ ((__nothrow__)); extern double __nearbyint (double __x) __attribute__ ((__nothrow__));



extern double round (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __round (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern double trunc (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern double __trunc (double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));




extern double remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__)); extern double __remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__));






extern long int lrint (double __x) __attribute__ ((__nothrow__)); extern long int __lrint (double __x) __attribute__ ((__nothrow__));
extern long long int llrint (double __x) __attribute__ ((__nothrow__)); extern long long int __llrint (double __x) __attribute__ ((__nothrow__));



extern long int lround (double __x) __attribute__ ((__nothrow__)); extern long int __lround (double __x) __attribute__ ((__nothrow__));
extern long long int llround (double __x) __attribute__ ((__nothrow__)); extern long long int __llround (double __x) __attribute__ ((__nothrow__));



extern double fdim (double __x, double __y) __attribute__ ((__nothrow__)); extern double __fdim (double __x, double __y) __attribute__ ((__nothrow__));


extern double fmax (double __x, double __y) __attribute__ ((__nothrow__)); extern double __fmax (double __x, double __y) __attribute__ ((__nothrow__));


extern double fmin (double __x, double __y) __attribute__ ((__nothrow__)); extern double __fmin (double __x, double __y) __attribute__ ((__nothrow__));



extern int __fpclassify (double __value) __attribute__ ((__nothrow__))
     __attribute__ ((__const__));


extern int __signbit (double __value) __attribute__ ((__nothrow__))
     __attribute__ ((__const__));



extern double fma (double __x, double __y, double __z) __attribute__ ((__nothrow__)); extern double __fma (double __x, double __y, double __z) __attribute__ ((__nothrow__));








extern double scalb (double __x, double __n) __attribute__ ((__nothrow__)); extern double __scalb (double __x, double __n) __attribute__ ((__nothrow__));
# 72 "/usr/include/math.h" 2 3 4
# 94 "/usr/include/math.h" 3 4
# 1 "/usr/include/bits/mathcalls.h" 1 3 4
# 53 "/usr/include/bits/mathcalls.h" 3 4


extern float acosf (float __x) __attribute__ ((__nothrow__)); extern float __acosf (float __x) __attribute__ ((__nothrow__));

extern float asinf (float __x) __attribute__ ((__nothrow__)); extern float __asinf (float __x) __attribute__ ((__nothrow__));

extern float atanf (float __x) __attribute__ ((__nothrow__)); extern float __atanf (float __x) __attribute__ ((__nothrow__));

extern float atan2f (float __y, float __x) __attribute__ ((__nothrow__)); extern float __atan2f (float __y, float __x) __attribute__ ((__nothrow__));


extern float cosf (float __x) __attribute__ ((__nothrow__)); extern float __cosf (float __x) __attribute__ ((__nothrow__));

extern float sinf (float __x) __attribute__ ((__nothrow__)); extern float __sinf (float __x) __attribute__ ((__nothrow__));

extern float tanf (float __x) __attribute__ ((__nothrow__)); extern float __tanf (float __x) __attribute__ ((__nothrow__));




extern float coshf (float __x) __attribute__ ((__nothrow__)); extern float __coshf (float __x) __attribute__ ((__nothrow__));

extern float sinhf (float __x) __attribute__ ((__nothrow__)); extern float __sinhf (float __x) __attribute__ ((__nothrow__));

extern float tanhf (float __x) __attribute__ ((__nothrow__)); extern float __tanhf (float __x) __attribute__ ((__nothrow__));

# 87 "/usr/include/bits/mathcalls.h" 3 4


extern float acoshf (float __x) __attribute__ ((__nothrow__)); extern float __acoshf (float __x) __attribute__ ((__nothrow__));

extern float asinhf (float __x) __attribute__ ((__nothrow__)); extern float __asinhf (float __x) __attribute__ ((__nothrow__));

extern float atanhf (float __x) __attribute__ ((__nothrow__)); extern float __atanhf (float __x) __attribute__ ((__nothrow__));







extern float expf (float __x) __attribute__ ((__nothrow__)); extern float __expf (float __x) __attribute__ ((__nothrow__));


extern float frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__)); extern float __frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__));


extern float ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__)); extern float __ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__));


extern float logf (float __x) __attribute__ ((__nothrow__)); extern float __logf (float __x) __attribute__ ((__nothrow__));


extern float log10f (float __x) __attribute__ ((__nothrow__)); extern float __log10f (float __x) __attribute__ ((__nothrow__));


extern float modff (float __x, float *__iptr) __attribute__ ((__nothrow__)); extern float __modff (float __x, float *__iptr) __attribute__ ((__nothrow__));

# 127 "/usr/include/bits/mathcalls.h" 3 4


extern float expm1f (float __x) __attribute__ ((__nothrow__)); extern float __expm1f (float __x) __attribute__ ((__nothrow__));


extern float log1pf (float __x) __attribute__ ((__nothrow__)); extern float __log1pf (float __x) __attribute__ ((__nothrow__));


extern float logbf (float __x) __attribute__ ((__nothrow__)); extern float __logbf (float __x) __attribute__ ((__nothrow__));






extern float exp2f (float __x) __attribute__ ((__nothrow__)); extern float __exp2f (float __x) __attribute__ ((__nothrow__));


extern float log2f (float __x) __attribute__ ((__nothrow__)); extern float __log2f (float __x) __attribute__ ((__nothrow__));








extern float powf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __powf (float __x, float __y) __attribute__ ((__nothrow__));


extern float sqrtf (float __x) __attribute__ ((__nothrow__)); extern float __sqrtf (float __x) __attribute__ ((__nothrow__));





extern float hypotf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __hypotf (float __x, float __y) __attribute__ ((__nothrow__));






extern float cbrtf (float __x) __attribute__ ((__nothrow__)); extern float __cbrtf (float __x) __attribute__ ((__nothrow__));








extern float ceilf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __ceilf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern float fabsf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __fabsf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern float floorf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __floorf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern float fmodf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __fmodf (float __x, float __y) __attribute__ ((__nothrow__));




extern int __isinff (float __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern int __finitef (float __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));





extern int isinff (float __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern int finitef (float __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern float dremf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __dremf (float __x, float __y) __attribute__ ((__nothrow__));



extern float significandf (float __x) __attribute__ ((__nothrow__)); extern float __significandf (float __x) __attribute__ ((__nothrow__));





extern float copysignf (float __x, float __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __copysignf (float __x, float __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));






extern float nanf (__const char *__tagb) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __nanf (__const char *__tagb) __attribute__ ((__nothrow__)) __attribute__ ((__const__));





extern int __isnanf (float __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern int isnanf (float __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern float j0f (float) __attribute__ ((__nothrow__)); extern float __j0f (float) __attribute__ ((__nothrow__));
extern float j1f (float) __attribute__ ((__nothrow__)); extern float __j1f (float) __attribute__ ((__nothrow__));
extern float jnf (int, float) __attribute__ ((__nothrow__)); extern float __jnf (int, float) __attribute__ ((__nothrow__));
extern float y0f (float) __attribute__ ((__nothrow__)); extern float __y0f (float) __attribute__ ((__nothrow__));
extern float y1f (float) __attribute__ ((__nothrow__)); extern float __y1f (float) __attribute__ ((__nothrow__));
extern float ynf (int, float) __attribute__ ((__nothrow__)); extern float __ynf (int, float) __attribute__ ((__nothrow__));






extern float erff (float) __attribute__ ((__nothrow__)); extern float __erff (float) __attribute__ ((__nothrow__));
extern float erfcf (float) __attribute__ ((__nothrow__)); extern float __erfcf (float) __attribute__ ((__nothrow__));
extern float lgammaf (float) __attribute__ ((__nothrow__)); extern float __lgammaf (float) __attribute__ ((__nothrow__));






extern float tgammaf (float) __attribute__ ((__nothrow__)); extern float __tgammaf (float) __attribute__ ((__nothrow__));





extern float gammaf (float) __attribute__ ((__nothrow__)); extern float __gammaf (float) __attribute__ ((__nothrow__));






extern float lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__)); extern float __lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__));







extern float rintf (float __x) __attribute__ ((__nothrow__)); extern float __rintf (float __x) __attribute__ ((__nothrow__));


extern float nextafterf (float __x, float __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __nextafterf (float __x, float __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));

extern float nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern float remainderf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __remainderf (float __x, float __y) __attribute__ ((__nothrow__));



extern float scalbnf (float __x, int __n) __attribute__ ((__nothrow__)); extern float __scalbnf (float __x, int __n) __attribute__ ((__nothrow__));



extern int ilogbf (float __x) __attribute__ ((__nothrow__)); extern int __ilogbf (float __x) __attribute__ ((__nothrow__));




extern float scalblnf (float __x, long int __n) __attribute__ ((__nothrow__)); extern float __scalblnf (float __x, long int __n) __attribute__ ((__nothrow__));



extern float nearbyintf (float __x) __attribute__ ((__nothrow__)); extern float __nearbyintf (float __x) __attribute__ ((__nothrow__));



extern float roundf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __roundf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern float truncf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern float __truncf (float __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));




extern float remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__)); extern float __remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__));






extern long int lrintf (float __x) __attribute__ ((__nothrow__)); extern long int __lrintf (float __x) __attribute__ ((__nothrow__));
extern long long int llrintf (float __x) __attribute__ ((__nothrow__)); extern long long int __llrintf (float __x) __attribute__ ((__nothrow__));



extern long int lroundf (float __x) __attribute__ ((__nothrow__)); extern long int __lroundf (float __x) __attribute__ ((__nothrow__));
extern long long int llroundf (float __x) __attribute__ ((__nothrow__)); extern long long int __llroundf (float __x) __attribute__ ((__nothrow__));



extern float fdimf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __fdimf (float __x, float __y) __attribute__ ((__nothrow__));


extern float fmaxf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __fmaxf (float __x, float __y) __attribute__ ((__nothrow__));


extern float fminf (float __x, float __y) __attribute__ ((__nothrow__)); extern float __fminf (float __x, float __y) __attribute__ ((__nothrow__));



extern int __fpclassifyf (float __value) __attribute__ ((__nothrow__))
     __attribute__ ((__const__));


extern int __signbitf (float __value) __attribute__ ((__nothrow__))
     __attribute__ ((__const__));



extern float fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__)); extern float __fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__));








extern float scalbf (float __x, float __n) __attribute__ ((__nothrow__)); extern float __scalbf (float __x, float __n) __attribute__ ((__nothrow__));
# 95 "/usr/include/math.h" 2 3 4
# 141 "/usr/include/math.h" 3 4
# 1 "/usr/include/bits/mathcalls.h" 1 3 4
# 53 "/usr/include/bits/mathcalls.h" 3 4


extern long double acosl (long double __x) __attribute__ ((__nothrow__)); extern long double __acosl (long double __x) __attribute__ ((__nothrow__));

extern long double asinl (long double __x) __attribute__ ((__nothrow__)); extern long double __asinl (long double __x) __attribute__ ((__nothrow__));

extern long double atanl (long double __x) __attribute__ ((__nothrow__)); extern long double __atanl (long double __x) __attribute__ ((__nothrow__));

extern long double atan2l (long double __y, long double __x) __attribute__ ((__nothrow__)); extern long double __atan2l (long double __y, long double __x) __attribute__ ((__nothrow__));


extern long double cosl (long double __x) __attribute__ ((__nothrow__)); extern long double __cosl (long double __x) __attribute__ ((__nothrow__));

extern long double sinl (long double __x) __attribute__ ((__nothrow__)); extern long double __sinl (long double __x) __attribute__ ((__nothrow__));

extern long double tanl (long double __x) __attribute__ ((__nothrow__)); extern long double __tanl (long double __x) __attribute__ ((__nothrow__));




extern long double coshl (long double __x) __attribute__ ((__nothrow__)); extern long double __coshl (long double __x) __attribute__ ((__nothrow__));

extern long double sinhl (long double __x) __attribute__ ((__nothrow__)); extern long double __sinhl (long double __x) __attribute__ ((__nothrow__));

extern long double tanhl (long double __x) __attribute__ ((__nothrow__)); extern long double __tanhl (long double __x) __attribute__ ((__nothrow__));

# 87 "/usr/include/bits/mathcalls.h" 3 4


extern long double acoshl (long double __x) __attribute__ ((__nothrow__)); extern long double __acoshl (long double __x) __attribute__ ((__nothrow__));

extern long double asinhl (long double __x) __attribute__ ((__nothrow__)); extern long double __asinhl (long double __x) __attribute__ ((__nothrow__));

extern long double atanhl (long double __x) __attribute__ ((__nothrow__)); extern long double __atanhl (long double __x) __attribute__ ((__nothrow__));







extern long double expl (long double __x) __attribute__ ((__nothrow__)); extern long double __expl (long double __x) __attribute__ ((__nothrow__));


extern long double frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__)); extern long double __frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__));


extern long double ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__)); extern long double __ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__));


extern long double logl (long double __x) __attribute__ ((__nothrow__)); extern long double __logl (long double __x) __attribute__ ((__nothrow__));


extern long double log10l (long double __x) __attribute__ ((__nothrow__)); extern long double __log10l (long double __x) __attribute__ ((__nothrow__));


extern long double modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__)); extern long double __modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__));

# 127 "/usr/include/bits/mathcalls.h" 3 4


extern long double expm1l (long double __x) __attribute__ ((__nothrow__)); extern long double __expm1l (long double __x) __attribute__ ((__nothrow__));


extern long double log1pl (long double __x) __attribute__ ((__nothrow__)); extern long double __log1pl (long double __x) __attribute__ ((__nothrow__));


extern long double logbl (long double __x) __attribute__ ((__nothrow__)); extern long double __logbl (long double __x) __attribute__ ((__nothrow__));






extern long double exp2l (long double __x) __attribute__ ((__nothrow__)); extern long double __exp2l (long double __x) __attribute__ ((__nothrow__));


extern long double log2l (long double __x) __attribute__ ((__nothrow__)); extern long double __log2l (long double __x) __attribute__ ((__nothrow__));








extern long double powl (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __powl (long double __x, long double __y) __attribute__ ((__nothrow__));


extern long double sqrtl (long double __x) __attribute__ ((__nothrow__)); extern long double __sqrtl (long double __x) __attribute__ ((__nothrow__));





extern long double hypotl (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __hypotl (long double __x, long double __y) __attribute__ ((__nothrow__));






extern long double cbrtl (long double __x) __attribute__ ((__nothrow__)); extern long double __cbrtl (long double __x) __attribute__ ((__nothrow__));








extern long double ceill (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __ceill (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern long double fabsl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __fabsl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern long double floorl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __floorl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern long double fmodl (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __fmodl (long double __x, long double __y) __attribute__ ((__nothrow__));




extern int __isinfl (long double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern int __finitel (long double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));





extern int isinfl (long double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern int finitel (long double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern long double dreml (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __dreml (long double __x, long double __y) __attribute__ ((__nothrow__));



extern long double significandl (long double __x) __attribute__ ((__nothrow__)); extern long double __significandl (long double __x) __attribute__ ((__nothrow__));





extern long double copysignl (long double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __copysignl (long double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));






extern long double nanl (__const char *__tagb) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __nanl (__const char *__tagb) __attribute__ ((__nothrow__)) __attribute__ ((__const__));





extern int __isnanl (long double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern int isnanl (long double __value) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern long double j0l (long double) __attribute__ ((__nothrow__)); extern long double __j0l (long double) __attribute__ ((__nothrow__));
extern long double j1l (long double) __attribute__ ((__nothrow__)); extern long double __j1l (long double) __attribute__ ((__nothrow__));
extern long double jnl (int, long double) __attribute__ ((__nothrow__)); extern long double __jnl (int, long double) __attribute__ ((__nothrow__));
extern long double y0l (long double) __attribute__ ((__nothrow__)); extern long double __y0l (long double) __attribute__ ((__nothrow__));
extern long double y1l (long double) __attribute__ ((__nothrow__)); extern long double __y1l (long double) __attribute__ ((__nothrow__));
extern long double ynl (int, long double) __attribute__ ((__nothrow__)); extern long double __ynl (int, long double) __attribute__ ((__nothrow__));






extern long double erfl (long double) __attribute__ ((__nothrow__)); extern long double __erfl (long double) __attribute__ ((__nothrow__));
extern long double erfcl (long double) __attribute__ ((__nothrow__)); extern long double __erfcl (long double) __attribute__ ((__nothrow__));
extern long double lgammal (long double) __attribute__ ((__nothrow__)); extern long double __lgammal (long double) __attribute__ ((__nothrow__));






extern long double tgammal (long double) __attribute__ ((__nothrow__)); extern long double __tgammal (long double) __attribute__ ((__nothrow__));





extern long double gammal (long double) __attribute__ ((__nothrow__)); extern long double __gammal (long double) __attribute__ ((__nothrow__));






extern long double lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__)); extern long double __lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__));







extern long double rintl (long double __x) __attribute__ ((__nothrow__)); extern long double __rintl (long double __x) __attribute__ ((__nothrow__));


extern long double nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));

extern long double nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern long double remainderl (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __remainderl (long double __x, long double __y) __attribute__ ((__nothrow__));



extern long double scalbnl (long double __x, int __n) __attribute__ ((__nothrow__)); extern long double __scalbnl (long double __x, int __n) __attribute__ ((__nothrow__));



extern int ilogbl (long double __x) __attribute__ ((__nothrow__)); extern int __ilogbl (long double __x) __attribute__ ((__nothrow__));




extern long double scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__)); extern long double __scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__));



extern long double nearbyintl (long double __x) __attribute__ ((__nothrow__)); extern long double __nearbyintl (long double __x) __attribute__ ((__nothrow__));



extern long double roundl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __roundl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));



extern long double truncl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)); extern long double __truncl (long double __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__));




extern long double remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__)); extern long double __remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__));






extern long int lrintl (long double __x) __attribute__ ((__nothrow__)); extern long int __lrintl (long double __x) __attribute__ ((__nothrow__));
extern long long int llrintl (long double __x) __attribute__ ((__nothrow__)); extern long long int __llrintl (long double __x) __attribute__ ((__nothrow__));



extern long int lroundl (long double __x) __attribute__ ((__nothrow__)); extern long int __lroundl (long double __x) __attribute__ ((__nothrow__));
extern long long int llroundl (long double __x) __attribute__ ((__nothrow__)); extern long long int __llroundl (long double __x) __attribute__ ((__nothrow__));



extern long double fdiml (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __fdiml (long double __x, long double __y) __attribute__ ((__nothrow__));


extern long double fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__));


extern long double fminl (long double __x, long double __y) __attribute__ ((__nothrow__)); extern long double __fminl (long double __x, long double __y) __attribute__ ((__nothrow__));



extern int __fpclassifyl (long double __value) __attribute__ ((__nothrow__))
     __attribute__ ((__const__));


extern int __signbitl (long double __value) __attribute__ ((__nothrow__))
     __attribute__ ((__const__));



extern long double fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__)); extern long double __fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__));








extern long double scalbl (long double __x, long double __n) __attribute__ ((__nothrow__)); extern long double __scalbl (long double __x, long double __n) __attribute__ ((__nothrow__));
# 142 "/usr/include/math.h" 2 3 4
# 157 "/usr/include/math.h" 3 4
extern int signgam;
# 198 "/usr/include/math.h" 3 4
enum
  {
    FP_NAN,

    FP_INFINITE,

    FP_ZERO,

    FP_SUBNORMAL,

    FP_NORMAL

  };
# 291 "/usr/include/math.h" 3 4
typedef enum
{
  _IEEE_ = -1,
  _SVID_,
  _XOPEN_,
  _POSIX_,
  _ISOC_
} _LIB_VERSION_TYPE;




extern _LIB_VERSION_TYPE _LIB_VERSION;
# 316 "/usr/include/math.h" 3 4
struct exception

  {
    int type;
    char *name;
    double arg1;
    double arg2;
    double retval;
  };




extern int matherr (struct exception *__exc);
# 472 "/usr/include/math.h" 3 4

# 14 "/home/dpnkarthik/petsc-rnet/include/petscmath.h" 2


extern MPI_Datatype MPIU_2SCALAR;
extern MPI_Datatype MPIU_2INT;
# 35 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
typedef double PetscReal;
# 126 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
typedef double PetscScalar;







static inline PetscReal PetscAbsScalar(PetscScalar a) {return a < 0.0 ? -a : a;}
# 164 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
typedef enum { PETSC_SCALAR_DOUBLE,PETSC_SCALAR_SINGLE, PETSC_SCALAR_LONG_DOUBLE } PetscScalarPrecision;


extern PetscScalar PETSC_i;
# 317 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
extern PetscErrorCode PetscIsInfOrNanScalar(PetscScalar);
extern PetscErrorCode PetscIsInfOrNanReal(PetscReal);







typedef double PetscLogDouble;
# 337 "/home/dpnkarthik/petsc-rnet/include/petscmath.h"
typedef PetscScalar MatScalar;
typedef PetscReal MatReal;



# 387 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2






# 407 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum { PETSC_FALSE,PETSC_TRUE } PetscBool;
extern const char *PetscBools[];
# 422 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum { PETSC_COPY_VALUES, PETSC_OWN_POINTER, PETSC_USE_POINTER} PetscCopyMode;
extern const char *PetscCopyModes[];
# 533 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern MPI_Comm PETSC_COMM_WORLD;
# 545 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscBool PetscInitializeCalled;
extern PetscBool PetscFinalizeCalled;

extern PetscErrorCode PetscSetHelpVersionFunctions(PetscErrorCode (*)(MPI_Comm),PetscErrorCode (*)(MPI_Comm));
extern PetscErrorCode PetscCommDuplicate(MPI_Comm,MPI_Comm*,int*);
extern PetscErrorCode PetscCommDestroy(MPI_Comm*);
# 1115 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode (*PetscTrMalloc)(size_t,int,const char[],const char[],const char[],void**);
extern PetscErrorCode (*PetscTrFree)(void*,int,const char[],const char[],const char[]);
extern PetscErrorCode PetscMallocSet(PetscErrorCode (*)(size_t,int,const char[],const char[],const char[],void**),PetscErrorCode (*)(void*,int,const char[],const char[],const char[]));
extern PetscErrorCode PetscMallocClear(void);




extern PetscErrorCode PetscMallocDump(FILE *);
extern PetscErrorCode PetscMallocDumpLog(FILE *);
extern PetscErrorCode PetscMallocGetCurrentUsage(PetscLogDouble *);
extern PetscErrorCode PetscMallocGetMaximumUsage(PetscLogDouble *);
extern PetscErrorCode PetscMallocDebug(PetscBool);
extern PetscErrorCode PetscMallocValidate(int,const char[],const char[],const char[]);
extern PetscErrorCode PetscMallocSetDumpLog(void);
# 1143 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum {PETSC_INT = 0,PETSC_DOUBLE = 1,PETSC_COMPLEX = 2, PETSC_LONG = 3 ,PETSC_SHORT = 4,PETSC_FLOAT = 5,
              PETSC_CHAR = 6,PETSC_BIT_LOGICAL = 7,PETSC_ENUM = 8,PETSC_BOOL=9, PETSC_LONG_DOUBLE = 10} PetscDataType;
extern const char *PetscDataTypes[];
# 1163 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscDataTypeToMPIDataType(PetscDataType,MPI_Datatype*);
extern PetscErrorCode PetscMPIDataTypeToPetscDataType(MPI_Datatype,PetscDataType*);
extern PetscErrorCode PetscDataTypeGetSize(PetscDataType,size_t*);






extern PetscErrorCode PetscBitMemcpy(void*,PetscInt,const void*,PetscInt,PetscInt,PetscDataType);
extern PetscErrorCode PetscMemmove(void*,void *,size_t);
extern PetscErrorCode PetscMemcmp(const void*,const void*,size_t,PetscBool *);
extern PetscErrorCode PetscStrlen(const char[],size_t*);
extern PetscErrorCode PetscStrToArray(const char[],int*,char ***);
extern PetscErrorCode PetscStrToArrayDestroy(int,char **);
extern PetscErrorCode PetscStrcmp(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscStrgrt(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscStrcasecmp(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscStrncmp(const char[],const char[],size_t,PetscBool *);
extern PetscErrorCode PetscStrcpy(char[],const char[]);
extern PetscErrorCode PetscStrcat(char[],const char[]);
extern PetscErrorCode PetscStrncat(char[],const char[],size_t);
extern PetscErrorCode PetscStrncpy(char[],const char[],size_t);
extern PetscErrorCode PetscStrchr(const char[],char,char *[]);
extern PetscErrorCode PetscStrtolower(char[]);
extern PetscErrorCode PetscStrrchr(const char[],char,char *[]);
extern PetscErrorCode PetscStrstr(const char[],const char[],char *[]);
extern PetscErrorCode PetscStrrstr(const char[],const char[],char *[]);
extern PetscErrorCode PetscStrendswith(const char[],const char[],PetscBool*);
extern PetscErrorCode PetscStrendswithwhich(const char[],const char *const*,PetscInt*);
extern PetscErrorCode PetscStrallocpy(const char[],char *[]);
extern PetscErrorCode PetscStrreplace(MPI_Comm,const char[],char[],size_t);
# 1203 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _p_PetscToken* PetscToken;

extern PetscErrorCode PetscTokenCreate(const char[],const char,PetscToken*);
extern PetscErrorCode PetscTokenFind(PetscToken,char *[]);
extern PetscErrorCode PetscTokenDestroy(PetscToken*);




extern MPI_Op PetscMaxSum_Op;
# 1225 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscMaxSum(MPI_Comm,const PetscInt[],PetscInt*,PetscInt*);

extern PetscErrorCode MPILong_Send(void*,PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm);
extern PetscErrorCode MPILong_Recv(void*,PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm);
# 1239 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _p_PetscObject* PetscObject;
# 1249 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _n_PetscFList *PetscFList;
# 1268 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum {FILE_MODE_READ, FILE_MODE_WRITE, FILE_MODE_APPEND, FILE_MODE_UPDATE, FILE_MODE_APPEND_UPDATE} PetscFileMode;

# 1 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h" 1
# 22 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
typedef struct _p_PetscViewer* PetscViewer;





# 1 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 1
# 29 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h" 2





extern PetscClassId PETSC_VIEWER_CLASSID;
# 45 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"

# 66 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
extern PetscFList PetscViewerList;
extern PetscErrorCode PetscViewerRegisterAll(const char *);
extern PetscErrorCode PetscViewerRegisterDestroy(void);
extern PetscErrorCode PetscViewerInitializePackage(const char[]);

extern PetscErrorCode PetscViewerRegister(const char*,const char*,const char*,PetscErrorCode (*)(PetscViewer));
# 116 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
extern PetscErrorCode PetscViewerCreate(MPI_Comm,PetscViewer*);

extern PetscErrorCode PetscViewerSetFromOptions(PetscViewer);
extern PetscErrorCode PetscViewerASCIIOpenWithFILE(MPI_Comm,FILE*,PetscViewer*);

extern PetscErrorCode PetscViewerASCIIOpen(MPI_Comm,const char[],PetscViewer*);
extern PetscErrorCode PetscViewerASCIISetFILE(PetscViewer,FILE*);
extern PetscErrorCode PetscViewerBinaryOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
extern PetscErrorCode PetscViewerBinaryGetFlowControl(PetscViewer,PetscInt*);
extern PetscErrorCode PetscViewerBinarySetFlowControl(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerBinarySetMPIIO(PetscViewer);
extern PetscErrorCode PetscViewerBinaryGetMPIIO(PetscViewer,PetscBool *);

extern PetscErrorCode PetscViewerBinaryGetMPIIODescriptor(PetscViewer,MPI_File*);
extern PetscErrorCode PetscViewerBinaryGetMPIIOOffset(PetscViewer,MPI_Offset*);
extern PetscErrorCode PetscViewerBinaryAddMPIIOOffset(PetscViewer,MPI_Offset);


extern PetscErrorCode PetscViewerSocketOpen(MPI_Comm,const char[],int,PetscViewer*);
extern PetscErrorCode PetscViewerStringOpen(MPI_Comm,char[],PetscInt,PetscViewer*);
extern PetscErrorCode PetscViewerDrawOpen(MPI_Comm,const char[],const char[],int,int,int,int,PetscViewer*);
extern PetscErrorCode PetscViewerMathematicaOpen(MPI_Comm, int, const char[], const char[], PetscViewer *);
extern PetscErrorCode PetscViewerSiloOpen(MPI_Comm, const char[], PetscViewer *);
extern PetscErrorCode PetscViewerMatlabOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);

extern PetscErrorCode PetscViewerGetType(PetscViewer,const char**);
extern PetscErrorCode PetscViewerSetType(PetscViewer,const char*);
extern PetscErrorCode PetscViewerDestroy(PetscViewer*);
extern PetscErrorCode PetscViewerGetSingleton(PetscViewer,PetscViewer*);
extern PetscErrorCode PetscViewerRestoreSingleton(PetscViewer,PetscViewer*);
extern PetscErrorCode PetscViewerGetSubcomm(PetscViewer,MPI_Comm,PetscViewer*);
extern PetscErrorCode PetscViewerRestoreSubcomm(PetscViewer,MPI_Comm,PetscViewer*);

extern PetscErrorCode PetscViewerSetUp(PetscViewer);
extern PetscErrorCode PetscViewerView(PetscViewer,PetscViewer);

extern PetscErrorCode PetscViewerSetOptionsPrefix(PetscViewer,const char[]);
extern PetscErrorCode PetscViewerAppendOptionsPrefix(PetscViewer,const char[]);
extern PetscErrorCode PetscViewerGetOptionsPrefix(PetscViewer,const char*[]);
# 166 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
typedef enum {
  PETSC_VIEWER_DEFAULT,
  PETSC_VIEWER_ASCII_MATLAB,
  PETSC_VIEWER_ASCII_MATHEMATICA,
  PETSC_VIEWER_ASCII_IMPL,
  PETSC_VIEWER_ASCII_INFO,
  PETSC_VIEWER_ASCII_INFO_DETAIL,
  PETSC_VIEWER_ASCII_COMMON,
  PETSC_VIEWER_ASCII_SYMMODU,
  PETSC_VIEWER_ASCII_INDEX,
  PETSC_VIEWER_ASCII_DENSE,
  PETSC_VIEWER_ASCII_MATRIXMARKET,
  PETSC_VIEWER_ASCII_VTK,
  PETSC_VIEWER_ASCII_VTK_CELL,
  PETSC_VIEWER_ASCII_VTK_COORDS,
  PETSC_VIEWER_ASCII_PCICE,
  PETSC_VIEWER_ASCII_PYTHON,
  PETSC_VIEWER_ASCII_FACTOR_INFO,
  PETSC_VIEWER_DRAW_BASIC,
  PETSC_VIEWER_DRAW_LG,
  PETSC_VIEWER_DRAW_CONTOUR,
  PETSC_VIEWER_DRAW_PORTS,
  PETSC_VIEWER_NATIVE,
  PETSC_VIEWER_NOFORMAT
  } PetscViewerFormat;

extern PetscErrorCode PetscViewerSetFormat(PetscViewer,PetscViewerFormat);
extern PetscErrorCode PetscViewerPushFormat(PetscViewer,PetscViewerFormat);
extern PetscErrorCode PetscViewerPopFormat(PetscViewer);
extern PetscErrorCode PetscViewerGetFormat(PetscViewer,PetscViewerFormat*);
extern PetscErrorCode PetscViewerFlush(PetscViewer);





extern PetscErrorCode PetscViewerASCIIGetPointer(PetscViewer,FILE**);
extern PetscErrorCode PetscViewerFileGetMode(PetscViewer,PetscFileMode*);
extern PetscErrorCode PetscViewerFileSetMode(PetscViewer,PetscFileMode);
extern PetscErrorCode PetscViewerASCIIPrintf(PetscViewer,const char[],...);
extern PetscErrorCode PetscViewerASCIISynchronizedPrintf(PetscViewer,const char[],...);
extern PetscErrorCode PetscViewerASCIISynchronizedAllow(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerASCIIPushTab(PetscViewer);
extern PetscErrorCode PetscViewerASCIIPopTab(PetscViewer);
extern PetscErrorCode PetscViewerASCIIUseTabs(PetscViewer,PetscBool );
extern PetscErrorCode PetscViewerASCIISetTab(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerASCIIAddTab(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerASCIISubtractTab(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerBinaryGetDescriptor(PetscViewer,int*);
extern PetscErrorCode PetscViewerBinaryGetInfoPointer(PetscViewer,FILE **);
extern PetscErrorCode PetscViewerBinaryRead(PetscViewer,void*,PetscInt,PetscDataType);
extern PetscErrorCode PetscViewerBinaryWrite(PetscViewer,void*,PetscInt,PetscDataType,PetscBool );
extern PetscErrorCode PetscViewerStringSPrintf(PetscViewer,const char[],...);
extern PetscErrorCode PetscViewerStringSetString(PetscViewer,char[],PetscInt);
extern PetscErrorCode PetscViewerDrawClear(PetscViewer);
extern PetscErrorCode PetscViewerDrawSetHold(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerDrawGetHold(PetscViewer,PetscBool*);
extern PetscErrorCode PetscViewerDrawSetPause(PetscViewer,PetscReal);
extern PetscErrorCode PetscViewerDrawGetPause(PetscViewer,PetscReal*);
extern PetscErrorCode PetscViewerDrawSetInfo(PetscViewer,const char[],const char[],int,int,int,int);
extern PetscErrorCode PetscViewerDrawResize(PetscViewer,int,int);
extern PetscErrorCode PetscViewerDrawSetBounds(PetscViewer,PetscInt,const PetscReal*);
extern PetscErrorCode PetscViewerDrawGetBounds(PetscViewer,PetscInt*,const PetscReal**);
extern PetscErrorCode PetscViewerSocketSetConnection(PetscViewer,const char[],int);
extern PetscErrorCode PetscViewerBinarySkipInfo(PetscViewer);
extern PetscErrorCode PetscViewerBinarySetSkipOptions(PetscViewer,PetscBool );
extern PetscErrorCode PetscViewerBinaryGetSkipOptions(PetscViewer,PetscBool *);
extern PetscErrorCode PetscViewerBinarySetSkipHeader(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerBinaryGetSkipHeader(PetscViewer,PetscBool*);
extern PetscErrorCode PetscViewerBinaryReadStringArray(PetscViewer,char***);
extern PetscErrorCode PetscViewerBinaryWriteStringArray(PetscViewer,char**);

extern PetscErrorCode PetscViewerFileSetName(PetscViewer,const char[]);
extern PetscErrorCode PetscViewerFileGetName(PetscViewer,const char**);

extern PetscErrorCode PetscPLAPACKInitializePackage(MPI_Comm com);
extern PetscErrorCode PetscPLAPACKFinalizePackage(void);

extern PetscErrorCode PetscViewerVUGetPointer(PetscViewer, FILE**);
extern PetscErrorCode PetscViewerVUSetVecSeen(PetscViewer, PetscBool );
extern PetscErrorCode PetscViewerVUGetVecSeen(PetscViewer, PetscBool *);
extern PetscErrorCode PetscViewerVUPrintDeferred(PetscViewer, const char [], ...);
extern PetscErrorCode PetscViewerVUFlushDeferred(PetscViewer);

extern PetscErrorCode PetscViewerMathematicaInitializePackage(const char[]);
extern PetscErrorCode PetscViewerMathematicaFinalizePackage(void);
extern PetscErrorCode PetscViewerMathematicaGetName(PetscViewer, const char **);
extern PetscErrorCode PetscViewerMathematicaSetName(PetscViewer, const char []);
extern PetscErrorCode PetscViewerMathematicaClearName(PetscViewer);
extern PetscErrorCode PetscViewerMathematicaSkipPackets(PetscViewer, int);

extern PetscErrorCode PetscViewerSiloGetName(PetscViewer, char **);
extern PetscErrorCode PetscViewerSiloSetName(PetscViewer, const char []);
extern PetscErrorCode PetscViewerSiloClearName(PetscViewer);
extern PetscErrorCode PetscViewerSiloGetMeshName(PetscViewer, char **);
extern PetscErrorCode PetscViewerSiloSetMeshName(PetscViewer, const char []);
extern PetscErrorCode PetscViewerSiloClearMeshName(PetscViewer);

extern PetscErrorCode PetscViewerNetcdfOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
extern PetscErrorCode PetscViewerNetcdfGetID(PetscViewer, int *);

extern PetscErrorCode PetscViewerHDF5WriteSDS(PetscViewer,float *,int,int *,int);

extern PetscErrorCode PetscViewerHDF5Open(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
extern PetscErrorCode PetscViewerHDF5PushGroup(PetscViewer,const char *);
extern PetscErrorCode PetscViewerHDF5PopGroup(PetscViewer);
extern PetscErrorCode PetscViewerHDF5GetGroup(PetscViewer, const char **);
extern PetscErrorCode PetscViewerHDF5IncrementTimestep(PetscViewer);
extern PetscErrorCode PetscViewerHDF5SetTimestep(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerHDF5GetTimestep(PetscViewer,PetscInt*);

# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 1
# 24 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 1
# 31 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5pubconf.h" 1
# 32 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2


# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5version.h" 1
# 35 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2





# 1 "/usr/include/sys/types.h" 1 3 4
# 28 "/usr/include/sys/types.h" 3 4






typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;




typedef __loff_t loff_t;



typedef __ino_t ino_t;
# 61 "/usr/include/sys/types.h" 3 4
typedef __dev_t dev_t;




typedef __gid_t gid_t;




typedef __mode_t mode_t;




typedef __nlink_t nlink_t;




typedef __uid_t uid_t;
# 99 "/usr/include/sys/types.h" 3 4
typedef __pid_t pid_t;





typedef __id_t id_t;
# 116 "/usr/include/sys/types.h" 3 4
typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;





typedef __key_t key_t;
# 133 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/time.h" 1 3 4
# 58 "/usr/include/time.h" 3 4


typedef __clock_t clock_t;



# 74 "/usr/include/time.h" 3 4


typedef __time_t time_t;



# 92 "/usr/include/time.h" 3 4
typedef __clockid_t clockid_t;
# 104 "/usr/include/time.h" 3 4
typedef __timer_t timer_t;
# 134 "/usr/include/sys/types.h" 2 3 4
# 147 "/usr/include/sys/types.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 148 "/usr/include/sys/types.h" 2 3 4



typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
# 195 "/usr/include/sys/types.h" 3 4
typedef int int8_t __attribute__ ((__mode__ (__QI__)));
typedef int int16_t __attribute__ ((__mode__ (__HI__)));
typedef int int32_t __attribute__ ((__mode__ (__SI__)));
typedef int int64_t __attribute__ ((__mode__ (__DI__)));


typedef unsigned int u_int8_t __attribute__ ((__mode__ (__QI__)));
typedef unsigned int u_int16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int u_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int u_int64_t __attribute__ ((__mode__ (__DI__)));

typedef int register_t __attribute__ ((__mode__ (__word__)));
# 217 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/endian.h" 1 3 4
# 37 "/usr/include/endian.h" 3 4
# 1 "/usr/include/bits/endian.h" 1 3 4
# 38 "/usr/include/endian.h" 2 3 4
# 61 "/usr/include/endian.h" 3 4
# 1 "/usr/include/bits/byteswap.h" 1 3 4
# 28 "/usr/include/bits/byteswap.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 29 "/usr/include/bits/byteswap.h" 2 3 4
# 62 "/usr/include/endian.h" 2 3 4
# 218 "/usr/include/sys/types.h" 2 3 4


# 1 "/usr/include/sys/select.h" 1 3 4
# 31 "/usr/include/sys/select.h" 3 4
# 1 "/usr/include/bits/select.h" 1 3 4
# 23 "/usr/include/bits/select.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 24 "/usr/include/bits/select.h" 2 3 4
# 32 "/usr/include/sys/select.h" 2 3 4


# 1 "/usr/include/bits/sigset.h" 1 3 4
# 24 "/usr/include/bits/sigset.h" 3 4
typedef int __sig_atomic_t;




typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;
# 35 "/usr/include/sys/select.h" 2 3 4



typedef __sigset_t sigset_t;





# 1 "/usr/include/time.h" 1 3 4
# 120 "/usr/include/time.h" 3 4
struct timespec
  {
    __time_t tv_sec;
    long int tv_nsec;
  };
# 45 "/usr/include/sys/select.h" 2 3 4

# 1 "/usr/include/bits/time.h" 1 3 4
# 75 "/usr/include/bits/time.h" 3 4
struct timeval
  {
    __time_t tv_sec;
    __suseconds_t tv_usec;
  };
# 47 "/usr/include/sys/select.h" 2 3 4


typedef __suseconds_t suseconds_t;





typedef long int __fd_mask;
# 67 "/usr/include/sys/select.h" 3 4
typedef struct
  {






    __fd_mask __fds_bits[1024 / (8 * (int) sizeof (__fd_mask))];


  } fd_set;






typedef __fd_mask fd_mask;
# 99 "/usr/include/sys/select.h" 3 4

# 109 "/usr/include/sys/select.h" 3 4
extern int select (int __nfds, fd_set *__restrict __readfds,
     fd_set *__restrict __writefds,
     fd_set *__restrict __exceptfds,
     struct timeval *__restrict __timeout);
# 121 "/usr/include/sys/select.h" 3 4
extern int pselect (int __nfds, fd_set *__restrict __readfds,
      fd_set *__restrict __writefds,
      fd_set *__restrict __exceptfds,
      const struct timespec *__restrict __timeout,
      const __sigset_t *__restrict __sigmask);



# 221 "/usr/include/sys/types.h" 2 3 4


# 1 "/usr/include/sys/sysmacros.h" 1 3 4
# 30 "/usr/include/sys/sysmacros.h" 3 4
__extension__
extern unsigned int gnu_dev_major (unsigned long long int __dev)
     __attribute__ ((__nothrow__));
__extension__
extern unsigned int gnu_dev_minor (unsigned long long int __dev)
     __attribute__ ((__nothrow__));
__extension__
extern unsigned long long int gnu_dev_makedev (unsigned int __major,
            unsigned int __minor)
     __attribute__ ((__nothrow__));
# 224 "/usr/include/sys/types.h" 2 3 4





typedef __blksize_t blksize_t;






typedef __blkcnt_t blkcnt_t;



typedef __fsblkcnt_t fsblkcnt_t;



typedef __fsfilcnt_t fsfilcnt_t;
# 271 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/bits/pthreadtypes.h" 1 3 4
# 23 "/usr/include/bits/pthreadtypes.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 24 "/usr/include/bits/pthreadtypes.h" 2 3 4
# 50 "/usr/include/bits/pthreadtypes.h" 3 4
typedef unsigned long int pthread_t;


typedef union
{
  char __size[56];
  long int __align;
} pthread_attr_t;



typedef struct __pthread_internal_list
{
  struct __pthread_internal_list *__prev;
  struct __pthread_internal_list *__next;
} __pthread_list_t;
# 76 "/usr/include/bits/pthreadtypes.h" 3 4
typedef union
{
  struct __pthread_mutex_s
  {
    int __lock;
    unsigned int __count;
    int __owner;

    unsigned int __nusers;



    int __kind;

    int __spins;
    __pthread_list_t __list;
# 101 "/usr/include/bits/pthreadtypes.h" 3 4
  } __data;
  char __size[40];
  long int __align;
} pthread_mutex_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_mutexattr_t;




typedef union
{
  struct
  {
    int __lock;
    unsigned int __futex;
    __extension__ unsigned long long int __total_seq;
    __extension__ unsigned long long int __wakeup_seq;
    __extension__ unsigned long long int __woken_seq;
    void *__mutex;
    unsigned int __nwaiters;
    unsigned int __broadcast_seq;
  } __data;
  char __size[48];
  __extension__ long long int __align;
} pthread_cond_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_condattr_t;



typedef unsigned int pthread_key_t;



typedef int pthread_once_t;





typedef union
{

  struct
  {
    int __lock;
    unsigned int __nr_readers;
    unsigned int __readers_wakeup;
    unsigned int __writer_wakeup;
    unsigned int __nr_readers_queued;
    unsigned int __nr_writers_queued;
    int __writer;
    int __shared;
    unsigned long int __pad1;
    unsigned long int __pad2;


    unsigned int __flags;
  } __data;
# 187 "/usr/include/bits/pthreadtypes.h" 3 4
  char __size[56];
  long int __align;
} pthread_rwlock_t;

typedef union
{
  char __size[8];
  long int __align;
} pthread_rwlockattr_t;





typedef volatile int pthread_spinlock_t;




typedef union
{
  char __size[32];
  long int __align;
} pthread_barrier_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_barrierattr_t;
# 272 "/usr/include/sys/types.h" 2 3 4



# 41 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2


# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 1 3 4
# 11 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/syslimits.h" 1 3 4






# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 1 3 4
# 122 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 3 4
# 1 "/usr/include/limits.h" 1 3 4
# 145 "/usr/include/limits.h" 3 4
# 1 "/usr/include/bits/posix1_lim.h" 1 3 4
# 157 "/usr/include/bits/posix1_lim.h" 3 4
# 1 "/usr/include/bits/local_lim.h" 1 3 4
# 39 "/usr/include/bits/local_lim.h" 3 4
# 1 "/usr/include/linux/limits.h" 1 3 4
# 40 "/usr/include/bits/local_lim.h" 2 3 4
# 158 "/usr/include/bits/posix1_lim.h" 2 3 4
# 146 "/usr/include/limits.h" 2 3 4



# 1 "/usr/include/bits/posix2_lim.h" 1 3 4
# 150 "/usr/include/limits.h" 2 3 4
# 123 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 2 3 4
# 8 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/syslimits.h" 2 3 4
# 12 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 2 3 4
# 44 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2



# 1 "/usr/include/stdint.h" 1 3 4
# 27 "/usr/include/stdint.h" 3 4
# 1 "/usr/include/bits/wchar.h" 1 3 4
# 28 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 29 "/usr/include/stdint.h" 2 3 4
# 49 "/usr/include/stdint.h" 3 4
typedef unsigned char uint8_t;
typedef unsigned short int uint16_t;

typedef unsigned int uint32_t;



typedef unsigned long int uint64_t;
# 66 "/usr/include/stdint.h" 3 4
typedef signed char int_least8_t;
typedef short int int_least16_t;
typedef int int_least32_t;

typedef long int int_least64_t;






typedef unsigned char uint_least8_t;
typedef unsigned short int uint_least16_t;
typedef unsigned int uint_least32_t;

typedef unsigned long int uint_least64_t;
# 91 "/usr/include/stdint.h" 3 4
typedef signed char int_fast8_t;

typedef long int int_fast16_t;
typedef long int int_fast32_t;
typedef long int int_fast64_t;
# 104 "/usr/include/stdint.h" 3 4
typedef unsigned char uint_fast8_t;

typedef unsigned long int uint_fast16_t;
typedef unsigned long int uint_fast32_t;
typedef unsigned long int uint_fast64_t;
# 120 "/usr/include/stdint.h" 3 4
typedef long int intptr_t;


typedef unsigned long int uintptr_t;
# 135 "/usr/include/stdint.h" 3 4
typedef long int intmax_t;
typedef unsigned long int uintmax_t;
# 48 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2







# 1 "/usr/include/inttypes.h" 1 3 4
# 35 "/usr/include/inttypes.h" 3 4
typedef int __gwchar_t;
# 274 "/usr/include/inttypes.h" 3 4





typedef struct
  {
    long int quot;
    long int rem;
  } imaxdiv_t;
# 298 "/usr/include/inttypes.h" 3 4
extern intmax_t imaxabs (intmax_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern imaxdiv_t imaxdiv (intmax_t __numer, intmax_t __denom)
      __attribute__ ((__nothrow__)) __attribute__ ((__const__));


extern intmax_t strtoimax (__const char *__restrict __nptr,
      char **__restrict __endptr, int __base) __attribute__ ((__nothrow__));


extern uintmax_t strtoumax (__const char *__restrict __nptr,
       char ** __restrict __endptr, int __base) __attribute__ ((__nothrow__));


extern intmax_t wcstoimax (__const __gwchar_t *__restrict __nptr,
      __gwchar_t **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__));


extern uintmax_t wcstoumax (__const __gwchar_t *__restrict __nptr,
       __gwchar_t ** __restrict __endptr, int __base)
     __attribute__ ((__nothrow__));
# 442 "/usr/include/inttypes.h" 3 4

# 56 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2


# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 149 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 3 4
typedef long int ptrdiff_t;
# 323 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 3 4
typedef int wchar_t;
# 59 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2
# 69 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5api_adpt.h" 1
# 70 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h" 2
# 96 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
typedef int herr_t;
# 114 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
typedef unsigned int hbool_t;
typedef int htri_t;
# 140 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
typedef unsigned long long hsize_t;
typedef signed long long hssize_t;
# 152 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
    typedef uint64_t haddr_t;
# 257 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
typedef enum {
    H5_ITER_UNKNOWN = -1,
    H5_ITER_INC,
    H5_ITER_DEC,
    H5_ITER_NATIVE,
    H5_ITER_N
} H5_iter_order_t;
# 278 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5public.h"
typedef enum H5_index_t {
    H5_INDEX_UNKNOWN = -1,
    H5_INDEX_NAME,
    H5_INDEX_CRT_ORDER,
    H5_INDEX_N
} H5_index_t;




typedef struct H5_ih_info_t {
    hsize_t index_size;
    hsize_t heap_size;
} H5_ih_info_t;


 herr_t H5open(void);
 herr_t H5close(void);
 herr_t H5dont_atexit(void);
 herr_t H5garbage_collect(void);
 herr_t H5set_free_list_limits (int reg_global_lim, int reg_list_lim,
                int arr_global_lim, int arr_list_lim, int blk_global_lim,
                int blk_list_lim);
 herr_t H5get_libversion(unsigned *majnum, unsigned *minnum,
    unsigned *relnum);
 herr_t H5check_version(unsigned majnum, unsigned minnum,
          unsigned relnum);
# 25 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Apublic.h" 1
# 23 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Apublic.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ipublic.h" 1
# 36 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ipublic.h"
typedef enum H5I_type_t {
    H5I_UNINIT = (-2),
    H5I_BADID = (-1),
    H5I_FILE = 1,
    H5I_GROUP,
    H5I_DATATYPE,
    H5I_DATASPACE,
    H5I_DATASET,
    H5I_ATTR,
    H5I_REFERENCE,
    H5I_VFL,
    H5I_GENPROP_CLS,
    H5I_GENPROP_LST,
    H5I_ERROR_CLASS,
    H5I_ERROR_MSG,
    H5I_ERROR_STACK,
    H5I_NTYPES
} H5I_type_t;


typedef int hid_t;
# 69 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ipublic.h"
typedef herr_t (*H5I_free_t)(void*);


typedef int (*H5I_search_func_t)(void *obj, hid_t id, void *key);







 hid_t H5Iregister(H5I_type_t type, const void *object);
 void *H5Iobject_verify(hid_t id, H5I_type_t id_type);
 void *H5Iremove_verify(hid_t id, H5I_type_t id_type);
 H5I_type_t H5Iget_type(hid_t id);
 hid_t H5Iget_file_id(hid_t id);
 ssize_t H5Iget_name(hid_t id, char *name , size_t size);
 int H5Iinc_ref(hid_t id);
 int H5Idec_ref(hid_t id);
 int H5Iget_ref(hid_t id);
 H5I_type_t H5Iregister_type(size_t hash_size, unsigned reserved, H5I_free_t free_func);
 herr_t H5Iclear_type(H5I_type_t type, hbool_t force);
 herr_t H5Idestroy_type(H5I_type_t type);
 int H5Iinc_type_ref(H5I_type_t type);
 int H5Idec_type_ref(H5I_type_t type);
 int H5Iget_type_ref(H5I_type_t type);
 void *H5Isearch(H5I_type_t type, H5I_search_func_t func, void *key);
 herr_t H5Inmembers(H5I_type_t type, hsize_t *num_members);
 htri_t H5Itype_exists(H5I_type_t type);
 htri_t H5Iis_valid(hid_t id);
# 24 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Apublic.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Opublic.h" 1
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Opublic.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Lpublic.h" 1
# 32 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Lpublic.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h" 1
# 30 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
typedef enum H5T_class_t {
    H5T_NO_CLASS = -1,
    H5T_INTEGER = 0,
    H5T_FLOAT = 1,
    H5T_TIME = 2,
    H5T_STRING = 3,
    H5T_BITFIELD = 4,
    H5T_OPAQUE = 5,
    H5T_COMPOUND = 6,
    H5T_REFERENCE = 7,
    H5T_ENUM = 8,
    H5T_VLEN = 9,
    H5T_ARRAY = 10,

    H5T_NCLASSES
} H5T_class_t;


typedef enum H5T_order_t {
    H5T_ORDER_ERROR = -1,
    H5T_ORDER_LE = 0,
    H5T_ORDER_BE = 1,
    H5T_ORDER_VAX = 2,
    H5T_ORDER_MIXED = 3,
    H5T_ORDER_NONE = 4

} H5T_order_t;


typedef enum H5T_sign_t {
    H5T_SGN_ERROR = -1,
    H5T_SGN_NONE = 0,
    H5T_SGN_2 = 1,

    H5T_NSGN = 2
} H5T_sign_t;


typedef enum H5T_norm_t {
    H5T_NORM_ERROR = -1,
    H5T_NORM_IMPLIED = 0,
    H5T_NORM_MSBSET = 1,
    H5T_NORM_NONE = 2

} H5T_norm_t;





typedef enum H5T_cset_t {
    H5T_CSET_ERROR = -1,
    H5T_CSET_ASCII = 0,
    H5T_CSET_UTF8 = 1,
    H5T_CSET_RESERVED_2 = 2,
    H5T_CSET_RESERVED_3 = 3,
    H5T_CSET_RESERVED_4 = 4,
    H5T_CSET_RESERVED_5 = 5,
    H5T_CSET_RESERVED_6 = 6,
    H5T_CSET_RESERVED_7 = 7,
    H5T_CSET_RESERVED_8 = 8,
    H5T_CSET_RESERVED_9 = 9,
    H5T_CSET_RESERVED_10 = 10,
    H5T_CSET_RESERVED_11 = 11,
    H5T_CSET_RESERVED_12 = 12,
    H5T_CSET_RESERVED_13 = 13,
    H5T_CSET_RESERVED_14 = 14,
    H5T_CSET_RESERVED_15 = 15
} H5T_cset_t;






typedef enum H5T_str_t {
    H5T_STR_ERROR = -1,
    H5T_STR_NULLTERM = 0,
    H5T_STR_NULLPAD = 1,
    H5T_STR_SPACEPAD = 2,
    H5T_STR_RESERVED_3 = 3,
    H5T_STR_RESERVED_4 = 4,
    H5T_STR_RESERVED_5 = 5,
    H5T_STR_RESERVED_6 = 6,
    H5T_STR_RESERVED_7 = 7,
    H5T_STR_RESERVED_8 = 8,
    H5T_STR_RESERVED_9 = 9,
    H5T_STR_RESERVED_10 = 10,
    H5T_STR_RESERVED_11 = 11,
    H5T_STR_RESERVED_12 = 12,
    H5T_STR_RESERVED_13 = 13,
    H5T_STR_RESERVED_14 = 14,
    H5T_STR_RESERVED_15 = 15
} H5T_str_t;



typedef enum H5T_pad_t {
    H5T_PAD_ERROR = -1,
    H5T_PAD_ZERO = 0,
    H5T_PAD_ONE = 1,
    H5T_PAD_BACKGROUND = 2,

    H5T_NPAD = 3
} H5T_pad_t;


typedef enum H5T_cmd_t {
    H5T_CONV_INIT = 0,
    H5T_CONV_CONV = 1,
    H5T_CONV_FREE = 2
} H5T_cmd_t;


typedef enum H5T_bkg_t {
    H5T_BKG_NO = 0,
    H5T_BKG_TEMP = 1,
    H5T_BKG_YES = 2
} H5T_bkg_t;


typedef struct H5T_cdata_t {
    H5T_cmd_t command;
    H5T_bkg_t need_bkg;
    hbool_t recalc;
    void *priv;
} H5T_cdata_t;


typedef enum H5T_pers_t {
    H5T_PERS_DONTCARE = -1,
    H5T_PERS_HARD = 0,
    H5T_PERS_SOFT = 1
} H5T_pers_t;


typedef enum H5T_direction_t {
    H5T_DIR_DEFAULT = 0,
    H5T_DIR_ASCEND = 1,
    H5T_DIR_DESCEND = 2
} H5T_direction_t;


typedef enum H5T_conv_except_t {
    H5T_CONV_EXCEPT_RANGE_HI = 0,
    H5T_CONV_EXCEPT_RANGE_LOW = 1,
    H5T_CONV_EXCEPT_PRECISION = 2,
    H5T_CONV_EXCEPT_TRUNCATE = 3,
    H5T_CONV_EXCEPT_PINF = 4,
    H5T_CONV_EXCEPT_NINF = 5,
    H5T_CONV_EXCEPT_NAN = 6
} H5T_conv_except_t;


typedef enum H5T_conv_ret_t {
    H5T_CONV_ABORT = -1,
    H5T_CONV_UNHANDLED = 0,
    H5T_CONV_HANDLED = 1
} H5T_conv_ret_t;



typedef struct {
    size_t len;
    void *p;
} hvl_t;
# 209 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
typedef herr_t (*H5T_conv_t) (hid_t src_id, hid_t dst_id, H5T_cdata_t *cdata,
      size_t nelmts, size_t buf_stride, size_t bkg_stride, void *buf,
      void *bkg, hid_t dset_xfer_plist);




typedef H5T_conv_ret_t (*H5T_conv_except_func_t)(H5T_conv_except_t except_type,
    hid_t src_id, hid_t dst_id, void *src_buf, void *dst_buf, void *user_data);
# 234 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
extern hid_t H5T_IEEE_F32BE_g;
extern hid_t H5T_IEEE_F32LE_g;
extern hid_t H5T_IEEE_F64BE_g;
extern hid_t H5T_IEEE_F64LE_g;
# 269 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
extern hid_t H5T_STD_I8BE_g;
extern hid_t H5T_STD_I8LE_g;
extern hid_t H5T_STD_I16BE_g;
extern hid_t H5T_STD_I16LE_g;
extern hid_t H5T_STD_I32BE_g;
extern hid_t H5T_STD_I32LE_g;
extern hid_t H5T_STD_I64BE_g;
extern hid_t H5T_STD_I64LE_g;
extern hid_t H5T_STD_U8BE_g;
extern hid_t H5T_STD_U8LE_g;
extern hid_t H5T_STD_U16BE_g;
extern hid_t H5T_STD_U16LE_g;
extern hid_t H5T_STD_U32BE_g;
extern hid_t H5T_STD_U32LE_g;
extern hid_t H5T_STD_U64BE_g;
extern hid_t H5T_STD_U64LE_g;
extern hid_t H5T_STD_B8BE_g;
extern hid_t H5T_STD_B8LE_g;
extern hid_t H5T_STD_B16BE_g;
extern hid_t H5T_STD_B16LE_g;
extern hid_t H5T_STD_B32BE_g;
extern hid_t H5T_STD_B32LE_g;
extern hid_t H5T_STD_B64BE_g;
extern hid_t H5T_STD_B64LE_g;
extern hid_t H5T_STD_REF_OBJ_g;
extern hid_t H5T_STD_REF_DSETREG_g;
# 303 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
extern hid_t H5T_UNIX_D32BE_g;
extern hid_t H5T_UNIX_D32LE_g;
extern hid_t H5T_UNIX_D64BE_g;
extern hid_t H5T_UNIX_D64LE_g;






extern hid_t H5T_C_S1_g;





extern hid_t H5T_FORTRAN_S1_g;
# 383 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
extern hid_t H5T_VAX_F32_g;
extern hid_t H5T_VAX_F64_g;
# 421 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
extern hid_t H5T_NATIVE_SCHAR_g;
extern hid_t H5T_NATIVE_UCHAR_g;
extern hid_t H5T_NATIVE_SHORT_g;
extern hid_t H5T_NATIVE_USHORT_g;
extern hid_t H5T_NATIVE_INT_g;
extern hid_t H5T_NATIVE_UINT_g;
extern hid_t H5T_NATIVE_LONG_g;
extern hid_t H5T_NATIVE_ULONG_g;
extern hid_t H5T_NATIVE_LLONG_g;
extern hid_t H5T_NATIVE_ULLONG_g;
extern hid_t H5T_NATIVE_FLOAT_g;
extern hid_t H5T_NATIVE_DOUBLE_g;

extern hid_t H5T_NATIVE_LDOUBLE_g;

extern hid_t H5T_NATIVE_B8_g;
extern hid_t H5T_NATIVE_B16_g;
extern hid_t H5T_NATIVE_B32_g;
extern hid_t H5T_NATIVE_B64_g;
extern hid_t H5T_NATIVE_OPAQUE_g;
extern hid_t H5T_NATIVE_HADDR_g;
extern hid_t H5T_NATIVE_HSIZE_g;
extern hid_t H5T_NATIVE_HSSIZE_g;
extern hid_t H5T_NATIVE_HERR_g;
extern hid_t H5T_NATIVE_HBOOL_g;
# 454 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
extern hid_t H5T_NATIVE_INT8_g;
extern hid_t H5T_NATIVE_UINT8_g;
extern hid_t H5T_NATIVE_INT_LEAST8_g;
extern hid_t H5T_NATIVE_UINT_LEAST8_g;
extern hid_t H5T_NATIVE_INT_FAST8_g;
extern hid_t H5T_NATIVE_UINT_FAST8_g;







extern hid_t H5T_NATIVE_INT16_g;
extern hid_t H5T_NATIVE_UINT16_g;
extern hid_t H5T_NATIVE_INT_LEAST16_g;
extern hid_t H5T_NATIVE_UINT_LEAST16_g;
extern hid_t H5T_NATIVE_INT_FAST16_g;
extern hid_t H5T_NATIVE_UINT_FAST16_g;







extern hid_t H5T_NATIVE_INT32_g;
extern hid_t H5T_NATIVE_UINT32_g;
extern hid_t H5T_NATIVE_INT_LEAST32_g;
extern hid_t H5T_NATIVE_UINT_LEAST32_g;
extern hid_t H5T_NATIVE_INT_FAST32_g;
extern hid_t H5T_NATIVE_UINT_FAST32_g;







extern hid_t H5T_NATIVE_INT64_g;
extern hid_t H5T_NATIVE_UINT64_g;
extern hid_t H5T_NATIVE_INT_LEAST64_g;
extern hid_t H5T_NATIVE_UINT_LEAST64_g;
extern hid_t H5T_NATIVE_INT_FAST64_g;
extern hid_t H5T_NATIVE_UINT_FAST64_g;


 hid_t H5Tcreate(H5T_class_t type, size_t size);
 hid_t H5Tcopy(hid_t type_id);
 herr_t H5Tclose(hid_t type_id);
 htri_t H5Tequal(hid_t type1_id, hid_t type2_id);
 herr_t H5Tlock(hid_t type_id);
 herr_t H5Tcommit2(hid_t loc_id, const char *name, hid_t type_id,
    hid_t lcpl_id, hid_t tcpl_id, hid_t tapl_id);
 hid_t H5Topen2(hid_t loc_id, const char *name, hid_t tapl_id);
 herr_t H5Tcommit_anon(hid_t loc_id, hid_t type_id, hid_t tcpl_id, hid_t tapl_id);
 hid_t H5Tget_create_plist(hid_t type_id);
 htri_t H5Tcommitted(hid_t type_id);
 herr_t H5Tencode(hid_t obj_id, void *buf, size_t *nalloc);
 hid_t H5Tdecode(const void *buf);


 herr_t H5Tinsert(hid_t parent_id, const char *name, size_t offset,
    hid_t member_id);
 herr_t H5Tpack(hid_t type_id);


 hid_t H5Tenum_create(hid_t base_id);
 herr_t H5Tenum_insert(hid_t type, const char *name, const void *value);
 herr_t H5Tenum_nameof(hid_t type, const void *value, char *name ,
        size_t size);
 herr_t H5Tenum_valueof(hid_t type, const char *name,
         void *value );


 hid_t H5Tvlen_create(hid_t base_id);


 hid_t H5Tarray_create2(hid_t base_id, unsigned ndims,
            const hsize_t dim[ ]);
 int H5Tget_array_ndims(hid_t type_id);
 int H5Tget_array_dims2(hid_t type_id, hsize_t dims[]);


 herr_t H5Tset_tag(hid_t type, const char *tag);
 char *H5Tget_tag(hid_t type);


 hid_t H5Tget_super(hid_t type);
 H5T_class_t H5Tget_class(hid_t type_id);
 htri_t H5Tdetect_class(hid_t type_id, H5T_class_t cls);
 size_t H5Tget_size(hid_t type_id);
 H5T_order_t H5Tget_order(hid_t type_id);
 size_t H5Tget_precision(hid_t type_id);
 int H5Tget_offset(hid_t type_id);
 herr_t H5Tget_pad(hid_t type_id, H5T_pad_t *lsb ,
     H5T_pad_t *msb );
 H5T_sign_t H5Tget_sign(hid_t type_id);
 herr_t H5Tget_fields(hid_t type_id, size_t *spos ,
        size_t *epos , size_t *esize ,
        size_t *mpos , size_t *msize );
 size_t H5Tget_ebias(hid_t type_id);
 H5T_norm_t H5Tget_norm(hid_t type_id);
 H5T_pad_t H5Tget_inpad(hid_t type_id);
 H5T_str_t H5Tget_strpad(hid_t type_id);
 int H5Tget_nmembers(hid_t type_id);
 char *H5Tget_member_name(hid_t type_id, unsigned membno);
 int H5Tget_member_index(hid_t type_id, const char *name);
 size_t H5Tget_member_offset(hid_t type_id, unsigned membno);
 H5T_class_t H5Tget_member_class(hid_t type_id, unsigned membno);
 hid_t H5Tget_member_type(hid_t type_id, unsigned membno);
 herr_t H5Tget_member_value(hid_t type_id, unsigned membno, void *value );
 H5T_cset_t H5Tget_cset(hid_t type_id);
 htri_t H5Tis_variable_str(hid_t type_id);
 hid_t H5Tget_native_type(hid_t type_id, H5T_direction_t direction);


 herr_t H5Tset_size(hid_t type_id, size_t size);
 herr_t H5Tset_order(hid_t type_id, H5T_order_t order);
 herr_t H5Tset_precision(hid_t type_id, size_t prec);
 herr_t H5Tset_offset(hid_t type_id, size_t offset);
 herr_t H5Tset_pad(hid_t type_id, H5T_pad_t lsb, H5T_pad_t msb);
 herr_t H5Tset_sign(hid_t type_id, H5T_sign_t sign);
 herr_t H5Tset_fields(hid_t type_id, size_t spos, size_t epos,
        size_t esize, size_t mpos, size_t msize);
 herr_t H5Tset_ebias(hid_t type_id, size_t ebias);
 herr_t H5Tset_norm(hid_t type_id, H5T_norm_t norm);
 herr_t H5Tset_inpad(hid_t type_id, H5T_pad_t pad);
 herr_t H5Tset_cset(hid_t type_id, H5T_cset_t cset);
 herr_t H5Tset_strpad(hid_t type_id, H5T_str_t strpad);


 herr_t H5Tregister(H5T_pers_t pers, const char *name, hid_t src_id,
      hid_t dst_id, H5T_conv_t func);
 herr_t H5Tunregister(H5T_pers_t pers, const char *name, hid_t src_id,
        hid_t dst_id, H5T_conv_t func);
 H5T_conv_t H5Tfind(hid_t src_id, hid_t dst_id, H5T_cdata_t **pcdata);
 htri_t H5Tcompiler_conv(hid_t src_id, hid_t dst_id);
 herr_t H5Tconvert(hid_t src_id, hid_t dst_id, size_t nelmts,
     void *buf, void *background, hid_t plist_id);
# 608 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Tpublic.h"
 herr_t H5Tcommit1(hid_t loc_id, const char *name, hid_t type_id);
 hid_t H5Topen1(hid_t loc_id, const char *name);
 hid_t H5Tarray_create1(hid_t base_id, int ndims,
            const hsize_t dim[ ],
            const int perm[ ]);
 int H5Tget_array_dims1(hid_t type_id, hsize_t dims[], int perm[]);
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Lpublic.h" 2
# 64 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Lpublic.h"
typedef enum {
    H5L_TYPE_ERROR = (-1),
    H5L_TYPE_HARD = 0,
    H5L_TYPE_SOFT = 1,
    H5L_TYPE_EXTERNAL = 64,
    H5L_TYPE_MAX = 255
} H5L_type_t;




typedef struct {
    H5L_type_t type;
    hbool_t corder_valid;
    int64_t corder;
    H5T_cset_t cset;
    union {
        haddr_t address;
        size_t val_size;
    } u;
} H5L_info_t;







typedef herr_t (*H5L_create_func_t)(const char *link_name, hid_t loc_group,
    const void *lnkdata, size_t lnkdata_size, hid_t lcpl_id);


typedef herr_t (*H5L_move_func_t)(const char *new_name, hid_t new_loc,
    const void *lnkdata, size_t lnkdata_size);


typedef herr_t (*H5L_copy_func_t)(const char *new_name, hid_t new_loc,
    const void *lnkdata, size_t lnkdata_size);


typedef herr_t (*H5L_traverse_func_t)(const char *link_name, hid_t cur_group,
    const void *lnkdata, size_t lnkdata_size, hid_t lapl_id);


typedef herr_t (*H5L_delete_func_t)(const char *link_name, hid_t file,
    const void *lnkdata, size_t lnkdata_size);



typedef ssize_t (*H5L_query_func_t)(const char *link_name, const void *lnkdata,
    size_t lnkdata_size, void *buf , size_t buf_size);


typedef struct {
    int version;
    H5L_type_t id;
    const char *comment;
    H5L_create_func_t create_func;
    H5L_move_func_t move_func;
    H5L_copy_func_t copy_func;
    H5L_traverse_func_t trav_func;
    H5L_delete_func_t del_func;
    H5L_query_func_t query_func;
} H5L_class_t;


typedef herr_t (*H5L_iterate_t)(hid_t group, const char *name, const H5L_info_t *info,
    void *op_data);


typedef herr_t (*H5L_elink_traverse_t)(const char *parent_file_name,
    const char *parent_group_name, const char *child_file_name,
    const char *child_object_name, unsigned *acc_flags, hid_t fapl_id,
    void *op_data);
# 148 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Lpublic.h"
 herr_t H5Lmove(hid_t src_loc, const char *src_name, hid_t dst_loc,
    const char *dst_name, hid_t lcpl_id, hid_t lapl_id);
 herr_t H5Lcopy(hid_t src_loc, const char *src_name, hid_t dst_loc,
    const char *dst_name, hid_t lcpl_id, hid_t lapl_id);
 herr_t H5Lcreate_hard(hid_t cur_loc, const char *cur_name,
    hid_t dst_loc, const char *dst_name, hid_t lcpl_id, hid_t lapl_id);
 herr_t H5Lcreate_soft(const char *link_target, hid_t link_loc_id,
    const char *link_name, hid_t lcpl_id, hid_t lapl_id);
 herr_t H5Ldelete(hid_t loc_id, const char *name, hid_t lapl_id);
 herr_t H5Ldelete_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n, hid_t lapl_id);
 herr_t H5Lget_val(hid_t loc_id, const char *name, void *buf ,
    size_t size, hid_t lapl_id);
 herr_t H5Lget_val_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n,
    void *buf , size_t size, hid_t lapl_id);
 htri_t H5Lexists(hid_t loc_id, const char *name, hid_t lapl_id);
 herr_t H5Lget_info(hid_t loc_id, const char *name,
    H5L_info_t *linfo , hid_t lapl_id);
 herr_t H5Lget_info_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n,
    H5L_info_t *linfo , hid_t lapl_id);
 ssize_t H5Lget_name_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n,
    char *name , size_t size, hid_t lapl_id);
 herr_t H5Literate(hid_t grp_id, H5_index_t idx_type,
    H5_iter_order_t order, hsize_t *idx, H5L_iterate_t op, void *op_data);
 herr_t H5Literate_by_name(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t *idx,
    H5L_iterate_t op, void *op_data, hid_t lapl_id);
 herr_t H5Lvisit(hid_t grp_id, H5_index_t idx_type, H5_iter_order_t order,
    H5L_iterate_t op, void *op_data);
 herr_t H5Lvisit_by_name(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, H5L_iterate_t op,
    void *op_data, hid_t lapl_id);


 herr_t H5Lcreate_ud(hid_t link_loc_id, const char *link_name,
    H5L_type_t link_type, const void *udata, size_t udata_size, hid_t lcpl_id,
    hid_t lapl_id);
 herr_t H5Lregister(const H5L_class_t *cls);
 herr_t H5Lunregister(H5L_type_t id);
 htri_t H5Lis_registered(H5L_type_t id);


 herr_t H5Lunpack_elink_val(const void *ext_linkval , size_t link_size,
   unsigned *flags, const char **filename , const char **obj_path );
 herr_t H5Lcreate_external(const char *file_name, const char *obj_name,
    hid_t link_loc_id, const char *link_name, hid_t lcpl_id, hid_t lapl_id);
# 34 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Opublic.h" 2
# 82 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Opublic.h"
typedef enum H5O_type_t {
    H5O_TYPE_UNKNOWN = -1,
    H5O_TYPE_GROUP,
    H5O_TYPE_DATASET,
    H5O_TYPE_NAMED_DATATYPE,
    H5O_TYPE_NTYPES
} H5O_type_t;


typedef struct H5O_hdr_info_t {
    unsigned version;
    unsigned nmesgs;
    unsigned nchunks;
    unsigned flags;
    struct {
        hsize_t total;
        hsize_t meta;
        hsize_t mesg;
        hsize_t free;
    } space;
    struct {
        uint64_t present;
        uint64_t shared;
    } mesg;
} H5O_hdr_info_t;


typedef struct H5O_info_t {
    unsigned long fileno;
    haddr_t addr;
    H5O_type_t type;
    unsigned rc;
    time_t atime;
    time_t mtime;
    time_t ctime;
    time_t btime;
    hsize_t num_attrs;
    H5O_hdr_info_t hdr;

    struct {
        H5_ih_info_t obj;
        H5_ih_info_t attr;
    } meta_size;
} H5O_info_t;


typedef uint32_t H5O_msg_crt_idx_t;


typedef herr_t (*H5O_iterate_t)(hid_t obj, const char *name, const H5O_info_t *info,
    void *op_data);
# 147 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Opublic.h"
 hid_t H5Oopen(hid_t loc_id, const char *name, hid_t lapl_id);
 hid_t H5Oopen_by_addr(hid_t loc_id, haddr_t addr);
 hid_t H5Oopen_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n, hid_t lapl_id);
 htri_t H5Oexists_by_name(hid_t loc_id, const char *name, hid_t lapl_id);
 herr_t H5Oget_info(hid_t loc_id, H5O_info_t *oinfo);
 herr_t H5Oget_info_by_name(hid_t loc_id, const char *name, H5O_info_t *oinfo,
    hid_t lapl_id);
 herr_t H5Oget_info_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n, H5O_info_t *oinfo,
    hid_t lapl_id);
 herr_t H5Olink(hid_t obj_id, hid_t new_loc_id, const char *new_name,
    hid_t lcpl_id, hid_t lapl_id);
 herr_t H5Oincr_refcount(hid_t object_id);
 herr_t H5Odecr_refcount(hid_t object_id);
 herr_t H5Ocopy(hid_t src_loc_id, const char *src_name, hid_t dst_loc_id,
    const char *dst_name, hid_t ocpypl_id, hid_t lcpl_id);
 herr_t H5Oset_comment(hid_t obj_id, const char *comment);
 herr_t H5Oset_comment_by_name(hid_t loc_id, const char *name,
    const char *comment, hid_t lapl_id);
 ssize_t H5Oget_comment(hid_t obj_id, char *comment, size_t bufsize);
 ssize_t H5Oget_comment_by_name(hid_t loc_id, const char *name,
    char *comment, size_t bufsize, hid_t lapl_id);
 herr_t H5Ovisit(hid_t obj_id, H5_index_t idx_type, H5_iter_order_t order,
    H5O_iterate_t op, void *op_data);
 herr_t H5Ovisit_by_name(hid_t loc_id, const char *obj_name,
    H5_index_t idx_type, H5_iter_order_t order, H5O_iterate_t op,
    void *op_data, hid_t lapl_id);
 herr_t H5Oclose(hid_t object_id);
# 188 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Opublic.h"
typedef struct H5O_stat_t {
    hsize_t size;
    hsize_t free;
    unsigned nmesgs;
    unsigned nchunks;
} H5O_stat_t;
# 25 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Apublic.h" 2







typedef struct {
    hbool_t corder_valid;
    H5O_msg_crt_idx_t corder;
    H5T_cset_t cset;
    hsize_t data_size;
} H5A_info_t;


typedef herr_t (*H5A_operator2_t)(hid_t location_id ,
    const char *attr_name , const H5A_info_t *ainfo , void *op_data );


 hid_t H5Acreate2(hid_t loc_id, const char *attr_name, hid_t type_id,
    hid_t space_id, hid_t acpl_id, hid_t aapl_id);
 hid_t H5Acreate_by_name(hid_t loc_id, const char *obj_name, const char *attr_name,
    hid_t type_id, hid_t space_id, hid_t acpl_id, hid_t aapl_id, hid_t lapl_id);
 hid_t H5Aopen(hid_t obj_id, const char *attr_name, hid_t aapl_id);
 hid_t H5Aopen_by_name(hid_t loc_id, const char *obj_name,
    const char *attr_name, hid_t aapl_id, hid_t lapl_id);
 hid_t H5Aopen_by_idx(hid_t loc_id, const char *obj_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n, hid_t aapl_id,
    hid_t lapl_id);
 herr_t H5Awrite(hid_t attr_id, hid_t type_id, const void *buf);
 herr_t H5Aread(hid_t attr_id, hid_t type_id, void *buf);
 herr_t H5Aclose(hid_t attr_id);
 hid_t H5Aget_space(hid_t attr_id);
 hid_t H5Aget_type(hid_t attr_id);
 hid_t H5Aget_create_plist(hid_t attr_id);
 ssize_t H5Aget_name(hid_t attr_id, size_t buf_size, char *buf);
 ssize_t H5Aget_name_by_idx(hid_t loc_id, const char *obj_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n,
    char *name , size_t size, hid_t lapl_id);
 hsize_t H5Aget_storage_size(hid_t attr_id);
 herr_t H5Aget_info(hid_t attr_id, H5A_info_t *ainfo );
 herr_t H5Aget_info_by_name(hid_t loc_id, const char *obj_name,
    const char *attr_name, H5A_info_t *ainfo , hid_t lapl_id);
 herr_t H5Aget_info_by_idx(hid_t loc_id, const char *obj_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n,
    H5A_info_t *ainfo , hid_t lapl_id);
 herr_t H5Arename(hid_t loc_id, const char *old_name, const char *new_name);
 herr_t H5Arename_by_name(hid_t loc_id, const char *obj_name,
    const char *old_attr_name, const char *new_attr_name, hid_t lapl_id);
 herr_t H5Aiterate2(hid_t loc_id, H5_index_t idx_type,
    H5_iter_order_t order, hsize_t *idx, H5A_operator2_t op, void *op_data);
 herr_t H5Aiterate_by_name(hid_t loc_id, const char *obj_name, H5_index_t idx_type,
    H5_iter_order_t order, hsize_t *idx, H5A_operator2_t op, void *op_data,
    hid_t lapd_id);
 herr_t H5Adelete(hid_t loc_id, const char *name);
 herr_t H5Adelete_by_name(hid_t loc_id, const char *obj_name,
    const char *attr_name, hid_t lapl_id);
 herr_t H5Adelete_by_idx(hid_t loc_id, const char *obj_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n, hid_t lapl_id);
 htri_t H5Aexists(hid_t obj_id, const char *attr_name);
 htri_t H5Aexists_by_name(hid_t obj_id, const char *obj_name,
    const char *attr_name, hid_t lapl_id);
# 100 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Apublic.h"
typedef herr_t (*H5A_operator1_t)(hid_t location_id ,
    const char *attr_name , void *operator_data );



 hid_t H5Acreate1(hid_t loc_id, const char *name, hid_t type_id,
    hid_t space_id, hid_t acpl_id);
 hid_t H5Aopen_name(hid_t loc_id, const char *name);
 hid_t H5Aopen_idx(hid_t loc_id, unsigned idx);
 int H5Aget_num_attrs(hid_t loc_id);
 herr_t H5Aiterate1(hid_t loc_id, unsigned *attr_num, H5A_operator1_t op,
    void *op_data);
# 26 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5ACpublic.h" 1
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5ACpublic.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Cpublic.h" 1
# 38 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Cpublic.h"
enum H5C_cache_incr_mode
{
    H5C_incr__off,
    H5C_incr__threshold
};

enum H5C_cache_flash_incr_mode
{
     H5C_flash_incr__off,
     H5C_flash_incr__add_space
};

enum H5C_cache_decr_mode
{
    H5C_decr__off,
    H5C_decr__threshold,
    H5C_decr__age_out,
    H5C_decr__age_out_with_threshold
};
# 34 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5ACpublic.h" 2
# 443 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5ACpublic.h"
typedef struct H5AC_cache_config_t
{

    int version;

    hbool_t rpt_fcn_enabled;

    hbool_t open_trace_file;
    hbool_t close_trace_file;
    char trace_file_name[1024 + 1];

    hbool_t evictions_enabled;

    hbool_t set_initial_size;
    size_t initial_size;

    double min_clean_fraction;

    size_t max_size;
    size_t min_size;

    long int epoch_length;



    enum H5C_cache_incr_mode incr_mode;

    double lower_hr_threshold;

    double increment;

    hbool_t apply_max_increment;
    size_t max_increment;

    enum H5C_cache_flash_incr_mode flash_incr_mode;
    double flash_multiple;
    double flash_threshold;



    enum H5C_cache_decr_mode decr_mode;

    double upper_hr_threshold;

    double decrement;

    hbool_t apply_max_decrement;
    size_t max_decrement;

    int epochs_before_eviction;

    hbool_t apply_empty_reserve;
    double empty_reserve;



    int dirty_bytes_threshold;
    int metadata_write_strategy;

} H5AC_cache_config_t;
# 27 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Dpublic.h" 1
# 42 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Dpublic.h"
typedef enum H5D_layout_t {
    H5D_LAYOUT_ERROR = -1,

    H5D_COMPACT = 0,
    H5D_CONTIGUOUS = 1,
    H5D_CHUNKED = 2,
    H5D_NLAYOUTS = 3
} H5D_layout_t;


typedef enum H5D_chunk_index_t {
    H5D_CHUNK_BTREE = 0
} H5D_chunk_index_t;


typedef enum H5D_alloc_time_t {
    H5D_ALLOC_TIME_ERROR = -1,
    H5D_ALLOC_TIME_DEFAULT = 0,
    H5D_ALLOC_TIME_EARLY = 1,
    H5D_ALLOC_TIME_LATE = 2,
    H5D_ALLOC_TIME_INCR = 3
} H5D_alloc_time_t;


typedef enum H5D_space_status_t {
    H5D_SPACE_STATUS_ERROR = -1,
    H5D_SPACE_STATUS_NOT_ALLOCATED = 0,
    H5D_SPACE_STATUS_PART_ALLOCATED = 1,
    H5D_SPACE_STATUS_ALLOCATED = 2
} H5D_space_status_t;


typedef enum H5D_fill_time_t {
    H5D_FILL_TIME_ERROR = -1,
    H5D_FILL_TIME_ALLOC = 0,
    H5D_FILL_TIME_NEVER = 1,
    H5D_FILL_TIME_IFSET = 2
} H5D_fill_time_t;


typedef enum H5D_fill_value_t {
    H5D_FILL_VALUE_ERROR =-1,
    H5D_FILL_VALUE_UNDEFINED =0,
    H5D_FILL_VALUE_DEFAULT =1,
    H5D_FILL_VALUE_USER_DEFINED =2
} H5D_fill_value_t;
# 101 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Dpublic.h"
typedef herr_t (*H5D_operator_t)(void *elem, hid_t type_id, unsigned ndim,
     const hsize_t *point, void *operator_data);

 hid_t H5Dcreate2(hid_t loc_id, const char *name, hid_t type_id,
    hid_t space_id, hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id);
 hid_t H5Dcreate_anon(hid_t file_id, hid_t type_id, hid_t space_id,
    hid_t plist_id, hid_t dapl_id);
 hid_t H5Dopen2(hid_t file_id, const char *name, hid_t dapl_id);
 herr_t H5Dclose(hid_t dset_id);
 hid_t H5Dget_space(hid_t dset_id);
 herr_t H5Dget_space_status(hid_t dset_id, H5D_space_status_t *allocation);
 hid_t H5Dget_type(hid_t dset_id);
 hid_t H5Dget_create_plist(hid_t dset_id);
 hid_t H5Dget_access_plist(hid_t dset_id);
 hsize_t H5Dget_storage_size(hid_t dset_id);
 haddr_t H5Dget_offset(hid_t dset_id);
 herr_t H5Dread(hid_t dset_id, hid_t mem_type_id, hid_t mem_space_id,
   hid_t file_space_id, hid_t plist_id, void *buf );
 herr_t H5Dwrite(hid_t dset_id, hid_t mem_type_id, hid_t mem_space_id,
    hid_t file_space_id, hid_t plist_id, const void *buf);
 herr_t H5Diterate(void *buf, hid_t type_id, hid_t space_id,
            H5D_operator_t op, void *operator_data);
 herr_t H5Dvlen_reclaim(hid_t type_id, hid_t space_id, hid_t plist_id, void *buf);
 herr_t H5Dvlen_get_buf_size(hid_t dataset_id, hid_t type_id, hid_t space_id, hsize_t *size);
 herr_t H5Dfill(const void *fill, hid_t fill_type, void *buf,
        hid_t buf_type, hid_t space);
 herr_t H5Dset_extent(hid_t dset_id, const hsize_t size[]);
 herr_t H5Ddebug(hid_t dset_id);
# 143 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Dpublic.h"
 hid_t H5Dcreate1(hid_t file_id, const char *name, hid_t type_id,
    hid_t space_id, hid_t dcpl_id);
 hid_t H5Dopen1(hid_t file_id, const char *name);
 herr_t H5Dextend(hid_t dset_id, const hsize_t size[]);
# 28 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epublic.h" 1
# 32 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epublic.h"
typedef enum H5E_type_t {
    H5E_MAJOR,
    H5E_MINOR
} H5E_type_t;


typedef struct H5E_error2_t {
    hid_t cls_id;
    hid_t maj_num;
    hid_t min_num;
    unsigned line;
    const char *func_name;
    const char *file_name;
    const char *desc;
} H5E_error2_t;
# 58 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epublic.h"
extern hid_t H5E_ERR_CLS_g;



# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h" 1
# 57 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_DATASET_g;
extern hid_t H5E_FUNC_g;
extern hid_t H5E_STORAGE_g;
extern hid_t H5E_FILE_g;
extern hid_t H5E_SOHM_g;
extern hid_t H5E_SYM_g;
extern hid_t H5E_VFL_g;
extern hid_t H5E_INTERNAL_g;
extern hid_t H5E_BTREE_g;
extern hid_t H5E_REFERENCE_g;
extern hid_t H5E_DATASPACE_g;
extern hid_t H5E_RESOURCE_g;
extern hid_t H5E_PLIST_g;
extern hid_t H5E_LINK_g;
extern hid_t H5E_DATATYPE_g;
extern hid_t H5E_RS_g;
extern hid_t H5E_HEAP_g;
extern hid_t H5E_OHDR_g;
extern hid_t H5E_ATOM_g;
extern hid_t H5E_ATTR_g;
extern hid_t H5E_NONE_MAJOR_g;
extern hid_t H5E_IO_g;
extern hid_t H5E_SLIST_g;
extern hid_t H5E_EFL_g;
extern hid_t H5E_TST_g;
extern hid_t H5E_ARGS_g;
extern hid_t H5E_ERROR_g;
extern hid_t H5E_PLINE_g;
extern hid_t H5E_FSPACE_g;
extern hid_t H5E_CACHE_g;
# 99 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_SEEKERROR_g;
extern hid_t H5E_READERROR_g;
extern hid_t H5E_WRITEERROR_g;
extern hid_t H5E_CLOSEERROR_g;
extern hid_t H5E_OVERFLOW_g;
extern hid_t H5E_FCNTL_g;
# 117 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_NOSPACE_g;
extern hid_t H5E_CANTALLOC_g;
extern hid_t H5E_CANTCOPY_g;
extern hid_t H5E_CANTFREE_g;
extern hid_t H5E_ALREADYEXISTS_g;
extern hid_t H5E_CANTLOCK_g;
extern hid_t H5E_CANTUNLOCK_g;
extern hid_t H5E_CANTGC_g;
extern hid_t H5E_CANTGETSIZE_g;
extern hid_t H5E_OBJOPEN_g;
# 135 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_CANTRESTORE_g;
extern hid_t H5E_CANTCOMPUTE_g;
extern hid_t H5E_CANTEXTEND_g;
extern hid_t H5E_CANTATTACH_g;
extern hid_t H5E_CANTUPDATE_g;
extern hid_t H5E_CANTOPERATE_g;





extern hid_t H5E_CANTINIT_g;
extern hid_t H5E_ALREADYINIT_g;
extern hid_t H5E_CANTRELEASE_g;





extern hid_t H5E_CANTGET_g;
extern hid_t H5E_CANTSET_g;
extern hid_t H5E_DUPCLASS_g;





extern hid_t H5E_CANTMERGE_g;
extern hid_t H5E_CANTREVIVE_g;
extern hid_t H5E_CANTSHRINK_g;
# 176 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_LINKCOUNT_g;
extern hid_t H5E_VERSION_g;
extern hid_t H5E_ALIGNMENT_g;
extern hid_t H5E_BADMESG_g;
extern hid_t H5E_CANTDELETE_g;
extern hid_t H5E_BADITER_g;
extern hid_t H5E_CANTPACK_g;
extern hid_t H5E_CANTRESET_g;
extern hid_t H5E_CANTRENAME_g;



extern hid_t H5E_SYSERRSTR_g;
# 197 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_NOFILTER_g;
extern hid_t H5E_CALLBACK_g;
extern hid_t H5E_CANAPPLY_g;
extern hid_t H5E_SETLOCAL_g;
extern hid_t H5E_NOENCODER_g;
extern hid_t H5E_CANTFILTER_g;






extern hid_t H5E_CANTOPENOBJ_g;
extern hid_t H5E_CANTCLOSEOBJ_g;
extern hid_t H5E_COMPLEN_g;
extern hid_t H5E_PATH_g;



extern hid_t H5E_NONE_MINOR_g;
# 228 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_FILEEXISTS_g;
extern hid_t H5E_FILEOPEN_g;
extern hid_t H5E_CANTCREATE_g;
extern hid_t H5E_CANTOPENFILE_g;
extern hid_t H5E_CANTCLOSEFILE_g;
extern hid_t H5E_NOTHDF5_g;
extern hid_t H5E_BADFILE_g;
extern hid_t H5E_TRUNCATED_g;
extern hid_t H5E_MOUNT_g;
# 245 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_BADATOM_g;
extern hid_t H5E_BADGROUP_g;
extern hid_t H5E_CANTREGISTER_g;
extern hid_t H5E_CANTINC_g;
extern hid_t H5E_CANTDEC_g;
extern hid_t H5E_NOIDS_g;
# 268 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_CANTFLUSH_g;
extern hid_t H5E_CANTSERIALIZE_g;
extern hid_t H5E_CANTLOAD_g;
extern hid_t H5E_PROTECT_g;
extern hid_t H5E_NOTCACHED_g;
extern hid_t H5E_SYSTEM_g;
extern hid_t H5E_CANTINS_g;
extern hid_t H5E_CANTPROTECT_g;
extern hid_t H5E_CANTUNPROTECT_g;
extern hid_t H5E_CANTPIN_g;
extern hid_t H5E_CANTUNPIN_g;
extern hid_t H5E_CANTMARKDIRTY_g;
extern hid_t H5E_CANTDIRTY_g;
extern hid_t H5E_CANTEXPUNGE_g;
extern hid_t H5E_CANTRESIZE_g;







extern hid_t H5E_TRAVERSE_g;
extern hid_t H5E_NLINKS_g;
extern hid_t H5E_NOTREGISTERED_g;
extern hid_t H5E_CANTMOVE_g;
extern hid_t H5E_CANTSORT_g;





extern hid_t H5E_MPI_g;
extern hid_t H5E_MPIERRSTR_g;
extern hid_t H5E_CANTRECV_g;
# 311 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_CANTCLIP_g;
extern hid_t H5E_CANTCOUNT_g;
extern hid_t H5E_CANTSELECT_g;
extern hid_t H5E_CANTNEXT_g;
extern hid_t H5E_BADSELECT_g;
extern hid_t H5E_CANTCOMPARE_g;







extern hid_t H5E_UNINITIALIZED_g;
extern hid_t H5E_UNSUPPORTED_g;
extern hid_t H5E_BADTYPE_g;
extern hid_t H5E_BADRANGE_g;
extern hid_t H5E_BADVALUE_g;
# 342 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epubgen.h"
extern hid_t H5E_NOTFOUND_g;
extern hid_t H5E_EXISTS_g;
extern hid_t H5E_CANTENCODE_g;
extern hid_t H5E_CANTDECODE_g;
extern hid_t H5E_CANTSPLIT_g;
extern hid_t H5E_CANTREDISTRIBUTE_g;
extern hid_t H5E_CANTSWAP_g;
extern hid_t H5E_CANTINSERT_g;
extern hid_t H5E_CANTLIST_g;
extern hid_t H5E_CANTMODIFY_g;
extern hid_t H5E_CANTREMOVE_g;




extern hid_t H5E_CANTCONVERT_g;
extern hid_t H5E_BADSIZE_g;
# 63 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epublic.h" 2
# 141 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epublic.h"
typedef enum H5E_direction_t {
    H5E_WALK_UPWARD = 0,
    H5E_WALK_DOWNWARD = 1
} H5E_direction_t;







typedef herr_t (*H5E_walk2_t)(unsigned n, const H5E_error2_t *err_desc,
    void *client_data);
typedef herr_t (*H5E_auto2_t)(hid_t estack, void *client_data);


 hid_t H5Eregister_class(const char *cls_name, const char *lib_name,
    const char *version);
 herr_t H5Eunregister_class(hid_t class_id);
 herr_t H5Eclose_msg(hid_t err_id);
 hid_t H5Ecreate_msg(hid_t cls, H5E_type_t msg_type, const char *msg);
 hid_t H5Ecreate_stack(void);
 hid_t H5Eget_current_stack(void);
 herr_t H5Eclose_stack(hid_t stack_id);
 ssize_t H5Eget_class_name(hid_t class_id, char *name, size_t size);
 herr_t H5Eset_current_stack(hid_t err_stack_id);
 herr_t H5Epush2(hid_t err_stack, const char *file, const char *func, unsigned line,
    hid_t cls_id, hid_t maj_id, hid_t min_id, const char *msg, ...);
 herr_t H5Epop(hid_t err_stack, size_t count);
 herr_t H5Eprint2(hid_t err_stack, FILE *stream);
 herr_t H5Ewalk2(hid_t err_stack, H5E_direction_t direction, H5E_walk2_t func,
    void *client_data);
 herr_t H5Eget_auto2(hid_t estack_id, H5E_auto2_t *func, void **client_data);
 herr_t H5Eset_auto2(hid_t estack_id, H5E_auto2_t func, void *client_data);
 herr_t H5Eclear2(hid_t err_stack);
 herr_t H5Eauto_is_v2(hid_t err_stack, unsigned *is_stack);
 ssize_t H5Eget_msg(hid_t msg_id, H5E_type_t *type, char *msg,
    size_t size);
 ssize_t H5Eget_num(hid_t error_stack_id);
# 193 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Epublic.h"
typedef hid_t H5E_major_t;
typedef hid_t H5E_minor_t;


typedef struct H5E_error1_t {
    H5E_major_t maj_num;
    H5E_minor_t min_num;
    const char *func_name;
    const char *file_name;
    unsigned line;
    const char *desc;
} H5E_error1_t;


typedef herr_t (*H5E_walk1_t)(int n, H5E_error1_t *err_desc, void *client_data);
typedef herr_t (*H5E_auto1_t)(void *client_data);


 herr_t H5Eclear1(void);
 herr_t H5Eget_auto1(H5E_auto1_t *func, void **client_data);
 herr_t H5Epush1(const char *file, const char *func, unsigned line,
    H5E_major_t maj, H5E_minor_t min, const char *str);
 herr_t H5Eprint1(FILE *stream);
 herr_t H5Eset_auto1(H5E_auto1_t func, void *client_data);
 herr_t H5Ewalk1(H5E_direction_t direction, H5E_walk1_t func,
    void *client_data);
 char *H5Eget_major(H5E_major_t maj);
 char *H5Eget_minor(H5E_minor_t min);
# 29 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Fpublic.h" 1
# 79 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Fpublic.h"
typedef enum H5F_scope_t {
    H5F_SCOPE_LOCAL = 0,
    H5F_SCOPE_GLOBAL = 1
} H5F_scope_t;
# 95 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Fpublic.h"
typedef enum H5F_close_degree_t {
    H5F_CLOSE_DEFAULT = 0,
    H5F_CLOSE_WEAK = 1,
    H5F_CLOSE_SEMI = 2,
    H5F_CLOSE_STRONG = 3
} H5F_close_degree_t;



typedef struct H5F_info_t {
    hsize_t super_ext_size;
    struct {
 hsize_t hdr_size;
 H5_ih_info_t msgs_info;
    } sohm;
} H5F_info_t;






typedef enum H5F_mem_t {
    H5FD_MEM_NOLIST = -1,
    H5FD_MEM_DEFAULT = 0,
    H5FD_MEM_SUPER = 1,
    H5FD_MEM_BTREE = 2,
    H5FD_MEM_DRAW = 3,
    H5FD_MEM_GHEAP = 4,
    H5FD_MEM_LHEAP = 5,
    H5FD_MEM_OHDR = 6,

    H5FD_MEM_NTYPES
} H5F_mem_t;


typedef enum H5F_libver_t {
    H5F_LIBVER_EARLIEST,
    H5F_LIBVER_LATEST
} H5F_libver_t;
# 145 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Fpublic.h"
 htri_t H5Fis_hdf5(const char *filename);
 hid_t H5Fcreate(const char *filename, unsigned flags,
       hid_t create_plist, hid_t access_plist);
 hid_t H5Fopen(const char *filename, unsigned flags,
          hid_t access_plist);
 hid_t H5Freopen(hid_t file_id);
 herr_t H5Fflush(hid_t object_id, H5F_scope_t scope);
 herr_t H5Fclose(hid_t file_id);
 hid_t H5Fget_create_plist(hid_t file_id);
 hid_t H5Fget_access_plist(hid_t file_id);
 herr_t H5Fget_intent(hid_t file_id, unsigned * intent);
 ssize_t H5Fget_obj_count(hid_t file_id, unsigned types);
 ssize_t H5Fget_obj_ids(hid_t file_id, unsigned types, size_t max_objs, hid_t *obj_id_list);
 herr_t H5Fget_vfd_handle(hid_t file_id, hid_t fapl, void **file_handle);
 herr_t H5Fmount(hid_t loc, const char *name, hid_t child, hid_t plist);
 herr_t H5Funmount(hid_t loc, const char *name);
 hssize_t H5Fget_freespace(hid_t file_id);
 herr_t H5Fget_filesize(hid_t file_id, hsize_t *size);
 herr_t H5Fget_mdc_config(hid_t file_id,
    H5AC_cache_config_t * config_ptr);
 herr_t H5Fset_mdc_config(hid_t file_id,
    H5AC_cache_config_t * config_ptr);
 herr_t H5Fget_mdc_hit_rate(hid_t file_id, double * hit_rate_ptr);
 herr_t H5Fget_mdc_size(hid_t file_id,
                              size_t * max_size_ptr,
                              size_t * min_clean_size_ptr,
                              size_t * cur_size_ptr,
                              int * cur_num_entries_ptr);
 herr_t H5Freset_mdc_hit_rate_stats(hid_t file_id);
 ssize_t H5Fget_name(hid_t obj_id, char *name, size_t size);
 herr_t H5Fget_info(hid_t obj_id, H5F_info_t *bh_info);
# 30 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDpublic.h" 1
# 30 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDpublic.h"
typedef enum H5F_mem_t H5FD_mem_t;
# 183 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDpublic.h"
typedef struct H5FD_t H5FD_t;


typedef struct H5FD_class_t {
    const char *name;
    haddr_t maxaddr;
    H5F_close_degree_t fc_degree;
    hsize_t (*sb_size)(H5FD_t *file);
    herr_t (*sb_encode)(H5FD_t *file, char *name ,
                         unsigned char *p );
    herr_t (*sb_decode)(H5FD_t *f, const char *name, const unsigned char *p);
    size_t fapl_size;
    void * (*fapl_get)(H5FD_t *file);
    void * (*fapl_copy)(const void *fapl);
    herr_t (*fapl_free)(void *fapl);
    size_t dxpl_size;
    void * (*dxpl_copy)(const void *dxpl);
    herr_t (*dxpl_free)(void *dxpl);
    H5FD_t *(*open)(const char *name, unsigned flags, hid_t fapl,
                    haddr_t maxaddr);
    herr_t (*close)(H5FD_t *file);
    int (*cmp)(const H5FD_t *f1, const H5FD_t *f2);
    herr_t (*query)(const H5FD_t *f1, unsigned long *flags);
    herr_t (*get_type_map)(const H5FD_t *file, H5FD_mem_t *type_map);
    haddr_t (*alloc)(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, hsize_t size);
    herr_t (*free)(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id,
                    haddr_t addr, hsize_t size);
    haddr_t (*get_eoa)(const H5FD_t *file, H5FD_mem_t type);
    herr_t (*set_eoa)(H5FD_t *file, H5FD_mem_t type, haddr_t addr);
    haddr_t (*get_eof)(const H5FD_t *file);
    herr_t (*get_handle)(H5FD_t *file, hid_t fapl, void**file_handle);
    herr_t (*read)(H5FD_t *file, H5FD_mem_t type, hid_t dxpl,
                    haddr_t addr, size_t size, void *buffer);
    herr_t (*write)(H5FD_t *file, H5FD_mem_t type, hid_t dxpl,
                     haddr_t addr, size_t size, const void *buffer);
    herr_t (*flush)(H5FD_t *file, hid_t dxpl_id, unsigned closing);
    herr_t (*truncate)(H5FD_t *file, hid_t dxpl_id, hbool_t closing);
    herr_t (*lock)(H5FD_t *file, unsigned char *oid, unsigned lock_type, hbool_t last);
    herr_t (*unlock)(H5FD_t *file, unsigned char *oid, hbool_t last);
    H5FD_mem_t fl_map[H5FD_MEM_NTYPES];
} H5FD_class_t;


typedef struct H5FD_free_t {
    haddr_t addr;
    hsize_t size;
    struct H5FD_free_t *next;
} H5FD_free_t;





struct H5FD_t {
    hid_t driver_id;
    const H5FD_class_t *cls;
    unsigned long fileno;
    unsigned long feature_flags;
    haddr_t maxaddr;
    haddr_t base_addr;


    hsize_t threshold;
    hsize_t alignment;
};






 hid_t H5FDregister(const H5FD_class_t *cls);
 herr_t H5FDunregister(hid_t driver_id);
 H5FD_t *H5FDopen(const char *name, unsigned flags, hid_t fapl_id,
                        haddr_t maxaddr);
 herr_t H5FDclose(H5FD_t *file);
 int H5FDcmp(const H5FD_t *f1, const H5FD_t *f2);
 int H5FDquery(const H5FD_t *f, unsigned long *flags);
 haddr_t H5FDalloc(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, hsize_t size);
 herr_t H5FDfree(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id,
                       haddr_t addr, hsize_t size);
 haddr_t H5FDget_eoa(H5FD_t *file, H5FD_mem_t type);
 herr_t H5FDset_eoa(H5FD_t *file, H5FD_mem_t type, haddr_t eoa);
 haddr_t H5FDget_eof(H5FD_t *file);
 herr_t H5FDget_vfd_handle(H5FD_t *file, hid_t fapl, void**file_handle);
 herr_t H5FDread(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id,
                       haddr_t addr, size_t size, void *buf );
 herr_t H5FDwrite(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id,
                        haddr_t addr, size_t size, const void *buf);
 herr_t H5FDflush(H5FD_t *file, hid_t dxpl_id, unsigned closing);
 herr_t H5FDtruncate(H5FD_t *file, hid_t dxpl_id, hbool_t closing);
# 31 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Gpublic.h" 1
# 51 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Gpublic.h"
typedef enum H5G_storage_type_t {
    H5G_STORAGE_TYPE_UNKNOWN = -1,
    H5G_STORAGE_TYPE_SYMBOL_TABLE,

    H5G_STORAGE_TYPE_COMPACT,
    H5G_STORAGE_TYPE_DENSE
} H5G_storage_type_t;


typedef struct H5G_info_t {
    H5G_storage_type_t storage_type;
    hsize_t nlinks;
    int64_t max_corder;
    hbool_t mounted;
} H5G_info_t;
# 75 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Gpublic.h"
 hid_t H5Gcreate2(hid_t loc_id, const char *name, hid_t lcpl_id,
    hid_t gcpl_id, hid_t gapl_id);
 hid_t H5Gcreate_anon(hid_t loc_id, hid_t gcpl_id, hid_t gapl_id);
 hid_t H5Gopen2(hid_t loc_id, const char *name, hid_t gapl_id);
 hid_t H5Gget_create_plist(hid_t group_id);
 herr_t H5Gget_info(hid_t loc_id, H5G_info_t *ginfo);
 herr_t H5Gget_info_by_name(hid_t loc_id, const char *name, H5G_info_t *ginfo,
    hid_t lapl_id);
 herr_t H5Gget_info_by_idx(hid_t loc_id, const char *group_name,
    H5_index_t idx_type, H5_iter_order_t order, hsize_t n, H5G_info_t *ginfo,
    hid_t lapl_id);
 herr_t H5Gclose(hid_t group_id);
# 119 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Gpublic.h"
typedef enum H5G_obj_t {
    H5G_UNKNOWN = -1,
    H5G_GROUP,
    H5G_DATASET,
    H5G_TYPE,
    H5G_LINK,
    H5G_UDLINK,
    H5G_RESERVED_5,
    H5G_RESERVED_6,
    H5G_RESERVED_7
} H5G_obj_t;


typedef herr_t (*H5G_iterate_t)(hid_t group, const char *name, void *op_data);


typedef struct H5G_stat_t {
    unsigned long fileno[2];
    unsigned long objno[2];
    unsigned nlink;
    H5G_obj_t type;
    time_t mtime;
    size_t linklen;
    H5O_stat_t ohdr;
} H5G_stat_t;



 hid_t H5Gcreate1(hid_t loc_id, const char *name, size_t size_hint);
 hid_t H5Gopen1(hid_t loc_id, const char *name);
 herr_t H5Glink(hid_t cur_loc_id, H5L_type_t type, const char *cur_name,
    const char *new_name);
 herr_t H5Glink2(hid_t cur_loc_id, const char *cur_name, H5L_type_t type,
    hid_t new_loc_id, const char *new_name);
 herr_t H5Gmove(hid_t src_loc_id, const char *src_name,
    const char *dst_name);
 herr_t H5Gmove2(hid_t src_loc_id, const char *src_name, hid_t dst_loc_id,
    const char *dst_name);
 herr_t H5Gunlink(hid_t loc_id, const char *name);
 herr_t H5Gget_linkval(hid_t loc_id, const char *name, size_t size,
    char *buf );
 herr_t H5Gset_comment(hid_t loc_id, const char *name, const char *comment);
 int H5Gget_comment(hid_t loc_id, const char *name, size_t bufsize,
    char *buf);
 herr_t H5Giterate(hid_t loc_id, const char *name, int *idx,
        H5G_iterate_t op, void *op_data);
 herr_t H5Gget_num_objs(hid_t loc_id, hsize_t *num_objs);
 herr_t H5Gget_objinfo(hid_t loc_id, const char *name,
    hbool_t follow_link, H5G_stat_t *statbuf );
 ssize_t H5Gget_objname_by_idx(hid_t loc_id, hsize_t idx, char* name,
    size_t size);
 H5G_obj_t H5Gget_objtype_by_idx(hid_t loc_id, hsize_t idx);
# 32 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2


# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5MMpublic.h" 1
# 36 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5MMpublic.h"
typedef void *(*H5MM_allocate_t)(size_t size, void *alloc_info);
typedef void (*H5MM_free_t)(void *mem, void *free_info);
# 35 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2

# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ppublic.h" 1
# 35 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ppublic.h"
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h" 1
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h"
typedef int H5Z_filter_t;
# 98 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h"
typedef enum H5Z_SO_scale_type_t {
    H5Z_SO_FLOAT_DSCALE = 0,
    H5Z_SO_FLOAT_ESCALE = 1,
    H5Z_SO_INT = 2
} H5Z_SO_scale_type_t;





typedef enum H5Z_EDC_t {
    H5Z_ERROR_EDC = -1,
    H5Z_DISABLE_EDC = 0,
    H5Z_ENABLE_EDC = 1,
    H5Z_NO_EDC = 2
} H5Z_EDC_t;






typedef enum H5Z_cb_return_t {
    H5Z_CB_ERROR = -1,
    H5Z_CB_FAIL = 0,
    H5Z_CB_CONT = 1,
    H5Z_CB_NO = 2
} H5Z_cb_return_t;


typedef H5Z_cb_return_t (*H5Z_filter_func_t)(H5Z_filter_t filter, void* buf,
                                size_t buf_size, void* op_data);


typedef struct H5Z_cb_t {
    H5Z_filter_func_t func;
    void* op_data;
} H5Z_cb_t;
# 161 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h"
typedef htri_t (*H5Z_can_apply_func_t)(hid_t dcpl_id, hid_t type_id, hid_t space_id);
# 184 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h"
typedef herr_t (*H5Z_set_local_func_t)(hid_t dcpl_id, hid_t type_id, hid_t space_id);
# 201 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h"
typedef size_t (*H5Z_func_t)(unsigned int flags, size_t cd_nelmts,
        const unsigned int cd_values[], size_t nbytes,
        size_t *buf_size, void **buf);





typedef struct H5Z_class2_t {
    int version;
    H5Z_filter_t id;
    unsigned encoder_present;
    unsigned decoder_present;
    const char *name;
    H5Z_can_apply_func_t can_apply;
    H5Z_set_local_func_t set_local;
    H5Z_func_t filter;
} H5Z_class2_t;

 herr_t H5Zregister(const void *cls);
 herr_t H5Zunregister(H5Z_filter_t id);
 htri_t H5Zfilter_avail(H5Z_filter_t id);
 herr_t H5Zget_filter_info(H5Z_filter_t filter, unsigned int *filter_config_flags);
# 235 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Zpublic.h"
typedef struct H5Z_class1_t {
    H5Z_filter_t id;
    const char *name;
    H5Z_can_apply_func_t can_apply;
    H5Z_set_local_func_t set_local;
    H5Z_func_t filter;
} H5Z_class1_t;
# 36 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ppublic.h" 2
# 104 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ppublic.h"
typedef herr_t (*H5P_cls_create_func_t)(hid_t prop_id, void *create_data);
typedef herr_t (*H5P_cls_copy_func_t)(hid_t new_prop_id, hid_t old_prop_id,
                                      void *copy_data);
typedef herr_t (*H5P_cls_close_func_t)(hid_t prop_id, void *close_data);


typedef herr_t (*H5P_prp_cb1_t)(const char *name, size_t size, void *value);
typedef herr_t (*H5P_prp_cb2_t)(hid_t prop_id, const char *name, size_t size, void *value);
typedef H5P_prp_cb1_t H5P_prp_create_func_t;
typedef H5P_prp_cb2_t H5P_prp_set_func_t;
typedef H5P_prp_cb2_t H5P_prp_get_func_t;
typedef H5P_prp_cb2_t H5P_prp_delete_func_t;
typedef H5P_prp_cb1_t H5P_prp_copy_func_t;
typedef int (*H5P_prp_compare_func_t)(const void *value1, const void *value2, size_t size);
typedef H5P_prp_cb1_t H5P_prp_close_func_t;


typedef herr_t (*H5P_iterate_t)(hid_t id, const char *name, void *iter_data);







extern hid_t H5P_CLS_ROOT_g;
extern hid_t H5P_CLS_OBJECT_CREATE_g;
extern hid_t H5P_CLS_FILE_CREATE_g;
extern hid_t H5P_CLS_FILE_ACCESS_g;
extern hid_t H5P_CLS_DATASET_CREATE_g;
extern hid_t H5P_CLS_DATASET_ACCESS_g;
extern hid_t H5P_CLS_DATASET_XFER_g;
extern hid_t H5P_CLS_FILE_MOUNT_g;
extern hid_t H5P_CLS_GROUP_CREATE_g;
extern hid_t H5P_CLS_GROUP_ACCESS_g;
extern hid_t H5P_CLS_DATATYPE_CREATE_g;
extern hid_t H5P_CLS_DATATYPE_ACCESS_g;
extern hid_t H5P_CLS_STRING_CREATE_g;
extern hid_t H5P_CLS_ATTRIBUTE_CREATE_g;
extern hid_t H5P_CLS_OBJECT_COPY_g;
extern hid_t H5P_CLS_LINK_CREATE_g;
extern hid_t H5P_CLS_LINK_ACCESS_g;



extern hid_t H5P_LST_FILE_CREATE_g;
extern hid_t H5P_LST_FILE_ACCESS_g;
extern hid_t H5P_LST_DATASET_CREATE_g;
extern hid_t H5P_LST_DATASET_ACCESS_g;
extern hid_t H5P_LST_DATASET_XFER_g;
extern hid_t H5P_LST_FILE_MOUNT_g;
extern hid_t H5P_LST_GROUP_CREATE_g;
extern hid_t H5P_LST_GROUP_ACCESS_g;
extern hid_t H5P_LST_DATATYPE_CREATE_g;
extern hid_t H5P_LST_DATATYPE_ACCESS_g;
extern hid_t H5P_LST_ATTRIBUTE_CREATE_g;
extern hid_t H5P_LST_OBJECT_COPY_g;
extern hid_t H5P_LST_LINK_CREATE_g;
extern hid_t H5P_LST_LINK_ACCESS_g;






 hid_t H5Pcreate_class(hid_t parent, const char *name,
    H5P_cls_create_func_t cls_create, void *create_data,
    H5P_cls_copy_func_t cls_copy, void *copy_data,
    H5P_cls_close_func_t cls_close, void *close_data);
 char *H5Pget_class_name(hid_t pclass_id);
 hid_t H5Pcreate(hid_t cls_id);
 herr_t H5Pregister2(hid_t cls_id, const char *name, size_t size,
    void *def_value, H5P_prp_create_func_t prp_create,
    H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
    H5P_prp_delete_func_t prp_del, H5P_prp_copy_func_t prp_copy,
    H5P_prp_compare_func_t prp_cmp, H5P_prp_close_func_t prp_close);
 herr_t H5Pinsert2(hid_t plist_id, const char *name, size_t size,
    void *value, H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
    H5P_prp_delete_func_t prp_delete, H5P_prp_copy_func_t prp_copy,
    H5P_prp_compare_func_t prp_cmp, H5P_prp_close_func_t prp_close);
 herr_t H5Pset(hid_t plist_id, const char *name, void *value);
 htri_t H5Pexist(hid_t plist_id, const char *name);
 herr_t H5Pget_size(hid_t id, const char *name, size_t *size);
 herr_t H5Pget_nprops(hid_t id, size_t *nprops);
 hid_t H5Pget_class(hid_t plist_id);
 hid_t H5Pget_class_parent(hid_t pclass_id);
 herr_t H5Pget(hid_t plist_id, const char *name, void * value);
 htri_t H5Pequal(hid_t id1, hid_t id2);
 htri_t H5Pisa_class(hid_t plist_id, hid_t pclass_id);
 int H5Piterate(hid_t id, int *idx, H5P_iterate_t iter_func,
            void *iter_data);
 herr_t H5Pcopy_prop(hid_t dst_id, hid_t src_id, const char *name);
 herr_t H5Premove(hid_t plist_id, const char *name);
 herr_t H5Punregister(hid_t pclass_id, const char *name);
 herr_t H5Pclose_class(hid_t plist_id);
 herr_t H5Pclose(hid_t plist_id);
 hid_t H5Pcopy(hid_t plist_id);


 herr_t H5Pset_attr_phase_change(hid_t plist_id, unsigned max_compact, unsigned min_dense);
 herr_t H5Pget_attr_phase_change(hid_t plist_id, unsigned *max_compact, unsigned *min_dense);
 herr_t H5Pset_attr_creation_order(hid_t plist_id, unsigned crt_order_flags);
 herr_t H5Pget_attr_creation_order(hid_t plist_id, unsigned *crt_order_flags);
 herr_t H5Pset_obj_track_times(hid_t plist_id, hbool_t track_times);
 herr_t H5Pget_obj_track_times(hid_t plist_id, hbool_t *track_times);
 herr_t H5Pmodify_filter(hid_t plist_id, H5Z_filter_t filter,
        unsigned int flags, size_t cd_nelmts,
        const unsigned int cd_values[ ]);
 herr_t H5Pset_filter(hid_t plist_id, H5Z_filter_t filter,
        unsigned int flags, size_t cd_nelmts,
        const unsigned int c_values[]);
 int H5Pget_nfilters(hid_t plist_id);
 H5Z_filter_t H5Pget_filter2(hid_t plist_id, unsigned filter,
       unsigned int *flags ,
       size_t *cd_nelmts ,
       unsigned cd_values[] ,
       size_t namelen, char name[],
       unsigned *filter_config );
 herr_t H5Pget_filter_by_id2(hid_t plist_id, H5Z_filter_t id,
       unsigned int *flags , size_t *cd_nelmts ,
       unsigned cd_values[] , size_t namelen, char name[] ,
       unsigned *filter_config );
 htri_t H5Pall_filters_avail(hid_t plist_id);
 herr_t H5Premove_filter(hid_t plist_id, H5Z_filter_t filter);
 herr_t H5Pset_deflate(hid_t plist_id, unsigned aggression);
 herr_t H5Pset_fletcher32(hid_t plist_id);


 herr_t H5Pget_version(hid_t plist_id, unsigned *boot ,
         unsigned *freelist , unsigned *stab ,
         unsigned *shhdr );
 herr_t H5Pset_userblock(hid_t plist_id, hsize_t size);
 herr_t H5Pget_userblock(hid_t plist_id, hsize_t *size);
 herr_t H5Pset_sizes(hid_t plist_id, size_t sizeof_addr,
       size_t sizeof_size);
 herr_t H5Pget_sizes(hid_t plist_id, size_t *sizeof_addr ,
       size_t *sizeof_size );
 herr_t H5Pset_sym_k(hid_t plist_id, unsigned ik, unsigned lk);
 herr_t H5Pget_sym_k(hid_t plist_id, unsigned *ik , unsigned *lk );
 herr_t H5Pset_istore_k(hid_t plist_id, unsigned ik);
 herr_t H5Pget_istore_k(hid_t plist_id, unsigned *ik );
 herr_t H5Pset_shared_mesg_nindexes(hid_t plist_id, unsigned nindexes);
 herr_t H5Pget_shared_mesg_nindexes(hid_t plist_id, unsigned *nindexes);
 herr_t H5Pset_shared_mesg_index(hid_t plist_id, unsigned index_num, unsigned mesg_type_flags, unsigned min_mesg_size);
 herr_t H5Pget_shared_mesg_index(hid_t plist_id, unsigned index_num, unsigned *mesg_type_flags, unsigned *min_mesg_size);
 herr_t H5Pset_shared_mesg_phase_change(hid_t plist_id, unsigned max_list, unsigned min_btree);
 herr_t H5Pget_shared_mesg_phase_change(hid_t plist_id, unsigned *max_list, unsigned *min_btree);



 herr_t H5Pset_alignment(hid_t fapl_id, hsize_t threshold,
    hsize_t alignment);
 herr_t H5Pget_alignment(hid_t fapl_id, hsize_t *threshold ,
    hsize_t *alignment );
 herr_t H5Pset_driver(hid_t plist_id, hid_t driver_id,
        const void *driver_info);
 hid_t H5Pget_driver(hid_t plist_id);
 void *H5Pget_driver_info(hid_t plist_id);
 herr_t H5Pset_family_offset(hid_t fapl_id, hsize_t offset);
 herr_t H5Pget_family_offset(hid_t fapl_id, hsize_t *offset);
 herr_t H5Pset_multi_type(hid_t fapl_id, H5FD_mem_t type);
 herr_t H5Pget_multi_type(hid_t fapl_id, H5FD_mem_t *type);
 herr_t H5Pset_cache(hid_t plist_id, int mdc_nelmts,
       size_t rdcc_nslots, size_t rdcc_nbytes,
       double rdcc_w0);
 herr_t H5Pget_cache(hid_t plist_id,
       int *mdc_nelmts,
       size_t *rdcc_nslots ,
       size_t *rdcc_nbytes , double *rdcc_w0);
 herr_t H5Pset_mdc_config(hid_t plist_id,
       H5AC_cache_config_t * config_ptr);
 herr_t H5Pget_mdc_config(hid_t plist_id,
       H5AC_cache_config_t * config_ptr);
 herr_t H5Pset_gc_references(hid_t fapl_id, unsigned gc_ref);
 herr_t H5Pget_gc_references(hid_t fapl_id, unsigned *gc_ref );
 herr_t H5Pset_fclose_degree(hid_t fapl_id, H5F_close_degree_t degree);
 herr_t H5Pget_fclose_degree(hid_t fapl_id, H5F_close_degree_t *degree);
 herr_t H5Pset_meta_block_size(hid_t fapl_id, hsize_t size);
 herr_t H5Pget_meta_block_size(hid_t fapl_id, hsize_t *size );
 herr_t H5Pset_sieve_buf_size(hid_t fapl_id, size_t size);
 herr_t H5Pget_sieve_buf_size(hid_t fapl_id, size_t *size );
 herr_t H5Pset_small_data_block_size(hid_t fapl_id, hsize_t size);
 herr_t H5Pget_small_data_block_size(hid_t fapl_id, hsize_t *size );
 herr_t H5Pset_libver_bounds(hid_t plist_id, H5F_libver_t low,
    H5F_libver_t high);
 herr_t H5Pget_libver_bounds(hid_t plist_id, H5F_libver_t *low,
    H5F_libver_t *high);


 herr_t H5Pset_layout(hid_t plist_id, H5D_layout_t layout);
 H5D_layout_t H5Pget_layout(hid_t plist_id);
 herr_t H5Pset_chunk(hid_t plist_id, int ndims, const hsize_t dim[ ]);
 int H5Pget_chunk(hid_t plist_id, int max_ndims, hsize_t dim[] );
 herr_t H5Pset_external(hid_t plist_id, const char *name, off_t offset,
          hsize_t size);
 int H5Pget_external_count(hid_t plist_id);
 herr_t H5Pget_external(hid_t plist_id, unsigned idx, size_t name_size,
          char *name , off_t *offset ,
          hsize_t *size );
 herr_t H5Pset_szip(hid_t plist_id, unsigned options_mask, unsigned pixels_per_block);
 herr_t H5Pset_shuffle(hid_t plist_id);
 herr_t H5Pset_nbit(hid_t plist_id);
 herr_t H5Pset_scaleoffset(hid_t plist_id, H5Z_SO_scale_type_t scale_type, int scale_factor);
 herr_t H5Pset_fill_value(hid_t plist_id, hid_t type_id,
     const void *value);
 herr_t H5Pget_fill_value(hid_t plist_id, hid_t type_id,
     void *value );
 herr_t H5Pfill_value_defined(hid_t plist, H5D_fill_value_t *status);
 herr_t H5Pset_alloc_time(hid_t plist_id, H5D_alloc_time_t
 alloc_time);
 herr_t H5Pget_alloc_time(hid_t plist_id, H5D_alloc_time_t
 *alloc_time );
 herr_t H5Pset_fill_time(hid_t plist_id, H5D_fill_time_t fill_time);
 herr_t H5Pget_fill_time(hid_t plist_id, H5D_fill_time_t
 *fill_time );


 herr_t H5Pset_chunk_cache(hid_t dapl_id, size_t rdcc_nslots,
       size_t rdcc_nbytes, double rdcc_w0);
 herr_t H5Pget_chunk_cache(hid_t dapl_id,
       size_t *rdcc_nslots ,
       size_t *rdcc_nbytes ,
       double *rdcc_w0 );


 herr_t H5Pset_data_transform(hid_t plist_id, const char* expression);
 ssize_t H5Pget_data_transform(hid_t plist_id, char* expression , size_t size);
 herr_t H5Pset_buffer(hid_t plist_id, size_t size, void *tconv,
        void *bkg);
 size_t H5Pget_buffer(hid_t plist_id, void **tconv ,
        void **bkg );
 herr_t H5Pset_preserve(hid_t plist_id, hbool_t status);
 int H5Pget_preserve(hid_t plist_id);
 herr_t H5Pset_edc_check(hid_t plist_id, H5Z_EDC_t check);
 H5Z_EDC_t H5Pget_edc_check(hid_t plist_id);
 herr_t H5Pset_filter_callback(hid_t plist_id, H5Z_filter_func_t func,
                                     void* op_data);
 herr_t H5Pset_btree_ratios(hid_t plist_id, double left, double middle,
       double right);
 herr_t H5Pget_btree_ratios(hid_t plist_id, double *left ,
       double *middle ,
       double *right );
 herr_t H5Pset_vlen_mem_manager(hid_t plist_id,
                                       H5MM_allocate_t alloc_func,
                                       void *alloc_info, H5MM_free_t free_func,
                                       void *free_info);
 herr_t H5Pget_vlen_mem_manager(hid_t plist_id,
                                       H5MM_allocate_t *alloc_func,
                                       void **alloc_info,
                                       H5MM_free_t *free_func,
                                       void **free_info);
 herr_t H5Pset_hyper_vector_size(hid_t fapl_id, size_t size);
 herr_t H5Pget_hyper_vector_size(hid_t fapl_id, size_t *size );
 herr_t H5Pset_type_conv_cb(hid_t dxpl_id, H5T_conv_except_func_t op, void* operate_data);
 herr_t H5Pget_type_conv_cb(hid_t dxpl_id, H5T_conv_except_func_t *op, void** operate_data);


 herr_t H5Pset_create_intermediate_group(hid_t plist_id, unsigned crt_intmd);
 herr_t H5Pget_create_intermediate_group(hid_t plist_id, unsigned *crt_intmd );


 herr_t H5Pset_local_heap_size_hint(hid_t plist_id, size_t size_hint);
 herr_t H5Pget_local_heap_size_hint(hid_t plist_id, size_t *size_hint );
 herr_t H5Pset_link_phase_change(hid_t plist_id, unsigned max_compact, unsigned min_dense);
 herr_t H5Pget_link_phase_change(hid_t plist_id, unsigned *max_compact , unsigned *min_dense );
 herr_t H5Pset_est_link_info(hid_t plist_id, unsigned est_num_entries, unsigned est_name_len);
 herr_t H5Pget_est_link_info(hid_t plist_id, unsigned *est_num_entries , unsigned *est_name_len );
 herr_t H5Pset_link_creation_order(hid_t plist_id, unsigned crt_order_flags);
 herr_t H5Pget_link_creation_order(hid_t plist_id, unsigned *crt_order_flags );


 herr_t H5Pset_char_encoding(hid_t plist_id, H5T_cset_t encoding);
 herr_t H5Pget_char_encoding(hid_t plist_id, H5T_cset_t *encoding );


 herr_t H5Pset_nlinks(hid_t plist_id, size_t nlinks);
 herr_t H5Pget_nlinks(hid_t plist_id, size_t *nlinks);
 herr_t H5Pset_elink_prefix(hid_t plist_id, const char *prefix);
 ssize_t H5Pget_elink_prefix(hid_t plist_id, char *prefix, size_t size);
 hid_t H5Pget_elink_fapl(hid_t lapl_id);
 herr_t H5Pset_elink_fapl(hid_t lapl_id, hid_t fapl_id);
 herr_t H5Pset_elink_acc_flags(hid_t lapl_id, unsigned flags);
 herr_t H5Pget_elink_acc_flags(hid_t lapl_id, unsigned *flags);
 herr_t H5Pset_elink_cb(hid_t lapl_id, H5L_elink_traverse_t func, void *op_data);
 herr_t H5Pget_elink_cb(hid_t lapl_id, H5L_elink_traverse_t *func, void **op_data);


 herr_t H5Pset_copy_object(hid_t plist_id, unsigned crt_intmd);
 herr_t H5Pget_copy_object(hid_t plist_id, unsigned *crt_intmd );
# 410 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Ppublic.h"
 herr_t H5Pregister1(hid_t cls_id, const char *name, size_t size,
    void *def_value, H5P_prp_create_func_t prp_create,
    H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
    H5P_prp_delete_func_t prp_del, H5P_prp_copy_func_t prp_copy,
    H5P_prp_close_func_t prp_close);
 herr_t H5Pinsert1(hid_t plist_id, const char *name, size_t size,
    void *value, H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
    H5P_prp_delete_func_t prp_delete, H5P_prp_copy_func_t prp_copy,
    H5P_prp_close_func_t prp_close);
 H5Z_filter_t H5Pget_filter1(hid_t plist_id, unsigned filter,
    unsigned int *flags , size_t *cd_nelmts ,
    unsigned cd_values[] , size_t namelen, char name[]);
 herr_t H5Pget_filter_by_id1(hid_t plist_id, H5Z_filter_t id,
    unsigned int *flags , size_t *cd_nelmts ,
    unsigned cd_values[] , size_t namelen, char name[] );
# 37 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Rpublic.h" 1
# 30 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Rpublic.h"
typedef enum {
    H5R_BADTYPE = (-1),
    H5R_OBJECT,
    H5R_DATASET_REGION,
    H5R_MAXTYPE
} H5R_type_t;
# 44 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Rpublic.h"
typedef haddr_t hobj_ref_t;






typedef unsigned char hdset_reg_ref_t[(sizeof(haddr_t)+4)];
# 61 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Rpublic.h"
 herr_t H5Rcreate(void *ref, hid_t loc_id, const char *name,
    H5R_type_t ref_type, hid_t space_id);
 hid_t H5Rdereference(hid_t dataset, H5R_type_t ref_type, const void *ref);
 hid_t H5Rget_region(hid_t dataset, H5R_type_t ref_type, const void *ref);
 herr_t H5Rget_obj_type2(hid_t id, H5R_type_t ref_type, const void *_ref,
    H5O_type_t *obj_type);
 ssize_t H5Rget_name(hid_t loc_id, H5R_type_t ref_type, const void *ref,
    char *name , size_t size);
# 83 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Rpublic.h"
 H5G_obj_t H5Rget_obj_type1(hid_t id, H5R_type_t ref_type, const void *_ref);
# 38 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Spublic.h" 1
# 34 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Spublic.h"
typedef enum H5S_class_t {
    H5S_NO_CLASS = -1,
    H5S_SCALAR = 0,
    H5S_SIMPLE = 1,
    H5S_NULL = 2
} H5S_class_t;


typedef enum H5S_seloper_t {
    H5S_SELECT_NOOP = -1,
    H5S_SELECT_SET = 0,
    H5S_SELECT_OR,





    H5S_SELECT_AND,





    H5S_SELECT_XOR,





    H5S_SELECT_NOTB,





    H5S_SELECT_NOTA,





    H5S_SELECT_APPEND,
    H5S_SELECT_PREPEND,
    H5S_SELECT_INVALID
} H5S_seloper_t;


typedef enum {
    H5S_SEL_ERROR = -1,
    H5S_SEL_NONE = 0,
    H5S_SEL_POINTS = 1,
    H5S_SEL_HYPERSLABS = 2,
    H5S_SEL_ALL = 3,
    H5S_SEL_N
}H5S_sel_type;






 hid_t H5Screate(H5S_class_t type);
 hid_t H5Screate_simple(int rank, const hsize_t dims[],
          const hsize_t maxdims[]);
 herr_t H5Sset_extent_simple(hid_t space_id, int rank,
        const hsize_t dims[],
        const hsize_t max[]);
 hid_t H5Scopy(hid_t space_id);
 herr_t H5Sclose(hid_t space_id);
 herr_t H5Sencode(hid_t obj_id, void *buf, size_t *nalloc);
 hid_t H5Sdecode(const void *buf);
 hssize_t H5Sget_simple_extent_npoints(hid_t space_id);
 int H5Sget_simple_extent_ndims(hid_t space_id);
 int H5Sget_simple_extent_dims(hid_t space_id, hsize_t dims[],
          hsize_t maxdims[]);
 htri_t H5Sis_simple(hid_t space_id);
 hssize_t H5Sget_select_npoints(hid_t spaceid);
 herr_t H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op,
       const hsize_t start[],
       const hsize_t _stride[],
       const hsize_t count[],
       const hsize_t _block[]);
# 128 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5Spublic.h"
 herr_t H5Sselect_elements(hid_t space_id, H5S_seloper_t op,
    size_t num_elem, const hsize_t *coord);
 H5S_class_t H5Sget_simple_extent_type(hid_t space_id);
 herr_t H5Sset_extent_none(hid_t space_id);
 herr_t H5Sextent_copy(hid_t dst_id,hid_t src_id);
 htri_t H5Sextent_equal(hid_t sid1, hid_t sid2);
 herr_t H5Sselect_all(hid_t spaceid);
 herr_t H5Sselect_none(hid_t spaceid);
 herr_t H5Soffset_simple(hid_t space_id, const hssize_t *offset);
 htri_t H5Sselect_valid(hid_t spaceid);
 hssize_t H5Sget_select_hyper_nblocks(hid_t spaceid);
 hssize_t H5Sget_select_elem_npoints(hid_t spaceid);
 herr_t H5Sget_select_hyper_blocklist(hid_t spaceid, hsize_t startblock,
    hsize_t numblocks, hsize_t buf[ ]);
 herr_t H5Sget_select_elem_pointlist(hid_t spaceid, hsize_t startpoint,
    hsize_t numpoints, hsize_t buf[ ]);
 herr_t H5Sget_select_bounds(hid_t spaceid, hsize_t start[],
    hsize_t end[]);
 H5S_sel_type H5Sget_select_type(hid_t spaceid);
# 39 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2




# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDcore.h" 1
# 32 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDcore.h"
 hid_t H5FD_core_init(void);
 void H5FD_core_term(void);
 herr_t H5Pset_fapl_core(hid_t fapl_id, size_t increment,
    hbool_t backing_store);
 herr_t H5Pget_fapl_core(hid_t fapl_id, size_t *increment ,
    hbool_t *backing_store );
# 44 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDfamily.h" 1
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDfamily.h"
 hid_t H5FD_family_init(void);
 void H5FD_family_term(void);
 herr_t H5Pset_fapl_family(hid_t fapl_id, hsize_t memb_size,
     hid_t memb_fapl_id);
 herr_t H5Pget_fapl_family(hid_t fapl_id, hsize_t *memb_size ,
     hid_t *memb_fapl_id );
# 45 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDlog.h" 1
# 61 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDlog.h"
 hid_t H5FD_log_init(void);
 void H5FD_log_term(void);
 herr_t H5Pset_fapl_log(hid_t fapl_id, const char *logfile, unsigned flags, size_t buf_size);
# 46 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpi.h" 1
# 40 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpi.h"
typedef enum H5FD_mpio_xfer_t {
    H5FD_MPIO_INDEPENDENT = 0,
    H5FD_MPIO_COLLECTIVE
} H5FD_mpio_xfer_t;


typedef enum H5FD_mpio_chunk_opt_t {
    H5FD_MPIO_CHUNK_DEFAULT = 0,
    H5FD_MPIO_CHUNK_ONE_IO,
    H5FD_MPIO_CHUNK_MULTI_IO
} H5FD_mpio_chunk_opt_t;


typedef enum H5FD_mpio_collective_opt_t {
    H5FD_MPIO_COLLECTIVE_IO = 0,
    H5FD_MPIO_INDIVIDUAL_IO
} H5FD_mpio_collective_opt_t;





typedef struct H5FD_class_mpi_t {
    H5FD_class_t super;
    int (*get_rank)(const H5FD_t *file);
    int (*get_size)(const H5FD_t *file);
    MPI_Comm (*get_comm)(const H5FD_t *file);
} H5FD_class_mpi_t;



# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpio.h" 1
# 48 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpio.h"
 hid_t H5FD_mpio_init(void);
 void H5FD_mpio_term(void);
 herr_t H5Pset_fapl_mpio(hid_t fapl_id, MPI_Comm comm, MPI_Info info);
 herr_t H5Pget_fapl_mpio(hid_t fapl_id, MPI_Comm *comm ,
   MPI_Info *info );
 herr_t H5Pset_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode);
 herr_t H5Pget_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t *xfer_mode );
 herr_t H5Pset_dxpl_mpio_collective_opt(hid_t dxpl_id, H5FD_mpio_collective_opt_t opt_mode);
 herr_t H5Pset_dxpl_mpio_chunk_opt(hid_t dxpl_id, H5FD_mpio_chunk_opt_t opt_mode);
 herr_t H5Pset_dxpl_mpio_chunk_opt_num(hid_t dxpl_id, unsigned num_chunk_per_proc);
 herr_t H5Pset_dxpl_mpio_chunk_opt_ratio(hid_t dxpl_id, unsigned percent_num_proc_per_chunk);
# 72 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpi.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpiposix.h" 1
# 44 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpiposix.h"
 hid_t H5FD_mpiposix_init(void);
 void H5FD_mpiposix_term(void);
 herr_t H5Pset_fapl_mpiposix(hid_t fapl_id, MPI_Comm comm, hbool_t use_gpfs);
 herr_t H5Pget_fapl_mpiposix(hid_t fapl_id, MPI_Comm *comm , hbool_t *use_gpfs );
# 73 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpi.h" 2
# 92 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmpi.h"
extern char H5FD_mpi_native_g[];






 haddr_t H5FD_mpi_MPIOff_to_haddr(MPI_Offset mpi_off);
 herr_t H5FD_mpi_haddr_to_MPIOff(haddr_t addr, MPI_Offset *mpi_off );
 herr_t H5FD_mpi_comm_info_dup(MPI_Comm comm, MPI_Info info,
    MPI_Comm *comm_new, MPI_Info *info_new);
 herr_t H5FD_mpi_comm_info_free(MPI_Comm *comm, MPI_Info *info);




 herr_t H5FD_mpi_setup_collective(hid_t dxpl_id, MPI_Datatype btype,
    MPI_Datatype ftype);
 herr_t H5FD_mpi_teardown_collective(hid_t dxpl_id);


 int H5FD_mpi_get_rank(const H5FD_t *file);
 int H5FD_mpi_get_size(const H5FD_t *file);
 MPI_Comm H5FD_mpi_get_comm(const H5FD_t *_file);
# 47 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmulti.h" 1
# 34 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDmulti.h"
 hid_t H5FD_multi_init(void);
 void H5FD_multi_term(void);
 herr_t H5Pset_fapl_multi(hid_t fapl_id, const H5FD_mem_t *memb_map,
    const hid_t *memb_fapl, const char * const *memb_name,
    const haddr_t *memb_addr, hbool_t relax);
 herr_t H5Pget_fapl_multi(hid_t fapl_id, H5FD_mem_t *memb_map ,
    hid_t *memb_fapl , char **memb_name ,
    haddr_t *memb_addr , hbool_t *relax );
 herr_t H5Pset_dxpl_multi(hid_t dxpl_id, const hid_t *memb_dxpl);
 herr_t H5Pget_dxpl_multi(hid_t dxpl_id, hid_t *memb_dxpl );

 herr_t H5Pset_fapl_split(hid_t fapl, const char *meta_ext,
    hid_t meta_plist_id, const char *raw_ext,
    hid_t raw_plist_id);
# 48 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDsec2.h" 1
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDsec2.h"
 hid_t H5FD_sec2_init(void);
 void H5FD_sec2_term(void);
 herr_t H5Pset_fapl_sec2(hid_t fapl_id);
# 49 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDstdio.h" 1
# 33 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDstdio.h"
 hid_t H5FD_stdio_init(void);
 void H5FD_stdio_term(void);
 herr_t H5Pset_fapl_stdio(hid_t fapl_id);
# 50 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2



# 1 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/H5FDdirect.h" 1
# 54 "/home/dpnkarthik/petsc-rnet/PETSC_RNET/include/hdf5.h" 2
# 278 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h" 2
extern PetscErrorCode PetscViewerHDF5GetFileId(PetscViewer,hid_t*);
extern PetscErrorCode PetscViewerHDF5OpenGroup(PetscViewer, hid_t *, hid_t *);






extern PetscViewer PETSC_VIEWER_STDOUT_(MPI_Comm);
extern PetscErrorCode PetscViewerASCIIGetStdout(MPI_Comm,PetscViewer*);
extern PetscViewer PETSC_VIEWER_STDERR_(MPI_Comm);
extern PetscErrorCode PetscViewerASCIIGetStderr(MPI_Comm,PetscViewer*);
extern PetscViewer PETSC_VIEWER_DRAW_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_SOCKET_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_BINARY_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_MATLAB_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE;
# 386 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
extern PetscErrorCode PetscViewerMatlabPutArray(PetscViewer,int,int,const PetscScalar*,const char*);
extern PetscErrorCode PetscViewerMatlabGetArray(PetscViewer,int,int,PetscScalar*,const char*);
extern PetscErrorCode PetscViewerMatlabPutVariable(PetscViewer,const char*,void*);
# 400 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"
typedef struct _n_PetscViewers* PetscViewers;
extern PetscErrorCode PetscViewersCreate(MPI_Comm,PetscViewers*);
extern PetscErrorCode PetscViewersDestroy(PetscViewers*);
extern PetscErrorCode PetscViewersGetViewer(PetscViewers,PetscInt,PetscViewer*);
# 417 "/home/dpnkarthik/petsc-rnet/include/petscviewer.h"

# 1271 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/include/petscoptions.h" 1








extern PetscErrorCode PetscOptionsHasName(const char[],const char[],PetscBool *);

extern PetscErrorCode PetscOptionsGetInt(const char[],const char [],PetscInt *,PetscBool *);

extern PetscErrorCode PetscOptionsGetBool(const char[],const char [],PetscBool *,PetscBool *);

extern PetscErrorCode PetscOptionsGetReal(const char[],const char[],PetscReal *,PetscBool *);

extern PetscErrorCode PetscOptionsGetScalar(const char[],const char[],PetscScalar *,PetscBool *);

extern PetscErrorCode PetscOptionsGetIntArray(const char[],const char[],PetscInt[],PetscInt *,PetscBool *);

extern PetscErrorCode PetscOptionsGetRealArray(const char[],const char[],PetscReal[],PetscInt *,PetscBool *);

extern PetscErrorCode PetscOptionsGetBoolArray(const char[],const char[],PetscBool [],PetscInt *,PetscBool *);

extern PetscErrorCode PetscOptionsGetString(const char[],const char[],char[],size_t,PetscBool *);

extern PetscErrorCode PetscOptionsGetStringArray(const char[],const char[],char*[],PetscInt*,PetscBool *);

extern PetscErrorCode PetscOptionsGetEList(const char[],const char[],const char*const*,PetscInt,PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsGetEnum(const char[],const char[],const char*const*,PetscEnum*,PetscBool *);
extern PetscErrorCode PetscOptionsValidKey(const char[],PetscBool *);

extern PetscErrorCode PetscOptionsSetAlias(const char[],const char[]);
extern PetscErrorCode PetscOptionsSetValue(const char[],const char[]);
extern PetscErrorCode PetscOptionsClearValue(const char[]);

extern PetscErrorCode PetscOptionsAllUsed(PetscInt*);
extern PetscErrorCode PetscOptionsLeft(void);
extern PetscErrorCode PetscOptionsView(PetscViewer);

extern PetscErrorCode PetscOptionsCreate(void);
extern PetscErrorCode PetscOptionsInsert(int*,char ***,const char[]);
extern PetscErrorCode PetscOptionsInsertFile(MPI_Comm,const char[],PetscBool );



extern PetscErrorCode PetscOptionsInsertString(const char[]);
extern PetscErrorCode PetscOptionsDestroy(void);
extern PetscErrorCode PetscOptionsClear(void);
extern PetscErrorCode PetscOptionsPrefixPush(const char[]);
extern PetscErrorCode PetscOptionsPrefixPop(void);

extern PetscErrorCode PetscOptionsReject(const char[],const char[]);
extern PetscErrorCode PetscOptionsGetAll(char*[]);

extern PetscErrorCode PetscOptionsGetenv(MPI_Comm,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOptionsStringToInt(const char[],PetscInt*);
extern PetscErrorCode PetscOptionsStringToReal(const char[],PetscReal*);
extern PetscErrorCode PetscOptionsStringToBool(const char[],PetscBool*);

extern PetscErrorCode PetscOptionsMonitorSet(PetscErrorCode (*)(const char[], const char[], void*), void *, PetscErrorCode (*)(void**));
extern PetscErrorCode PetscOptionsMonitorCancel(void);
extern PetscErrorCode PetscOptionsMonitorDefault(const char[], const char[], void *);

extern PetscBool PetscOptionsPublish;
extern PetscInt PetscOptionsPublishCount;
# 166 "/home/dpnkarthik/petsc-rnet/include/petscoptions.h"
extern PetscErrorCode PetscOptionsBegin_Private(MPI_Comm,const char[],const char[],const char[]);
extern PetscErrorCode PetscObjectOptionsBegin_Private(PetscObject);
extern PetscErrorCode PetscOptionsEnd_Private(void);
extern PetscErrorCode PetscOptionsHead(const char[]);
# 202 "/home/dpnkarthik/petsc-rnet/include/petscoptions.h"
extern PetscErrorCode PetscOptionsEnum(const char[],const char[],const char[],const char *const*,PetscEnum,PetscEnum*,PetscBool *);
extern PetscErrorCode PetscOptionsInt(const char[],const char[],const char[],PetscInt,PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsReal(const char[],const char[],const char[],PetscReal,PetscReal*,PetscBool *);
extern PetscErrorCode PetscOptionsScalar(const char[],const char[],const char[],PetscScalar,PetscScalar*,PetscBool *);
extern PetscErrorCode PetscOptionsName(const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsString(const char[],const char[],const char[],const char[],char*,size_t,PetscBool *);
extern PetscErrorCode PetscOptionsBool(const char[],const char[],const char[],PetscBool ,PetscBool *,PetscBool *);
extern PetscErrorCode PetscOptionsBoolGroupBegin(const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsBoolGroup(const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsBoolGroupEnd(const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsList(const char[],const char[],const char[],PetscFList,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOptionsEList(const char[],const char[],const char[],const char*const*,PetscInt,const char[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsRealArray(const char[],const char[],const char[],PetscReal[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsIntArray(const char[],const char[],const char[],PetscInt[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsStringArray(const char[],const char[],const char[],char*[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsBoolArray(const char[],const char[],const char[],PetscBool [],PetscInt*,PetscBool *);

extern PetscErrorCode PetscOptionsSetFromOptions(void);
extern PetscErrorCode PetscOptionsAMSDestroy(void);





typedef enum {OPTION_INT,OPTION_LOGICAL,OPTION_REAL,OPTION_LIST,OPTION_STRING,OPTION_REAL_ARRAY,OPTION_HEAD,OPTION_INT_ARRAY,OPTION_ELIST,OPTION_LOGICAL_ARRAY,OPTION_STRING_ARRAY} PetscOptionType;
typedef struct _n_PetscOptions* PetscOptions;
struct _n_PetscOptions {
  char *option;
  char *text;
  void *data;
  PetscFList flist;
  const char *const *list;
  char nlist;
  char *man;
  size_t arraylength;
  PetscBool set;
  PetscOptionType type;
  PetscOptions next;
  char *pman;
  void *edata;
};

typedef struct {
  PetscOptions next;
  char *prefix,*pprefix;
  char *title;
  MPI_Comm comm;
  PetscBool printhelp,changedmethod,alreadyprinted;
  PetscObject object;
} PetscOptionsObjectType;
# 1272 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2


extern PetscClassId PETSC_LARGEST_CLASSID;
extern PetscClassId PETSC_OBJECT_CLASSID;
extern PetscErrorCode PetscClassIdRegister(const char[],PetscClassId *);




extern PetscErrorCode PetscMemoryGetCurrentUsage(PetscLogDouble *);
extern PetscErrorCode PetscMemoryGetMaximumUsage(PetscLogDouble *);
extern PetscErrorCode PetscMemorySetGetMaximumUsage(void);
extern PetscErrorCode PetscMemoryShowUsage(PetscViewer,const char[]);

extern PetscErrorCode PetscInfoAllow(PetscBool ,const char []);
extern PetscErrorCode PetscGetTime(PetscLogDouble*);
extern PetscErrorCode PetscGetCPUTime(PetscLogDouble*);
extern PetscErrorCode PetscSleep(PetscReal);




extern PetscErrorCode PetscInitialize(int*,char***,const char[],const char[]);

extern PetscErrorCode PetscInitializeNoArguments(void);
extern PetscErrorCode PetscInitialized(PetscBool *);
extern PetscErrorCode PetscFinalized(PetscBool *);
extern PetscErrorCode PetscFinalize(void);
extern PetscErrorCode PetscInitializeFortran(void);
extern PetscErrorCode PetscGetArgs(int*,char ***);
extern PetscErrorCode PetscGetArguments(char ***);
extern PetscErrorCode PetscFreeArguments(char **);

extern PetscErrorCode PetscEnd(void);
extern PetscErrorCode PetscSysInitializePackage(const char[]);

extern MPI_Comm PETSC_COMM_LOCAL_WORLD;
extern PetscErrorCode PetscHMPIMerge(PetscMPIInt,PetscErrorCode (*)(void*),void*);
extern PetscErrorCode PetscHMPISpawn(PetscMPIInt);
extern PetscErrorCode PetscHMPIFinalize(void);
extern PetscErrorCode PetscHMPIRun(MPI_Comm,PetscErrorCode (*)(MPI_Comm,void *),void*);
extern PetscErrorCode PetscHMPIRunCtx(MPI_Comm,PetscErrorCode (*)(MPI_Comm,void*,void *),void*);
extern PetscErrorCode PetscHMPIFree(MPI_Comm,void*);
extern PetscErrorCode PetscHMPIMalloc(MPI_Comm,size_t,void**);

extern PetscErrorCode PetscPythonInitialize(const char[],const char[]);
extern PetscErrorCode PetscPythonFinalize(void);
extern PetscErrorCode PetscPythonPrintError(void);
extern PetscErrorCode PetscPythonMonitorSet(PetscObject,const char[]);






typedef void (**PetscVoidStarFunction)(void);
typedef void (*PetscVoidFunction)(void);
typedef PetscErrorCode (*PetscErrorCodeFunction)(void);
# 1363 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscObjectDestroy(PetscObject*);
extern PetscErrorCode PetscObjectGetComm(PetscObject,MPI_Comm *);
extern PetscErrorCode PetscObjectGetClassId(PetscObject,PetscClassId *);
extern PetscErrorCode PetscObjectGetClassName(PetscObject,const char *[]);
extern PetscErrorCode PetscObjectSetType(PetscObject,const char []);
extern PetscErrorCode PetscObjectSetPrecision(PetscObject,PetscPrecision);
extern PetscErrorCode PetscObjectGetType(PetscObject,const char *[]);
extern PetscErrorCode PetscObjectSetName(PetscObject,const char[]);
extern PetscErrorCode PetscObjectGetName(PetscObject,const char*[]);
extern PetscErrorCode PetscObjectPrintClassNamePrefixType(PetscObject,PetscViewer,const char[]);
extern PetscErrorCode PetscObjectSetTabLevel(PetscObject,PetscInt);
extern PetscErrorCode PetscObjectGetTabLevel(PetscObject,PetscInt*);
extern PetscErrorCode PetscObjectIncrementTabLevel(PetscObject,PetscObject,PetscInt);
extern PetscErrorCode PetscObjectReference(PetscObject);
extern PetscErrorCode PetscObjectGetReference(PetscObject,PetscInt*);
extern PetscErrorCode PetscObjectDereference(PetscObject);
extern PetscErrorCode PetscObjectGetNewTag(PetscObject,PetscMPIInt *);
extern PetscErrorCode PetscObjectView(PetscObject,PetscViewer);
extern PetscErrorCode PetscObjectCompose(PetscObject,const char[],PetscObject);
extern PetscErrorCode PetscObjectRemoveReference(PetscObject,const char[]);
extern PetscErrorCode PetscObjectQuery(PetscObject,const char[],PetscObject *);
extern PetscErrorCode PetscObjectComposeFunction(PetscObject,const char[],const char[],void (*)(void));
extern PetscErrorCode PetscObjectSetFromOptions(PetscObject);
extern PetscErrorCode PetscObjectSetUp(PetscObject);
extern PetscErrorCode PetscCommGetNewTag(MPI_Comm,PetscMPIInt *);
extern PetscErrorCode PetscObjectAddOptionsHandler(PetscObject,PetscErrorCode (*)(PetscObject,void*),PetscErrorCode (*)(PetscObject,void*),void*);
extern PetscErrorCode PetscObjectProcessOptionsHandlers(PetscObject);
extern PetscErrorCode PetscObjectDestroyOptionsHandlers(PetscObject);
# 1434 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscObjectQueryFunction(PetscObject,const char[],void (**)(void));
extern PetscErrorCode PetscObjectSetOptionsPrefix(PetscObject,const char[]);
extern PetscErrorCode PetscObjectAppendOptionsPrefix(PetscObject,const char[]);
extern PetscErrorCode PetscObjectPrependOptionsPrefix(PetscObject,const char[]);
extern PetscErrorCode PetscObjectGetOptionsPrefix(PetscObject,const char*[]);
extern PetscErrorCode PetscObjectAMSPublish(PetscObject);
extern PetscErrorCode PetscObjectUnPublish(PetscObject);
extern PetscErrorCode PetscObjectChangeTypeName(PetscObject,const char[]);
extern PetscErrorCode PetscObjectRegisterDestroy(PetscObject);
extern PetscErrorCode PetscObjectRegisterDestroyAll(void);
extern PetscErrorCode PetscObjectName(PetscObject);
extern PetscErrorCode PetscTypeCompare(PetscObject,const char[],PetscBool *);
extern PetscErrorCode PetscTypeCompareAny(PetscObject,PetscBool*,const char[],...);
extern PetscErrorCode PetscRegisterFinalize(PetscErrorCode (*)(void));
extern PetscErrorCode PetscRegisterFinalizeAll(void);




# 1 "/home/dpnkarthik/petsc-rnet/include/petscerror.h" 1







# 1 "/usr/include/string.h" 1 3 4
# 29 "/usr/include/string.h" 3 4





# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 35 "/usr/include/string.h" 2 3 4









extern void *memcpy (void *__restrict __dest,
       __const void *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern void *memmove (void *__dest, __const void *__src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));






extern void *memccpy (void *__restrict __dest, __const void *__restrict __src,
        int __c, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));





extern void *memset (void *__s, int __c, size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int memcmp (__const void *__s1, __const void *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 95 "/usr/include/string.h" 3 4
extern void *memchr (__const void *__s, int __c, size_t __n)
      __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));


# 126 "/usr/include/string.h" 3 4


extern char *strcpy (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern char *strncpy (char *__restrict __dest,
        __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern char *strcat (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern char *strncat (char *__restrict __dest, __const char *__restrict __src,
        size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strcmp (__const char *__s1, __const char *__s2)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));

extern int strncmp (__const char *__s1, __const char *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strcoll (__const char *__s1, __const char *__s2)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));

extern size_t strxfrm (char *__restrict __dest,
         __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));






# 1 "/usr/include/xlocale.h" 1 3 4
# 28 "/usr/include/xlocale.h" 3 4
typedef struct __locale_struct
{

  struct __locale_data *__locales[13];


  const unsigned short int *__ctype_b;
  const int *__ctype_tolower;
  const int *__ctype_toupper;


  const char *__names[13];
} *__locale_t;


typedef __locale_t locale_t;
# 163 "/usr/include/string.h" 2 3 4


extern int strcoll_l (__const char *__s1, __const char *__s2, __locale_t __l)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2, 3)));

extern size_t strxfrm_l (char *__dest, __const char *__src, size_t __n,
    __locale_t __l) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 4)));





extern char *strdup (__const char *__s)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) __attribute__ ((__nonnull__ (1)));






extern char *strndup (__const char *__string, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) __attribute__ ((__nonnull__ (1)));
# 210 "/usr/include/string.h" 3 4

# 235 "/usr/include/string.h" 3 4
extern char *strchr (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));
# 262 "/usr/include/string.h" 3 4
extern char *strrchr (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));


# 281 "/usr/include/string.h" 3 4



extern size_t strcspn (__const char *__s, __const char *__reject)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern size_t strspn (__const char *__s, __const char *__accept)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 314 "/usr/include/string.h" 3 4
extern char *strpbrk (__const char *__s, __const char *__accept)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 342 "/usr/include/string.h" 3 4
extern char *strstr (__const char *__haystack, __const char *__needle)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));




extern char *strtok (char *__restrict __s, __const char *__restrict __delim)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));




extern char *__strtok_r (char *__restrict __s,
    __const char *__restrict __delim,
    char **__restrict __save_ptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 3)));

extern char *strtok_r (char *__restrict __s, __const char *__restrict __delim,
         char **__restrict __save_ptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 3)));
# 397 "/usr/include/string.h" 3 4


extern size_t strlen (__const char *__s)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));





extern size_t strnlen (__const char *__string, size_t __maxlen)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));





extern char *strerror (int __errnum) __attribute__ ((__nothrow__));

# 427 "/usr/include/string.h" 3 4
extern int strerror_r (int __errnum, char *__buf, size_t __buflen) __asm__ ("" "__xpg_strerror_r") __attribute__ ((__nothrow__))

                        __attribute__ ((__nonnull__ (2)));
# 445 "/usr/include/string.h" 3 4
extern char *strerror_l (int __errnum, __locale_t __l) __attribute__ ((__nothrow__));





extern void __bzero (void *__s, size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));



extern void bcopy (__const void *__src, void *__dest, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern void bzero (void *__s, size_t __n) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern int bcmp (__const void *__s1, __const void *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 489 "/usr/include/string.h" 3 4
extern char *index (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));
# 517 "/usr/include/string.h" 3 4
extern char *rindex (__const char *__s, int __c)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));




extern int ffs (int __i) __attribute__ ((__nothrow__)) __attribute__ ((__const__));
# 536 "/usr/include/string.h" 3 4
extern int strcasecmp (__const char *__s1, __const char *__s2)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strncasecmp (__const char *__s1, __const char *__s2, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 559 "/usr/include/string.h" 3 4
extern char *strsep (char **__restrict __stringp,
       __const char *__restrict __delim)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));




extern char *strsignal (int __sig) __attribute__ ((__nothrow__));


extern char *__stpcpy (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern char *stpcpy (char *__restrict __dest, __const char *__restrict __src)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));



extern char *__stpncpy (char *__restrict __dest,
   __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern char *stpncpy (char *__restrict __dest,
        __const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
# 646 "/usr/include/string.h" 3 4

# 9 "/home/dpnkarthik/petsc-rnet/include/petscerror.h" 2



# 356 "/home/dpnkarthik/petsc-rnet/include/petscerror.h"
typedef enum {PETSC_ERROR_INITIAL=0,PETSC_ERROR_REPEAT=1,PETSC_ERROR_IN_CXX = 2} PetscErrorType;

extern PetscErrorCode PetscErrorPrintfInitialize(void);
extern PetscErrorCode PetscErrorMessage(int,const char*[],char **);
extern PetscErrorCode PetscTraceBackErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);




extern PetscErrorCode PetscIgnoreErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscEmacsClientErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscMPIAbortErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscAbortErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscAttachDebuggerErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscReturnErrorHandler(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscError(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,...);
extern PetscErrorCode PetscPushErrorHandler(PetscErrorCode (*handler)(MPI_Comm,int,const char*,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*),void*);
extern PetscErrorCode PetscPopErrorHandler(void);
extern PetscErrorCode PetscDefaultSignalHandler(int,void*);
extern PetscErrorCode PetscPushSignalHandler(PetscErrorCode (*)(int,void *),void*);
extern PetscErrorCode PetscPopSignalHandler(void);

typedef enum {PETSC_FP_TRAP_OFF=0,PETSC_FP_TRAP_ON=1} PetscFPTrap;
extern PetscErrorCode PetscSetFPTrap(PetscFPTrap);
# 389 "/home/dpnkarthik/petsc-rnet/include/petscerror.h"
typedef struct {
  const char *function[64];
  const char *file[64];
  const char *directory[64];
        int line[64];
        int currentsize;
} PetscStack;

extern PetscStack *petscstack;
extern PetscErrorCode PetscStackCopy(PetscStack*,PetscStack*);
extern PetscErrorCode PetscStackPrint(PetscStack*,FILE* fp);
# 552 "/home/dpnkarthik/petsc-rnet/include/petscerror.h"
extern PetscErrorCode PetscStackCreate(void);
extern PetscErrorCode PetscStackView(PetscViewer);
extern PetscErrorCode PetscStackDestroy(void);
extern PetscErrorCode PetscStackPublish(void);
extern PetscErrorCode PetscStackDepublish(void);



# 1454 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 1464 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _n_PetscOList *PetscOList;

extern PetscErrorCode PetscOListDestroy(PetscOList*);
extern PetscErrorCode PetscOListFind(PetscOList,const char[],PetscObject*);
extern PetscErrorCode PetscOListReverseFind(PetscOList,PetscObject,char**,PetscBool*);
extern PetscErrorCode PetscOListAdd(PetscOList *,const char[],PetscObject);
extern PetscErrorCode PetscOListRemoveReference(PetscOList *,const char[]);
extern PetscErrorCode PetscOListDuplicate(PetscOList,PetscOList *);





extern PetscErrorCode PetscFListAdd(PetscFList*,const char[],const char[],void (*)(void));
extern PetscErrorCode PetscFListDestroy(PetscFList*);
extern PetscErrorCode PetscFListFind(PetscFList,MPI_Comm,const char[],PetscBool,void (**)(void));
extern PetscErrorCode PetscFListPrintTypes(MPI_Comm,FILE*,const char[],const char[],const char[],const char[],PetscFList,const char[]);





extern PetscErrorCode PetscFListDuplicate(PetscFList,PetscFList *);
extern PetscErrorCode PetscFListView(PetscFList,PetscViewer);
extern PetscErrorCode PetscFListConcat(const char [],const char [],char []);
extern PetscErrorCode PetscFListGet(PetscFList,char ***,int*);
# 1500 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _n_PetscDLLibrary *PetscDLLibrary;
extern PetscDLLibrary DLLibrariesLoaded;
extern PetscErrorCode PetscDLLibraryAppend(MPI_Comm,PetscDLLibrary *,const char[]);
extern PetscErrorCode PetscDLLibraryPrepend(MPI_Comm,PetscDLLibrary *,const char[]);
extern PetscErrorCode PetscDLLibrarySym(MPI_Comm,PetscDLLibrary *,const char[],const char[],void **);
extern PetscErrorCode PetscDLLibraryPrintPath(PetscDLLibrary);
extern PetscErrorCode PetscDLLibraryRetrieve(MPI_Comm,const char[],char *,size_t,PetscBool *);
extern PetscErrorCode PetscDLLibraryOpen(MPI_Comm,const char[],PetscDLLibrary *);
extern PetscErrorCode PetscDLLibraryClose(PetscDLLibrary);
extern PetscErrorCode PetscDLLibraryCCAAppend(MPI_Comm,PetscDLLibrary *,const char[]);





# 1 "/home/dpnkarthik/petsc-rnet/include/petscfwk.h" 1





extern PetscClassId PETSC_FWK_CLASSID;
# 16 "/home/dpnkarthik/petsc-rnet/include/petscfwk.h"
struct _p_PetscFwk;
typedef struct _p_PetscFwk *PetscFwk;

extern PetscFwk PETSC_FWK_DEFAULT_(MPI_Comm);





extern PetscErrorCode PetscFwkInitializePackage(const char path[]);
extern PetscErrorCode PetscFwkFinalizePackage(void);


extern PetscErrorCode PetscFwkCall(PetscFwk component, const char *message);
extern PetscErrorCode PetscFwkGetURL(PetscFwk component, const char **outurl);
extern PetscErrorCode PetscFwkSetURL(PetscFwk component, const char *inurl);

extern PetscErrorCode PetscFwkCreate(MPI_Comm comm, PetscFwk *fwk);
extern PetscErrorCode PetscFwkView(PetscFwk fwk, PetscViewer viewerASCII);
extern PetscErrorCode PetscFwkRegisterComponent(PetscFwk fwk, const char key[]);
extern PetscErrorCode PetscFwkRegisterDependence(PetscFwk fwk, const char server_key[], const char client_key[]);
extern PetscErrorCode PetscFwkRegisterComponentURL(PetscFwk fwk, const char key[], const char url[]);
extern PetscErrorCode PetscFwkGetComponent(PetscFwk fwk, const char key[], PetscFwk *component, PetscBool *found);
extern PetscErrorCode PetscFwkGetParent(PetscFwk fwk, PetscFwk *parent);
extern PetscErrorCode PetscFwkVisit(PetscFwk fwk, const char *message);
extern PetscErrorCode PetscFwkDestroy(PetscFwk* fwk);
# 1516 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2




extern PetscErrorCode PetscSplitOwnership(MPI_Comm,PetscInt*,PetscInt*);
extern PetscErrorCode PetscSplitOwnershipBlock(MPI_Comm,PetscInt,PetscInt*,PetscInt*);
extern PetscErrorCode PetscSequentialPhaseBegin(MPI_Comm,PetscMPIInt);


extern PetscErrorCode PetscSequentialPhaseEnd(MPI_Comm,PetscMPIInt);


extern PetscErrorCode PetscBarrier(PetscObject);
extern PetscErrorCode PetscMPIDump(FILE*);
# 1544 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h" 1








extern PetscClassId PETSC_DRAW_CLASSID;
# 32 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef struct _p_PetscDraw* PetscDraw;

extern PetscFList PetscDrawList;
extern PetscErrorCode PetscDrawRegisterAll(const char[]);
extern PetscErrorCode PetscDrawInitializePackage(const char[]);
extern PetscErrorCode PetscDrawRegisterDestroy(void);

extern PetscErrorCode PetscDrawRegister(const char*,const char*,const char*,PetscErrorCode(*)(PetscDraw));
# 85 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
extern PetscErrorCode PetscDrawGetType(PetscDraw,const char**);
extern PetscErrorCode PetscDrawSetType(PetscDraw,const char*);
extern PetscErrorCode PetscDrawCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDraw*);
extern PetscErrorCode PetscDrawSetFromOptions(PetscDraw);
extern PetscErrorCode PetscDrawSetSave(PetscDraw,const char*);
# 133 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
extern PetscErrorCode PetscDrawOpenX(MPI_Comm,const char[],const char[],int,int,int,int,PetscDraw*);
extern PetscErrorCode PetscDrawOpenPS(MPI_Comm,char *,PetscDraw *);





extern PetscErrorCode PetscDrawOpenNull(MPI_Comm,PetscDraw *);
extern PetscErrorCode PetscDrawDestroy(PetscDraw*);
extern PetscErrorCode PetscDrawIsNull(PetscDraw,PetscBool *);

extern PetscErrorCode PetscDrawGetPopup(PetscDraw,PetscDraw*);
extern PetscErrorCode PetscDrawCheckResizedWindow(PetscDraw);
extern PetscErrorCode PetscDrawResizeWindow(PetscDraw,int,int);

extern PetscErrorCode PetscDrawScalePopup(PetscDraw,PetscReal,PetscReal);

extern PetscErrorCode PetscDrawLine(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
extern PetscErrorCode PetscDrawArrow(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
extern PetscErrorCode PetscDrawLineSetWidth(PetscDraw,PetscReal);
extern PetscErrorCode PetscDrawLineGetWidth(PetscDraw,PetscReal*);

extern PetscErrorCode PetscDrawPoint(PetscDraw,PetscReal,PetscReal,int);
extern PetscErrorCode PetscDrawPointSetSize(PetscDraw,PetscReal);

extern PetscErrorCode PetscDrawRectangle(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int,int,int,int);
extern PetscErrorCode PetscDrawTriangle(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,int,int,int);
extern PetscErrorCode PetscDrawEllipse(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal,int);
extern PetscErrorCode PetscDrawTensorContourPatch(PetscDraw,int,int,PetscReal*,PetscReal*,PetscReal,PetscReal,PetscReal*);
extern PetscErrorCode PetscDrawTensorContour(PetscDraw,int,int,const PetscReal[],const PetscReal[],PetscReal *);

extern PetscErrorCode PetscDrawString(PetscDraw,PetscReal,PetscReal,int,const char[]);
extern PetscErrorCode PetscDrawStringVertical(PetscDraw,PetscReal,PetscReal,int,const char[]);
extern PetscErrorCode PetscDrawStringSetSize(PetscDraw,PetscReal,PetscReal);
extern PetscErrorCode PetscDrawStringGetSize(PetscDraw,PetscReal*,PetscReal*);

extern PetscErrorCode PetscDrawSetViewPort(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode PetscDrawSplitViewPort(PetscDraw);

extern PetscErrorCode PetscDrawSetCoordinates(PetscDraw,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode PetscDrawGetCoordinates(PetscDraw,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

extern PetscErrorCode PetscDrawSetTitle(PetscDraw,const char[]);
extern PetscErrorCode PetscDrawAppendTitle(PetscDraw,const char[]);
extern PetscErrorCode PetscDrawGetTitle(PetscDraw,char **);

extern PetscErrorCode PetscDrawSetPause(PetscDraw,PetscReal);
extern PetscErrorCode PetscDrawGetPause(PetscDraw,PetscReal*);
extern PetscErrorCode PetscDrawPause(PetscDraw);
extern PetscErrorCode PetscDrawSetDoubleBuffer(PetscDraw);
extern PetscErrorCode PetscDrawFlush(PetscDraw);
extern PetscErrorCode PetscDrawSynchronizedFlush(PetscDraw);
extern PetscErrorCode PetscDrawClear(PetscDraw);
extern PetscErrorCode PetscDrawSave(PetscDraw);
extern PetscErrorCode PetscDrawSynchronizedClear(PetscDraw);
extern PetscErrorCode PetscDrawBOP(PetscDraw);
extern PetscErrorCode PetscDrawEOP(PetscDraw);

extern PetscErrorCode PetscDrawSetDisplay(PetscDraw,char*);


extern PetscErrorCode PetscDrawGetSingleton(PetscDraw,PetscDraw*);
extern PetscErrorCode PetscDrawRestoreSingleton(PetscDraw,PetscDraw*);
# 204 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef enum {BUTTON_NONE,BUTTON_LEFT,BUTTON_CENTER,BUTTON_RIGHT,BUTTON_LEFT_SHIFT,BUTTON_CENTER_SHIFT,BUTTON_RIGHT_SHIFT} PetscDrawButton;

extern PetscErrorCode PetscDrawGetMouseButton(PetscDraw,PetscDrawButton *,PetscReal*,PetscReal *,PetscReal *,PetscReal *);
extern PetscErrorCode PetscDrawSynchronizedGetMouseButton(PetscDraw,PetscDrawButton *,PetscReal*,PetscReal *,PetscReal *,PetscReal *);

extern PetscErrorCode PetscDrawZoom(PetscDraw,PetscErrorCode (*)(PetscDraw,void *),void *);
# 220 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef struct {
  PetscInt nports;
  PetscReal *xl;
  PetscReal *xr;
  PetscReal *yl;
  PetscReal *yr;
  PetscDraw draw;
} PetscDrawViewPorts;
extern PetscErrorCode PetscDrawViewPortsCreate(PetscDraw,PetscInt,PetscDrawViewPorts**);
extern PetscErrorCode PetscDrawViewPortsCreateRect(PetscDraw,PetscInt,PetscInt,PetscDrawViewPorts**);
extern PetscErrorCode PetscDrawViewPortsDestroy(PetscDrawViewPorts*);
extern PetscErrorCode PetscDrawViewPortsSet(PetscDrawViewPorts*,PetscInt);
# 242 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef struct _p_DrawAxis* PetscDrawAxis;

extern PetscClassId DRAWAXIS_CLASSID;

extern PetscErrorCode PetscDrawAxisCreate(PetscDraw,PetscDrawAxis *);
extern PetscErrorCode PetscDrawAxisDestroy(PetscDrawAxis*);
extern PetscErrorCode PetscDrawAxisDraw(PetscDrawAxis);
extern PetscErrorCode PetscDrawAxisSetLimits(PetscDrawAxis,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode PetscDrawAxisSetHoldLimits(PetscDrawAxis,PetscBool );
extern PetscErrorCode PetscDrawAxisSetColors(PetscDrawAxis,int,int,int);
extern PetscErrorCode PetscDrawAxisSetLabels(PetscDrawAxis,const char[],const char[],const char[]);
# 263 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef struct _p_PetscDrawLG* PetscDrawLG;

extern PetscClassId DRAWLG_CLASSID;

extern PetscErrorCode PetscDrawLGCreate(PetscDraw,int,PetscDrawLG *);
extern PetscErrorCode PetscDrawLGDestroy(PetscDrawLG*);
extern PetscErrorCode PetscDrawLGAddPoint(PetscDrawLG,PetscReal*,PetscReal*);
extern PetscErrorCode PetscDrawLGAddPoints(PetscDrawLG,int,PetscReal**,PetscReal**);
extern PetscErrorCode PetscDrawLGDraw(PetscDrawLG);
extern PetscErrorCode PetscDrawLGPrint(PetscDrawLG);
extern PetscErrorCode PetscDrawLGReset(PetscDrawLG);
extern PetscErrorCode PetscDrawLGSetDimension(PetscDrawLG,PetscInt);
extern PetscErrorCode PetscDrawLGSetLegend(PetscDrawLG,const char *const*);
extern PetscErrorCode PetscDrawLGGetAxis(PetscDrawLG,PetscDrawAxis *);
extern PetscErrorCode PetscDrawLGGetDraw(PetscDrawLG,PetscDraw *);
extern PetscErrorCode PetscDrawLGIndicateDataPoints(PetscDrawLG);
extern PetscErrorCode PetscDrawLGSetLimits(PetscDrawLG,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode PetscDrawLGSetColors(PetscDrawLG,const int*);
# 291 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef struct _p_DrawSP* PetscDrawSP;

extern PetscClassId DRAWSP_CLASSID;

extern PetscErrorCode PetscDrawSPCreate(PetscDraw,int,PetscDrawSP *);
extern PetscErrorCode PetscDrawSPDestroy(PetscDrawSP*);
extern PetscErrorCode PetscDrawSPAddPoint(PetscDrawSP,PetscReal*,PetscReal*);
extern PetscErrorCode PetscDrawSPAddPoints(PetscDrawSP,int,PetscReal**,PetscReal**);
extern PetscErrorCode PetscDrawSPDraw(PetscDrawSP);
extern PetscErrorCode PetscDrawSPReset(PetscDrawSP);
extern PetscErrorCode PetscDrawSPSetDimension(PetscDrawSP,int);
extern PetscErrorCode PetscDrawSPGetAxis(PetscDrawSP,PetscDrawAxis *);
extern PetscErrorCode PetscDrawSPGetDraw(PetscDrawSP,PetscDraw *);
extern PetscErrorCode PetscDrawSPSetLimits(PetscDrawSP,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode PetscDrawLGSPDraw(PetscDrawLG,PetscDrawSP);
# 316 "/home/dpnkarthik/petsc-rnet/include/petscdraw.h"
typedef struct _p_DrawHG* PetscDrawHG;

extern PetscClassId DRAWHG_CLASSID;

extern PetscErrorCode PetscDrawHGCreate(PetscDraw,int,PetscDrawHG *);
extern PetscErrorCode PetscDrawHGDestroy(PetscDrawHG*);
extern PetscErrorCode PetscDrawHGAddValue(PetscDrawHG,PetscReal);
extern PetscErrorCode PetscDrawHGDraw(PetscDrawHG);
extern PetscErrorCode PetscDrawHGPrint(PetscDrawHG);
extern PetscErrorCode PetscDrawHGReset(PetscDrawHG);
extern PetscErrorCode PetscDrawHGGetAxis(PetscDrawHG,PetscDrawAxis *);
extern PetscErrorCode PetscDrawHGGetDraw(PetscDrawHG,PetscDraw *);
extern PetscErrorCode PetscDrawHGSetLimits(PetscDrawHG,PetscReal,PetscReal,int,int);
extern PetscErrorCode PetscDrawHGSetNumberBins(PetscDrawHG,int);
extern PetscErrorCode PetscDrawHGSetColor(PetscDrawHG,int);
extern PetscErrorCode PetscDrawHGCalcStats(PetscDrawHG, PetscBool );
extern PetscErrorCode PetscDrawHGIntegerBins(PetscDrawHG, PetscBool );




extern PetscErrorCode PetscViewerDrawGetDraw(PetscViewer,PetscInt,PetscDraw*);
extern PetscErrorCode PetscViewerDrawBaseAdd(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerDrawBaseSet(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerDrawGetDrawLG(PetscViewer,PetscInt,PetscDrawLG*);
extern PetscErrorCode PetscViewerDrawGetDrawAxis(PetscViewer,PetscInt,PetscDrawAxis*);

extern PetscErrorCode PetscDrawUtilitySetCmapHue(unsigned char *,unsigned char *,unsigned char *,int);
extern PetscErrorCode PetscDrawUtilitySetGamma(PetscReal);


# 1545 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2




# 1 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h" 1







# 1 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 1
# 9 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h" 2

# 35 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
typedef struct {
   PetscErrorCode (*getcomm)(PetscObject,MPI_Comm *);
   PetscErrorCode (*view)(PetscObject,PetscViewer);
   PetscErrorCode (*destroy)(PetscObject*);
   PetscErrorCode (*compose)(PetscObject,const char[],PetscObject);
   PetscErrorCode (*query)(PetscObject,const char[],PetscObject *);
   PetscErrorCode (*composefunction)(PetscObject,const char[],const char[],void (*)(void));
   PetscErrorCode (*queryfunction)(PetscObject,const char[],void (**)(void));
   PetscErrorCode (*publish)(PetscObject);
} PetscOps;
# 53 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
typedef struct _p_PetscObject {
  PetscClassId classid;
  PetscOps *bops;
  MPI_Comm comm;
  PetscInt type;
  PetscLogDouble flops,time,mem;
  PetscInt id;
  PetscInt refct;
  PetscMPIInt tag;
  PetscFList qlist;
  PetscOList olist;
  char *class_name;
  char *description;
  char *mansec;
  char *type_name;
  PetscObject parent;
  PetscInt parentid;
  char* name;
  char *prefix;
  PetscInt tablevel;
  void *cpp;
  PetscInt amem;
  PetscInt state;
  PetscInt int_idmax, intstar_idmax;
  PetscInt *intcomposedstate,*intstarcomposedstate;
  PetscInt *intcomposeddata, **intstarcomposeddata;
  PetscInt real_idmax, realstar_idmax;
  PetscInt *realcomposedstate,*realstarcomposedstate;
  PetscReal *realcomposeddata, **realstarcomposeddata;
  PetscInt scalar_idmax, scalarstar_idmax;
  PetscInt *scalarcomposedstate,*scalarstarcomposedstate;
  PetscScalar *scalarcomposeddata, **scalarstarcomposeddata;
  void (**fortran_func_pointers)(void);
  void *python_context;
  PetscErrorCode (*python_destroy)(void*);

  PetscInt noptionhandler;
  PetscErrorCode (*optionhandler[5])(PetscObject,void*);
  PetscErrorCode (*optiondestroy[5])(PetscObject,void*);
  void *optionctx[5];
  PetscPrecision precision;
  PetscBool optionsprinted;
} _p_PetscObject;







typedef PetscErrorCode (*PetscObjectFunction)(PetscObject*);
typedef PetscErrorCode (*PetscObjectViewerFunction)(PetscObject,PetscViewer);
# 135 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
extern PetscErrorCode PetscComposedQuantitiesDestroy(PetscObject obj);
extern PetscErrorCode PetscHeaderCreate_Private(PetscObject,PetscClassId,PetscInt,const char[],const char[],const char[],MPI_Comm,PetscErrorCode (*)(PetscObject*),PetscErrorCode (*)(PetscObject,PetscViewer));
# 155 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
extern PetscErrorCode PetscHeaderDestroy_Private(PetscObject);
# 430 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
extern PetscErrorCode PetscObjectStateQuery(PetscObject,PetscInt*);
extern PetscErrorCode PetscObjectSetState(PetscObject,PetscInt);
extern PetscErrorCode PetscObjectComposedDataRegister(PetscInt*);
extern PetscErrorCode PetscObjectComposedDataIncreaseInt(PetscObject);
extern PetscErrorCode PetscObjectComposedDataIncreaseIntstar(PetscObject);
extern PetscErrorCode PetscObjectComposedDataIncreaseReal(PetscObject);
extern PetscErrorCode PetscObjectComposedDataIncreaseRealstar(PetscObject);
extern PetscErrorCode PetscObjectComposedDataIncreaseScalar(PetscObject);
extern PetscErrorCode PetscObjectComposedDataIncreaseScalarstar(PetscObject);
extern PetscInt globalcurrentstate;
extern PetscInt globalmaxstate;
# 745 "/home/dpnkarthik/petsc-rnet/include/private/petscimpl.h"
extern PetscBool PetscPreLoadingUsed;
extern PetscBool PetscPreLoadingOn;

extern PetscMPIInt Petsc_Counter_keyval;
extern PetscMPIInt Petsc_InnerComm_keyval;
extern PetscMPIInt Petsc_OuterComm_keyval;





typedef struct {
  PetscMPIInt tag;
  PetscInt refcount;
  PetscInt namecount;
} PetscCommCounter;


# 1550 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 1592 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode (*PetscErrorPrintf)(const char[],...);
# 1615 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode (*PetscHelpPrintf)(MPI_Comm,const char[],...);




# 1 "/home/dpnkarthik/petsc-rnet/include/petsclog.h" 1








# 18 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
typedef int PetscLogEvent;
# 27 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
typedef int PetscLogStage;


extern PetscLogEvent PETSC_LARGEST_EVENT;


extern PetscLogDouble _TotalFlops;
extern PetscLogDouble petsc_tmp_flops;


extern PetscErrorCode PetscInfo_Private(const char[],void*,const char[],...);
# 57 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
extern PetscErrorCode PetscInfoDeactivateClass(PetscClassId);
extern PetscErrorCode PetscInfoActivateClass(PetscClassId);
extern PetscBool PetscLogPrintInfo;







typedef struct _n_PetscIntStack *PetscIntStack;
# 76 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
typedef struct {
  char *name;
  PetscClassId classid;
} PetscClassRegInfo;

typedef struct {
  PetscClassId id;
  int creations;
  int destructions;
  PetscLogDouble mem;
  PetscLogDouble descMem;
} PetscClassPerfInfo;

typedef struct _n_PetscClassRegLog *PetscClassRegLog;
struct _n_PetscClassRegLog {
  int numClasses;
  int maxClasses;
  PetscClassRegInfo *classInfo;
};

typedef struct _n_PetscClassPerfLog *PetscClassPerfLog;
struct _n_PetscClassPerfLog {
  int numClasses;
  int maxClasses;
  PetscClassPerfInfo *classInfo;
};
# 112 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
typedef struct {
  char *name;
  PetscClassId classid;




} PetscEventRegInfo;

typedef struct {
  int id;
  PetscBool active;
  PetscBool visible;
  int depth;
  int count;
  PetscLogDouble flops;
  PetscLogDouble time;
  PetscLogDouble numMessages;
  PetscLogDouble messageLength;
  PetscLogDouble numReductions;
} PetscEventPerfInfo;

typedef struct _n_PetscEventRegLog *PetscEventRegLog;
struct _n_PetscEventRegLog {
  int numEvents;
  int maxEvents;
  PetscEventRegInfo *eventInfo;
};

typedef struct _n_PetscEventPerfLog *PetscEventPerfLog;
struct _n_PetscEventPerfLog {
  int numEvents;
  int maxEvents;
  PetscEventPerfInfo *eventInfo;
};






typedef struct _PetscStageInfo {
  char *name;
  PetscBool used;
  PetscEventPerfInfo perfInfo;
  PetscEventPerfLog eventLog;
  PetscClassPerfLog classLog;
} PetscStageInfo;

typedef struct _n_PetscStageLog *PetscStageLog;
extern PetscStageLog _stageLog;
struct _n_PetscStageLog {
  int numStages;
  int maxStages;
  PetscIntStack stack;
  int curStage;
  PetscStageInfo *stageInfo;
  PetscEventRegLog eventLog;
  PetscClassRegLog classLog;
};
# 220 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
extern PetscErrorCode (*_PetscLogPLB)(PetscLogEvent,int,PetscObject,PetscObject,PetscObject,PetscObject);
extern PetscErrorCode (*_PetscLogPLE)(PetscLogEvent,int,PetscObject,PetscObject,PetscObject,PetscObject);
extern PetscErrorCode (*_PetscLogPHC)(PetscObject);
extern PetscErrorCode (*_PetscLogPHD)(PetscObject);
# 233 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
extern PetscErrorCode PetscLogBegin(void);
extern PetscErrorCode PetscLogAllBegin(void);
extern PetscErrorCode PetscLogTraceBegin(FILE *);
extern PetscErrorCode PetscLogActions(PetscBool);
extern PetscErrorCode PetscLogObjects(PetscBool);

extern PetscErrorCode PetscLogGetRGBColor(const char*[]);
extern PetscErrorCode PetscLogDestroy(void);
extern PetscErrorCode PetscLogSet(PetscErrorCode (*)(int, int, PetscObject, PetscObject, PetscObject, PetscObject),
                   PetscErrorCode (*)(int, int, PetscObject, PetscObject, PetscObject, PetscObject));
extern PetscErrorCode PetscLogObjectState(PetscObject, const char[], ...);

extern PetscErrorCode PetscLogView(PetscViewer);
extern PetscErrorCode PetscLogViewPython(PetscViewer);
extern PetscErrorCode PetscLogPrintDetailed(MPI_Comm, const char[]);
extern PetscErrorCode PetscLogDump(const char[]);

extern PetscErrorCode PetscGetFlops(PetscLogDouble *);

extern PetscErrorCode PetscLogStageRegister(const char[],PetscLogStage*);
extern PetscErrorCode PetscLogStagePush(PetscLogStage);
extern PetscErrorCode PetscLogStagePop(void);
extern PetscErrorCode PetscLogStageSetActive(PetscLogStage, PetscBool );
extern PetscErrorCode PetscLogStageGetActive(PetscLogStage, PetscBool *);
extern PetscErrorCode PetscLogStageSetVisible(PetscLogStage, PetscBool );
extern PetscErrorCode PetscLogStageGetVisible(PetscLogStage, PetscBool *);
extern PetscErrorCode PetscLogStageGetId(const char [], PetscLogStage *);

extern PetscErrorCode PetscLogEventRegister(const char[], PetscClassId,PetscLogEvent*);
extern PetscErrorCode PetscLogEventActivate(PetscLogEvent);
extern PetscErrorCode PetscLogEventDeactivate(PetscLogEvent);
extern PetscErrorCode PetscLogEventSetActiveAll(PetscLogEvent, PetscBool );
extern PetscErrorCode PetscLogEventActivateClass(PetscClassId);
extern PetscErrorCode PetscLogEventDeactivateClass(PetscClassId);



extern PetscLogDouble petsc_irecv_ct;
extern PetscLogDouble petsc_isend_ct;
extern PetscLogDouble petsc_recv_ct;
extern PetscLogDouble petsc_send_ct;
extern PetscLogDouble petsc_irecv_len;
extern PetscLogDouble petsc_isend_len;
extern PetscLogDouble petsc_recv_len;
extern PetscLogDouble petsc_send_len;
extern PetscLogDouble petsc_allreduce_ct;
extern PetscLogDouble petsc_gather_ct;
extern PetscLogDouble petsc_scatter_ct;
extern PetscLogDouble petsc_wait_ct;
extern PetscLogDouble petsc_wait_any_ct;
extern PetscLogDouble petsc_wait_all_ct;
extern PetscLogDouble petsc_sum_of_waits_ct;
# 303 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
extern PetscErrorCode PetscLogEventGetFlops(PetscLogEvent, PetscLogDouble*);
extern PetscErrorCode PetscLogEventZeroFlops(PetscLogEvent);
# 323 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
static inline PetscErrorCode PetscMPITypeSize(PetscLogDouble *buff,PetscMPIInt count,MPI_Datatype type)
{
  PetscMPIInt mysize; return (MPI_Type_size(type,&mysize) || ((*buff += (PetscLogDouble) (count*mysize)),0));
}

static inline PetscErrorCode PetscMPITypeSizeComm(MPI_Comm comm, PetscLogDouble *buff,PetscMPIInt *counts,MPI_Datatype type)
{
  PetscMPIInt mysize, commsize, p;
  PetscErrorCode _myierr;

  _myierr = MPI_Comm_size(comm,&commsize);do {if (__builtin_expect(!!(_myierr),0)) return PetscError(((MPI_Comm)0x44000001),333,__func__,"/home/dpnkarthik/petsc-rnet/include/petsclog.h","src/mat/impls/baij/seq/",_myierr,PETSC_ERROR_REPEAT," ");} while (0);
  _myierr = MPI_Type_size(type,&mysize);do {if (__builtin_expect(!!(_myierr),0)) return PetscError(((MPI_Comm)0x44000001),334,__func__,"/home/dpnkarthik/petsc-rnet/include/petsclog.h","src/mat/impls/baij/seq/",_myierr,PETSC_ERROR_REPEAT," ");} while (0);
  for(p = 0; p < commsize; ++p) {
    *buff += (PetscLogDouble) (counts[p]*mysize);
  }
  return 0;
}
# 473 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
extern PetscErrorCode PetscIntStackCreate(PetscIntStack *);
extern PetscErrorCode PetscIntStackDestroy(PetscIntStack);
extern PetscErrorCode PetscIntStackPush(PetscIntStack, int);
extern PetscErrorCode PetscIntStackPop(PetscIntStack, int *);
extern PetscErrorCode PetscIntStackTop(PetscIntStack, int *);
extern PetscErrorCode PetscIntStackEmpty(PetscIntStack, PetscBool *);
# 497 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
static inline PetscErrorCode PetscLogGetStageLog(PetscStageLog *stageLog)
{
  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petsclog.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 499; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscLogGetStageLog") && strcmp("PetscLogGetStageLog","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petsclog.h",499,"PetscLogGetStageLog","__func__",__func__); } } while (0); } while (0);
  do {} while (0);
  if (!_stageLog) {
    fprintf(stderr, "PETSC ERROR: Logging has not been enabled.\nYou might have forgotten to call PetscInitialize().\n");
    MPI_Abort(((MPI_Comm)0x44000000), 56);
  }
  *stageLog = _stageLog;
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 532 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
static inline PetscErrorCode PetscStageLogGetCurrent(PetscStageLog stageLog, int *stage)
{
  PetscBool empty;
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petsclog.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 537; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscStageLogGetCurrent") && strcmp("PetscStageLogGetCurrent","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petsclog.h",537,"PetscStageLogGetCurrent","__func__",__func__); } } while (0); } while (0);
  ierr = PetscIntStackEmpty(stageLog->stack, &empty);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),538,__func__,"/home/dpnkarthik/petsc-rnet/include/petsclog.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (empty) {
    *stage = -1;
  } else {
    ierr = PetscIntStackTop(stageLog->stack, stage);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),542,__func__,"/home/dpnkarthik/petsc-rnet/include/petsclog.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }





  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 573 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"
static inline PetscErrorCode PetscStageLogGetEventPerfLog(PetscStageLog stageLog, int stage, PetscEventPerfLog *eventLog)
{
  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petsclog.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 575; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscStageLogGetEventPerfLog") && strcmp("PetscStageLogGetEventPerfLog","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petsclog.h",575,"PetscStageLogGetEventPerfLog","__func__",__func__); } } while (0); } while (0);
  do {} while (0);
  if ((stage < 0) || (stage >= stageLog->numStages)) {
    return PetscError(((MPI_Comm)0x44000001),578,__func__,"/home/dpnkarthik/petsc-rnet/include/petsclog.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"Invalid stage %d should be in [0,%d)",stage,stageLog->numStages);
  }
  *eventLog = stageLog->stageInfo[stage].eventLog;
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}


# 1 "/home/dpnkarthik/petsc-rnet/include/petsclog.hh" 1
# 586 "/home/dpnkarthik/petsc-rnet/include/petsclog.h" 2
# 624 "/home/dpnkarthik/petsc-rnet/include/petsclog.h"

# 1621 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 1639 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscFixFilename(const char[],char[]);
extern PetscErrorCode PetscFOpen(MPI_Comm,const char[],const char[],FILE**);
extern PetscErrorCode PetscFClose(MPI_Comm,FILE*);
extern PetscErrorCode PetscFPrintf(MPI_Comm,FILE*,const char[],...);
extern PetscErrorCode PetscPrintf(MPI_Comm,const char[],...);
extern PetscErrorCode PetscSNPrintf(char*,size_t,const char [],...);
extern PetscErrorCode PetscSNPrintfCount(char*,size_t,const char [],size_t*,...);




# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stdarg.h" 1 3 4
# 1651 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
extern PetscErrorCode PetscVSNPrintf(char*,size_t,const char[],size_t*,va_list);
extern PetscErrorCode (*PetscVFPrintf)(FILE*,const char[],va_list);
extern PetscErrorCode PetscVFPrintfDefault(FILE*,const char[],va_list);
extern PetscErrorCode PetscVFPrintfRegress(FILE*,const char *,va_list);
extern PetscErrorCode PetscVFPrintfRegressSetUp(MPI_Comm,const char *);





extern PetscErrorCode PetscErrorPrintfDefault(const char [],...);
extern PetscErrorCode PetscErrorPrintfNone(const char [],...);
extern PetscErrorCode PetscHelpPrintfDefault(MPI_Comm,const char [],...);


extern PetscErrorCode PetscPOpen(MPI_Comm,const char[],const char[],const char[],FILE **);
extern PetscErrorCode PetscPClose(MPI_Comm,FILE*);


extern PetscErrorCode PetscSynchronizedPrintf(MPI_Comm,const char[],...);
extern PetscErrorCode PetscSynchronizedFPrintf(MPI_Comm,FILE*,const char[],...);
extern PetscErrorCode PetscSynchronizedFlush(MPI_Comm);
extern PetscErrorCode PetscSynchronizedFGets(MPI_Comm,FILE*,size_t,char[]);
extern PetscErrorCode PetscStartMatlab(MPI_Comm,const char[],const char[],FILE**);
extern PetscErrorCode PetscStartJava(MPI_Comm,const char[],const char[],FILE**);
extern PetscErrorCode PetscGetPetscDir(const char*[]);

extern PetscErrorCode PetscPopUpSelect(MPI_Comm,const char*,const char*,int,const char**,int*);
# 1687 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscClassId PETSC_CONTAINER_CLASSID;
typedef struct _p_PetscContainer* PetscContainer;
extern PetscErrorCode PetscContainerGetPointer(PetscContainer,void **);
extern PetscErrorCode PetscContainerSetPointer(PetscContainer,void *);
extern PetscErrorCode PetscContainerDestroy(PetscContainer*);
extern PetscErrorCode PetscContainerCreate(MPI_Comm,PetscContainer *);
extern PetscErrorCode PetscContainerSetUserDestroy(PetscContainer, PetscErrorCode (*)(void*));




extern PetscMPIInt PetscGlobalRank;
extern PetscMPIInt PetscGlobalSize;
extern PetscErrorCode PetscIntView(PetscInt,const PetscInt[],PetscViewer);
extern PetscErrorCode PetscRealView(PetscInt,const PetscReal[],PetscViewer);
extern PetscErrorCode PetscScalarView(PetscInt,const PetscScalar[],PetscViewer);


# 1 "/usr/include/memory.h" 1 3 4
# 1706 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2


# 1 "/usr/include/stdlib.h" 1 3 4
# 33 "/usr/include/stdlib.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 34 "/usr/include/stdlib.h" 2 3 4








# 1 "/usr/include/bits/waitflags.h" 1 3 4
# 43 "/usr/include/stdlib.h" 2 3 4
# 1 "/usr/include/bits/waitstatus.h" 1 3 4
# 67 "/usr/include/bits/waitstatus.h" 3 4
union wait
  {
    int w_status;
    struct
      {

 unsigned int __w_termsig:7;
 unsigned int __w_coredump:1;
 unsigned int __w_retcode:8;
 unsigned int:16;







      } __wait_terminated;
    struct
      {

 unsigned int __w_stopval:8;
 unsigned int __w_stopsig:8;
 unsigned int:16;






      } __wait_stopped;
  };
# 44 "/usr/include/stdlib.h" 2 3 4
# 68 "/usr/include/stdlib.h" 3 4
typedef union
  {
    union wait *__uptr;
    int *__iptr;
  } __WAIT_STATUS __attribute__ ((__transparent_union__));
# 96 "/usr/include/stdlib.h" 3 4


typedef struct
  {
    int quot;
    int rem;
  } div_t;



typedef struct
  {
    long int quot;
    long int rem;
  } ldiv_t;







__extension__ typedef struct
  {
    long long int quot;
    long long int rem;
  } lldiv_t;


# 140 "/usr/include/stdlib.h" 3 4
extern size_t __ctype_get_mb_cur_max (void) __attribute__ ((__nothrow__)) ;




extern double atof (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

extern int atoi (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

extern long int atol (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;





__extension__ extern long long int atoll (__const char *__nptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;





extern double strtod (__const char *__restrict __nptr,
        char **__restrict __endptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern float strtof (__const char *__restrict __nptr,
       char **__restrict __endptr) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

extern long double strtold (__const char *__restrict __nptr,
       char **__restrict __endptr)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern long int strtol (__const char *__restrict __nptr,
   char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

extern unsigned long int strtoul (__const char *__restrict __nptr,
      char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;




__extension__
extern long long int strtoq (__const char *__restrict __nptr,
        char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

__extension__
extern unsigned long long int strtouq (__const char *__restrict __nptr,
           char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





__extension__
extern long long int strtoll (__const char *__restrict __nptr,
         char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

__extension__
extern unsigned long long int strtoull (__const char *__restrict __nptr,
     char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;

# 311 "/usr/include/stdlib.h" 3 4
extern char *l64a (long int __n) __attribute__ ((__nothrow__)) ;


extern long int a64l (__const char *__s)
     __attribute__ ((__nothrow__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;
# 327 "/usr/include/stdlib.h" 3 4
extern long int random (void) __attribute__ ((__nothrow__));


extern void srandom (unsigned int __seed) __attribute__ ((__nothrow__));





extern char *initstate (unsigned int __seed, char *__statebuf,
   size_t __statelen) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));



extern char *setstate (char *__statebuf) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));







struct random_data
  {
    int32_t *fptr;
    int32_t *rptr;
    int32_t *state;
    int rand_type;
    int rand_deg;
    int rand_sep;
    int32_t *end_ptr;
  };

extern int random_r (struct random_data *__restrict __buf,
       int32_t *__restrict __result) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern int srandom_r (unsigned int __seed, struct random_data *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));

extern int initstate_r (unsigned int __seed, char *__restrict __statebuf,
   size_t __statelen,
   struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2, 4)));

extern int setstate_r (char *__restrict __statebuf,
         struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));






extern int rand (void) __attribute__ ((__nothrow__));

extern void srand (unsigned int __seed) __attribute__ ((__nothrow__));




extern int rand_r (unsigned int *__seed) __attribute__ ((__nothrow__));







extern double drand48 (void) __attribute__ ((__nothrow__));
extern double erand48 (unsigned short int __xsubi[3]) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern long int lrand48 (void) __attribute__ ((__nothrow__));
extern long int nrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern long int mrand48 (void) __attribute__ ((__nothrow__));
extern long int jrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));


extern void srand48 (long int __seedval) __attribute__ ((__nothrow__));
extern unsigned short int *seed48 (unsigned short int __seed16v[3])
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
extern void lcong48 (unsigned short int __param[7]) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





struct drand48_data
  {
    unsigned short int __x[3];
    unsigned short int __old_x[3];
    unsigned short int __c;
    unsigned short int __init;
    unsigned long long int __a;
  };


extern int drand48_r (struct drand48_data *__restrict __buffer,
        double *__restrict __result) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern int erand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        double *__restrict __result) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int lrand48_r (struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern int nrand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int mrand48_r (struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));
extern int jrand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));


extern int srand48_r (long int __seedval, struct drand48_data *__buffer)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));

extern int seed48_r (unsigned short int __seed16v[3],
       struct drand48_data *__buffer) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));

extern int lcong48_r (unsigned short int __param[7],
        struct drand48_data *__buffer)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2)));









extern void *malloc (size_t __size) __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;

extern void *calloc (size_t __nmemb, size_t __size)
     __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;










extern void *realloc (void *__ptr, size_t __size)
     __attribute__ ((__nothrow__)) __attribute__ ((__warn_unused_result__));

extern void free (void *__ptr) __attribute__ ((__nothrow__));




extern void cfree (void *__ptr) __attribute__ ((__nothrow__));



# 1 "/usr/include/alloca.h" 1 3 4
# 25 "/usr/include/alloca.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 26 "/usr/include/alloca.h" 2 3 4







extern void *alloca (size_t __size) __attribute__ ((__nothrow__));






# 498 "/usr/include/stdlib.h" 2 3 4





extern void *valloc (size_t __size) __attribute__ ((__nothrow__)) __attribute__ ((__malloc__)) ;




extern int posix_memalign (void **__memptr, size_t __alignment, size_t __size)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;




extern void abort (void) __attribute__ ((__nothrow__)) __attribute__ ((__noreturn__));



extern int atexit (void (*__func) (void)) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 531 "/usr/include/stdlib.h" 3 4





extern int on_exit (void (*__func) (int __status, void *__arg), void *__arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));






extern void exit (int __status) __attribute__ ((__nothrow__)) __attribute__ ((__noreturn__));
# 554 "/usr/include/stdlib.h" 3 4






extern void _Exit (int __status) __attribute__ ((__nothrow__)) __attribute__ ((__noreturn__));






extern char *getenv (__const char *__name) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;




extern char *__secure_getenv (__const char *__name)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;





extern int putenv (char *__string) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));





extern int setenv (__const char *__name, __const char *__value, int __replace)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (2)));


extern int unsetenv (__const char *__name) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));






extern int clearenv (void) __attribute__ ((__nothrow__));
# 606 "/usr/include/stdlib.h" 3 4
extern char *mktemp (char *__template) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
# 620 "/usr/include/stdlib.h" 3 4
extern int mkstemp (char *__template) __attribute__ ((__nonnull__ (1))) ;
# 642 "/usr/include/stdlib.h" 3 4
extern int mkstemps (char *__template, int __suffixlen) __attribute__ ((__nonnull__ (1))) ;
# 663 "/usr/include/stdlib.h" 3 4
extern char *mkdtemp (char *__template) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
# 712 "/usr/include/stdlib.h" 3 4





extern int system (__const char *__command) ;

# 734 "/usr/include/stdlib.h" 3 4
extern char *realpath (__const char *__restrict __name,
         char *__restrict __resolved) __attribute__ ((__nothrow__)) ;






typedef int (*__compar_fn_t) (__const void *, __const void *);
# 752 "/usr/include/stdlib.h" 3 4



extern void *bsearch (__const void *__key, __const void *__base,
        size_t __nmemb, size_t __size, __compar_fn_t __compar)
     __attribute__ ((__nonnull__ (1, 2, 5))) ;



extern void qsort (void *__base, size_t __nmemb, size_t __size,
     __compar_fn_t __compar) __attribute__ ((__nonnull__ (1, 4)));
# 771 "/usr/include/stdlib.h" 3 4
extern int abs (int __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;
extern long int labs (long int __x) __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;



__extension__ extern long long int llabs (long long int __x)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;







extern div_t div (int __numer, int __denom)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;
extern ldiv_t ldiv (long int __numer, long int __denom)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;




__extension__ extern lldiv_t lldiv (long long int __numer,
        long long int __denom)
     __attribute__ ((__nothrow__)) __attribute__ ((__const__)) ;

# 808 "/usr/include/stdlib.h" 3 4
extern char *ecvt (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;




extern char *fcvt (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;




extern char *gcvt (double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3))) ;




extern char *qecvt (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qfcvt (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qgcvt (long double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3))) ;




extern int ecvt_r (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign, char *__restrict __buf,
     size_t __len) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int fcvt_r (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign, char *__restrict __buf,
     size_t __len) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));

extern int qecvt_r (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int qfcvt_r (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (3, 4, 5)));







extern int mblen (__const char *__s, size_t __n) __attribute__ ((__nothrow__)) ;


extern int mbtowc (wchar_t *__restrict __pwc,
     __const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__)) ;


extern int wctomb (char *__s, wchar_t __wchar) __attribute__ ((__nothrow__)) ;



extern size_t mbstowcs (wchar_t *__restrict __pwcs,
   __const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__));

extern size_t wcstombs (char *__restrict __s,
   __const wchar_t *__restrict __pwcs, size_t __n)
     __attribute__ ((__nothrow__));








extern int rpmatch (__const char *__response) __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1))) ;
# 896 "/usr/include/stdlib.h" 3 4
extern int getsubopt (char **__restrict __optionp,
        char *__const *__restrict __tokens,
        char **__restrict __valuep)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1, 2, 3))) ;
# 948 "/usr/include/stdlib.h" 3 4
extern int getloadavg (double __loadavg[], int __nelem)
     __attribute__ ((__nothrow__)) __attribute__ ((__nonnull__ (1)));
# 964 "/usr/include/stdlib.h" 3 4

# 1709 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2


# 1 "/usr/include/strings.h" 1 3 4
# 1712 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 1762 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
static inline PetscErrorCode PetscMemcpy(void *a,const void *b,size_t n)
{







  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petscsys.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1771; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscMemcpy") && strcmp("PetscMemcpy","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petscsys.h",1771,"PetscMemcpy","__func__",__func__); } } while (0); } while (0);

  if (a != b) {
# 1798 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
    memcpy((char*)(a),(char*)(b),n);

  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 1826 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
static inline PetscErrorCode PetscMemzero(void *a,size_t n)
{
  if (n > 0) {
# 1847 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
      memset((char*)a,0,n);




  }
  return 0;
}
# 1895 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscmatlab.h" 1








extern PetscClassId MATLABENGINE_CLASSID;
# 21 "/home/dpnkarthik/petsc-rnet/include/petscmatlab.h"
typedef struct _p_PetscMatlabEngine* PetscMatlabEngine;

extern PetscErrorCode PetscMatlabEngineCreate(MPI_Comm,const char[],PetscMatlabEngine*);
extern PetscErrorCode PetscMatlabEngineDestroy(PetscMatlabEngine*);
extern PetscErrorCode PetscMatlabEngineEvaluate(PetscMatlabEngine,const char[],...);
extern PetscErrorCode PetscMatlabEngineGetOutput(PetscMatlabEngine,char **);
extern PetscErrorCode PetscMatlabEnginePrintOutput(PetscMatlabEngine,FILE*);
extern PetscErrorCode PetscMatlabEnginePut(PetscMatlabEngine,PetscObject);
extern PetscErrorCode PetscMatlabEngineGet(PetscMatlabEngine,PetscObject);
extern PetscErrorCode PetscMatlabEnginePutArray(PetscMatlabEngine,int,int,const PetscScalar*,const char[]);
extern PetscErrorCode PetscMatlabEngineGetArray(PetscMatlabEngine,int,int,PetscScalar*,const char[]);

extern PetscMatlabEngine PETSC_MATLAB_ENGINE_(MPI_Comm);
# 49 "/home/dpnkarthik/petsc-rnet/include/petscmatlab.h"

# 1896 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 2045 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode MPIU_File_write_all(MPI_File,void*,PetscMPIInt,MPI_Datatype,MPI_Status*);
extern PetscErrorCode MPIU_File_read_all(MPI_File,void*,PetscMPIInt,MPI_Datatype,MPI_Status*);
# 2101 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 1 3 4
# 2102 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2


# 1 "/usr/include/sys/param.h" 1 3 4
# 26 "/usr/include/sys/param.h" 3 4
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 1 3 4
# 27 "/usr/include/sys/param.h" 2 3 4

# 1 "/usr/include/linux/param.h" 1 3 4



# 1 "/usr/include/asm/param.h" 1 3 4
# 1 "/usr/include/asm-generic/param.h" 1 3 4
# 1 "/usr/include/asm/param.h" 2 3 4
# 5 "/usr/include/linux/param.h" 2 3 4
# 29 "/usr/include/sys/param.h" 2 3 4
# 2105 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 2120 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscsys.hh" 1
# 2121 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2
# 2202 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscGetArchType(char[],size_t);
extern PetscErrorCode PetscGetHostName(char[],size_t);
extern PetscErrorCode PetscGetUserName(char[],size_t);
extern PetscErrorCode PetscGetProgramName(char[],size_t);
extern PetscErrorCode PetscSetProgramName(const char[]);
extern PetscErrorCode PetscGetDate(char[],size_t);

extern PetscErrorCode PetscSortInt(PetscInt,PetscInt[]);
extern PetscErrorCode PetscSortRemoveDupsInt(PetscInt*,PetscInt[]);
extern PetscErrorCode PetscSortIntWithPermutation(PetscInt,const PetscInt[],PetscInt[]);
extern PetscErrorCode PetscSortStrWithPermutation(PetscInt,const char*[],PetscInt[]);
extern PetscErrorCode PetscSortIntWithArray(PetscInt,PetscInt[],PetscInt[]);
extern PetscErrorCode PetscSortIntWithArrayPair(PetscInt,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode PetscSortMPIIntWithArray(PetscMPIInt,PetscMPIInt[],PetscMPIInt[]);
extern PetscErrorCode PetscSortIntWithScalarArray(PetscInt,PetscInt[],PetscScalar[]);
extern PetscErrorCode PetscSortReal(PetscInt,PetscReal[]);
extern PetscErrorCode PetscSortRealWithPermutation(PetscInt,const PetscReal[],PetscInt[]);
extern PetscErrorCode PetscSortSplit(PetscInt,PetscInt,PetscScalar[],PetscInt[]);
extern PetscErrorCode PetscSortSplitReal(PetscInt,PetscInt,PetscReal[],PetscInt[]);
extern PetscErrorCode PetscProcessTree(PetscInt,const PetscBool [],const PetscInt[],PetscInt*,PetscInt**,PetscInt**,PetscInt**,PetscInt**);
extern PetscErrorCode PetscMergeIntArrayPair(PetscInt,const PetscInt*,const PetscInt*,PetscInt,const PetscInt*,const PetscInt*,PetscInt*,PetscInt**,PetscInt**);

extern PetscErrorCode PetscSetDisplay(void);
extern PetscErrorCode PetscGetDisplay(char[],size_t);
# 2245 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscClassId PETSC_RANDOM_CLASSID;

extern PetscErrorCode PetscRandomInitializePackage(const char[]);
# 2258 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _p_PetscRandom* PetscRandom;


extern PetscFList PetscRandomList;
extern PetscBool PetscRandomRegisterAllCalled;

extern PetscErrorCode PetscRandomRegisterAll(const char []);
extern PetscErrorCode PetscRandomRegister(const char[],const char[],const char[],PetscErrorCode (*)(PetscRandom));
extern PetscErrorCode PetscRandomRegisterDestroy(void);
extern PetscErrorCode PetscRandomSetType(PetscRandom, const char*);
extern PetscErrorCode PetscRandomSetFromOptions(PetscRandom);
extern PetscErrorCode PetscRandomGetType(PetscRandom, const char**);
extern PetscErrorCode PetscRandomViewFromOptions(PetscRandom,char*);
extern PetscErrorCode PetscRandomView(PetscRandom,PetscViewer);
# 2323 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
extern PetscErrorCode PetscRandomCreate(MPI_Comm,PetscRandom*);
extern PetscErrorCode PetscRandomGetValue(PetscRandom,PetscScalar*);
extern PetscErrorCode PetscRandomGetValueReal(PetscRandom,PetscReal*);
extern PetscErrorCode PetscRandomGetInterval(PetscRandom,PetscScalar*,PetscScalar*);
extern PetscErrorCode PetscRandomSetInterval(PetscRandom,PetscScalar,PetscScalar);
extern PetscErrorCode PetscRandomSetSeed(PetscRandom,unsigned long);
extern PetscErrorCode PetscRandomGetSeed(PetscRandom,unsigned long *);
extern PetscErrorCode PetscRandomSeed(PetscRandom);
extern PetscErrorCode PetscRandomDestroy(PetscRandom*);

extern PetscErrorCode PetscGetFullPath(const char[],char[],size_t);
extern PetscErrorCode PetscGetRelativePath(const char[],char[],size_t);
extern PetscErrorCode PetscGetWorkingDirectory(char[],size_t);
extern PetscErrorCode PetscGetRealPath(const char[],char[]);
extern PetscErrorCode PetscGetHomeDirectory(char[],size_t);
extern PetscErrorCode PetscTestFile(const char[],char,PetscBool *);
extern PetscErrorCode PetscTestDirectory(const char[],char,PetscBool *);

extern PetscErrorCode PetscBinaryRead(int,void*,PetscInt,PetscDataType);
extern PetscErrorCode PetscBinarySynchronizedRead(MPI_Comm,int,void*,PetscInt,PetscDataType);
extern PetscErrorCode PetscBinarySynchronizedWrite(MPI_Comm,int,void*,PetscInt,PetscDataType,PetscBool );
extern PetscErrorCode PetscBinaryWrite(int,void*,PetscInt,PetscDataType,PetscBool );
extern PetscErrorCode PetscBinaryOpen(const char[],PetscFileMode,int *);
extern PetscErrorCode PetscBinaryClose(int);
extern PetscErrorCode PetscSharedTmp(MPI_Comm,PetscBool *);
extern PetscErrorCode PetscSharedWorkingDirectory(MPI_Comm,PetscBool *);
extern PetscErrorCode PetscGetTmp(MPI_Comm,char[],size_t);
extern PetscErrorCode PetscFileRetrieve(MPI_Comm,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscLs(MPI_Comm,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOpenSocket(char*,int,int*);
extern PetscErrorCode PetscWebServe(MPI_Comm,int);
# 2375 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef enum {PETSC_BINARY_SEEK_SET = 0,PETSC_BINARY_SEEK_CUR = 1,PETSC_BINARY_SEEK_END = 2} PetscBinarySeekType;
extern PetscErrorCode PetscBinarySeek(int,off_t,PetscBinarySeekType,off_t*);
extern PetscErrorCode PetscBinarySynchronizedSeek(MPI_Comm,int,off_t,PetscBinarySeekType,off_t*);

extern PetscErrorCode PetscSetDebugTerminal(const char[]);
extern PetscErrorCode PetscSetDebugger(const char[],PetscBool );
extern PetscErrorCode PetscSetDefaultDebugger(void);
extern PetscErrorCode PetscSetDebuggerFromString(char*);
extern PetscErrorCode PetscAttachDebugger(void);
extern PetscErrorCode PetscStopForDebugger(void);

extern PetscErrorCode PetscGatherNumberOfMessages(MPI_Comm,const PetscMPIInt[],const PetscMPIInt[],PetscMPIInt*);
extern PetscErrorCode PetscGatherMessageLengths(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],PetscMPIInt**,PetscMPIInt**);
extern PetscErrorCode PetscGatherMessageLengths2(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],const PetscMPIInt[],PetscMPIInt**,PetscMPIInt**,PetscMPIInt**);
extern PetscErrorCode PetscPostIrecvInt(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],const PetscMPIInt[],PetscInt***,MPI_Request**);
extern PetscErrorCode PetscPostIrecvScalar(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],const PetscMPIInt[],PetscScalar***,MPI_Request**);

extern PetscErrorCode PetscSSEIsEnabled(MPI_Comm,PetscBool *,PetscBool *);
# 2403 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
 typedef enum {NOT_SET_VALUES, INSERT_VALUES, ADD_VALUES, MAX_VALUES, INSERT_ALL_VALUES, ADD_ALL_VALUES} InsertMode;
# 2444 "/home/dpnkarthik/petsc-rnet/include/petscsys.h"
typedef struct _n_PetscSubcomm* PetscSubcomm;

struct _n_PetscSubcomm {
  MPI_Comm parent;
  MPI_Comm dupparent;
  MPI_Comm comm;
  PetscInt n;
  PetscInt color;
};

typedef enum {PETSC_SUBCOMM_GENERAL=0,PETSC_SUBCOMM_CONTIGUOUS=1,PETSC_SUBCOMM_INTERLACED=2} PetscSubcommType;
extern const char *PetscSubcommTypes[];

extern PetscErrorCode PetscSubcommCreate(MPI_Comm,PetscSubcomm*);
extern PetscErrorCode PetscSubcommDestroy(PetscSubcomm*);
extern PetscErrorCode PetscSubcommSetNumber(PetscSubcomm,PetscInt);
extern PetscErrorCode PetscSubcommSetType(PetscSubcomm,const PetscSubcommType);
extern PetscErrorCode PetscSubcommSetTypeGeneral(PetscSubcomm,PetscMPIInt,PetscMPIInt,PetscMPIInt);

# 1 "/home/dpnkarthik/petsc-rnet/include/petscctable.h" 1





struct _n_PetscTable {
  PetscInt *keytable;
  PetscInt *table;
  PetscInt count;
  PetscInt tablesize;
  PetscInt head;
  PetscInt maxkey;
};

typedef struct _n_PetscTable* PetscTable;
typedef PetscInt* PetscTablePosition;




extern PetscErrorCode PetscTableCreate(const PetscInt,PetscInt,PetscTable*);
extern PetscErrorCode PetscTableCreateCopy(const PetscTable,PetscTable*);
extern PetscErrorCode PetscTableDestroy(PetscTable*);
extern PetscErrorCode PetscTableGetCount(const PetscTable,PetscInt*);
extern PetscErrorCode PetscTableIsEmpty(const PetscTable,PetscInt*);
extern PetscErrorCode PetscTableAddExpand(PetscTable,PetscInt,PetscInt);
extern PetscErrorCode PetscTableAddCountExpand(PetscTable,PetscInt);
extern PetscErrorCode PetscTableGetHeadPosition(PetscTable,PetscTablePosition*);
extern PetscErrorCode PetscTableGetNext(PetscTable,PetscTablePosition*,PetscInt*,PetscInt*);
extern PetscErrorCode PetscTableRemoveAll(PetscTable);



static inline PetscErrorCode PetscTableAdd(PetscTable ta,PetscInt key,PetscInt data)
{
  PetscErrorCode ierr;
  PetscInt ii = 0,hash = ((unsigned long)((79943*(unsigned long)key)%ta->tablesize));

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petscctable.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 39; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscTableAdd") && strcmp("PetscTableAdd","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petscctable.h",39,"PetscTableAdd","__func__",__func__); } } while (0); } while (0);
  if (key <= 0) return PetscError(((MPI_Comm)0x44000001),40,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"key <= 0");
  if (key > ta->maxkey) return PetscError(((MPI_Comm)0x44000001),41,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"key %D is greater than largest key allowed %D",key,ta->maxkey);
  if (!data) return PetscError(((MPI_Comm)0x44000001),42,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"Null data");

  if (ta->count < 5*(ta->tablesize/6) - 1) {
    while (ii++ < ta->tablesize){
      if (ta->keytable[hash] == key) {
 ta->table[hash] = data;
 do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
      } else if (!ta->keytable[hash]) {
 ta->count++;
 ta->keytable[hash] = key; ta->table[hash] = data;
 do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
      }
      hash = (hash == (ta->tablesize-1)) ? 0 : hash+1;
    }
    return PetscError(((MPI_Comm)0x44000001),56,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",74,PETSC_ERROR_INITIAL,"Full table");
  } else {
    ierr = PetscTableAddExpand(ta,key,data);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),58,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode PetscTableAddCount(PetscTable ta,PetscInt key)
{
  PetscErrorCode ierr;
  PetscInt ii = 0,hash = ((unsigned long)((79943*(unsigned long)key)%ta->tablesize));

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petscctable.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 70; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscTableAddCount") && strcmp("PetscTableAddCount","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petscctable.h",70,"PetscTableAddCount","__func__",__func__); } } while (0); } while (0);
  if (key <= 0) return PetscError(((MPI_Comm)0x44000001),71,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"key <= 0");
  if (key > ta->maxkey) return PetscError(((MPI_Comm)0x44000001),72,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"key %D is greater than largest key allowed %D",key,ta->maxkey);

  if (ta->count < 5*(ta->tablesize/6) - 1) {
    while (ii++ < ta->tablesize){
      if (ta->keytable[hash] == key) {
 do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
      } else if (!ta->keytable[hash]) {
 ta->count++;
 ta->keytable[hash] = key; ta->table[hash] = ta->count;
 do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
      }
      hash = (hash == (ta->tablesize-1)) ? 0 : hash+1;
    }
    return PetscError(((MPI_Comm)0x44000001),85,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",74,PETSC_ERROR_INITIAL,"Full table");
  } else {
    ierr = PetscTableAddCountExpand(ta,key);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),87,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 101 "/home/dpnkarthik/petsc-rnet/include/petscctable.h"
static inline PetscErrorCode PetscTableFind(PetscTable ta,PetscInt key,PetscInt *data)
{
  PetscInt hash,ii = 0;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petscctable.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 105; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscTableFind") && strcmp("PetscTableFind","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petscctable.h",105,"PetscTableFind","__func__",__func__); } } while (0); } while (0);
  if (key <= 0) return PetscError(((MPI_Comm)0x44000001),106,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"Key <= 0");
  if (key > ta->maxkey) return PetscError(((MPI_Comm)0x44000001),107,__func__,"/home/dpnkarthik/petsc-rnet/include/petscctable.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"key %D is greater than largest key allowed %D",key,ta->maxkey);

  hash = ((unsigned long)((79943*(unsigned long)key)%ta->tablesize));
  *data = 0;
  while (ii++ < ta->tablesize) {
    if (!ta->keytable[hash]) break;
    else if (ta->keytable[hash] == key) {
      *data = ta->table[hash];
      break;
    }
    hash = (hash == (ta->tablesize-1)) ? 0 : hash+1;
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}


# 2464 "/home/dpnkarthik/petsc-rnet/include/petscsys.h" 2


# 8 "/home/dpnkarthik/petsc-rnet/include/petscis.h" 2



extern PetscClassId IS_CLASSID;

extern PetscErrorCode ISInitializePackage(const char[]);
# 24 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
typedef struct _p_IS* IS;
# 41 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
extern PetscFList ISList;
extern PetscBool ISRegisterAllCalled;
extern PetscErrorCode ISSetType(IS, const char*);
extern PetscErrorCode ISGetType(IS, const char* *);
extern PetscErrorCode ISRegister(const char[],const char[],const char[],PetscErrorCode (*)(IS));
extern PetscErrorCode ISRegisterAll(const char []);
extern PetscErrorCode ISRegisterDestroy(void);
extern PetscErrorCode ISCreate(MPI_Comm,IS*);
# 101 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
extern PetscErrorCode ISCreateGeneral(MPI_Comm,PetscInt,const PetscInt[],PetscCopyMode,IS *);
extern PetscErrorCode ISGeneralSetIndices(IS,PetscInt,const PetscInt[],PetscCopyMode);
extern PetscErrorCode ISCreateBlock(MPI_Comm,PetscInt,PetscInt,const PetscInt[],PetscCopyMode,IS *);
extern PetscErrorCode ISBlockSetIndices(IS,PetscInt,PetscInt,const PetscInt[],PetscCopyMode);
extern PetscErrorCode ISCreateStride(MPI_Comm,PetscInt,PetscInt,PetscInt,IS *);
extern PetscErrorCode ISStrideSetStride(IS,PetscInt,PetscInt,PetscInt);

extern PetscErrorCode ISDestroy(IS*);
extern PetscErrorCode ISSetPermutation(IS);
extern PetscErrorCode ISPermutation(IS,PetscBool *);
extern PetscErrorCode ISSetIdentity(IS);
extern PetscErrorCode ISIdentity(IS,PetscBool *);
extern PetscErrorCode ISContiguousLocal(IS,PetscInt,PetscInt,PetscInt*,PetscBool*);

extern PetscErrorCode ISGetIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISRestoreIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISGetTotalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISRestoreTotalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISGetNonlocalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISRestoreNonlocalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISGetNonlocalIS(IS, IS *is);
extern PetscErrorCode ISRestoreNonlocalIS(IS, IS *is);
extern PetscErrorCode ISGetSize(IS,PetscInt *);
extern PetscErrorCode ISGetLocalSize(IS,PetscInt *);
extern PetscErrorCode ISInvertPermutation(IS,PetscInt,IS*);
extern PetscErrorCode ISView(IS,PetscViewer);
extern PetscErrorCode ISEqual(IS,IS,PetscBool *);
extern PetscErrorCode ISSort(IS);
extern PetscErrorCode ISSorted(IS,PetscBool *);
extern PetscErrorCode ISDifference(IS,IS,IS*);
extern PetscErrorCode ISSum(IS,IS,IS*);
extern PetscErrorCode ISExpand(IS,IS,IS*);

extern PetscErrorCode ISBlockGetIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISBlockRestoreIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISBlockGetLocalSize(IS,PetscInt *);
extern PetscErrorCode ISBlockGetSize(IS,PetscInt *);
extern PetscErrorCode ISGetBlockSize(IS,PetscInt*);
extern PetscErrorCode ISSetBlockSize(IS,PetscInt);

extern PetscErrorCode ISStrideGetInfo(IS,PetscInt *,PetscInt*);

extern PetscErrorCode ISToGeneral(IS);

extern PetscErrorCode ISDuplicate(IS,IS*);
extern PetscErrorCode ISCopy(IS,IS);
extern PetscErrorCode ISAllGather(IS,IS*);
extern PetscErrorCode ISComplement(IS,PetscInt,PetscInt,IS*);
extern PetscErrorCode ISConcatenate(MPI_Comm,PetscInt,const IS[],IS*);
extern PetscErrorCode ISListToColoring(MPI_Comm,PetscInt, IS[],IS*,IS*);
extern PetscErrorCode ISColoringToList(IS, IS, PetscInt*, IS *[]);
extern PetscErrorCode ISOnComm(IS,MPI_Comm,PetscCopyMode,IS*);


extern PetscClassId IS_LTOGM_CLASSID;
# 174 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
struct _p_ISLocalToGlobalMapping{
  _p_PetscObject hdr; int *ops;
  PetscInt n;
  PetscInt *indices;
  PetscInt globalstart;
  PetscInt globalend;
  PetscInt *globals;
};
typedef struct _p_ISLocalToGlobalMapping* ISLocalToGlobalMapping;
# 195 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
typedef enum {IS_GTOLM_MASK,IS_GTOLM_DROP} ISGlobalToLocalMappingType;

extern PetscErrorCode ISLocalToGlobalMappingCreate(MPI_Comm,PetscInt,const PetscInt[],PetscCopyMode,ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingCreateIS(IS,ISLocalToGlobalMapping *);
extern PetscErrorCode ISLocalToGlobalMappingView(ISLocalToGlobalMapping,PetscViewer);
extern PetscErrorCode ISLocalToGlobalMappingDestroy(ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingApplyIS(ISLocalToGlobalMapping,IS,IS*);
extern PetscErrorCode ISGlobalToLocalMappingApply(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,const PetscInt[],PetscInt*,PetscInt[]);
extern PetscErrorCode ISLocalToGlobalMappingGetSize(ISLocalToGlobalMapping,PetscInt*);
extern PetscErrorCode ISLocalToGlobalMappingGetInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
extern PetscErrorCode ISLocalToGlobalMappingRestoreInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
extern PetscErrorCode ISLocalToGlobalMappingGetIndices(ISLocalToGlobalMapping,const PetscInt**);
extern PetscErrorCode ISLocalToGlobalMappingRestoreIndices(ISLocalToGlobalMapping,const PetscInt**);
extern PetscErrorCode ISLocalToGlobalMappingBlock(ISLocalToGlobalMapping,PetscInt,ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingUnBlock(ISLocalToGlobalMapping,PetscInt,ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingConcatenate(MPI_Comm,PetscInt,const ISLocalToGlobalMapping[],ISLocalToGlobalMapping*);



static inline PetscErrorCode ISLocalToGlobalMappingApply(ISLocalToGlobalMapping mapping,PetscInt N,const PetscInt in[],PetscInt out[])
{
  PetscInt i,Nmax = mapping->n;
  const PetscInt *idx = mapping->indices;
  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/petscis.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 218; petscstack->currentsize++; } do { if (strcmp(__func__,"ISLocalToGlobalMappingApply") && strcmp("ISLocalToGlobalMappingApply","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/petscis.h",218,"ISLocalToGlobalMappingApply","__func__",__func__); } } while (0); } while (0);
  for (i=0; i<N; i++) {
    if (in[i] < 0) {out[i] = in[i]; continue;}
    if (in[i] >= Nmax) return PetscError(((MPI_Comm)0x44000001),221,__func__,"/home/dpnkarthik/petsc-rnet/include/petscis.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"Local index %D too large %D (max) at %D",in[i],Nmax,i);
    out[i] = idx[in[i]];
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 244 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
typedef enum {IS_COLORING_GLOBAL,IS_COLORING_GHOSTED} ISColoringType;
extern const char *ISColoringTypes[];
typedef unsigned short ISColoringValue;
extern PetscErrorCode ISAllGatherColors(MPI_Comm,PetscInt,ISColoringValue*,PetscInt*,ISColoringValue*[]);
# 262 "/home/dpnkarthik/petsc-rnet/include/petscis.h"
struct _n_ISColoring {
  PetscInt refct;
  PetscInt n;
  IS *is;
  MPI_Comm comm;
  ISColoringValue *colors;
  PetscInt N;
  ISColoringType ctype;
};
typedef struct _n_ISColoring* ISColoring;

extern PetscErrorCode ISColoringCreate(MPI_Comm,PetscInt,PetscInt,const ISColoringValue[],ISColoring*);
extern PetscErrorCode ISColoringDestroy(ISColoring*);
extern PetscErrorCode ISColoringView(ISColoring,PetscViewer);
extern PetscErrorCode ISColoringGetIS(ISColoring,PetscInt*,IS*[]);
extern PetscErrorCode ISColoringRestoreIS(ISColoring,IS*[]);





extern PetscErrorCode ISPartitioningToNumbering(IS,IS*);
extern PetscErrorCode ISPartitioningCount(IS,PetscInt,PetscInt[]);

extern PetscErrorCode ISCompressIndicesGeneral(PetscInt,PetscInt,PetscInt,PetscInt,const IS[],IS[]);
extern PetscErrorCode ISCompressIndicesSorted(PetscInt,PetscInt,PetscInt,const IS[],IS[]);
extern PetscErrorCode ISExpandIndicesGeneral(PetscInt,PetscInt,PetscInt,PetscInt,const IS[],IS[]);


# 10 "/home/dpnkarthik/petsc-rnet/include/petscvec.h" 2


# 22 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
typedef struct _p_Vec* Vec;
# 34 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
typedef struct _p_VecScatter* VecScatter;
# 43 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
typedef enum {SCATTER_FORWARD=0, SCATTER_REVERSE=1, SCATTER_FORWARD_LOCAL=2, SCATTER_REVERSE_LOCAL=3, SCATTER_LOCAL=2} ScatterMode;
# 114 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
extern PetscClassId VEC_CLASSID;
extern PetscClassId VEC_SCATTER_CLASSID;


extern PetscErrorCode VecInitializePackage(const char[]);
extern PetscErrorCode VecFinalizePackage(void);

extern PetscErrorCode VecCreate(MPI_Comm,Vec*);

extern PetscErrorCode VecCreateSeq(MPI_Comm,PetscInt,Vec*);

extern PetscErrorCode VecCreateMPI(MPI_Comm,PetscInt,PetscInt,Vec*);

extern PetscErrorCode VecCreateSeqWithArray(MPI_Comm,PetscInt,const PetscScalar[],Vec*);

extern PetscErrorCode VecCreateMPIWithArray(MPI_Comm,PetscInt,PetscInt,const PetscScalar[],Vec*);

extern PetscErrorCode VecCreateShared(MPI_Comm,PetscInt,PetscInt,Vec*);
extern PetscErrorCode VecSetFromOptions(Vec);
extern PetscErrorCode VecSetUp(Vec);
extern PetscErrorCode VecDestroy(Vec*);
extern PetscErrorCode VecZeroEntries(Vec);
extern PetscErrorCode VecSetOptionsPrefix(Vec,const char[]);
extern PetscErrorCode VecAppendOptionsPrefix(Vec,const char[]);
extern PetscErrorCode VecGetOptionsPrefix(Vec,const char*[]);

extern PetscErrorCode VecSetSizes(Vec,PetscInt,PetscInt);

extern PetscErrorCode VecDotNorm2(Vec,Vec,PetscScalar*,PetscScalar*);
extern PetscErrorCode VecDot(Vec,Vec,PetscScalar*);

extern PetscErrorCode VecTDot(Vec,Vec,PetscScalar*);

extern PetscErrorCode VecMDot(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMTDot(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecGetSubVector(Vec,IS,Vec*);
extern PetscErrorCode VecRestoreSubVector(Vec,IS,Vec*);
# 159 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
extern const char *NormTypes[];
# 220 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
extern PetscErrorCode VecNorm(Vec,NormType,PetscReal *);
extern PetscErrorCode VecNormAvailable(Vec,NormType,PetscBool *,PetscReal *);



extern PetscErrorCode VecNormalize(Vec,PetscReal *);
extern PetscErrorCode VecSum(Vec,PetscScalar*);
extern PetscErrorCode VecMax(Vec,PetscInt*,PetscReal *);

extern PetscErrorCode VecMin(Vec,PetscInt*,PetscReal *);

extern PetscErrorCode VecScale(Vec,PetscScalar);
extern PetscErrorCode VecCopy(Vec,Vec);
extern PetscErrorCode VecSetRandom(Vec,PetscRandom);
extern PetscErrorCode VecSet(Vec,PetscScalar);
extern PetscErrorCode VecSwap(Vec,Vec);
extern PetscErrorCode VecAXPY(Vec,PetscScalar,Vec);
extern PetscErrorCode VecAXPBY(Vec,PetscScalar,PetscScalar,Vec);
extern PetscErrorCode VecMAXPY(Vec,PetscInt,const PetscScalar[],Vec[]);
extern PetscErrorCode VecAYPX(Vec,PetscScalar,Vec);
extern PetscErrorCode VecWAXPY(Vec,PetscScalar,Vec,Vec);
extern PetscErrorCode VecAXPBYPCZ(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
extern PetscErrorCode VecPointwiseMax(Vec,Vec,Vec);

extern PetscErrorCode VecPointwiseMaxAbs(Vec,Vec,Vec);

extern PetscErrorCode VecPointwiseMin(Vec,Vec,Vec);

extern PetscErrorCode VecPointwiseMult(Vec,Vec,Vec);

extern PetscErrorCode VecPointwiseDivide(Vec,Vec,Vec);

extern PetscErrorCode VecMaxPointwiseDivide(Vec,Vec,PetscReal*);
extern PetscErrorCode VecShift(Vec,PetscScalar);
extern PetscErrorCode VecReciprocal(Vec);
extern PetscErrorCode VecPermute(Vec, IS, PetscBool );
extern PetscErrorCode VecSqrtAbs(Vec);
extern PetscErrorCode VecLog(Vec);
extern PetscErrorCode VecExp(Vec);
extern PetscErrorCode VecAbs(Vec);
extern PetscErrorCode VecDuplicate(Vec,Vec*);
extern PetscErrorCode VecDuplicateVecs(Vec,PetscInt,Vec*[]);
extern PetscErrorCode VecDestroyVecs(PetscInt, Vec*[]);
extern PetscErrorCode VecStrideNormAll(Vec,NormType,PetscReal[]);
extern PetscErrorCode VecStrideMaxAll(Vec,PetscInt [],PetscReal []);
extern PetscErrorCode VecStrideMinAll(Vec,PetscInt [],PetscReal []);
extern PetscErrorCode VecStrideScaleAll(Vec,const PetscScalar[]);

extern PetscErrorCode VecStrideNorm(Vec,PetscInt,NormType,PetscReal*);


extern PetscErrorCode VecStrideMax(Vec,PetscInt,PetscInt *,PetscReal *);

extern PetscErrorCode VecStrideMin(Vec,PetscInt,PetscInt *,PetscReal *);

extern PetscErrorCode VecStrideScale(Vec,PetscInt,PetscScalar);
extern PetscErrorCode VecStrideSet(Vec,PetscInt,PetscScalar);


extern PetscErrorCode VecStrideGather(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode VecStrideScatter(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode VecStrideGatherAll(Vec,Vec[],InsertMode);
extern PetscErrorCode VecStrideScatterAll(Vec[],Vec,InsertMode);

extern PetscErrorCode VecSetValues(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode VecGetValues(Vec,PetscInt,const PetscInt[],PetscScalar[]);
extern PetscErrorCode VecAssemblyBegin(Vec);
extern PetscErrorCode VecAssemblyEnd(Vec);
extern PetscErrorCode VecStashSetInitialSize(Vec,PetscInt,PetscInt);
extern PetscErrorCode VecStashView(Vec,PetscViewer);
extern PetscErrorCode VecStashGetInfo(Vec,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
# 319 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
static inline PetscErrorCode VecSetValue(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValues(v,1,&i,&va,mode);}


extern PetscErrorCode VecSetBlockSize(Vec,PetscInt);
extern PetscErrorCode VecGetBlockSize(Vec,PetscInt*);

extern PetscErrorCode VecSetValuesBlocked(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);


extern PetscFList VecList;
extern PetscBool VecRegisterAllCalled;
extern PetscErrorCode VecSetType(Vec, const char*);
extern PetscErrorCode VecGetType(Vec, const char* *);
extern PetscErrorCode VecRegister(const char[],const char[],const char[],PetscErrorCode (*)(Vec));
extern PetscErrorCode VecRegisterAll(const char []);
extern PetscErrorCode VecRegisterDestroy(void);
# 385 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
extern PetscErrorCode VecScatterCreate(Vec,IS,Vec,IS,VecScatter *);





extern PetscErrorCode VecScatterCreateEmpty(MPI_Comm,VecScatter *);
extern PetscErrorCode VecScatterCreateLocal(VecScatter,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],PetscInt);
extern PetscErrorCode VecScatterBegin(VecScatter,Vec,Vec,InsertMode,ScatterMode);
extern PetscErrorCode VecScatterEnd(VecScatter,Vec,Vec,InsertMode,ScatterMode);
extern PetscErrorCode VecScatterDestroy(VecScatter*);
extern PetscErrorCode VecScatterCopy(VecScatter,VecScatter *);
extern PetscErrorCode VecScatterView(VecScatter,PetscViewer);
extern PetscErrorCode VecScatterRemap(VecScatter,PetscInt *,PetscInt*);
extern PetscErrorCode VecScatterGetMerged(VecScatter,PetscBool *);

extern PetscErrorCode VecGetArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode VecRestoreArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode VecGetArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode VecRestoreArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode VecGetArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode VecRestoreArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode VecGetArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);
extern PetscErrorCode VecRestoreArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);

extern PetscErrorCode VecPlaceArray(Vec,const PetscScalar[]);
extern PetscErrorCode VecResetArray(Vec);
extern PetscErrorCode VecReplaceArray(Vec,const PetscScalar[]);
extern PetscErrorCode VecGetArrays(const Vec[],PetscInt,PetscScalar**[]);
extern PetscErrorCode VecRestoreArrays(const Vec[],PetscInt,PetscScalar**[]);

extern PetscErrorCode VecView(Vec,PetscViewer);
extern PetscErrorCode VecViewFromOptions(Vec, const char *);
extern PetscErrorCode VecEqual(Vec,Vec,PetscBool *);

extern PetscErrorCode VecLoad(Vec, PetscViewer);

extern PetscErrorCode VecGetSize(Vec,PetscInt*);

extern PetscErrorCode VecGetLocalSize(Vec,PetscInt*);

extern PetscErrorCode VecGetOwnershipRange(Vec,PetscInt*,PetscInt*);
extern PetscErrorCode VecGetOwnershipRanges(Vec,const PetscInt *[]);

extern PetscErrorCode VecSetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping);
extern PetscErrorCode VecSetValuesLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
# 459 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
static inline PetscErrorCode VecSetValueLocal(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValuesLocal(v,1,&i,&va,mode);}

extern PetscErrorCode VecSetLocalToGlobalMappingBlock(Vec,ISLocalToGlobalMapping);
extern PetscErrorCode VecSetValuesBlockedLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode VecGetLocalToGlobalMappingBlock(Vec,ISLocalToGlobalMapping*);
extern PetscErrorCode VecGetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping*);

extern PetscErrorCode VecDotBegin(Vec,Vec,PetscScalar *);

extern PetscErrorCode VecDotEnd(Vec,Vec,PetscScalar *);

extern PetscErrorCode VecTDotBegin(Vec,Vec,PetscScalar *);

extern PetscErrorCode VecTDotEnd(Vec,Vec,PetscScalar *);

extern PetscErrorCode VecNormBegin(Vec,NormType,PetscReal *);


extern PetscErrorCode VecNormEnd(Vec,NormType,PetscReal *);



extern PetscErrorCode VecMDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMTDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMTDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);


typedef enum {VEC_IGNORE_OFF_PROC_ENTRIES,VEC_IGNORE_NEGATIVE_INDICES} VecOption;
extern PetscErrorCode VecSetOption(Vec,VecOption,PetscBool );






# 1 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h" 1
# 11 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscvec.h" 1
# 12 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h" 2

# 22 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
typedef struct _n_PetscLayout* PetscLayout;
struct _n_PetscLayout{
  MPI_Comm comm;
  PetscInt n,N;
  PetscInt rstart,rend;
  PetscInt *range;
  PetscInt bs;
  PetscInt refcnt;
  ISLocalToGlobalMapping mapping;
  ISLocalToGlobalMapping bmapping;
};

extern PetscErrorCode PetscLayoutCreate(MPI_Comm,PetscLayout*);
extern PetscErrorCode PetscLayoutSetUp(PetscLayout);
extern PetscErrorCode PetscLayoutDestroy(PetscLayout*);
extern PetscErrorCode PetscLayoutCopy(PetscLayout,PetscLayout*);
extern PetscErrorCode PetscLayoutReference(PetscLayout,PetscLayout*);
extern PetscErrorCode PetscLayoutSetLocalSize(PetscLayout,PetscInt);
extern PetscErrorCode PetscLayoutGetLocalSize(PetscLayout,PetscInt *);

extern PetscErrorCode PetscLayoutSetSize(PetscLayout,PetscInt);
extern PetscErrorCode PetscLayoutGetSize(PetscLayout,PetscInt *);

extern PetscErrorCode PetscLayoutSetBlockSize(PetscLayout,PetscInt);
extern PetscErrorCode PetscLayoutGetBlockSize(PetscLayout,PetscInt*);
extern PetscErrorCode PetscLayoutGetRange(PetscLayout,PetscInt *,PetscInt *);
extern PetscErrorCode PetscLayoutGetRanges(PetscLayout,const PetscInt *[]);
extern PetscErrorCode PetscLayoutSetISLocalToGlobalMapping(PetscLayout,ISLocalToGlobalMapping);
extern PetscErrorCode PetscLayoutSetISLocalToGlobalMappingBlock(PetscLayout,ISLocalToGlobalMapping);
# 72 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
static inline PetscErrorCode PetscLayoutFindOwner(PetscLayout map,PetscInt idx,PetscInt *owner)
{
  PetscErrorCode ierr;
  PetscMPIInt lo = 0,hi,t;
  PetscInt bs = map->bs;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 78; petscstack->currentsize++; } do { if (strcmp(__func__,"PetscLayoutFindOwner") && strcmp("PetscLayoutFindOwner","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h",78,"PetscLayoutFindOwner","__func__",__func__); } } while (0); } while (0);
  if (!((map->n >= 0) && (map->N >= 0) && (map->range))) return PetscError(((MPI_Comm)0x44000001),79,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",73,PETSC_ERROR_INITIAL,"PetscLayoutSetUp() must be called first");
  if (idx < 0 || idx > map->N) return PetscError(((MPI_Comm)0x44000001),80,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"Index %D is out of range",idx);
  ierr = MPI_Comm_size(map->comm,&hi);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),81,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  while (hi - lo > 1) {
    t = lo + (hi - lo) / 2;
    if (idx < map->range[t]/bs) hi = t;
    else lo = t;
  }
  *owner = lo;
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}


typedef struct _n_PetscUniformSection *PetscUniformSection;
struct _n_PetscUniformSection {
  MPI_Comm comm;
  PetscInt pStart, pEnd;
  PetscInt numDof;
};
# 113 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
typedef struct _n_PetscSection *PetscSection;
struct _n_PetscSection {
  struct _n_PetscUniformSection atlasLayout;
  PetscInt *atlasDof;
  PetscInt *atlasOff;
  PetscSection bc;
  PetscInt *bcIndices;
  PetscInt refcnt;
};

extern PetscErrorCode PetscSectionCreate(MPI_Comm,PetscSection*);
extern PetscErrorCode PetscSectionGetChart(PetscSection, PetscInt *, PetscInt *);
extern PetscErrorCode PetscSectionSetChart(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetDof(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetDof(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetConstraintDof(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetConstraintDof(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetConstraintIndices(PetscSection, PetscInt, PetscInt**);
extern PetscErrorCode PetscSectionSetConstraintIndices(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetUp(PetscSection);
extern PetscErrorCode PetscSectionGetStorageSize(PetscSection, PetscInt*);
extern PetscErrorCode PetscSectionGetOffset(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionDestroy(PetscSection*);

extern PetscErrorCode VecGetValuesSection(Vec, PetscSection, PetscInt, PetscScalar **);
extern PetscErrorCode VecSetValuesSection(Vec, PetscSection, PetscInt, PetscScalar [], InsertMode);



typedef struct _VecOps *VecOps;
struct _VecOps {
  PetscErrorCode (*duplicate)(Vec,Vec*);
  PetscErrorCode (*duplicatevecs)(Vec,PetscInt,Vec**);
  PetscErrorCode (*destroyvecs)(PetscInt,Vec[]);
  PetscErrorCode (*dot)(Vec,Vec,PetscScalar*);
  PetscErrorCode (*mdot)(Vec,PetscInt,const Vec[],PetscScalar*);
  PetscErrorCode (*norm)(Vec,NormType,PetscReal*);
  PetscErrorCode (*tdot)(Vec,Vec,PetscScalar*);
  PetscErrorCode (*mtdot)(Vec,PetscInt,const Vec[],PetscScalar*);
  PetscErrorCode (*scale)(Vec,PetscScalar);
  PetscErrorCode (*copy)(Vec,Vec);
  PetscErrorCode (*set)(Vec,PetscScalar);
  PetscErrorCode (*swap)(Vec,Vec);
  PetscErrorCode (*axpy)(Vec,PetscScalar,Vec);
  PetscErrorCode (*axpby)(Vec,PetscScalar,PetscScalar,Vec);
  PetscErrorCode (*maxpy)(Vec,PetscInt,const PetscScalar*,Vec*);
  PetscErrorCode (*aypx)(Vec,PetscScalar,Vec);
  PetscErrorCode (*waxpy)(Vec,PetscScalar,Vec,Vec);
  PetscErrorCode (*axpbypcz)(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
  PetscErrorCode (*pointwisemult)(Vec,Vec,Vec);
  PetscErrorCode (*pointwisedivide)(Vec,Vec,Vec);
  PetscErrorCode (*setvalues)(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*assemblybegin)(Vec);
  PetscErrorCode (*assemblyend)(Vec);
  PetscErrorCode (*getarray)(Vec,PetscScalar**);
  PetscErrorCode (*getsize)(Vec,PetscInt*);
  PetscErrorCode (*getlocalsize)(Vec,PetscInt*);
  PetscErrorCode (*restorearray)(Vec,PetscScalar**);
  PetscErrorCode (*max)(Vec,PetscInt*,PetscReal*);
  PetscErrorCode (*min)(Vec,PetscInt*,PetscReal*);
  PetscErrorCode (*setrandom)(Vec,PetscRandom);
  PetscErrorCode (*setoption)(Vec,VecOption,PetscBool );
  PetscErrorCode (*setvaluesblocked)(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*destroy)(Vec);
  PetscErrorCode (*view)(Vec,PetscViewer);
  PetscErrorCode (*placearray)(Vec,const PetscScalar*);
  PetscErrorCode (*replacearray)(Vec,const PetscScalar*);
  PetscErrorCode (*dot_local)(Vec,Vec,PetscScalar*);
  PetscErrorCode (*tdot_local)(Vec,Vec,PetscScalar*);
  PetscErrorCode (*norm_local)(Vec,NormType,PetscReal*);
  PetscErrorCode (*mdot_local)(Vec,PetscInt,const Vec[],PetscScalar*);
  PetscErrorCode (*mtdot_local)(Vec,PetscInt,const Vec[],PetscScalar*);
  PetscErrorCode (*load)(Vec,PetscViewer);
  PetscErrorCode (*reciprocal)(Vec);
  PetscErrorCode (*conjugate)(Vec);
  PetscErrorCode (*setlocaltoglobalmapping)(Vec,ISLocalToGlobalMapping);
  PetscErrorCode (*setvalueslocal)(Vec,PetscInt,const PetscInt *,const PetscScalar *,InsertMode);
  PetscErrorCode (*resetarray)(Vec);
  PetscErrorCode (*setfromoptions)(Vec);
  PetscErrorCode (*maxpointwisedivide)(Vec,Vec,PetscReal*);
  PetscErrorCode (*pointwisemax)(Vec,Vec,Vec);
  PetscErrorCode (*pointwisemaxabs)(Vec,Vec,Vec);
  PetscErrorCode (*pointwisemin)(Vec,Vec,Vec);
  PetscErrorCode (*getvalues)(Vec,PetscInt,const PetscInt[],PetscScalar[]);
  PetscErrorCode (*sqrt)(Vec);
  PetscErrorCode (*abs)(Vec);
  PetscErrorCode (*exp)(Vec);
  PetscErrorCode (*log)(Vec);
  PetscErrorCode (*shift)(Vec);
  PetscErrorCode (*create)(Vec);
  PetscErrorCode (*stridegather)(Vec,PetscInt,Vec,InsertMode);
  PetscErrorCode (*stridescatter)(Vec,PetscInt,Vec,InsertMode);
  PetscErrorCode (*dotnorm2)(Vec,Vec,PetscScalar*,PetscScalar*);
  PetscErrorCode (*getsubvector)(Vec,IS,Vec*);
  PetscErrorCode (*restoresubvector)(Vec,IS,Vec*);
};







typedef struct {
  PetscInt nmax;
  PetscInt umax;
  PetscInt oldnmax;
  PetscInt n;
  PetscInt bs;
  PetscInt reallocs;
  PetscInt *idx;
  PetscScalar *array;

  MPI_Comm comm;
  PetscMPIInt size,rank;
  PetscMPIInt tag1,tag2;
  MPI_Request *send_waits;
  MPI_Request *recv_waits;
  MPI_Status *send_status;
  PetscInt nsends,nrecvs;
  PetscScalar *svalues,*rvalues;
  PetscInt *sindices,*rindices;
  PetscInt rmax;
  PetscInt *nprocs;
  PetscInt nprocessed;
  PetscBool donotstash;
  PetscBool ignorenegidx;
  InsertMode insertmode;
  PetscInt *bowners;
} VecStash;


typedef enum {PETSC_CUSP_UNALLOCATED,PETSC_CUSP_GPU,PETSC_CUSP_CPU,PETSC_CUSP_BOTH} PetscCUSPFlag;




struct _p_Vec {
  _p_PetscObject hdr; struct _VecOps *ops;
  PetscLayout map;
  void *data;
  PetscBool array_gotten;
  VecStash stash,bstash;
  PetscBool petscnative;

  PetscCUSPFlag valid_GPU_array;
  void *spptr;

};

extern PetscLogEvent VEC_View, VEC_Max, VEC_Min, VEC_DotBarrier, VEC_Dot, VEC_MDotBarrier, VEC_MDot, VEC_TDot, VEC_MTDot;
extern PetscLogEvent VEC_Norm, VEC_Normalize, VEC_Scale, VEC_Copy, VEC_Set, VEC_AXPY, VEC_AYPX, VEC_WAXPY, VEC_MAXPY;
extern PetscLogEvent VEC_AssemblyEnd, VEC_PointwiseMult, VEC_SetValues, VEC_Load, VEC_ScatterBarrier, VEC_ScatterBegin, VEC_ScatterEnd;
extern PetscLogEvent VEC_SetRandom, VEC_ReduceArithmetic, VEC_ReduceBarrier, VEC_ReduceCommunication;
extern PetscLogEvent VEC_Swap, VEC_AssemblyBegin, VEC_NormBarrier, VEC_DotNormBarrier, VEC_DotNorm, VEC_AXPBYPCZ, VEC_Ops;
extern PetscLogEvent VEC_CUSPCopyToGPU, VEC_CUSPCopyFromGPU;
extern PetscLogEvent VEC_CUSPCopyToGPUSome, VEC_CUSPCopyFromGPUSome;


extern PetscErrorCode VecCUSPCopyFromGPU(Vec v);




static inline PetscErrorCode VecGetArrayRead(Vec x,const PetscScalar *a[])
{
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 281; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGetArrayRead") && strcmp("VecGetArrayRead","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h",281,"VecGetArrayRead","__func__",__func__); } } while (0); } while (0);
  if (x->petscnative){

    if (x->valid_GPU_array == PETSC_CUSP_GPU || !*((PetscScalar**)x->data)){
      ierr = VecCUSPCopyFromGPU(x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),285,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }

    *a = *((PetscScalar **)x->data);
  } else {
    ierr = (*x->ops->getarray)(x,(PetscScalar**)a);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),290,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode VecRestoreArrayRead(Vec x,const PetscScalar *a[])
{
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 301; petscstack->currentsize++; } do { if (strcmp(__func__,"VecRestoreArrayRead") && strcmp("VecRestoreArrayRead","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h",301,"VecRestoreArrayRead","__func__",__func__); } } while (0); } while (0);
  if (x->petscnative){

    if (x->valid_GPU_array != PETSC_CUSP_UNALLOCATED) {
      x->valid_GPU_array = PETSC_CUSP_BOTH;
    }

  } else {
    ierr = (*x->ops->restorearray)(x,(PetscScalar**)a);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),309,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode VecGetArray(Vec x,PetscScalar *a[])
{
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 320; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGetArray") && strcmp("VecGetArray","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h",320,"VecGetArray","__func__",__func__); } } while (0); } while (0);
  if (x->petscnative){

    if (x->valid_GPU_array == PETSC_CUSP_GPU || !*((PetscScalar**)x->data)){
      ierr = VecCUSPCopyFromGPU(x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),324,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }

    *a = *((PetscScalar **)x->data);
  } else {
    ierr = (*x->ops->getarray)(x,a);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),329,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode VecRestoreArray(Vec x,PetscScalar *a[])
{
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 340; petscstack->currentsize++; } do { if (strcmp(__func__,"VecRestoreArray") && strcmp("VecRestoreArray","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h",340,"VecRestoreArray","__func__",__func__); } } while (0); } while (0);
  if (x->petscnative){

    if (x->valid_GPU_array != PETSC_CUSP_UNALLOCATED) {
      x->valid_GPU_array = PETSC_CUSP_CPU;
    }

  } else {
    ierr = (*x->ops->restorearray)(x,a);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),348,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (((PetscObject)x)->state++,0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),350,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 366 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
extern PetscErrorCode VecDuplicateVecs_Default(Vec,PetscInt,Vec *[]);
extern PetscErrorCode VecDestroyVecs_Default(PetscInt,Vec []);
extern PetscErrorCode VecLoad_Binary(Vec, PetscViewer);
extern PetscErrorCode VecLoad_Default(Vec, PetscViewer);

extern PetscInt NormIds[7];





typedef enum { VEC_SCATTER_SEQ_GENERAL,VEC_SCATTER_SEQ_STRIDE,
               VEC_SCATTER_MPI_GENERAL,VEC_SCATTER_MPI_TOALL,
               VEC_SCATTER_MPI_TOONE} VecScatterType;




typedef struct {
  VecScatterType type;
  PetscInt n;
  PetscInt *vslots;






  PetscBool nonmatching_computed;
  PetscInt n_nonmatching;
  PetscInt *slots_nonmatching;
  PetscBool is_copy;
  PetscInt copy_start;
  PetscInt copy_length;
} VecScatter_Seq_General;

typedef struct {
  VecScatterType type;
  PetscInt n;
  PetscInt first;
  PetscInt step;
} VecScatter_Seq_Stride;




typedef struct {
  VecScatterType type;
  PetscMPIInt *count;
  PetscMPIInt *displx;
  PetscScalar *work1;
  PetscScalar *work2;
} VecScatter_MPI_ToAll;




typedef struct {
  VecScatterType type;
  PetscInt n;
  PetscInt *starts;
  PetscInt *indices;
  PetscMPIInt *procs;
  MPI_Request *requests,*rev_requests;
  PetscScalar *values;
  VecScatter_Seq_General local;
  MPI_Status *sstatus,*rstatus;
  PetscBool use_readyreceiver;
  PetscInt bs;
  PetscBool sendfirst;
  PetscBool contiq;

  PetscBool use_alltoallv;
  PetscMPIInt *counts,*displs;

  PetscBool use_alltoallw;

  PetscMPIInt *wcounts,*wdispls;
  MPI_Datatype *types;

  PetscBool use_window;

  MPI_Win window;
  PetscInt *winstarts;

} VecScatter_MPI_General;

struct _p_VecScatter {
  _p_PetscObject hdr; int *ops;
  PetscInt to_n,from_n;
  PetscBool inuse;
  PetscBool beginandendtogether;

  PetscBool packtogether;
  PetscBool reproduce;
  PetscErrorCode (*begin)(VecScatter,Vec,Vec,InsertMode,ScatterMode);
  PetscErrorCode (*end)(VecScatter,Vec,Vec,InsertMode,ScatterMode);
  PetscErrorCode (*copy)(VecScatter,VecScatter);
  PetscErrorCode (*destroy)(VecScatter);
  PetscErrorCode (*view)(VecScatter,PetscViewer);
  void *fromdata,*todata;
  void *spptr;
};

extern PetscErrorCode VecStashCreate_Private(MPI_Comm,PetscInt,VecStash*);
extern PetscErrorCode VecStashDestroy_Private(VecStash*);
extern PetscErrorCode VecStashExpand_Private(VecStash*,PetscInt);
extern PetscErrorCode VecStashScatterEnd_Private(VecStash*);
extern PetscErrorCode VecStashSetInitialSize_Private(VecStash*,PetscInt);
extern PetscErrorCode VecStashGetInfo_Private(VecStash*,PetscInt*,PetscInt*);
extern PetscErrorCode VecStashScatterBegin_Private(VecStash*,PetscInt*);
extern PetscErrorCode VecStashScatterGetMesg_Private(VecStash*,PetscMPIInt*,PetscInt**,PetscScalar**,PetscInt*);
# 487 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
static inline PetscErrorCode VecStashValue_Private(VecStash *stash,PetscInt row,PetscScalar value)
{
  PetscErrorCode ierr;

  if (((stash)->n + 1) > (stash)->nmax) {
    ierr = VecStashExpand_Private(stash,1);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),492,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  (stash)->idx[(stash)->n] = row;
  (stash)->array[(stash)->n] = value;
  (stash)->n++;
  return 0;
}
# 508 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"
static inline PetscErrorCode VecStashValuesBlocked_Private(VecStash *stash,PetscInt row,PetscScalar *values)
{
  PetscInt jj,stash_bs=(stash)->bs;
  PetscScalar *array;
  PetscErrorCode ierr;
  if (((stash)->n+1) > (stash)->nmax) {
    ierr = VecStashExpand_Private(stash,1);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),514,__func__,"/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  array = (stash)->array + stash_bs*(stash)->n;
  (stash)->idx[(stash)->n] = row;
  for (jj=0; jj<stash_bs; jj++) { array[jj] = values[jj];}
  (stash)->n++;
  return 0;
}

extern PetscErrorCode VecStrideGather_Default(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode VecStrideScatter_Default(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode VecReciprocal_Default(Vec);
# 534 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h"

# 496 "/home/dpnkarthik/petsc-rnet/include/petscvec.h" 2

extern PetscErrorCode VecContourScale(Vec,PetscReal,PetscReal);





typedef enum { VECOP_VIEW = 33, VECOP_LOAD = 41, VECOP_DUPLICATE = 0} VecOperation;
extern PetscErrorCode VecSetOperation(Vec,VecOperation,void(*)(void));





extern PetscErrorCode VecCreateGhost(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);
extern PetscErrorCode VecCreateGhostWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);
extern PetscErrorCode VecCreateGhostBlock(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);
extern PetscErrorCode VecCreateGhostBlockWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);
extern PetscErrorCode VecGhostGetLocalForm(Vec,Vec*);

extern PetscErrorCode VecGhostRestoreLocalForm(Vec,Vec*);
extern PetscErrorCode VecGhostUpdateBegin(Vec,InsertMode,ScatterMode);
extern PetscErrorCode VecGhostUpdateEnd(Vec,InsertMode,ScatterMode);

extern PetscErrorCode VecConjugate(Vec);

extern PetscErrorCode VecScatterCreateToAll(Vec,VecScatter*,Vec*);
extern PetscErrorCode VecScatterCreateToZero(Vec,VecScatter*,Vec*);

extern PetscErrorCode PetscViewerMathematicaGetVector(PetscViewer, Vec);
extern PetscErrorCode PetscViewerMathematicaPutVector(PetscViewer, Vec);
# 543 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
        struct _n_Vecs {PetscInt n; Vec v;};
typedef struct _n_Vecs* Vecs;
# 553 "/home/dpnkarthik/petsc-rnet/include/petscvec.h"
typedef struct _p_PetscCUSPIndices* PetscCUSPIndices;
extern PetscErrorCode PetscCUSPIndicesCreate(PetscInt,const PetscInt*,PetscCUSPIndices*);
extern PetscErrorCode PetscCUSPIndicesDestroy(PetscCUSPIndices*);
extern PetscErrorCode VecCUSPCopyToGPUSome_Public(Vec,PetscCUSPIndices);
extern PetscErrorCode VecCUSPCopyFromGPUSome_Public(Vec,PetscCUSPIndices);

extern PetscErrorCode VecCreateSeqCUSP(MPI_Comm,PetscInt,Vec*);
extern PetscErrorCode VecCreateMPICUSP(MPI_Comm,PetscInt,PetscInt,Vec*);







extern PetscErrorCode VecNestGetSubVecs(Vec,PetscInt*,Vec**);
extern PetscErrorCode VecNestGetSubVec(Vec,PetscInt,Vec*);
extern PetscErrorCode VecCreateNest(MPI_Comm,PetscInt,IS*,Vec*,Vec*);
extern PetscErrorCode VecNestGetSize(Vec,PetscInt*);


# 7 "/home/dpnkarthik/petsc-rnet/include/petscmat.h" 2

# 18 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct _p_Mat* Mat;
# 135 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_FACTOR_NONE, MAT_FACTOR_LU, MAT_FACTOR_CHOLESKY, MAT_FACTOR_ILU, MAT_FACTOR_ICC,MAT_FACTOR_ILUDT} MatFactorType;
extern const char *const MatFactorTypes[];

extern PetscErrorCode MatGetFactor(Mat,const char*,MatFactorType,Mat*);
extern PetscErrorCode MatGetFactorAvailable(Mat,const char*,MatFactorType,PetscBool *);
extern PetscErrorCode MatFactorGetSolverPackage(Mat,const char**);
extern PetscErrorCode MatGetFactorType(Mat,MatFactorType*);



extern PetscClassId MAT_CLASSID;
extern PetscClassId MAT_FDCOLORING_CLASSID;
extern PetscClassId MAT_PARTITIONING_CLASSID;
extern PetscClassId MAT_NULLSPACE_CLASSID;
extern PetscClassId MATMFFD_CLASSID;
# 162 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_INITIAL_MATRIX,MAT_REUSE_MATRIX,MAT_IGNORE_MATRIX} MatReuse;
# 172 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_DO_NOT_GET_VALUES,MAT_GET_VALUES} MatGetSubMatrixOption;

extern PetscErrorCode MatInitializePackage(const char[]);

extern PetscErrorCode MatCreate(MPI_Comm,Mat*);


extern PetscErrorCode MatSetSizes(Mat,PetscInt,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode MatSetType(Mat,const char*);
extern PetscErrorCode MatSetFromOptions(Mat);
extern PetscErrorCode MatSetUpPreallocation(Mat);
extern PetscErrorCode MatRegisterAll(const char[]);
extern PetscErrorCode MatRegister(const char[],const char[],const char[],PetscErrorCode(*)(Mat));
extern PetscErrorCode MatRegisterBaseName(const char[],const char[],const char[]);
extern PetscErrorCode MatSetOptionsPrefix(Mat,const char[]);
extern PetscErrorCode MatAppendOptionsPrefix(Mat,const char[]);
extern PetscErrorCode MatGetOptionsPrefix(Mat,const char*[]);
# 237 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscBool MatRegisterAllCalled;
extern PetscFList MatList;
extern PetscFList MatColoringList;
extern PetscFList MatPartitioningList;
# 251 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {DIFFERENT_NONZERO_PATTERN,SUBSET_NONZERO_PATTERN,SAME_NONZERO_PATTERN,SAME_PRECONDITIONER} MatStructure;

extern PetscErrorCode MatCreateSeqDense(MPI_Comm,PetscInt,PetscInt,PetscScalar[],Mat*);
extern PetscErrorCode MatCreateMPIDense(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);







extern PetscErrorCode MatCreateMPIAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);














extern PetscErrorCode MatCreateMPIAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat *);
extern PetscErrorCode MatCreateMPIAIJWithSplitArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],PetscInt[],PetscInt[],PetscScalar[],Mat*);

extern PetscErrorCode MatCreateSeqBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);







extern PetscErrorCode MatCreateMPIBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);














extern PetscErrorCode MatCreateMPIBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat*);

extern PetscErrorCode MatCreateMPIAdj(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscInt[],Mat*);
extern PetscErrorCode MatCreateSeqSBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);








extern PetscErrorCode MatCreateMPISBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);














extern PetscErrorCode MatCreateMPISBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat *);
extern PetscErrorCode MatMPISBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);

extern PetscErrorCode MatCreateShell(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,void *,Mat*);


extern PetscErrorCode MatCreateAdic(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,void (*)(void),Mat*);
extern PetscErrorCode MatCreateNormal(Mat,Mat*);

extern PetscErrorCode MatCreateLRC(Mat,Mat,Mat,Mat*);
extern PetscErrorCode MatCreateIS(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,ISLocalToGlobalMapping,Mat*);
extern PetscErrorCode MatCreateSeqAIJCRL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIAIJCRL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateSeqBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateSeqSBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPISBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateScatter(MPI_Comm,VecScatter,Mat*);
extern PetscErrorCode MatScatterSetVecScatter(Mat,VecScatter);
extern PetscErrorCode MatScatterGetVecScatter(Mat,VecScatter*);
extern PetscErrorCode MatCreateBlockMat(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,Mat*);
extern PetscErrorCode MatCompositeAddMat(Mat,Mat);
extern PetscErrorCode MatCompositeMerge(Mat);
extern PetscErrorCode MatCreateComposite(MPI_Comm,PetscInt,const Mat*,Mat*);
typedef enum {MAT_COMPOSITE_ADDITIVE,MAT_COMPOSITE_MULTIPLICATIVE} MatCompositeType;
extern PetscErrorCode MatCompositeSetType(Mat,MatCompositeType);

extern PetscErrorCode MatCreateFFT(MPI_Comm,PetscInt,const PetscInt[],const char*,Mat*);
extern PetscErrorCode MatCreateSeqCUFFT(MPI_Comm,PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateTranspose(Mat,Mat*);
extern PetscErrorCode MatCreateSubMatrix(Mat,IS,IS,Mat*);
extern PetscErrorCode MatSubMatrixUpdate(Mat,Mat,IS,IS);
extern PetscErrorCode MatCreateLocalRef(Mat,IS,IS,Mat*);

extern PetscErrorCode MatPythonSetType(Mat,const char[]);

extern PetscErrorCode MatSetUp(Mat);
extern PetscErrorCode MatDestroy(Mat*);

extern PetscErrorCode MatConjugate(Mat);
extern PetscErrorCode MatRealPart(Mat);
extern PetscErrorCode MatImaginaryPart(Mat);
extern PetscErrorCode MatGetDiagonalBlock(Mat,Mat*);
extern PetscErrorCode MatGetTrace(Mat,PetscScalar*);
extern PetscErrorCode MatInvertBlockDiagonal(Mat,PetscScalar **);


extern PetscErrorCode MatSetValues(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesBlocked(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesRow(Mat,PetscInt,const PetscScalar[]);
extern PetscErrorCode MatSetValuesRowLocal(Mat,PetscInt,const PetscScalar[]);
extern PetscErrorCode MatSetValuesBatch(Mat,PetscInt,PetscInt,PetscInt[],const PetscScalar[]);
# 397 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct {
  PetscInt k,j,i,c;
} MatStencil;

extern PetscErrorCode MatSetValuesStencil(Mat,PetscInt,const MatStencil[],PetscInt,const MatStencil[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesBlockedStencil(Mat,PetscInt,const MatStencil[],PetscInt,const MatStencil[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetStencil(Mat,PetscInt,const PetscInt[],const PetscInt[],PetscInt);

extern PetscErrorCode MatSetColoring(Mat,ISColoring);
extern PetscErrorCode MatSetValuesAdic(Mat,void*);
extern PetscErrorCode MatSetValuesAdifor(Mat,PetscInt,void*);
# 417 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_FLUSH_ASSEMBLY=1,MAT_FINAL_ASSEMBLY=0} MatAssemblyType;
extern PetscErrorCode MatAssemblyBegin(Mat,MatAssemblyType);
extern PetscErrorCode MatAssemblyEnd(Mat,MatAssemblyType);
extern PetscErrorCode MatAssembled(Mat,PetscBool *);
# 433 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_ROW_ORIENTED,MAT_NEW_NONZERO_LOCATIONS,
              MAT_SYMMETRIC,
              MAT_STRUCTURALLY_SYMMETRIC,
              MAT_NEW_DIAGONALS,MAT_IGNORE_OFF_PROC_ENTRIES,
              MAT_NEW_NONZERO_LOCATION_ERR,
              MAT_NEW_NONZERO_ALLOCATION_ERR,MAT_USE_HASH_TABLE,
              MAT_KEEP_NONZERO_PATTERN,MAT_IGNORE_ZERO_ENTRIES,
              MAT_USE_INODES,
              MAT_HERMITIAN,
              MAT_SYMMETRY_ETERNAL,
              MAT_CHECK_COMPRESSED_ROW,
              MAT_IGNORE_LOWER_TRIANGULAR,MAT_ERROR_LOWER_TRIANGULAR,
              MAT_GETROW_UPPERTRIANGULAR,MAT_UNUSED_NONZERO_LOCATION_ERR,
              MAT_SPD,MAT_NO_OFF_PROC_ENTRIES,MAT_NO_OFF_PROC_ZERO_ROWS,
              NUM_MAT_OPTIONS} MatOption;
extern const char *MatOptions[];
extern PetscErrorCode MatSetOption(Mat,MatOption,PetscBool );
extern PetscErrorCode MatGetType(Mat,const char**);


extern PetscErrorCode MatGetValues(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscScalar[]);
extern PetscErrorCode MatGetRow(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatRestoreRow(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatGetRowUpperTriangular(Mat);
extern PetscErrorCode MatRestoreRowUpperTriangular(Mat);
extern PetscErrorCode MatGetColumn(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatRestoreColumn(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatGetColumnVector(Mat,Vec,PetscInt);
extern PetscErrorCode MatGetArray(Mat,PetscScalar *[]);

extern PetscErrorCode MatRestoreArray(Mat,PetscScalar *[]);
extern PetscErrorCode MatGetBlockSize(Mat,PetscInt *);

extern PetscErrorCode MatSetBlockSize(Mat,PetscInt);


extern PetscErrorCode MatMult(Mat,Vec,Vec);
extern PetscErrorCode MatMultDiagonalBlock(Mat,Vec,Vec);
extern PetscErrorCode MatMultAdd(Mat,Vec,Vec,Vec);

extern PetscErrorCode MatMultTranspose(Mat,Vec,Vec);
extern PetscErrorCode MatMultHermitianTranspose(Mat,Vec,Vec);
extern PetscErrorCode MatIsTranspose(Mat,Mat,PetscReal,PetscBool *);


extern PetscErrorCode MatIsHermitianTranspose(Mat,Mat,PetscReal,PetscBool *);
extern PetscErrorCode MatMultTransposeAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultHermitianTransposeAdd(Mat,Vec,Vec,Vec);

extern PetscErrorCode MatMultConstrained(Mat,Vec,Vec);
extern PetscErrorCode MatMultTransposeConstrained(Mat,Vec,Vec);
extern PetscErrorCode MatMatSolve(Mat,Mat,Mat);
# 500 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_DO_NOT_COPY_VALUES,MAT_COPY_VALUES,MAT_SHARE_NONZERO_PATTERN} MatDuplicateOption;

extern PetscErrorCode MatConvert(Mat,const char*,MatReuse,Mat*);

extern PetscErrorCode MatDuplicate(Mat,MatDuplicateOption,Mat*);




extern PetscErrorCode MatCopy(Mat,Mat,MatStructure);
extern PetscErrorCode MatView(Mat,PetscViewer);
extern PetscErrorCode MatIsSymmetric(Mat,PetscReal,PetscBool *);


extern PetscErrorCode MatIsStructurallySymmetric(Mat,PetscBool *);

extern PetscErrorCode MatIsHermitian(Mat,PetscReal,PetscBool *);

extern PetscErrorCode MatIsSymmetricKnown(Mat,PetscBool *,PetscBool *);
extern PetscErrorCode MatIsHermitianKnown(Mat,PetscBool *,PetscBool *);
extern PetscErrorCode MatMissingDiagonal(Mat,PetscBool *,PetscInt *);
extern PetscErrorCode MatLoad(Mat, PetscViewer);

extern PetscErrorCode MatGetRowIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatRestoreRowIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatGetColumnIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatRestoreColumnIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,PetscInt *[],PetscInt *[],PetscBool *);
# 539 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct {
  PetscLogDouble block_size;
  PetscLogDouble nz_allocated,nz_used,nz_unneeded;
  PetscLogDouble memory;
  PetscLogDouble assemblies;
  PetscLogDouble mallocs;
  PetscLogDouble fill_ratio_given,fill_ratio_needed;
  PetscLogDouble factor_mallocs;
} MatInfo;
# 559 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MAT_LOCAL=1,MAT_GLOBAL_MAX=2,MAT_GLOBAL_SUM=3} MatInfoType;
extern PetscErrorCode MatGetInfo(Mat,MatInfoType,MatInfo*);
extern PetscErrorCode MatGetDiagonal(Mat,Vec);
extern PetscErrorCode MatGetRowMax(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowMin(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowMaxAbs(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowMinAbs(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowSum(Mat,Vec);
extern PetscErrorCode MatTranspose(Mat,MatReuse,Mat*);

extern PetscErrorCode MatHermitianTranspose(Mat,MatReuse,Mat*);
extern PetscErrorCode MatPermute(Mat,IS,IS,Mat *);

extern PetscErrorCode MatDiagonalScale(Mat,Vec,Vec);
extern PetscErrorCode MatDiagonalSet(Mat,Vec,InsertMode);
extern PetscErrorCode MatEqual(Mat,Mat,PetscBool *);

extern PetscErrorCode MatMultEqual(Mat,Mat,PetscInt,PetscBool *);
extern PetscErrorCode MatMultAddEqual(Mat,Mat,PetscInt,PetscBool *);
extern PetscErrorCode MatMultTransposeEqual(Mat,Mat,PetscInt,PetscBool *);
extern PetscErrorCode MatMultTransposeAddEqual(Mat,Mat,PetscInt,PetscBool *);

extern PetscErrorCode MatNorm(Mat,NormType,PetscReal *);

extern PetscErrorCode MatGetColumnNorms(Mat,NormType,PetscReal *);
extern PetscErrorCode MatZeroEntries(Mat);
extern PetscErrorCode MatZeroRows(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsIS(Mat,IS,PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsStencil(Mat,PetscInt,const MatStencil [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsStencil(Mat,PetscInt,const MatStencil[],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumns(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsIS(Mat,IS,PetscScalar,Vec,Vec);

extern PetscErrorCode MatUseScaledForm(Mat,PetscBool );
extern PetscErrorCode MatScaleSystem(Mat,Vec,Vec);
extern PetscErrorCode MatUnScaleSystem(Mat,Vec,Vec);

extern PetscErrorCode MatGetSize(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetLocalSize(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetOwnershipRange(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetOwnershipRanges(Mat,const PetscInt**);
extern PetscErrorCode MatGetOwnershipRangeColumn(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetOwnershipRangesColumn(Mat,const PetscInt**);

extern PetscErrorCode MatGetSubMatrices(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
extern PetscErrorCode MatGetSubMatricesParallel(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
extern PetscErrorCode MatDestroyMatrices(PetscInt,Mat *[]);
extern PetscErrorCode MatGetSubMatrix(Mat,IS,IS,MatReuse,Mat *);
extern PetscErrorCode MatGetLocalSubMatrix(Mat,IS,IS,Mat*);
extern PetscErrorCode MatRestoreLocalSubMatrix(Mat,IS,IS,Mat*);
extern PetscErrorCode MatGetSeqNonzeroStructure(Mat,Mat*);
extern PetscErrorCode MatDestroySeqNonzeroStructure(Mat*);

extern PetscErrorCode MatMerge(MPI_Comm,Mat,PetscInt,MatReuse,Mat*);
extern PetscErrorCode MatMerge_SeqsToMPI(MPI_Comm,Mat,PetscInt,PetscInt,MatReuse,Mat*);
extern PetscErrorCode MatMerge_SeqsToMPISymbolic(MPI_Comm,Mat,PetscInt,PetscInt,Mat*);
extern PetscErrorCode MatMerge_SeqsToMPINumeric(Mat,Mat);
extern PetscErrorCode MatMPIAIJGetLocalMat(Mat,MatReuse,Mat*);
extern PetscErrorCode MatMPIAIJGetLocalMatCondensed(Mat,MatReuse,IS*,IS*,Mat*);
extern PetscErrorCode MatGetBrowsOfAcols(Mat,Mat,MatReuse,IS*,IS*,PetscInt*,Mat*);
extern PetscErrorCode MatGetBrowsOfAoCols(Mat,Mat,MatReuse,PetscInt**,PetscInt**,MatScalar**,Mat*);

extern PetscErrorCode MatGetCommunicationStructs(Mat, Vec *, PetscTable *, VecScatter *);



extern PetscErrorCode MatGetGhosts(Mat, PetscInt *,const PetscInt *[]);

extern PetscErrorCode MatIncreaseOverlap(Mat,PetscInt,IS[],PetscInt);

extern PetscErrorCode MatMatMult(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatMatMultSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatMultNumeric(Mat,Mat,Mat);

extern PetscErrorCode MatPtAP(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatPtAPSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatPtAPNumeric(Mat,Mat,Mat);

extern PetscErrorCode MatMatMultTranspose(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatMatMultTransposeSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatMultTransposeNumeric(Mat,Mat,Mat);

extern PetscErrorCode MatAXPY(Mat,PetscScalar,Mat,MatStructure);
extern PetscErrorCode MatAYPX(Mat,PetscScalar,Mat,MatStructure);

extern PetscErrorCode MatScale(Mat,PetscScalar);
extern PetscErrorCode MatShift(Mat,PetscScalar);

extern PetscErrorCode MatSetLocalToGlobalMapping(Mat,ISLocalToGlobalMapping,ISLocalToGlobalMapping);
extern PetscErrorCode MatSetLocalToGlobalMappingBlock(Mat,ISLocalToGlobalMapping,ISLocalToGlobalMapping);
extern PetscErrorCode MatGetLocalToGlobalMapping(Mat,ISLocalToGlobalMapping*,ISLocalToGlobalMapping*);
extern PetscErrorCode MatGetLocalToGlobalMappingBlock(Mat,ISLocalToGlobalMapping*,ISLocalToGlobalMapping*);
extern PetscErrorCode MatZeroRowsLocal(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsLocalIS(Mat,IS,PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsLocal(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsLocalIS(Mat,IS,PetscScalar,Vec,Vec);
extern PetscErrorCode MatSetValuesLocal(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesBlockedLocal(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

extern PetscErrorCode MatStashSetInitialSize(Mat,PetscInt,PetscInt);
extern PetscErrorCode MatStashGetInfo(Mat,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

extern PetscErrorCode MatInterpolate(Mat,Vec,Vec);
extern PetscErrorCode MatInterpolateAdd(Mat,Vec,Vec,Vec);

extern PetscErrorCode MatRestrict(Mat,Vec,Vec);
extern PetscErrorCode MatGetVecs(Mat,Vec*,Vec*);
extern PetscErrorCode MatGetRedundantMatrix(Mat,PetscInt,MPI_Comm,PetscInt,MatReuse,Mat*);
extern PetscErrorCode MatGetMultiProcBlock(Mat,MPI_Comm,Mat*);
extern PetscErrorCode MatFindZeroDiagonals(Mat,IS*);
# 690 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
static inline PetscErrorCode MatSetValue(Mat v,PetscInt i,PetscInt j,PetscScalar va,InsertMode mode) {return MatSetValues(v,1,&i,1,&j,&va,mode);}

static inline PetscErrorCode MatGetValue(Mat v,PetscInt i,PetscInt j,PetscScalar *va) {return MatGetValues(v,1,&i,1,&j,va);}

static inline PetscErrorCode MatSetValueLocal(Mat v,PetscInt i,PetscInt j,PetscScalar va,InsertMode mode) {return MatSetValuesLocal(v,1,&i,1,&j,&va,mode);}
# 1015 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatShellGetContext(Mat,void *);


extern PetscErrorCode MatInodeAdjustForInodes(Mat,IS*,IS*);
extern PetscErrorCode MatInodeGetInodeSizes(Mat,PetscInt *,PetscInt *[],PetscInt *);

extern PetscErrorCode MatSeqAIJSetColumnIndices(Mat,PetscInt[]);
extern PetscErrorCode MatSeqBAIJSetColumnIndices(Mat,PetscInt[]);
extern PetscErrorCode MatCreateSeqAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqSBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqAIJFromTriple(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*,PetscInt,PetscBool);



extern PetscErrorCode MatSeqBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[]);

extern PetscErrorCode MatSeqSBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[]);

extern PetscErrorCode MatSeqAIJSetPreallocation(Mat,PetscInt,const PetscInt[]);


extern PetscErrorCode MatMPIBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);

extern PetscErrorCode MatMPISBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
extern PetscErrorCode MatMPIAIJSetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
extern PetscErrorCode MatSeqAIJSetPreallocationCSR(Mat,const PetscInt [],const PetscInt [],const PetscScalar []);
extern PetscErrorCode MatSeqBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPIAIJSetPreallocationCSR(Mat,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPIBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPIAdjSetPreallocation(Mat,PetscInt[],PetscInt[],PetscInt[]);
extern PetscErrorCode MatMPIDenseSetPreallocation(Mat,PetscScalar[]);
extern PetscErrorCode MatSeqDenseSetPreallocation(Mat,PetscScalar[]);
extern PetscErrorCode MatMPIAIJGetSeqAIJ(Mat,Mat*,Mat*,PetscInt*[]);
extern PetscErrorCode MatMPIBAIJGetSeqBAIJ(Mat,Mat*,Mat*,PetscInt*[]);
extern PetscErrorCode MatAdicSetLocalFunction(Mat,void (*)(void));

extern PetscErrorCode MatSeqDenseSetLDA(Mat,PetscInt);
extern PetscErrorCode MatDenseGetLocalMatrix(Mat,Mat*);

extern PetscErrorCode MatStoreValues(Mat);
extern PetscErrorCode MatRetrieveValues(Mat);

extern PetscErrorCode MatDAADSetCtx(Mat,void*);

extern PetscErrorCode MatFindNonzeroRows(Mat,IS*);
# 1086 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatGetOrdering(Mat,const char*,IS*,IS*);
extern PetscErrorCode MatGetOrderingList(PetscFList *list);
extern PetscErrorCode MatOrderingRegister(const char[],const char[],const char[],PetscErrorCode(*)(Mat,const char*,IS*,IS*));
# 1132 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatOrderingRegisterDestroy(void);
extern PetscErrorCode MatOrderingRegisterAll(const char[]);
extern PetscBool MatOrderingRegisterAllCalled;
extern PetscFList MatOrderingList;

extern PetscErrorCode MatReorderForNonzeroDiagonal(Mat,PetscReal,IS,IS);







typedef enum {MAT_SHIFT_NONE,MAT_SHIFT_NONZERO,MAT_SHIFT_POSITIVE_DEFINITE,MAT_SHIFT_INBLOCKS} MatFactorShiftType;
extern const char *MatFactorShiftTypes[];
# 1164 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct {
  PetscReal diagonal_fill;
  PetscReal usedt;
  PetscReal dt;
  PetscReal dtcol;
  PetscReal dtcount;
  PetscReal fill;
  PetscReal levels;
  PetscReal pivotinblocks;

  PetscReal zeropivot;
  PetscReal shifttype;
  PetscReal shiftamount;
} MatFactorInfo;

extern PetscErrorCode MatFactorInfoInitialize(MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactor(Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorSymbolic(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorNumeric(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactor(Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactor(Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorSymbolic(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactorSymbolic(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatICCFactorSymbolic(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatICCFactor(Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatGetInertia(Mat,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode MatSolve(Mat,Vec,Vec);
extern PetscErrorCode MatForwardSolve(Mat,Vec,Vec);
extern PetscErrorCode MatBackwardSolve(Mat,Vec,Vec);
extern PetscErrorCode MatSolveAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolveTranspose(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTransposeAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolves(Mat,Vecs,Vecs);

extern PetscErrorCode MatSetUnfactored(Mat);
# 1214 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {SOR_FORWARD_SWEEP=1,SOR_BACKWARD_SWEEP=2,SOR_SYMMETRIC_SWEEP=3,
              SOR_LOCAL_FORWARD_SWEEP=4,SOR_LOCAL_BACKWARD_SWEEP=8,
              SOR_LOCAL_SYMMETRIC_SWEEP=12,SOR_ZERO_INITIAL_GUESS=16,
              SOR_EISENSTAT=32,SOR_APPLY_UPPER=64,SOR_APPLY_LOWER=128} MatSORType;
extern PetscErrorCode MatSOR(Mat,Vec,PetscReal,MatSORType,PetscReal,PetscInt,PetscInt,Vec);
# 1239 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatGetColoring(Mat,const char*,ISColoring*);
extern PetscErrorCode MatColoringRegister(const char[],const char[],const char[],PetscErrorCode(*)(Mat,char*,ISColoring *));
# 1285 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscBool MatColoringRegisterAllCalled;

extern PetscErrorCode MatColoringRegisterAll(const char[]);
extern PetscErrorCode MatColoringRegisterDestroy(void);
extern PetscErrorCode MatColoringPatch(Mat,PetscInt,PetscInt,ISColoringValue[],ISColoring*);
# 1301 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct _p_MatFDColoring* MatFDColoring;

extern PetscErrorCode MatFDColoringCreate(Mat,ISColoring,MatFDColoring *);
extern PetscErrorCode MatFDColoringDestroy(MatFDColoring*);
extern PetscErrorCode MatFDColoringView(MatFDColoring,PetscViewer);
extern PetscErrorCode MatFDColoringSetFunction(MatFDColoring,PetscErrorCode (*)(void),void*);
extern PetscErrorCode MatFDColoringGetFunction(MatFDColoring,PetscErrorCode (**)(void),void**);
extern PetscErrorCode MatFDColoringSetParameters(MatFDColoring,PetscReal,PetscReal);
extern PetscErrorCode MatFDColoringSetFromOptions(MatFDColoring);
extern PetscErrorCode MatFDColoringApply(Mat,MatFDColoring,Vec,MatStructure*,void *);
extern PetscErrorCode MatFDColoringApplyTS(Mat,MatFDColoring,PetscReal,Vec,MatStructure*,void *);
extern PetscErrorCode MatFDColoringSetF(MatFDColoring,Vec);
extern PetscErrorCode MatFDColoringGetPerturbedColumns(MatFDColoring,PetscInt*,PetscInt*[]);
# 1328 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct _p_MatPartitioning* MatPartitioning;
# 1348 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatPartitioningCreate(MPI_Comm,MatPartitioning*);
extern PetscErrorCode MatPartitioningSetType(MatPartitioning,const char*);
extern PetscErrorCode MatPartitioningSetNParts(MatPartitioning,PetscInt);
extern PetscErrorCode MatPartitioningSetAdjacency(MatPartitioning,Mat);
extern PetscErrorCode MatPartitioningSetVertexWeights(MatPartitioning,const PetscInt[]);
extern PetscErrorCode MatPartitioningSetPartitionWeights(MatPartitioning,const PetscReal []);
extern PetscErrorCode MatPartitioningApply(MatPartitioning,IS*);
extern PetscErrorCode MatPartitioningDestroy(MatPartitioning*);

extern PetscErrorCode MatPartitioningRegister(const char[],const char[],const char[],PetscErrorCode (*)(MatPartitioning));
# 1402 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscBool MatPartitioningRegisterAllCalled;

extern PetscErrorCode MatPartitioningRegisterAll(const char[]);
extern PetscErrorCode MatPartitioningRegisterDestroy(void);

extern PetscErrorCode MatPartitioningView(MatPartitioning,PetscViewer);
extern PetscErrorCode MatPartitioningSetFromOptions(MatPartitioning);
extern PetscErrorCode MatPartitioningGetType(MatPartitioning,const char**);

extern PetscErrorCode MatPartitioningParmetisSetCoarseSequential(MatPartitioning);
extern PetscErrorCode MatPartitioningParmetisGetEdgeCut(MatPartitioning, PetscInt *);

typedef enum { MP_CHACO_MULTILEVEL=1,MP_CHACO_SPECTRAL=2,MP_CHACO_LINEAR=4,MP_CHACO_RANDOM=5,MP_CHACO_SCATTERED=6 } MPChacoGlobalType;
extern const char *MPChacoGlobalTypes[];
typedef enum { MP_CHACO_KERNIGHAN=1,MP_CHACO_NONE=2 } MPChacoLocalType;
extern const char *MPChacoLocalTypes[];
typedef enum { MP_CHACO_LANCZOS=0,MP_CHACO_RQI=1 } MPChacoEigenType;
extern const char *MPChacoEigenTypes[];

extern PetscErrorCode MatPartitioningChacoSetGlobal(MatPartitioning,MPChacoGlobalType);
extern PetscErrorCode MatPartitioningChacoGetGlobal(MatPartitioning,MPChacoGlobalType*);
extern PetscErrorCode MatPartitioningChacoSetLocal(MatPartitioning,MPChacoLocalType);
extern PetscErrorCode MatPartitioningChacoGetLocal(MatPartitioning,MPChacoLocalType*);
extern PetscErrorCode MatPartitioningChacoSetCoarseLevel(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningChacoSetEigenSolver(MatPartitioning,MPChacoEigenType);
extern PetscErrorCode MatPartitioningChacoGetEigenSolver(MatPartitioning,MPChacoEigenType*);
extern PetscErrorCode MatPartitioningChacoSetEigenTol(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningChacoGetEigenTol(MatPartitioning,PetscReal*);
extern PetscErrorCode MatPartitioningChacoSetEigenNumber(MatPartitioning,PetscInt);
extern PetscErrorCode MatPartitioningChacoGetEigenNumber(MatPartitioning,PetscInt*);
# 1441 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatPartitioningPartySetGlobal(MatPartitioning,const char*);



extern PetscErrorCode MatPartitioningPartySetLocal(MatPartitioning,const char*);
extern PetscErrorCode MatPartitioningPartySetCoarseLevel(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningPartySetBipart(MatPartitioning,PetscBool);
extern PetscErrorCode MatPartitioningPartySetMatchOptimization(MatPartitioning,PetscBool);

typedef enum { MP_PTSCOTCH_QUALITY,MP_PTSCOTCH_SPEED,MP_PTSCOTCH_BALANCE,MP_PTSCOTCH_SAFETY,MP_PTSCOTCH_SCALABILITY } MPPTScotchStrategyType;
extern const char *MPPTScotchStrategyTypes[];

extern PetscErrorCode MatPartitioningPTScotchSetImbalance(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningPTScotchGetImbalance(MatPartitioning,PetscReal*);
extern PetscErrorCode MatPartitioningPTScotchSetStrategy(MatPartitioning,MPPTScotchStrategyType);
extern PetscErrorCode MatPartitioningPTScotchGetStrategy(MatPartitioning,MPPTScotchStrategyType*);

extern PetscErrorCode MatMeshToVertexGraph(Mat,PetscInt,Mat*);
extern PetscErrorCode MatMeshToCellGraph(Mat,PetscInt,Mat*);




typedef enum { MATOP_SET_VALUES=0,
               MATOP_GET_ROW=1,
               MATOP_RESTORE_ROW=2,
               MATOP_MULT=3,
               MATOP_MULT_ADD=4,
               MATOP_MULT_TRANSPOSE=5,
               MATOP_MULT_TRANSPOSE_ADD=6,
               MATOP_SOLVE=7,
               MATOP_SOLVE_ADD=8,
               MATOP_SOLVE_TRANSPOSE=9,
               MATOP_SOLVE_TRANSPOSE_ADD=10,
               MATOP_LUFACTOR=11,
               MATOP_CHOLESKYFACTOR=12,
               MATOP_SOR=13,
               MATOP_TRANSPOSE=14,
               MATOP_GETINFO=15,
               MATOP_EQUAL=16,
               MATOP_GET_DIAGONAL=17,
               MATOP_DIAGONAL_SCALE=18,
               MATOP_NORM=19,
               MATOP_ASSEMBLY_BEGIN=20,
               MATOP_ASSEMBLY_END=21,
               MATOP_SET_OPTION=22,
               MATOP_ZERO_ENTRIES=23,
               MATOP_ZERO_ROWS=24,
               MATOP_LUFACTOR_SYMBOLIC=25,
               MATOP_LUFACTOR_NUMERIC=26,
               MATOP_CHOLESKY_FACTOR_SYMBOLIC=27,
               MATOP_CHOLESKY_FACTOR_NUMERIC=28,
               MATOP_SETUP_PREALLOCATION=29,
               MATOP_ILUFACTOR_SYMBOLIC=30,
               MATOP_ICCFACTOR_SYMBOLIC=31,
               MATOP_GET_ARRAY=32,
               MATOP_RESTORE_ARRAY=33,
               MATOP_DUPLICATE=34,
               MATOP_FORWARD_SOLVE=35,
               MATOP_BACKWARD_SOLVE=36,
               MATOP_ILUFACTOR=37,
               MATOP_ICCFACTOR=38,
               MATOP_AXPY=39,
               MATOP_GET_SUBMATRICES=40,
               MATOP_INCREASE_OVERLAP=41,
               MATOP_GET_VALUES=42,
               MATOP_COPY=43,
               MATOP_GET_ROW_MAX=44,
               MATOP_SCALE=45,
               MATOP_SHIFT=46,
               MATOP_DIAGONAL_SET=47,
               MATOP_ILUDT_FACTOR=48,
               MATOP_SET_BLOCK_SIZE=49,
               MATOP_GET_ROW_IJ=50,
               MATOP_RESTORE_ROW_IJ=51,
               MATOP_GET_COLUMN_IJ=52,
               MATOP_RESTORE_COLUMN_IJ=53,
               MATOP_FDCOLORING_CREATE=54,
               MATOP_COLORING_PATCH=55,
               MATOP_SET_UNFACTORED=56,
               MATOP_PERMUTE=57,
               MATOP_SET_VALUES_BLOCKED=58,
               MATOP_GET_SUBMATRIX=59,
               MATOP_DESTROY=60,
               MATOP_VIEW=61,
               MATOP_CONVERT_FROM=62,
               MATOP_USE_SCALED_FORM=63,
               MATOP_SCALE_SYSTEM=64,
               MATOP_UNSCALE_SYSTEM=65,
               MATOP_SET_LOCAL_TO_GLOBAL_MAP=66,
               MATOP_SET_VALUES_LOCAL=67,
               MATOP_ZERO_ROWS_LOCAL=68,
               MATOP_GET_ROW_MAX_ABS=69,
               MATOP_GET_ROW_MIN_ABS=70,
               MATOP_CONVERT=71,
               MATOP_SET_COLORING=72,
               MATOP_SET_VALUES_ADIC=73,
               MATOP_SET_VALUES_ADIFOR=74,
               MATOP_FD_COLORING_APPLY=75,
               MATOP_SET_FROM_OPTIONS=76,
               MATOP_MULT_CON=77,
               MATOP_MULT_TRANSPOSE_CON=78,
               MATOP_PERMUTE_SPARSIFY=79,
               MATOP_MULT_MULTIPLE=80,
               MATOP_SOLVE_MULTIPLE=81,
               MATOP_GET_INERTIA=82,
               MATOP_LOAD=83,
               MATOP_IS_SYMMETRIC=84,
               MATOP_IS_HERMITIAN=85,
               MATOP_IS_STRUCTURALLY_SYMMETRIC=86,
               MATOP_DUMMY=87,
               MATOP_GET_VECS=88,
               MATOP_MAT_MULT=89,
               MATOP_MAT_MULT_SYMBOLIC=90,
               MATOP_MAT_MULT_NUMERIC=91,
               MATOP_PTAP=92,
               MATOP_PTAP_SYMBOLIC=93,
               MATOP_PTAP_NUMERIC=94,
               MATOP_MAT_MULTTRANSPOSE=95,
               MATOP_MAT_MULTTRANSPOSE_SYM=96,
               MATOP_MAT_MULTTRANSPOSE_NUM=97,
               MATOP_PTAP_SYMBOLIC_SEQAIJ=98,
               MATOP_PTAP_NUMERIC_SEQAIJ=99,
               MATOP_PTAP_SYMBOLIC_MPIAIJ=100,
               MATOP_PTAP_NUMERIC_MPIAIJ=101,
               MATOP_CONJUGATE=102,
               MATOP_SET_SIZES=103,
               MATOP_SET_VALUES_ROW=104,
               MATOP_REAL_PART=105,
               MATOP_IMAG_PART=106,
               MATOP_GET_ROW_UTRIANGULAR=107,
               MATOP_RESTORE_ROW_UTRIANGULAR=108,
               MATOP_MATSOLVE=109,
               MATOP_GET_REDUNDANTMATRIX=110,
               MATOP_GET_ROW_MIN=111,
               MATOP_GET_COLUMN_VEC=112,
               MATOP_MISSING_DIAGONAL=113,
               MATOP_MATGETSEQNONZEROSTRUCTURE=114,
               MATOP_CREATE=115,
               MATOP_GET_GHOSTS=116,
               MATOP_GET_LOCALSUBMATRIX=117,
               MATOP_RESTORE_LOCALSUBMATRIX=118,
               MATOP_MULT_DIAGONAL_BLOCK=119,
               MATOP_HERMITIANTRANSPOSE=120,
               MATOP_MULTHERMITIANTRANSPOSE=121,
               MATOP_MULTHERMITIANTRANSPOSEADD=122,
               MATOP_GETMULTIPROCBLOCK=123,
               MATOP_GETCOLUMNNORMS=125,

        MATOP_GET_SUBMATRICES_PARALLEL=128,
               MATOP_SET_VALUES_BATCH=129,


        MATOP_SET_STENCIL=130,

               MATOP_SET_GRID=131

             } MatOperation;
extern PetscErrorCode MatHasOperation(Mat,MatOperation,PetscBool *);
extern PetscErrorCode MatShellSetOperation(Mat,MatOperation,void(*)(void));
extern PetscErrorCode MatShellGetOperation(Mat,MatOperation,void(**)(void));
extern PetscErrorCode MatShellSetContext(Mat,void*);
# 1614 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatMPIBAIJSetHashTableFactor(Mat,PetscReal);
extern PetscErrorCode MatISGetLocalMat(Mat,Mat*);
# 1630 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct _p_MatNullSpace* MatNullSpace;

extern PetscErrorCode MatNullSpaceCreate(MPI_Comm,PetscBool ,PetscInt,const Vec[],MatNullSpace*);
extern PetscErrorCode MatNullSpaceSetFunction(MatNullSpace,PetscErrorCode (*)(MatNullSpace,Vec,void*),void*);
extern PetscErrorCode MatNullSpaceDestroy(MatNullSpace*);
extern PetscErrorCode MatNullSpaceRemove(MatNullSpace,Vec,Vec*);
extern PetscErrorCode MatSetNullSpace(Mat,MatNullSpace);
extern PetscErrorCode MatSetNearNullSpace(Mat,MatNullSpace);
extern PetscErrorCode MatNullSpaceTest(MatNullSpace,Mat,PetscBool *);
extern PetscErrorCode MatNullSpaceView(MatNullSpace,PetscViewer);

extern PetscErrorCode MatReorderingSeqSBAIJ(Mat,IS);
extern PetscErrorCode MatMPISBAIJSetHashTableFactor(Mat,PetscReal);
extern PetscErrorCode MatSeqSBAIJSetColumnIndices(Mat,PetscInt *);
extern PetscErrorCode MatSeqBAIJInvertBlockDiagonal(Mat);

extern PetscErrorCode MatCreateMAIJ(Mat,PetscInt,Mat*);
extern PetscErrorCode MatMAIJRedimension(Mat,PetscInt,Mat*);
extern PetscErrorCode MatMAIJGetAIJ(Mat,Mat*);

extern PetscErrorCode MatComputeExplicitOperator(Mat,Mat*);

extern PetscErrorCode MatDiagonalScaleLocal(Mat,Vec);

extern PetscErrorCode MatCreateMFFD(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Mat*);
extern PetscErrorCode MatMFFDSetBase(Mat,Vec,Vec);
extern PetscErrorCode MatMFFDSetFunction(Mat,PetscErrorCode(*)(void*,Vec,Vec),void*);
extern PetscErrorCode MatMFFDSetFunctioni(Mat,PetscErrorCode (*)(void*,PetscInt,Vec,PetscScalar*));
extern PetscErrorCode MatMFFDSetFunctioniBase(Mat,PetscErrorCode (*)(void*,Vec));
extern PetscErrorCode MatMFFDAddNullSpace(Mat,MatNullSpace);
extern PetscErrorCode MatMFFDSetHHistory(Mat,PetscScalar[],PetscInt);
extern PetscErrorCode MatMFFDResetHHistory(Mat);
extern PetscErrorCode MatMFFDSetFunctionError(Mat,PetscReal);
extern PetscErrorCode MatMFFDSetPeriod(Mat,PetscInt);
extern PetscErrorCode MatMFFDGetH(Mat,PetscScalar *);
extern PetscErrorCode MatMFFDSetOptionsPrefix(Mat,const char[]);
extern PetscErrorCode MatMFFDCheckPositivity(void*,Vec,Vec,PetscScalar*);
extern PetscErrorCode MatMFFDSetCheckh(Mat,PetscErrorCode (*)(void*,Vec,Vec,PetscScalar*),void*);
# 1681 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef struct _p_MatMFFD* MatMFFD;
# 1694 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatMFFDSetType(Mat,const char*);
extern PetscErrorCode MatMFFDRegister(const char[],const char[],const char[],PetscErrorCode (*)(MatMFFD));
# 1740 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatMFFDRegisterAll(const char[]);
extern PetscErrorCode MatMFFDRegisterDestroy(void);
extern PetscErrorCode MatMFFDDSSetUmin(Mat,PetscReal);
extern PetscErrorCode MatMFFDWPSetComputeNormU(Mat,PetscBool );


extern PetscErrorCode PetscViewerMathematicaPutMatrix(PetscViewer, PetscInt, PetscInt, PetscReal *);
extern PetscErrorCode PetscViewerMathematicaPutCSRMatrix(PetscViewer, PetscInt, PetscInt, PetscInt *, PetscInt *, PetscReal *);
# 1764 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatCreateSeqAIJCUSP(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIAIJCUSP(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
# 1777 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
extern PetscErrorCode MatCreateNest(MPI_Comm,PetscInt,const IS[],PetscInt,const IS[],const Mat[],Mat*);
extern PetscErrorCode MatNestGetSize(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatNestGetSubMats(Mat,PetscInt*,PetscInt*,Mat***);
extern PetscErrorCode MatNestGetSubMat(Mat,PetscInt,PetscInt,Mat*);
extern PetscErrorCode MatNestSetVecType(Mat,const char*);
extern PetscErrorCode MatNestSetSubMats(Mat,PetscInt,const IS[],PetscInt,const IS[],const Mat[]);
extern PetscErrorCode MatNestSetSubMat(Mat,PetscInt,PetscInt,Mat);
# 1792 "/home/dpnkarthik/petsc-rnet/include/petscmat.h"
typedef enum {MATIJ_LOCAL, MATIJ_GLOBAL} MatIJIndexType;
extern PetscErrorCode MatIJSetMultivalued(Mat, PetscBool);
extern PetscErrorCode MatIJGetMultivalued(Mat, PetscBool*);
extern PetscErrorCode MatIJSetEdges(Mat, PetscInt, const PetscInt*, const PetscInt*);
extern PetscErrorCode MatIJGetEdges(Mat, PetscInt *, PetscInt **, PetscInt **);
extern PetscErrorCode MatIJSetEdgesIS(Mat, IS, IS);
extern PetscErrorCode MatIJGetEdgesIS(Mat, IS*, IS*);
extern PetscErrorCode MatIJGetRowSizes(Mat, MatIJIndexType, PetscInt, const PetscInt *, PetscInt **);
extern PetscErrorCode MatIJGetMinRowSize(Mat, PetscInt *);
extern PetscErrorCode MatIJGetMaxRowSize(Mat, PetscInt *);
extern PetscErrorCode MatIJGetSupport(Mat, PetscInt *, PetscInt **);
extern PetscErrorCode MatIJGetSupportIS(Mat, IS *);
extern PetscErrorCode MatIJGetImage(Mat, PetscInt*, PetscInt**);
extern PetscErrorCode MatIJGetImageIS(Mat, IS *);
extern PetscErrorCode MatIJGetSupportSize(Mat, PetscInt *);
extern PetscErrorCode MatIJGetImageSize(Mat, PetscInt *);

extern PetscErrorCode MatIJBinRenumber(Mat, Mat*);

extern PetscErrorCode MatIJMap(Mat, MatIJIndexType, PetscInt,const PetscInt*,const PetscInt*,const PetscScalar*, MatIJIndexType,PetscInt*,PetscInt**,PetscInt**,PetscScalar**,PetscInt**);
extern PetscErrorCode MatIJBin(Mat, MatIJIndexType, PetscInt,const PetscInt*,const PetscInt*,const PetscScalar*,PetscInt*,PetscInt**,PetscInt**,PetscScalar**,PetscInt**);
extern PetscErrorCode MatIJBinMap(Mat,Mat, MatIJIndexType,PetscInt,const PetscInt*,const PetscInt*,const PetscScalar*,MatIJIndexType,PetscInt*,PetscInt**,PetscInt**,PetscScalar**,PetscInt**);


# 6 "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h" 2
# 16 "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"
typedef struct _MatOps *MatOps;
struct _MatOps {

  PetscErrorCode (*setvalues)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*getrow)(Mat,PetscInt,PetscInt *,PetscInt*[],PetscScalar*[]);
  PetscErrorCode (*restorerow)(Mat,PetscInt,PetscInt *,PetscInt *[],PetscScalar *[]);
  PetscErrorCode (*mult)(Mat,Vec,Vec);
  PetscErrorCode (*multadd)(Mat,Vec,Vec,Vec);

  PetscErrorCode (*multtranspose)(Mat,Vec,Vec);
  PetscErrorCode (*multtransposeadd)(Mat,Vec,Vec,Vec);
  PetscErrorCode (*solve)(Mat,Vec,Vec);
  PetscErrorCode (*solveadd)(Mat,Vec,Vec,Vec);
  PetscErrorCode (*solvetranspose)(Mat,Vec,Vec);

  PetscErrorCode (*solvetransposeadd)(Mat,Vec,Vec,Vec);
  PetscErrorCode (*lufactor)(Mat,IS,IS,const MatFactorInfo*);
  PetscErrorCode (*choleskyfactor)(Mat,IS,const MatFactorInfo*);
  PetscErrorCode (*sor)(Mat,Vec,PetscReal,MatSORType,PetscReal,PetscInt,PetscInt,Vec);
  PetscErrorCode (*transpose)(Mat,MatReuse,Mat *);

  PetscErrorCode (*getinfo)(Mat,MatInfoType,MatInfo*);
  PetscErrorCode (*equal)(Mat,Mat,PetscBool *);
  PetscErrorCode (*getdiagonal)(Mat,Vec);
  PetscErrorCode (*diagonalscale)(Mat,Vec,Vec);
  PetscErrorCode (*norm)(Mat,NormType,PetscReal*);

  PetscErrorCode (*assemblybegin)(Mat,MatAssemblyType);
  PetscErrorCode (*assemblyend)(Mat,MatAssemblyType);
  PetscErrorCode (*setoption)(Mat,MatOption,PetscBool );
  PetscErrorCode (*zeroentries)(Mat);

  PetscErrorCode (*zerorows)(Mat,PetscInt,const PetscInt[],PetscScalar,Vec,Vec);
  PetscErrorCode (*lufactorsymbolic)(Mat,Mat,IS,IS,const MatFactorInfo*);
  PetscErrorCode (*lufactornumeric)(Mat,Mat,const MatFactorInfo*);
  PetscErrorCode (*choleskyfactorsymbolic)(Mat,Mat,IS,const MatFactorInfo*);
  PetscErrorCode (*choleskyfactornumeric)(Mat,Mat,const MatFactorInfo*);

  PetscErrorCode (*setuppreallocation)(Mat);
  PetscErrorCode (*ilufactorsymbolic)(Mat,Mat,IS,IS,const MatFactorInfo*);
  PetscErrorCode (*iccfactorsymbolic)(Mat,Mat,IS,const MatFactorInfo*);
  PetscErrorCode (*getarray)(Mat,PetscScalar**);
  PetscErrorCode (*restorearray)(Mat,PetscScalar**);

  PetscErrorCode (*duplicate)(Mat,MatDuplicateOption,Mat*);
  PetscErrorCode (*forwardsolve)(Mat,Vec,Vec);
  PetscErrorCode (*backwardsolve)(Mat,Vec,Vec);
  PetscErrorCode (*ilufactor)(Mat,IS,IS,const MatFactorInfo*);
  PetscErrorCode (*iccfactor)(Mat,IS,const MatFactorInfo*);

  PetscErrorCode (*axpy)(Mat,PetscScalar,Mat,MatStructure);
  PetscErrorCode (*getsubmatrices)(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
  PetscErrorCode (*increaseoverlap)(Mat,PetscInt,IS[],PetscInt);
  PetscErrorCode (*getvalues)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscScalar []);
  PetscErrorCode (*copy)(Mat,Mat,MatStructure);

  PetscErrorCode (*getrowmax)(Mat,Vec,PetscInt[]);
  PetscErrorCode (*scale)(Mat,PetscScalar);
  PetscErrorCode (*shift)(Mat,PetscScalar);
  PetscErrorCode (*diagonalset)(Mat,Vec,InsertMode);
  PetscErrorCode (*zerorowscolumns)(Mat,PetscInt,const PetscInt[],PetscScalar,Vec,Vec);

  PetscErrorCode (*setblocksize)(Mat,PetscInt);
  PetscErrorCode (*getrowij)(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);
  PetscErrorCode (*restorerowij)(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,PetscInt *[],PetscInt *[],PetscBool *);
  PetscErrorCode (*getcolumnij)(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);
  PetscErrorCode (*restorecolumnij)(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);

  PetscErrorCode (*fdcoloringcreate)(Mat,ISColoring,MatFDColoring);
  PetscErrorCode (*coloringpatch)(Mat,PetscInt,PetscInt,ISColoringValue[],ISColoring*);
  PetscErrorCode (*setunfactored)(Mat);
  PetscErrorCode (*permute)(Mat,IS,IS,Mat*);
  PetscErrorCode (*setvaluesblocked)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

  PetscErrorCode (*getsubmatrix)(Mat,IS,IS,MatReuse,Mat*);
  PetscErrorCode (*destroy)(Mat);
  PetscErrorCode (*view)(Mat,PetscViewer);
  PetscErrorCode (*convertfrom)(Mat, const char*,MatReuse,Mat*);
  PetscErrorCode (*usescaledform)(Mat,PetscBool );

  PetscErrorCode (*scalesystem)(Mat,Vec,Vec);
  PetscErrorCode (*unscalesystem)(Mat,Vec,Vec);
  PetscErrorCode (*setlocaltoglobalmapping)(Mat,ISLocalToGlobalMapping,ISLocalToGlobalMapping);
  PetscErrorCode (*setvalueslocal)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*zerorowslocal)(Mat,PetscInt,const PetscInt[],PetscScalar,Vec,Vec);

  PetscErrorCode (*getrowmaxabs)(Mat,Vec,PetscInt[]);
  PetscErrorCode (*getrowminabs)(Mat,Vec,PetscInt[]);
  PetscErrorCode (*convert)(Mat, const char*,MatReuse,Mat*);
  PetscErrorCode (*setcoloring)(Mat,ISColoring);
  PetscErrorCode (*setvaluesadic)(Mat,void*);

  PetscErrorCode (*setvaluesadifor)(Mat,PetscInt,void*);
  PetscErrorCode (*fdcoloringapply)(Mat,MatFDColoring,Vec,MatStructure*,void*);
  PetscErrorCode (*setfromoptions)(Mat);
  PetscErrorCode (*multconstrained)(Mat,Vec,Vec);
  PetscErrorCode (*multtransposeconstrained)(Mat,Vec,Vec);

  PetscErrorCode (*findzerodiagonals)(Mat,IS*);
  PetscErrorCode (*mults)(Mat, Vecs, Vecs);
  PetscErrorCode (*solves)(Mat, Vecs, Vecs);
  PetscErrorCode (*getinertia)(Mat,PetscInt*,PetscInt*,PetscInt*);
  PetscErrorCode (*load)(Mat, PetscViewer);

  PetscErrorCode (*issymmetric)(Mat,PetscReal,PetscBool *);
  PetscErrorCode (*ishermitian)(Mat,PetscReal,PetscBool *);
  PetscErrorCode (*isstructurallysymmetric)(Mat,PetscBool *);
  PetscErrorCode (*setvaluesblockedlocal)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*getvecs)(Mat,Vec*,Vec*);

  PetscErrorCode (*matmult)(Mat,Mat,MatReuse,PetscReal,Mat*);
  PetscErrorCode (*matmultsymbolic)(Mat,Mat,PetscReal,Mat*);
  PetscErrorCode (*matmultnumeric)(Mat,Mat,Mat);
  PetscErrorCode (*ptap)(Mat,Mat,MatReuse,PetscReal,Mat*);
  PetscErrorCode (*ptapsymbolic)(Mat,Mat,PetscReal,Mat*);

  PetscErrorCode (*ptapnumeric)(Mat,Mat,Mat);
  PetscErrorCode (*matmulttranspose)(Mat,Mat,MatReuse,PetscReal,Mat*);
  PetscErrorCode (*matmulttransposesymbolic)(Mat,Mat,PetscReal,Mat*);
  PetscErrorCode (*matmulttransposenumeric)(Mat,Mat,Mat);
  PetscErrorCode (*ptapsymbolic_seqaij)(Mat,Mat,PetscReal,Mat*);

  PetscErrorCode (*ptapnumeric_seqaij)(Mat,Mat,Mat);
  PetscErrorCode (*ptapsymbolic_mpiaij)(Mat,Mat,PetscReal,Mat*);
  PetscErrorCode (*ptapnumeric_mpiaij)(Mat,Mat,Mat);
  PetscErrorCode (*conjugate)(Mat);
  PetscErrorCode (*setsizes)(Mat,PetscInt,PetscInt,PetscInt,PetscInt);

  PetscErrorCode (*setvaluesrow)(Mat,PetscInt,const PetscScalar[]);
  PetscErrorCode (*realpart)(Mat);
  PetscErrorCode (*imaginarypart)(Mat);
  PetscErrorCode (*getrowuppertriangular)(Mat);
  PetscErrorCode (*restorerowuppertriangular)(Mat);

  PetscErrorCode (*matsolve)(Mat,Mat,Mat);
  PetscErrorCode (*getredundantmatrix)(Mat,PetscInt,MPI_Comm,PetscInt,MatReuse,Mat*);
  PetscErrorCode (*getrowmin)(Mat,Vec,PetscInt[]);
  PetscErrorCode (*getcolumnvector)(Mat,Vec,PetscInt);
  PetscErrorCode (*missingdiagonal)(Mat,PetscBool *,PetscInt*);

  PetscErrorCode (*getseqnonzerostructure)(Mat,Mat *);
  PetscErrorCode (*create)(Mat);
  PetscErrorCode (*getghosts)(Mat,PetscInt*,const PetscInt *[]);
  PetscErrorCode (*getlocalsubmatrix)(Mat,IS,IS,Mat*);
  PetscErrorCode (*restorelocalsubmatrix)(Mat,IS,IS,Mat*);

  PetscErrorCode (*multdiagonalblock)(Mat,Vec,Vec);
  PetscErrorCode (*hermitiantranspose)(Mat,MatReuse,Mat*);
  PetscErrorCode (*multhermitiantranspose)(Mat,Vec,Vec);
  PetscErrorCode (*multhermitiantransposeadd)(Mat,Vec,Vec,Vec);
  PetscErrorCode (*getmultiprocblock)(Mat,MPI_Comm,Mat*);

  PetscErrorCode (*findnonzerorows)(Mat,IS*);
  PetscErrorCode (*getcolumnnorms)(Mat,NormType,PetscReal*);
  PetscErrorCode (*invertblockdiagonal)(Mat,PetscScalar**);
  PetscErrorCode (*dummy4)(Mat,Vec,Vec,Vec);
  PetscErrorCode (*getsubmatricesparallel)(Mat,PetscInt,const IS[], const IS[], MatReuse, Mat**);

  PetscErrorCode (*setvaluesbatch)(Mat,PetscInt,PetscInt,PetscInt*,const PetscScalar*);

  PetscErrorCode (*setstencil)(Mat, PetscInt, const PetscInt[],const PetscInt[], PetscInt);

  PetscErrorCode (*setgrid)(Mat,PetscInt, PetscInt, PetscInt);

};





typedef struct _p_MatBaseName* MatBaseName;
struct _p_MatBaseName {
  char *bname,*sname,*mname;
  MatBaseName next;
};

extern MatBaseName MatBaseNameList;




extern PetscErrorCode MatConvert_Basic(Mat, const char*,MatReuse,Mat*);
extern PetscErrorCode MatCopy_Basic(Mat,Mat,MatStructure);
extern PetscErrorCode MatView_Private(Mat);

extern PetscErrorCode MatHeaderMerge(Mat,Mat);
extern PetscErrorCode MatHeaderReplace(Mat,Mat);
extern PetscErrorCode MatAXPYGetxtoy_Private(PetscInt,PetscInt*,PetscInt*,PetscInt*, PetscInt*,PetscInt*,PetscInt*, PetscInt**);
extern PetscErrorCode MatPtAP_Basic(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatDiagonalSet_Default(Mat,Vec,InsertMode);







typedef struct _MatStashSpace *PetscMatStashSpace;

struct _MatStashSpace {
  PetscMatStashSpace next;
  PetscScalar *space_head,*val;
  PetscInt *idx,*idy;
  PetscInt total_space_size;
  PetscInt local_used;
  PetscInt local_remaining;
};

extern PetscErrorCode PetscMatStashSpaceGet(PetscInt,PetscInt,PetscMatStashSpace *);
extern PetscErrorCode PetscMatStashSpaceContiguous(PetscInt,PetscMatStashSpace *,PetscScalar *,PetscInt *,PetscInt *);
extern PetscErrorCode PetscMatStashSpaceDestroy(PetscMatStashSpace*);

typedef struct {
  PetscInt nmax;
  PetscInt umax;
  PetscInt oldnmax;
  PetscInt n;
  PetscInt bs;
  PetscInt reallocs;
  PetscMatStashSpace space_head,space;

  MPI_Comm comm;
  PetscMPIInt size,rank;
  PetscMPIInt tag1,tag2;
  MPI_Request *send_waits;
  MPI_Request *recv_waits;
  MPI_Status *send_status;
  PetscInt nsends,nrecvs;
  PetscScalar *svalues;
  PetscInt *sindices;
  PetscScalar **rvalues;
  PetscInt **rindices;
  PetscInt nprocessed;
  PetscMPIInt *flg_v;
  PetscBool reproduce;
  PetscInt reproduce_count;
} MatStash;

extern PetscErrorCode MatStashCreate_Private(MPI_Comm,PetscInt,MatStash*);
extern PetscErrorCode MatStashDestroy_Private(MatStash*);
extern PetscErrorCode MatStashScatterEnd_Private(MatStash*);
extern PetscErrorCode MatStashSetInitialSize_Private(MatStash*,PetscInt);
extern PetscErrorCode MatStashGetInfo_Private(MatStash*,PetscInt*,PetscInt*);
extern PetscErrorCode MatStashValuesRow_Private(MatStash*,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],PetscBool );
extern PetscErrorCode MatStashValuesCol_Private(MatStash*,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],PetscInt,PetscBool );
extern PetscErrorCode MatStashValuesRowBlocked_Private(MatStash*,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],PetscInt,PetscInt,PetscInt);
extern PetscErrorCode MatStashValuesColBlocked_Private(MatStash*,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],PetscInt,PetscInt,PetscInt);
extern PetscErrorCode MatStashScatterBegin_Private(Mat,MatStash*,PetscInt*);
extern PetscErrorCode MatStashScatterGetMesg_Private(MatStash*,PetscMPIInt*,PetscInt**,PetscInt**,PetscScalar**,PetscInt*);

typedef struct {
  PetscInt dim;
  PetscInt dims[4];
  PetscInt starts[4];
  PetscBool noc;
} MatStencilInfo;


typedef struct {
  PetscBool check;
  PetscBool use;
  PetscInt nrows;
  PetscInt *i;
  PetscInt *rindex;
} Mat_CompressedRow;
extern PetscErrorCode MatCheckCompressedRow(Mat,Mat_CompressedRow*,PetscInt*,PetscInt,PetscReal);

struct _p_Mat {
  _p_PetscObject hdr; struct _MatOps *ops;
  PetscLayout rmap,cmap;
  void *data;
  MatFactorType factortype;
  PetscBool assembled;
  PetscBool was_assembled;
  PetscInt num_ass;
  PetscBool same_nonzero;
  MatInfo info;
  InsertMode insertmode;
  MatStash stash,bstash;
  MatNullSpace nullsp;
  MatNullSpace nearnullsp;
  PetscBool preallocated;
  MatStencilInfo stencil;
  PetscBool symmetric,hermitian,structurally_symmetric,spd;
  PetscBool symmetric_set,hermitian_set,structurally_symmetric_set,spd_set;
  PetscBool symmetric_eternal;
  PetscBool nooffprocentries,nooffproczerorows;

  PetscCUSPFlag valid_GPU_matrix;

  void *spptr;
  char* solvertype;
  };


extern PetscErrorCode MatAXPY_Basic(Mat,PetscScalar,Mat,MatStructure);
extern PetscErrorCode MatAXPY_BasicWithPreallocation(Mat,Mat,PetscScalar,Mat,MatStructure);




typedef struct _MatPartitioningOps *MatPartitioningOps;
struct _MatPartitioningOps {
  PetscErrorCode (*apply)(MatPartitioning,IS*);
  PetscErrorCode (*setfromoptions)(MatPartitioning);
  PetscErrorCode (*destroy)(MatPartitioning);
  PetscErrorCode (*view)(MatPartitioning,PetscViewer);
};

struct _p_MatPartitioning {
  _p_PetscObject hdr; struct _MatPartitioningOps *ops;
  Mat adj;
  PetscInt *vertex_weights;
  PetscReal *part_weights;
  PetscInt n;
  void *data;
  PetscInt setupcalled;
};
# 376 "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"
struct _p_MatFDColoring{
  _p_PetscObject hdr; int *ops;
  PetscInt M,N,m;
  PetscInt rstart;
  PetscInt ncolors;
  PetscInt *ncolumns;
  PetscInt **columns;
  PetscInt *nrows;
  PetscInt **rows;
  PetscInt **columnsforrow;
  PetscReal error_rel;
  PetscReal umin;
  Vec w1,w2,w3;
  PetscErrorCode (*f)(void);
  void *fctx;
  PetscInt **vscaleforrow;
  Vec vscale;
  Vec F;
  PetscInt currentcolor;
  const char *htype;
  ISColoringType ctype;

  void *ftn_func_pointer,*ftn_func_cntx;
};




struct _p_MatNullSpace {
  _p_PetscObject hdr; int *ops;
  PetscBool has_cnst;
  PetscInt n;
  Vec* vecs;
  PetscScalar* alpha;
  Vec vec;
  PetscErrorCode (*remove)(MatNullSpace,Vec,void*);
  void* rmctx;
};




typedef struct {
  PetscInt nshift,nshift_max;
  PetscReal shift_amount,shift_lo,shift_hi,shift_top,shift_fraction;
  PetscBool newshift;
  PetscReal rs;
  PetscScalar pv;
} FactorShiftCtx;

extern PetscErrorCode MatFactorDumpMatrix(Mat);



static inline PetscErrorCode MatPivotCheck_nz(Mat mat,const MatFactorInfo *info,FactorShiftCtx *sctx,PetscInt row)
{
  PetscReal _rs = sctx->rs;
  PetscReal _zero = info->zeropivot*_rs;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 435; petscstack->currentsize++; } do { if (strcmp(__func__,"MatPivotCheck_nz") && strcmp("MatPivotCheck_nz","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/matimpl.h",435,"MatPivotCheck_nz","__func__",__func__); } } while (0); } while (0);
  if (PetscAbsScalar(sctx->pv) <= _zero){

    if (!sctx->nshift) sctx->shift_amount = info->shiftamount;
    else sctx->shift_amount *= 2.0;
    sctx->newshift = PETSC_TRUE;
    (sctx->nshift)++;
  } else {
    sctx->newshift = PETSC_FALSE;
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode MatPivotCheck_pd(Mat mat,const MatFactorInfo *info,FactorShiftCtx *sctx,PetscInt row)
{
  PetscReal _rs = sctx->rs;
  PetscReal _zero = info->zeropivot*_rs;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 455; petscstack->currentsize++; } do { if (strcmp(__func__,"MatPivotCheck_pd") && strcmp("MatPivotCheck_pd","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/matimpl.h",455,"MatPivotCheck_pd","__func__",__func__); } } while (0); } while (0);
  if ((sctx->pv) <= _zero){

    if (sctx->nshift == sctx->nshift_max) {
      sctx->shift_fraction = sctx->shift_hi;
    } else {
      sctx->shift_lo = sctx->shift_fraction;
      sctx->shift_fraction = (sctx->shift_hi+sctx->shift_lo)/2.;
    }
    sctx->shift_amount = sctx->shift_fraction * sctx->shift_top;
    sctx->nshift++;
    sctx->newshift = PETSC_TRUE;
  } else {
    sctx->newshift = PETSC_FALSE;
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode MatPivotCheck_inblocks(Mat mat,const MatFactorInfo *info,FactorShiftCtx *sctx,PetscInt row)
{
  PetscReal _zero = info->zeropivot;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 479; petscstack->currentsize++; } do { if (strcmp(__func__,"MatPivotCheck_inblocks") && strcmp("MatPivotCheck_inblocks","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/matimpl.h",479,"MatPivotCheck_inblocks","__func__",__func__); } } while (0); } while (0);
  if (PetscAbsScalar(sctx->pv) <= _zero){
    sctx->pv += info->shiftamount;
    sctx->shift_amount = 0.0;
    sctx->nshift++;
  }
  sctx->newshift = PETSC_FALSE;
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode MatPivotCheck_none(Mat mat,const MatFactorInfo *info,FactorShiftCtx *sctx,PetscInt row)
{
  PetscReal _zero = info->zeropivot;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 495; petscstack->currentsize++; } do { if (strcmp(__func__,"MatPivotCheck_none") && strcmp("MatPivotCheck_none","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/matimpl.h",495,"MatPivotCheck_none","__func__",__func__); } } while (0); } while (0);
  sctx->newshift = PETSC_FALSE;
  if (PetscAbsScalar(sctx->pv) <= _zero) {
    PetscErrorCode ierr;
    PetscBool flg = PETSC_FALSE;

    ierr = PetscOptionsGetBool(0,"-mat_dump",&flg,0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),501,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    if (flg) {
      ierr = MatView(mat,PETSC_VIEWER_BINARY_(((PetscObject)mat)->comm));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),503,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    return PetscError(((MPI_Comm)0x44000001),505,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",71,PETSC_ERROR_INITIAL,"Zero pivot row %D value %G tolerance %G",row,PetscAbsScalar(sctx->pv),_zero);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



static inline PetscErrorCode MatPivotCheck(Mat mat,const MatFactorInfo *info,FactorShiftCtx *sctx,PetscInt row)
{
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 516; petscstack->currentsize++; } do { if (strcmp(__func__,"MatPivotCheck") && strcmp("MatPivotCheck","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/private/matimpl.h",516,"MatPivotCheck","__func__",__func__); } } while (0); } while (0);
  if (info->shifttype == (PetscReal) MAT_SHIFT_NONZERO){
    ierr = MatPivotCheck_nz(mat,info,sctx,row);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),518,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else if (info->shifttype == (PetscReal) MAT_SHIFT_POSITIVE_DEFINITE){
    ierr = MatPivotCheck_pd(mat,info,sctx,row);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),520,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else if (info->shifttype == (PetscReal) MAT_SHIFT_INBLOCKS){
    ierr = MatPivotCheck_inblocks(mat,info,sctx,row);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),522,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    ierr = MatPivotCheck_none(mat,info,sctx,row);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),524,__func__,"/home/dpnkarthik/petsc-rnet/include/private/matimpl.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 995 "/home/dpnkarthik/petsc-rnet/include/private/matimpl.h"
extern PetscLogEvent MAT_Mult, MAT_MultMatrixFree, MAT_Mults, MAT_MultConstrained, MAT_MultAdd, MAT_MultTranspose;
extern PetscLogEvent MAT_MultTransposeConstrained, MAT_MultTransposeAdd, MAT_Solve, MAT_Solves, MAT_SolveAdd, MAT_SolveTranspose;
extern PetscLogEvent MAT_SolveTransposeAdd, MAT_SOR, MAT_ForwardSolve, MAT_BackwardSolve, MAT_LUFactor, MAT_LUFactorSymbolic;
extern PetscLogEvent MAT_LUFactorNumeric, MAT_CholeskyFactor, MAT_CholeskyFactorSymbolic, MAT_CholeskyFactorNumeric, MAT_ILUFactor;
extern PetscLogEvent MAT_ILUFactorSymbolic, MAT_ICCFactorSymbolic, MAT_Copy, MAT_Convert, MAT_Scale, MAT_AssemblyBegin;
extern PetscLogEvent MAT_AssemblyEnd, MAT_SetValues, MAT_GetValues, MAT_GetRow, MAT_GetRowIJ, MAT_GetSubMatrices, MAT_GetColoring, MAT_GetOrdering, MAT_GetRedundantMatrix;
extern PetscLogEvent MAT_IncreaseOverlap, MAT_Partitioning, MAT_ZeroEntries, MAT_Load, MAT_View, MAT_AXPY, MAT_FDColoringCreate;
extern PetscLogEvent MAT_FDColoringApply, MAT_Transpose, MAT_FDColoringFunction;
extern PetscLogEvent MAT_MatMult, MAT_MatSolve,MAT_MatMultSymbolic, MAT_MatMultNumeric,MAT_Getlocalmatcondensed,MAT_GetBrowsOfAcols,MAT_GetBrowsOfAocols;
extern PetscLogEvent MAT_PtAP, MAT_PtAPSymbolic, MAT_PtAPNumeric,MAT_Seqstompinum,MAT_Seqstompisym,MAT_Seqstompi,MAT_Getlocalmat;

extern PetscLogEvent MAT_MatMultTranspose, MAT_MatMultTransposeSymbolic, MAT_MatMultTransposeNumeric;
extern PetscLogEvent MAT_Applypapt, MAT_Applypapt_symbolic, MAT_Applypapt_numeric;
extern PetscLogEvent MAT_Getsymtranspose, MAT_Transpose_SeqAIJ, MAT_Getsymtransreduced,MAT_GetSequentialNonzeroStructure;

extern PetscLogEvent MATMFFD_Mult;
extern PetscLogEvent MAT_GetMultiProcBlock;
extern PetscLogEvent MAT_CUSPCopyToGPU, MAT_SetValuesBatch, MAT_SetValuesBatchI, MAT_SetValuesBatchII, MAT_SetValuesBatchIII, MAT_SetValuesBatchIV;
# 5 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/baij.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h" 1





# 1 "/home/dpnkarthik/softwares/cuda/include/cuda.h" 1
# 149 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
typedef unsigned long long CUdeviceptr;






typedef int CUdevice;
typedef struct CUctx_st *CUcontext;
typedef struct CUmod_st *CUmodule;
typedef struct CUfunc_st *CUfunction;
typedef struct CUarray_st *CUarray;
typedef struct CUtexref_st *CUtexref;
typedef struct CUsurfref_st *CUsurfref;
typedef struct CUevent_st *CUevent;
typedef struct CUstream_st *CUstream;
typedef struct CUgraphicsResource_st *CUgraphicsResource;

typedef struct CUuuid_st {
    char bytes[16];
} CUuuid;




typedef enum CUctx_flags_enum {
    CU_CTX_SCHED_AUTO = 0x00,
    CU_CTX_SCHED_SPIN = 0x01,
    CU_CTX_SCHED_YIELD = 0x02,
    CU_CTX_SCHED_BLOCKING_SYNC = 0x04,
    CU_CTX_BLOCKING_SYNC = 0x04,
    CU_CTX_SCHED_MASK = 0x07,
    CU_CTX_MAP_HOST = 0x08,
    CU_CTX_LMEM_RESIZE_TO_MAX = 0x10,
    CU_CTX_FLAGS_MASK = 0x1f
} CUctx_flags;




typedef enum CUevent_flags_enum {
    CU_EVENT_DEFAULT = 0,
    CU_EVENT_BLOCKING_SYNC = 1,
    CU_EVENT_DISABLE_TIMING = 2
} CUevent_flags;




typedef enum CUarray_format_enum {
    CU_AD_FORMAT_UNSIGNED_INT8 = 0x01,
    CU_AD_FORMAT_UNSIGNED_INT16 = 0x02,
    CU_AD_FORMAT_UNSIGNED_INT32 = 0x03,
    CU_AD_FORMAT_SIGNED_INT8 = 0x08,
    CU_AD_FORMAT_SIGNED_INT16 = 0x09,
    CU_AD_FORMAT_SIGNED_INT32 = 0x0a,
    CU_AD_FORMAT_HALF = 0x10,
    CU_AD_FORMAT_FLOAT = 0x20
} CUarray_format;




typedef enum CUaddress_mode_enum {
    CU_TR_ADDRESS_MODE_WRAP = 0,
    CU_TR_ADDRESS_MODE_CLAMP = 1,
    CU_TR_ADDRESS_MODE_MIRROR = 2,
    CU_TR_ADDRESS_MODE_BORDER = 3
} CUaddress_mode;




typedef enum CUfilter_mode_enum {
    CU_TR_FILTER_MODE_POINT = 0,
    CU_TR_FILTER_MODE_LINEAR = 1
} CUfilter_mode;




typedef enum CUdevice_attribute_enum {
    CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 1,
    CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X = 2,
    CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y = 3,
    CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z = 4,
    CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X = 5,
    CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y = 6,
    CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z = 7,
    CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK = 8,
    CU_DEVICE_ATTRIBUTE_SHARED_MEMORY_PER_BLOCK = 8,
    CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY = 9,
    CU_DEVICE_ATTRIBUTE_WARP_SIZE = 10,
    CU_DEVICE_ATTRIBUTE_MAX_PITCH = 11,
    CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK = 12,
    CU_DEVICE_ATTRIBUTE_REGISTERS_PER_BLOCK = 12,
    CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13,
    CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT = 14,
    CU_DEVICE_ATTRIBUTE_GPU_OVERLAP = 15,
    CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16,
    CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT = 17,
    CU_DEVICE_ATTRIBUTE_INTEGRATED = 18,
    CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY = 19,
    CU_DEVICE_ATTRIBUTE_COMPUTE_MODE = 20,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_WIDTH = 21,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_WIDTH = 22,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_HEIGHT = 23,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_WIDTH = 24,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_HEIGHT = 25,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE3D_DEPTH = 26,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_WIDTH = 27,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_HEIGHT = 28,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_LAYERED_LAYERS = 29,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_WIDTH = 27,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_HEIGHT = 28,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE2D_ARRAY_NUMSLICES = 29,
    CU_DEVICE_ATTRIBUTE_SURFACE_ALIGNMENT = 30,
    CU_DEVICE_ATTRIBUTE_CONCURRENT_KERNELS = 31,
    CU_DEVICE_ATTRIBUTE_ECC_ENABLED = 32,
    CU_DEVICE_ATTRIBUTE_PCI_BUS_ID = 33,
    CU_DEVICE_ATTRIBUTE_PCI_DEVICE_ID = 34,
    CU_DEVICE_ATTRIBUTE_TCC_DRIVER = 35,
    CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE = 36,
    CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH = 37,
    CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE = 38,
    CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR = 39,
    CU_DEVICE_ATTRIBUTE_ASYNC_ENGINE_COUNT = 40,
    CU_DEVICE_ATTRIBUTE_UNIFIED_ADDRESSING = 41,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_LAYERED_WIDTH = 42,
    CU_DEVICE_ATTRIBUTE_MAXIMUM_TEXTURE1D_LAYERED_LAYERS = 43,
    CU_DEVICE_ATTRIBUTE_PCI_DOMAIN_ID = 50
} CUdevice_attribute;




typedef struct CUdevprop_st {
    int maxThreadsPerBlock;
    int maxThreadsDim[3];
    int maxGridSize[3];
    int sharedMemPerBlock;
    int totalConstantMemory;
    int SIMDWidth;
    int memPitch;
    int regsPerBlock;
    int clockRate;
    int textureAlign;
} CUdevprop;




typedef enum CUpointer_attribute_enum {
    CU_POINTER_ATTRIBUTE_CONTEXT = 1,
    CU_POINTER_ATTRIBUTE_MEMORY_TYPE = 2,
    CU_POINTER_ATTRIBUTE_DEVICE_POINTER = 3,
    CU_POINTER_ATTRIBUTE_HOST_POINTER = 4,
} CUpointer_attribute;




typedef enum CUfunction_attribute_enum {





    CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK = 0,






    CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES = 1,





    CU_FUNC_ATTRIBUTE_CONST_SIZE_BYTES = 2,




    CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES = 3,




    CU_FUNC_ATTRIBUTE_NUM_REGS = 4,
# 349 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CU_FUNC_ATTRIBUTE_PTX_VERSION = 5,
# 358 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CU_FUNC_ATTRIBUTE_BINARY_VERSION = 6,

    CU_FUNC_ATTRIBUTE_MAX
} CUfunction_attribute;




typedef enum CUfunc_cache_enum {
    CU_FUNC_CACHE_PREFER_NONE = 0x00,
    CU_FUNC_CACHE_PREFER_SHARED = 0x01,
    CU_FUNC_CACHE_PREFER_L1 = 0x02
} CUfunc_cache;




typedef enum CUmemorytype_enum {
    CU_MEMORYTYPE_HOST = 0x01,
    CU_MEMORYTYPE_DEVICE = 0x02,
    CU_MEMORYTYPE_ARRAY = 0x03,
    CU_MEMORYTYPE_UNIFIED = 0x04
} CUmemorytype;




typedef enum CUcomputemode_enum {
    CU_COMPUTEMODE_DEFAULT = 0,
    CU_COMPUTEMODE_EXCLUSIVE = 1,
    CU_COMPUTEMODE_PROHIBITED = 2,
    CU_COMPUTEMODE_EXCLUSIVE_PROCESS = 3
} CUcomputemode;




typedef enum CUjit_option_enum
{




    CU_JIT_MAX_REGISTERS = 0,
# 414 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CU_JIT_THREADS_PER_BLOCK,






    CU_JIT_WALL_TIME,







    CU_JIT_INFO_LOG_BUFFER,







    CU_JIT_INFO_LOG_BUFFER_SIZE_BYTES,







    CU_JIT_ERROR_LOG_BUFFER,







    CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES,






    CU_JIT_OPTIMIZATION_LEVEL,






    CU_JIT_TARGET_FROM_CUCONTEXT,





    CU_JIT_TARGET,






    CU_JIT_FALLBACK_STRATEGY

} CUjit_option;




typedef enum CUjit_target_enum
{
    CU_TARGET_COMPUTE_10 = 0,
    CU_TARGET_COMPUTE_11,
    CU_TARGET_COMPUTE_12,
    CU_TARGET_COMPUTE_13,
    CU_TARGET_COMPUTE_20,
    CU_TARGET_COMPUTE_21
} CUjit_target;




typedef enum CUjit_fallback_enum
{
    CU_PREFER_PTX = 0,

    CU_PREFER_BINARY

} CUjit_fallback;




typedef enum CUgraphicsRegisterFlags_enum {
    CU_GRAPHICS_REGISTER_FLAGS_NONE = 0x00,
    CU_GRAPHICS_REGISTER_FLAGS_READ_ONLY = 0x01,
    CU_GRAPHICS_REGISTER_FLAGS_WRITE_DISCARD = 0x02,
    CU_GRAPHICS_REGISTER_FLAGS_SURFACE_LDST = 0x04
} CUgraphicsRegisterFlags;




typedef enum CUgraphicsMapResourceFlags_enum {
    CU_GRAPHICS_MAP_RESOURCE_FLAGS_NONE = 0x00,
    CU_GRAPHICS_MAP_RESOURCE_FLAGS_READ_ONLY = 0x01,
    CU_GRAPHICS_MAP_RESOURCE_FLAGS_WRITE_DISCARD = 0x02
} CUgraphicsMapResourceFlags;




typedef enum CUarray_cubemap_face_enum {
    CU_CUBEMAP_FACE_POSITIVE_X = 0x00,
    CU_CUBEMAP_FACE_NEGATIVE_X = 0x01,
    CU_CUBEMAP_FACE_POSITIVE_Y = 0x02,
    CU_CUBEMAP_FACE_NEGATIVE_Y = 0x03,
    CU_CUBEMAP_FACE_POSITIVE_Z = 0x04,
    CU_CUBEMAP_FACE_NEGATIVE_Z = 0x05
} CUarray_cubemap_face;




typedef enum CUlimit_enum {
    CU_LIMIT_STACK_SIZE = 0x00,
    CU_LIMIT_PRINTF_FIFO_SIZE = 0x01,
    CU_LIMIT_MALLOC_HEAP_SIZE = 0x02
} CUlimit;




typedef enum cudaError_enum {





    CUDA_SUCCESS = 0,





    CUDA_ERROR_INVALID_VALUE = 1,





    CUDA_ERROR_OUT_OF_MEMORY = 2,





    CUDA_ERROR_NOT_INITIALIZED = 3,




    CUDA_ERROR_DEINITIALIZED = 4,





    CUDA_ERROR_PROFILER_DISABLED = 5,




    CUDA_ERROR_PROFILER_NOT_INITIALIZED = 6,




    CUDA_ERROR_PROFILER_ALREADY_STARTED = 7,




    CUDA_ERROR_PROFILER_ALREADY_STOPPED = 8,




    CUDA_ERROR_NO_DEVICE = 100,





    CUDA_ERROR_INVALID_DEVICE = 101,






    CUDA_ERROR_INVALID_IMAGE = 200,
# 629 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CUDA_ERROR_INVALID_CONTEXT = 201,
# 638 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CUDA_ERROR_CONTEXT_ALREADY_CURRENT = 202,




    CUDA_ERROR_MAP_FAILED = 205,




    CUDA_ERROR_UNMAP_FAILED = 206,





    CUDA_ERROR_ARRAY_IS_MAPPED = 207,




    CUDA_ERROR_ALREADY_MAPPED = 208,







    CUDA_ERROR_NO_BINARY_FOR_GPU = 209,




    CUDA_ERROR_ALREADY_ACQUIRED = 210,




    CUDA_ERROR_NOT_MAPPED = 211,





    CUDA_ERROR_NOT_MAPPED_AS_ARRAY = 212,





    CUDA_ERROR_NOT_MAPPED_AS_POINTER = 213,





    CUDA_ERROR_ECC_UNCORRECTABLE = 214,





    CUDA_ERROR_UNSUPPORTED_LIMIT = 215,






    CUDA_ERROR_CONTEXT_ALREADY_IN_USE = 216,




    CUDA_ERROR_INVALID_SOURCE = 300,




    CUDA_ERROR_FILE_NOT_FOUND = 301,




    CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND = 302,




    CUDA_ERROR_SHARED_OBJECT_INIT_FAILED = 303,




    CUDA_ERROR_OPERATING_SYSTEM = 304,






    CUDA_ERROR_INVALID_HANDLE = 400,






    CUDA_ERROR_NOT_FOUND = 500,
# 756 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CUDA_ERROR_NOT_READY = 600,
# 767 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CUDA_ERROR_LAUNCH_FAILED = 700,
# 778 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES = 701,
# 789 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
    CUDA_ERROR_LAUNCH_TIMEOUT = 702,





    CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING = 703,






    CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED = 704,






    CUDA_ERROR_PEER_ACCESS_NOT_ENABLED = 705,





    CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE = 708,






    CUDA_ERROR_CONTEXT_IS_DESTROYED = 709,




    CUDA_ERROR_UNKNOWN = 999
} CUresult;
# 869 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
typedef struct CUDA_MEMCPY2D_st {
    size_t srcXInBytes;
    size_t srcY;

    CUmemorytype srcMemoryType;
    const void *srcHost;
    CUdeviceptr srcDevice;
    CUarray srcArray;
    size_t srcPitch;

    size_t dstXInBytes;
    size_t dstY;

    CUmemorytype dstMemoryType;
    void *dstHost;
    CUdeviceptr dstDevice;
    CUarray dstArray;
    size_t dstPitch;

    size_t WidthInBytes;
    size_t Height;
} CUDA_MEMCPY2D;




typedef struct CUDA_MEMCPY3D_st {
    size_t srcXInBytes;
    size_t srcY;
    size_t srcZ;
    size_t srcLOD;
    CUmemorytype srcMemoryType;
    const void *srcHost;
    CUdeviceptr srcDevice;
    CUarray srcArray;
    void *reserved0;
    size_t srcPitch;
    size_t srcHeight;

    size_t dstXInBytes;
    size_t dstY;
    size_t dstZ;
    size_t dstLOD;
    CUmemorytype dstMemoryType;
    void *dstHost;
    CUdeviceptr dstDevice;
    CUarray dstArray;
    void *reserved1;
    size_t dstPitch;
    size_t dstHeight;

    size_t WidthInBytes;
    size_t Height;
    size_t Depth;
} CUDA_MEMCPY3D;




typedef struct CUDA_MEMCPY3D_PEER_st {
    size_t srcXInBytes;
    size_t srcY;
    size_t srcZ;
    size_t srcLOD;
    CUmemorytype srcMemoryType;
    const void *srcHost;
    CUdeviceptr srcDevice;
    CUarray srcArray;
    CUcontext srcContext;
    size_t srcPitch;
    size_t srcHeight;

    size_t dstXInBytes;
    size_t dstY;
    size_t dstZ;
    size_t dstLOD;
    CUmemorytype dstMemoryType;
    void *dstHost;
    CUdeviceptr dstDevice;
    CUarray dstArray;
    CUcontext dstContext;
    size_t dstPitch;
    size_t dstHeight;

    size_t WidthInBytes;
    size_t Height;
    size_t Depth;
} CUDA_MEMCPY3D_PEER;




typedef struct CUDA_ARRAY_DESCRIPTOR_st
{
    size_t Width;
    size_t Height;

    CUarray_format Format;
    unsigned int NumChannels;
} CUDA_ARRAY_DESCRIPTOR;




typedef struct CUDA_ARRAY3D_DESCRIPTOR_st
{
    size_t Width;
    size_t Height;
    size_t Depth;

    CUarray_format Format;
    unsigned int NumChannels;
    unsigned int Flags;
} CUDA_ARRAY3D_DESCRIPTOR;
# 1095 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuInit(unsigned int Flags);
# 1122 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDriverGetVersion(int *driverVersion);
# 1160 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceGet(CUdevice *device, int ordinal);
# 1186 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceGetCount(int *count);
# 1215 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceGetName(char *name, int len, CUdevice dev);
# 1244 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceComputeCapability(int *major, int *minor, CUdevice dev);
# 1272 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceTotalMem_v2(size_t *bytes, CUdevice dev);
# 1332 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceGetProperties(CUdevprop *prop, CUdevice dev);
# 1444 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceGetAttribute(int *pi, CUdevice_attribute attrib, CUdevice dev);
# 1544 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxCreate_v2(CUcontext *pctx, unsigned int flags, CUdevice dev);
# 1583 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxDestroy_v2(CUcontext ctx);
# 1633 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxAttach(CUcontext *pctx, unsigned int flags);
# 1668 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxDetach(CUcontext ctx);
# 1704 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxPushCurrent_v2(CUcontext ctx);
# 1737 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxPopCurrent_v2(CUcontext *pctx);
# 1763 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxSetCurrent(CUcontext ctx);
# 1782 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxGetCurrent(CUcontext *pctx);
# 1811 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxGetDevice(CUdevice *device);
# 1839 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxSynchronize(void);
# 1900 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxSetLimit(CUlimit limit, size_t value);
# 1933 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxGetLimit(size_t *pvalue, CUlimit limit);
# 1974 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxGetCacheConfig(CUfunc_cache *pconfig);
# 2022 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxSetCacheConfig(CUfunc_cache config);
# 2057 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxGetApiVersion(CUcontext ctx, unsigned int *version);
# 2106 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleLoad(CUmodule *module, const char *fname);
# 2140 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleLoadData(CUmodule *module, const void *image);
# 2219 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleLoadDataEx(CUmodule *module, const void *image, unsigned int numOptions, CUjit_option *options, void **optionValues);
# 2259 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleLoadFatBinary(CUmodule *module, const void *fatCubin);
# 2284 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleUnload(CUmodule hmod);
# 2314 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleGetFunction(CUfunction *hfunc, CUmodule hmod, const char *name);
# 2348 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleGetGlobal_v2(CUdeviceptr *dptr, size_t *bytes, CUmodule hmod, const char *name);
# 2382 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleGetTexRef(CUtexref *pTexRef, CUmodule hmod, const char *name);
# 2413 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuModuleGetSurfRef(CUsurfref *pSurfRef, CUmodule hmod, const char *name);
# 2456 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemGetInfo_v2(size_t *free, size_t *total);
# 2489 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemAlloc_v2(CUdeviceptr *dptr, size_t bytesize);
# 2550 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemAllocPitch_v2(CUdeviceptr *dptr, size_t *pPitch, size_t WidthInBytes, size_t Height, unsigned int ElementSizeBytes);
# 2579 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemFree_v2(CUdeviceptr dptr);
# 2612 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemGetAddressRange_v2(CUdeviceptr *pbase, size_t *psize, CUdeviceptr dptr);
# 2658 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemAllocHost_v2(void **pp, size_t bytesize);
# 2688 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemFreeHost(void *p);
# 2770 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemHostAlloc(void **pp, size_t bytesize, unsigned int Flags);
# 2808 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemHostGetDevicePointer_v2(CUdeviceptr *pdptr, void *p, unsigned int Flags);
# 2833 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemHostGetFlags(unsigned int *pFlags, void *p);
# 2896 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemHostRegister(void *p, size_t bytesize, unsigned int Flags);
# 2919 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemHostUnregister(void *p);
# 2955 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy(CUdeviceptr dst, CUdeviceptr src, size_t ByteCount);
# 2988 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyPeer(CUdeviceptr dstDevice, CUcontext dstContext, CUdeviceptr srcDevice, CUcontext srcContext, size_t ByteCount);
# 3026 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyHtoD_v2(CUdeviceptr dstDevice, const void *srcHost, size_t ByteCount);
# 3059 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyDtoH_v2(void *dstHost, CUdeviceptr srcDevice, size_t ByteCount);
# 3092 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyDtoD_v2(CUdeviceptr dstDevice, CUdeviceptr srcDevice, size_t ByteCount);
# 3126 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyDtoA_v2(CUarray dstArray, size_t dstOffset, CUdeviceptr srcDevice, size_t ByteCount);
# 3162 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyAtoD_v2(CUdeviceptr dstDevice, CUarray srcArray, size_t srcOffset, size_t ByteCount);
# 3196 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyHtoA_v2(CUarray dstArray, size_t dstOffset, const void *srcHost, size_t ByteCount);
# 3230 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyAtoH_v2(void *dstHost, CUarray srcArray, size_t srcOffset, size_t ByteCount);
# 3268 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyAtoA_v2(CUarray dstArray, size_t dstOffset, CUarray srcArray, size_t srcOffset, size_t ByteCount);
# 3428 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy2D_v2(const CUDA_MEMCPY2D *pCopy);
# 3586 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy2DUnaligned_v2(const CUDA_MEMCPY2D *pCopy);
# 3753 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy3D_v2(const CUDA_MEMCPY3D *pCopy);
# 3784 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy3DPeer(const CUDA_MEMCPY3D_PEER *pCopy);
# 3824 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyAsync(CUdeviceptr dst, CUdeviceptr src, size_t ByteCount, CUstream hStream);
# 3855 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyPeerAsync(CUdeviceptr dstDevice, CUcontext dstContext, CUdeviceptr srcDevice, CUcontext srcContext, size_t ByteCount, CUstream hStream);
# 3897 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyHtoDAsync_v2(CUdeviceptr dstDevice, const void *srcHost, size_t ByteCount, CUstream hStream);
# 3937 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyDtoHAsync_v2(void *dstHost, CUdeviceptr srcDevice, size_t ByteCount, CUstream hStream);
# 3974 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyDtoDAsync_v2(CUdeviceptr dstDevice, CUdeviceptr srcDevice, size_t ByteCount, CUstream hStream);
# 4016 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyHtoAAsync_v2(CUarray dstArray, size_t dstOffset, const void *srcHost, size_t ByteCount, CUstream hStream);
# 4058 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpyAtoHAsync_v2(void *dstHost, CUarray srcArray, size_t srcOffset, size_t ByteCount, CUstream hStream);
# 4229 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy2DAsync_v2(const CUDA_MEMCPY2D *pCopy, CUstream hStream);
# 4404 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy3DAsync_v2(const CUDA_MEMCPY3D *pCopy, CUstream hStream);
# 4429 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemcpy3DPeerAsync(const CUDA_MEMCPY3D_PEER *pCopy, CUstream hStream);
# 4464 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD8_v2(CUdeviceptr dstDevice, unsigned char uc, size_t N);
# 4497 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD16_v2(CUdeviceptr dstDevice, unsigned short us, size_t N);
# 4530 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD32_v2(CUdeviceptr dstDevice, unsigned int ui, size_t N);
# 4568 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD2D8_v2(CUdeviceptr dstDevice, size_t dstPitch, unsigned char uc, size_t Width, size_t Height);
# 4606 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD2D16_v2(CUdeviceptr dstDevice, size_t dstPitch, unsigned short us, size_t Width, size_t Height);
# 4644 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD2D32_v2(CUdeviceptr dstDevice, size_t dstPitch, unsigned int ui, size_t Width, size_t Height);
# 4681 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD8Async(CUdeviceptr dstDevice, unsigned char uc, size_t N, CUstream hStream);
# 4718 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD16Async(CUdeviceptr dstDevice, unsigned short us, size_t N, CUstream hStream);
# 4754 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD32Async(CUdeviceptr dstDevice, unsigned int ui, size_t N, CUstream hStream);
# 4796 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD2D8Async(CUdeviceptr dstDevice, size_t dstPitch, unsigned char uc, size_t Width, size_t Height, CUstream hStream);
# 4838 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD2D16Async(CUdeviceptr dstDevice, size_t dstPitch, unsigned short us, size_t Width, size_t Height, CUstream hStream);
# 4880 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuMemsetD2D32Async(CUdeviceptr dstDevice, size_t dstPitch, unsigned int ui, size_t Width, size_t Height, CUstream hStream);
# 4983 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuArrayCreate_v2(CUarray *pHandle, const CUDA_ARRAY_DESCRIPTOR *pAllocateArray);
# 5016 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuArrayGetDescriptor_v2(CUDA_ARRAY_DESCRIPTOR *pArrayDescriptor, CUarray hArray);
# 5047 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuArrayDestroy(CUarray hArray);
# 5155 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuArray3DCreate_v2(CUarray *pHandle, const CUDA_ARRAY3D_DESCRIPTOR *pAllocateArray);
# 5191 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuArray3DGetDescriptor_v2(CUDA_ARRAY3D_DESCRIPTOR *pArrayDescriptor, CUarray hArray);
# 5398 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuPointerGetAttribute(void *data, CUpointer_attribute attribute, CUdeviceptr ptr);
# 5435 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuStreamCreate(CUstream *phStream, unsigned int Flags);
# 5477 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuStreamWaitEvent(CUstream hStream, CUevent hEvent, unsigned int Flags);
# 5501 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuStreamQuery(CUstream hStream);
# 5526 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuStreamSynchronize(CUstream hStream);
# 5554 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuStreamDestroy_v2(CUstream hStream);
# 5603 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuEventCreate(CUevent *phEvent, unsigned int Flags);
# 5641 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuEventRecord(CUevent hEvent, CUstream hStream);
# 5672 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuEventQuery(CUevent hEvent);
# 5706 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuEventSynchronize(CUevent hEvent);
# 5735 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuEventDestroy_v2(CUevent hEvent);
# 5779 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuEventElapsedTime(float *pMilliseconds, CUevent hStart, CUevent hEnd);
# 5842 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuFuncGetAttribute(int *pi, CUfunction_attribute attrib, CUfunction hfunc);
# 5883 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuFuncSetCacheConfig(CUfunction hfunc, CUfunc_cache config);
# 5999 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuLaunchKernel(CUfunction f,
                                unsigned int gridDimX,
                                unsigned int gridDimY,
                                unsigned int gridDimZ,
                                unsigned int blockDimX,
                                unsigned int blockDimY,
                                unsigned int blockDimZ,
                                unsigned int sharedMemBytes,
                                CUstream hStream,
                                void **kernelParams,
                                void **extra);
# 6055 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuFuncSetBlockShape(CUfunction hfunc, int x, int y, int z);
# 6089 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuFuncSetSharedSize(CUfunction hfunc, unsigned int bytes);
# 6121 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuParamSetSize(CUfunction hfunc, unsigned int numbytes);
# 6154 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuParamSeti(CUfunction hfunc, int offset, unsigned int value);
# 6187 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuParamSetf(CUfunction hfunc, int offset, float value);
# 6222 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuParamSetv(CUfunction hfunc, int offset, void *ptr, unsigned int numbytes);
# 6259 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuLaunch(CUfunction f);
# 6298 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuLaunchGrid(CUfunction f, int grid_width, int grid_height);
# 6342 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuLaunchGridAsync(CUfunction f, int grid_width, int grid_height, CUstream hStream);
# 6367 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuParamSetTexRef(CUfunction hfunc, int texunit, CUtexref hTexRef);
# 6408 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetArray(CUtexref hTexRef, CUarray hArray, unsigned int Flags);
# 6446 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetAddress_v2(size_t *ByteOffset, CUtexref hTexRef, CUdeviceptr dptr, size_t bytes);
# 6487 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetAddress2D_v2(CUtexref hTexRef, const CUDA_ARRAY_DESCRIPTOR *desc, CUdeviceptr dptr, size_t Pitch);
# 6516 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetFormat(CUtexref hTexRef, CUarray_format fmt, int NumPackedComponents);
# 6556 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetAddressMode(CUtexref hTexRef, int dim, CUaddress_mode am);
# 6589 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetFilterMode(CUtexref hTexRef, CUfilter_mode fm);
# 6621 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefSetFlags(CUtexref hTexRef, unsigned int Flags);
# 6647 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefGetAddress_v2(CUdeviceptr *pdptr, CUtexref hTexRef);
# 6673 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefGetArray(CUarray *phArray, CUtexref hTexRef);
# 6699 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefGetAddressMode(CUaddress_mode *pam, CUtexref hTexRef, int dim);
# 6723 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefGetFilterMode(CUfilter_mode *pfm, CUtexref hTexRef);
# 6749 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefGetFormat(CUarray_format *pFormat, int *pNumChannels, CUtexref hTexRef);
# 6772 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefGetFlags(unsigned int *pFlags, CUtexref hTexRef);
# 6806 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefCreate(CUtexref *pTexRef);
# 6826 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuTexRefDestroy(CUtexref hTexRef);
# 6864 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuSurfRefSetArray(CUsurfref hSurfRef, CUarray hArray, unsigned int Flags);
# 6885 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuSurfRefGetArray(CUarray *phArray, CUsurfref hSurfRef);
# 6923 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuDeviceCanAccessPeer(int *canAccessPeer, CUdevice dev, CUdevice peerDev);
# 6966 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxEnablePeerAccess(CUcontext peerContext, unsigned int Flags);
# 6991 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuCtxDisablePeerAccess(CUcontext peerContext);
# 7032 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuGraphicsUnregisterResource(CUgraphicsResource resource);
# 7070 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuGraphicsSubResourceGetMappedArray(CUarray *pArray, CUgraphicsResource resource, unsigned int arrayIndex, unsigned int mipLevel);
# 7104 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuGraphicsResourceGetMappedPointer_v2(CUdeviceptr *pDevPtr, size_t *pSize, CUgraphicsResource resource);
# 7145 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuGraphicsResourceSetMapFlags(CUgraphicsResource resource, unsigned int flags);
# 7183 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuGraphicsMapResources(unsigned int count, CUgraphicsResource *resources, CUstream hStream);
# 7218 "/home/dpnkarthik/softwares/cuda/include/cuda.h"
CUresult cuGraphicsUnmapResources(unsigned int count, CUgraphicsResource *resources, CUstream hStream);



CUresult cuGetExportTable(const void **ppExportTable, const CUuuid *pExportTableId);
# 7 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h" 2
# 1 "/home/dpnkarthik/softwares/cuda/include/cusparse.h" 1
# 65 "/home/dpnkarthik/softwares/cuda/include/cusparse.h"
# 1 "/home/dpnkarthik/softwares/cuda/include/driver_types.h" 1
# 74 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/limits.h" 1 3 4
# 75 "/home/dpnkarthik/softwares/cuda/include/driver_types.h" 2
# 1 "/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include/stddef.h" 1 3 4
# 76 "/home/dpnkarthik/softwares/cuda/include/driver_types.h" 2
# 118 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
enum cudaError
{





  cudaSuccess = 0,





  cudaErrorMissingConfiguration = 1,





  cudaErrorMemoryAllocation = 2,





  cudaErrorInitializationError = 3,
# 153 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorLaunchFailure = 4,
# 162 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorPriorLaunchFailure = 5,
# 172 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorLaunchTimeout = 6,
# 181 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorLaunchOutOfResources = 7,





  cudaErrorInvalidDeviceFunction = 8,
# 196 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorInvalidConfiguration = 9,





  cudaErrorInvalidDevice = 10,





  cudaErrorInvalidValue = 11,





  cudaErrorInvalidPitchValue = 12,





  cudaErrorInvalidSymbol = 13,




  cudaErrorMapBufferObjectFailed = 14,




  cudaErrorUnmapBufferObjectFailed = 15,





  cudaErrorInvalidHostPointer = 16,





  cudaErrorInvalidDevicePointer = 17,





  cudaErrorInvalidTexture = 18,





  cudaErrorInvalidTextureBinding = 19,






  cudaErrorInvalidChannelDescriptor = 20,





  cudaErrorInvalidMemcpyDirection = 21,
# 277 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorAddressOfConstant = 22,
# 286 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorTextureFetchFailed = 23,
# 295 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorTextureNotBound = 24,
# 304 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorSynchronizationError = 25,





  cudaErrorInvalidFilterSetting = 26,





  cudaErrorInvalidNormSetting = 27,







  cudaErrorMixedDeviceExecution = 28,







  cudaErrorCudartUnloading = 29,




  cudaErrorUnknown = 30,





  cudaErrorNotYetImplemented = 31,
# 352 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorMemoryValueTooLarge = 32,






  cudaErrorInvalidResourceHandle = 33,







  cudaErrorNotReady = 34,






  cudaErrorInsufficientDriver = 35,
# 387 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorSetOnActiveProcess = 36,





  cudaErrorInvalidSurface = 37,





  cudaErrorNoDevice = 38,





  cudaErrorECCUncorrectable = 39,




  cudaErrorSharedObjectSymbolNotFound = 40,




  cudaErrorSharedObjectInitFailed = 41,





  cudaErrorUnsupportedLimit = 42,





  cudaErrorDuplicateVariableName = 43,





  cudaErrorDuplicateTextureName = 44,





  cudaErrorDuplicateSurfaceName = 45,
# 449 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorDevicesUnavailable = 46,




  cudaErrorInvalidKernelImage = 47,







  cudaErrorNoKernelImageForDevice = 48,
# 475 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
  cudaErrorIncompatibleDriverContext = 49,






  cudaErrorPeerAccessAlreadyEnabled = 50,






  cudaErrorPeerAccessNotEnabled = 51,





  cudaErrorDeviceAlreadyInUse = 54,







  cudaErrorProfilerDisabled = 55,






  cudaErrorProfilerNotInitialized = 56,






  cudaErrorProfilerAlreadyStarted = 57,





   cudaErrorProfilerAlreadyStopped = 58,




  cudaErrorStartupFailure = 0x7f,





  cudaErrorApiFailureBase = 10000
};





enum cudaChannelFormatKind
{
  cudaChannelFormatKindSigned = 0,
  cudaChannelFormatKindUnsigned = 1,
  cudaChannelFormatKindFloat = 2,
  cudaChannelFormatKindNone = 3
};





struct cudaChannelFormatDesc
{
  int x;
  int y;
  int z;
  int w;
  enum cudaChannelFormatKind f;
};





struct cudaArray;





enum cudaMemoryType
{
  cudaMemoryTypeHost = 1,
  cudaMemoryTypeDevice = 2
};





enum cudaMemcpyKind
{
  cudaMemcpyHostToHost = 0,
  cudaMemcpyHostToDevice = 1,
  cudaMemcpyDeviceToHost = 2,
  cudaMemcpyDeviceToDevice = 3,
  cudaMemcpyDefault = 4
};






struct cudaPitchedPtr
{
  void *ptr;
  size_t pitch;
  size_t xsize;
  size_t ysize;
};






struct cudaExtent
{
  size_t width;
  size_t height;
  size_t depth;
};






struct cudaPos
{
  size_t x;
  size_t y;
  size_t z;
};





struct cudaMemcpy3DParms
{
  struct cudaArray *srcArray;
  struct cudaPos srcPos;
  struct cudaPitchedPtr srcPtr;

  struct cudaArray *dstArray;
  struct cudaPos dstPos;
  struct cudaPitchedPtr dstPtr;

  struct cudaExtent extent;
  enum cudaMemcpyKind kind;
};





struct cudaMemcpy3DPeerParms
{
  struct cudaArray *srcArray;
  struct cudaPos srcPos;
  struct cudaPitchedPtr srcPtr;
  int srcDevice;

  struct cudaArray *dstArray;
  struct cudaPos dstPos;
  struct cudaPitchedPtr dstPtr;
  int dstDevice;

  struct cudaExtent extent;
};





struct cudaGraphicsResource;





enum cudaGraphicsRegisterFlags
{
  cudaGraphicsRegisterFlagsNone = 0,
  cudaGraphicsRegisterFlagsReadOnly = 1,
  cudaGraphicsRegisterFlagsWriteDiscard = 2,
  cudaGraphicsRegisterFlagsSurfaceLoadStore = 4
};





enum cudaGraphicsMapFlags
{
  cudaGraphicsMapFlagsNone = 0,
  cudaGraphicsMapFlagsReadOnly = 1,
  cudaGraphicsMapFlagsWriteDiscard = 2
};





enum cudaGraphicsCubeFace {
  cudaGraphicsCubeFacePositiveX = 0x00,
  cudaGraphicsCubeFaceNegativeX = 0x01,
  cudaGraphicsCubeFacePositiveY = 0x02,
  cudaGraphicsCubeFaceNegativeY = 0x03,
  cudaGraphicsCubeFacePositiveZ = 0x04,
  cudaGraphicsCubeFaceNegativeZ = 0x05
};





struct cudaPointerAttributes
{




    enum cudaMemoryType memoryType;
# 728 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
    int device;





    void *devicePointer;





    void *hostPointer;
};





struct cudaFuncAttributes
{





   size_t sharedSizeBytes;





   size_t constSizeBytes;




   size_t localSizeBytes;






   int maxThreadsPerBlock;




   int numRegs;






   int ptxVersion;






   int binaryVersion;
};





enum cudaFuncCache
{
  cudaFuncCachePreferNone = 0,
  cudaFuncCachePreferShared = 1,
  cudaFuncCachePreferL1 = 2
};





enum cudaComputeMode
{
  cudaComputeModeDefault = 0,
  cudaComputeModeExclusive = 1,
  cudaComputeModeProhibited = 2,
  cudaComputeModeExclusiveProcess = 3
};





enum cudaLimit
{
    cudaLimitStackSize = 0x00,
    cudaLimitPrintfFifoSize = 0x01,
    cudaLimitMallocHeapSize = 0x02
};





enum cudaOutputMode
{
    cudaKeyValuePair = 0x00,
    cudaCSV = 0x01
};





struct cudaDeviceProp
{
  char name[256];
  size_t totalGlobalMem;
  size_t sharedMemPerBlock;
  int regsPerBlock;
  int warpSize;
  size_t memPitch;
  int maxThreadsPerBlock;
  int maxThreadsDim[3];
  int maxGridSize[3];
  int clockRate;
  size_t totalConstMem;
  int major;
  int minor;
  size_t textureAlignment;
  int deviceOverlap;
  int multiProcessorCount;
  int kernelExecTimeoutEnabled;
  int integrated;
  int canMapHostMemory;
  int computeMode;
  int maxTexture1D;
  int maxTexture2D[2];
  int maxTexture3D[3];
  int maxTexture1DLayered[2];
  int maxTexture2DLayered[3];
  size_t surfaceAlignment;
  int concurrentKernels;
  int ECCEnabled;
  int pciBusID;
  int pciDeviceID;
  int pciDomainID;
  int tccDriver;
  int asyncEngineCount;
  int unifiedAddressing;
  int memoryClockRate;
  int memoryBusWidth;
  int l2CacheSize;
  int maxThreadsPerMultiProcessor;
};
# 936 "/home/dpnkarthik/softwares/cuda/include/driver_types.h"
typedef enum cudaError cudaError_t;





typedef struct CUstream_st *cudaStream_t;





typedef struct CUevent_st *cudaEvent_t;





typedef struct cudaGraphicsResource *cudaGraphicsResource_t;





typedef struct CUuuid_st cudaUUID_t;





typedef enum cudaOutputMode cudaOutputMode_t;
# 66 "/home/dpnkarthik/softwares/cuda/include/cusparse.h" 2
# 1 "/home/dpnkarthik/softwares/cuda/include/cuComplex.h" 1
# 58 "/home/dpnkarthik/softwares/cuda/include/cuComplex.h"
# 1 "/home/dpnkarthik/softwares/cuda/include/vector_types.h" 1
# 59 "/home/dpnkarthik/softwares/cuda/include/vector_types.h"
# 1 "/home/dpnkarthik/softwares/cuda/include/builtin_types.h" 1
# 56 "/home/dpnkarthik/softwares/cuda/include/builtin_types.h"
# 1 "/home/dpnkarthik/softwares/cuda/include/device_types.h" 1
# 60 "/home/dpnkarthik/softwares/cuda/include/device_types.h"
enum cudaRoundMode
{
  cudaRoundNearest,
  cudaRoundZero,
  cudaRoundPosInf,
  cudaRoundMinInf
};
# 57 "/home/dpnkarthik/softwares/cuda/include/builtin_types.h" 2

# 1 "/home/dpnkarthik/softwares/cuda/include/surface_types.h" 1
# 77 "/home/dpnkarthik/softwares/cuda/include/surface_types.h"
enum cudaSurfaceBoundaryMode
{
  cudaBoundaryModeZero = 0,
  cudaBoundaryModeClamp = 1,
  cudaBoundaryModeTrap = 2
};





enum cudaSurfaceFormatMode
{
  cudaFormatModeForced = 0,
  cudaFormatModeAuto = 1
};





struct surfaceReference
{



  struct cudaChannelFormatDesc channelDesc;
};
# 59 "/home/dpnkarthik/softwares/cuda/include/builtin_types.h" 2
# 1 "/home/dpnkarthik/softwares/cuda/include/texture_types.h" 1
# 83 "/home/dpnkarthik/softwares/cuda/include/texture_types.h"
enum cudaTextureAddressMode
{
  cudaAddressModeWrap = 0,
  cudaAddressModeClamp = 1,
  cudaAddressModeMirror = 2,
  cudaAddressModeBorder = 3
};





enum cudaTextureFilterMode
{
  cudaFilterModePoint = 0,
  cudaFilterModeLinear = 1
};





enum cudaTextureReadMode
{
  cudaReadModeElementType = 0,
  cudaReadModeNormalizedFloat = 1
};





struct textureReference
{



  int normalized;



  enum cudaTextureFilterMode filterMode;



  enum cudaTextureAddressMode addressMode[3];



  struct cudaChannelFormatDesc channelDesc;



  int sRGB;
  int __cudaReserved[15];
};
# 60 "/home/dpnkarthik/softwares/cuda/include/builtin_types.h" 2
# 1 "/home/dpnkarthik/softwares/cuda/include/vector_types.h" 1
# 60 "/home/dpnkarthik/softwares/cuda/include/builtin_types.h" 2
# 60 "/home/dpnkarthik/softwares/cuda/include/vector_types.h" 2
# 1 "/home/dpnkarthik/softwares/cuda/include/host_defines.h" 1
# 61 "/home/dpnkarthik/softwares/cuda/include/vector_types.h" 2
# 92 "/home/dpnkarthik/softwares/cuda/include/vector_types.h"
struct char1
{
  signed char x;
};


struct uchar1
{
  unsigned char x;
};


struct __attribute__((aligned(2))) char2
{
  signed char x, y;
};


struct __attribute__((aligned(2))) uchar2
{
  unsigned char x, y;
};


struct char3
{
  signed char x, y, z;
};


struct uchar3
{
  unsigned char x, y, z;
};


struct __attribute__((aligned(4))) char4
{
  signed char x, y, z, w;
};


struct __attribute__((aligned(4))) uchar4
{
  unsigned char x, y, z, w;
};


struct short1
{
  short x;
};


struct ushort1
{
  unsigned short x;
};


struct __attribute__((aligned(4))) short2
{
  short x, y;
};


struct __attribute__((aligned(4))) ushort2
{
  unsigned short x, y;
};


struct short3
{
  short x, y, z;
};


struct ushort3
{
  unsigned short x, y, z;
};


struct __attribute__((aligned(8))) short4 { short x; short y; short z; short w; };


struct __attribute__((aligned(8))) ushort4 { unsigned short x; unsigned short y; unsigned short z; unsigned short w; };


struct int1
{
  int x;
};


struct uint1
{
  unsigned int x;
};


struct __attribute__((aligned(8))) int2 { int x; int y; };


struct __attribute__((aligned(8))) uint2 { unsigned int x; unsigned int y; };


struct int3
{
  int x, y, z;
};


struct uint3
{
  unsigned int x, y, z;
};


struct __attribute__((aligned(16))) int4
{
  int x, y, z, w;
};


struct __attribute__((aligned(16))) uint4
{
  unsigned int x, y, z, w;
};


struct long1
{
  long int x;
};


struct ulong1
{
  unsigned long x;
};
# 246 "/home/dpnkarthik/softwares/cuda/include/vector_types.h"
struct __attribute__((aligned(2*sizeof(long int)))) long2
{
  long int x, y;
};


struct __attribute__((aligned(2*sizeof(unsigned long int)))) ulong2
{
  unsigned long int x, y;
};




struct long3
{
  long int x, y, z;
};


struct ulong3
{
  unsigned long int x, y, z;
};


struct __attribute__((aligned(16))) long4
{
  long int x, y, z, w;
};


struct __attribute__((aligned(16))) ulong4
{
  unsigned long int x, y, z, w;
};


struct float1
{
  float x;
};


struct __attribute__((aligned(8))) float2 { float x; float y; };


struct float3
{
  float x, y, z;
};


struct __attribute__((aligned(16))) float4
{
  float x, y, z, w;
};


struct longlong1
{
  long long int x;
};


struct ulonglong1
{
  unsigned long long int x;
};


struct __attribute__((aligned(16))) longlong2
{
  long long int x, y;
};


struct __attribute__((aligned(16))) ulonglong2
{
  unsigned long long int x, y;
};


struct longlong3
{
  long long int x, y, z;
};


struct ulonglong3
{
  unsigned long long int x, y, z;
};


struct __attribute__((aligned(16))) longlong4
{
  long long int x, y, z ,w;
};


struct __attribute__((aligned(16))) ulonglong4
{
  unsigned long long int x, y, z, w;
};


struct double1
{
  double x;
};


struct __attribute__((aligned(16))) double2
{
  double x, y;
};


struct double3
{
  double x, y, z;
};


struct __attribute__((aligned(16))) double4
{
  double x, y, z, w;
};
# 390 "/home/dpnkarthik/softwares/cuda/include/vector_types.h"
typedef struct char1 char1;

typedef struct uchar1 uchar1;

typedef struct char2 char2;

typedef struct uchar2 uchar2;

typedef struct char3 char3;

typedef struct uchar3 uchar3;

typedef struct char4 char4;

typedef struct uchar4 uchar4;

typedef struct short1 short1;

typedef struct ushort1 ushort1;

typedef struct short2 short2;

typedef struct ushort2 ushort2;

typedef struct short3 short3;

typedef struct ushort3 ushort3;

typedef struct short4 short4;

typedef struct ushort4 ushort4;

typedef struct int1 int1;

typedef struct uint1 uint1;

typedef struct int2 int2;

typedef struct uint2 uint2;

typedef struct int3 int3;

typedef struct uint3 uint3;

typedef struct int4 int4;

typedef struct uint4 uint4;

typedef struct long1 long1;

typedef struct ulong1 ulong1;

typedef struct long2 long2;

typedef struct ulong2 ulong2;

typedef struct long3 long3;

typedef struct ulong3 ulong3;

typedef struct long4 long4;

typedef struct ulong4 ulong4;

typedef struct float1 float1;

typedef struct float2 float2;

typedef struct float3 float3;

typedef struct float4 float4;

typedef struct longlong1 longlong1;

typedef struct ulonglong1 ulonglong1;

typedef struct longlong2 longlong2;

typedef struct ulonglong2 ulonglong2;

typedef struct longlong3 longlong3;

typedef struct ulonglong3 ulonglong3;

typedef struct longlong4 longlong4;

typedef struct ulonglong4 ulonglong4;

typedef struct double1 double1;

typedef struct double2 double2;

typedef struct double3 double3;

typedef struct double4 double4;
# 493 "/home/dpnkarthik/softwares/cuda/include/vector_types.h"
struct dim3
{
    unsigned int x, y, z;





};


typedef struct dim3 dim3;
# 59 "/home/dpnkarthik/softwares/cuda/include/cuComplex.h" 2

typedef float2 cuFloatComplex;

 static __inline__ float cuCrealf (cuFloatComplex x)
{
    return x.x;
}

 static __inline__ float cuCimagf (cuFloatComplex x)
{
    return x.y;
}

 static __inline__ cuFloatComplex make_cuFloatComplex
                                                             (float r, float i)
{
    cuFloatComplex res;
    res.x = r;
    res.y = i;
    return res;
}

 static __inline__ cuFloatComplex cuConjf (cuFloatComplex x)
{
    return make_cuFloatComplex (cuCrealf(x), -cuCimagf(x));
}
 static __inline__ cuFloatComplex cuCaddf (cuFloatComplex x,
                                                              cuFloatComplex y)
{
    return make_cuFloatComplex (cuCrealf(x) + cuCrealf(y),
                                cuCimagf(x) + cuCimagf(y));
}

 static __inline__ cuFloatComplex cuCsubf (cuFloatComplex x,
                                                              cuFloatComplex y)
{
        return make_cuFloatComplex (cuCrealf(x) - cuCrealf(y),
                                    cuCimagf(x) - cuCimagf(y));
}






 static __inline__ cuFloatComplex cuCmulf (cuFloatComplex x,
                                                              cuFloatComplex y)
{
    cuFloatComplex prod;
    prod = make_cuFloatComplex ((cuCrealf(x) * cuCrealf(y)) -
                                 (cuCimagf(x) * cuCimagf(y)),
                                 (cuCrealf(x) * cuCimagf(y)) +
                                 (cuCimagf(x) * cuCrealf(y)));
    return prod;
}






 static __inline__ cuFloatComplex cuCdivf (cuFloatComplex x,
                                                              cuFloatComplex y)
{
    cuFloatComplex quot;
    float s = fabsf(cuCrealf(y)) + fabsf(cuCimagf(y));
    float oos = 1.0f / s;
    float ars = cuCrealf(x) * oos;
    float ais = cuCimagf(x) * oos;
    float brs = cuCrealf(y) * oos;
    float bis = cuCimagf(y) * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0f / s;
    quot = make_cuFloatComplex (((ars * brs) + (ais * bis)) * oos,
                                ((ais * brs) - (ars * bis)) * oos);
    return quot;
}
# 145 "/home/dpnkarthik/softwares/cuda/include/cuComplex.h"
 static __inline__ float cuCabsf (cuFloatComplex x)
{
    float a = cuCrealf(x);
    float b = cuCimagf(x);
    float v, w, t;
    a = fabsf(a);
    b = fabsf(b);
    if (a > b) {
        v = a;
        w = b;
    } else {
        v = b;
        w = a;
    }
    t = w / v;
    t = 1.0f + t * t;
    t = v * sqrtf(t);
    if ((v == 0.0f) || (v > 3.402823466e38f) || (w > 3.402823466e38f)) {
        t = v + w;
    }
    return t;
}


typedef double2 cuDoubleComplex;

 static __inline__ double cuCreal (cuDoubleComplex x)
{
    return x.x;
}

 static __inline__ double cuCimag (cuDoubleComplex x)
{
    return x.y;
}

 static __inline__ cuDoubleComplex make_cuDoubleComplex
                                                           (double r, double i)
{
    cuDoubleComplex res;
    res.x = r;
    res.y = i;
    return res;
}

 static __inline__ cuDoubleComplex cuConj(cuDoubleComplex x)
{
    return make_cuDoubleComplex (cuCreal(x), -cuCimag(x));
}

 static __inline__ cuDoubleComplex cuCadd(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
    return make_cuDoubleComplex (cuCreal(x) + cuCreal(y),
                                 cuCimag(x) + cuCimag(y));
}

 static __inline__ cuDoubleComplex cuCsub(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
    return make_cuDoubleComplex (cuCreal(x) - cuCreal(y),
                                 cuCimag(x) - cuCimag(y));
}






 static __inline__ cuDoubleComplex cuCmul(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
    cuDoubleComplex prod;
    prod = make_cuDoubleComplex ((cuCreal(x) * cuCreal(y)) -
                                 (cuCimag(x) * cuCimag(y)),
                                 (cuCreal(x) * cuCimag(y)) +
                                 (cuCimag(x) * cuCreal(y)));
    return prod;
}






 static __inline__ cuDoubleComplex cuCdiv(cuDoubleComplex x,
                                                             cuDoubleComplex y)
{
    cuDoubleComplex quot;
    double s = (fabs(cuCreal(y))) + (fabs(cuCimag(y)));
    double oos = 1.0 / s;
    double ars = cuCreal(x) * oos;
    double ais = cuCimag(x) * oos;
    double brs = cuCreal(y) * oos;
    double bis = cuCimag(y) * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0 / s;
    quot = make_cuDoubleComplex (((ars * brs) + (ais * bis)) * oos,
                                 ((ais * brs) - (ars * bis)) * oos);
    return quot;
}







 static __inline__ double cuCabs (cuDoubleComplex x)
{
    double a = cuCreal(x);
    double b = cuCimag(x);
    double v, w, t;
    a = fabs(a);
    b = fabs(b);
    if (a > b) {
        v = a;
        w = b;
    } else {
        v = b;
        w = a;
    }
    t = w / v;
    t = 1.0 + t * t;
    t = v * sqrt(t);
    if ((v == 0.0) ||
        (v > 1.79769313486231570e+308) || (w > 1.79769313486231570e+308)) {
        t = v + w;
    }
    return t;
}






typedef cuFloatComplex cuComplex;
 static __inline__ cuComplex make_cuComplex (float x,
                                                                float y)
{
    return make_cuFloatComplex (x, y);
}


 static __inline__ cuDoubleComplex cuComplexFloatToDouble
                                                      (cuFloatComplex c)
{
    return make_cuDoubleComplex ((double)cuCrealf(c), (double)cuCimagf(c));
}

 static __inline__ cuFloatComplex cuComplexDoubleToFloat
(cuDoubleComplex c)
{
 return make_cuFloatComplex ((float)cuCreal(c), (float)cuCimag(c));
}


 static __inline__ cuComplex cuCfmaf( cuComplex x, cuComplex y, cuComplex d)
{
    float real_res;
    float imag_res;

    real_res = (cuCrealf(x) * cuCrealf(y)) + cuCrealf(d);
    imag_res = (cuCrealf(x) * cuCimagf(y)) + cuCimagf(d);

    real_res = -(cuCimagf(x) * cuCimagf(y)) + real_res;
    imag_res = (cuCimagf(x) * cuCrealf(y)) + imag_res;

    return make_cuComplex(real_res, imag_res);
}

 static __inline__ cuDoubleComplex cuCfma( cuDoubleComplex x, cuDoubleComplex y, cuDoubleComplex d)
{
    double real_res;
    double imag_res;

    real_res = (cuCreal(x) * cuCreal(y)) + cuCreal(d);
    imag_res = (cuCreal(x) * cuCimag(y)) + cuCimag(d);

    real_res = -(cuCimag(x) * cuCimag(y)) + real_res;
    imag_res = (cuCimag(x) * cuCreal(y)) + imag_res;

    return make_cuDoubleComplex(real_res, imag_res);
}
# 67 "/home/dpnkarthik/softwares/cuda/include/cusparse.h" 2


typedef enum{
    CUSPARSE_STATUS_SUCCESS=0,
    CUSPARSE_STATUS_NOT_INITIALIZED=1,
    CUSPARSE_STATUS_ALLOC_FAILED=2,
    CUSPARSE_STATUS_INVALID_VALUE=3,
    CUSPARSE_STATUS_ARCH_MISMATCH=4,
    CUSPARSE_STATUS_MAPPING_ERROR=5,
    CUSPARSE_STATUS_EXECUTION_FAILED=6,
    CUSPARSE_STATUS_INTERNAL_ERROR=7,
    CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED=8
} cusparseStatus_t;


struct cusparseContext;
typedef struct cusparseContext *cusparseHandle_t;

cusparseStatus_t cusparseCreate (cusparseHandle_t *handle);
cusparseStatus_t cusparseDestroy (cusparseHandle_t handle);
cusparseStatus_t cusparseGetVersion(cusparseHandle_t handle, int *version);
cusparseStatus_t cusparseSetKernelStream (cusparseHandle_t handle, cudaStream_t streamId);
# 103 "/home/dpnkarthik/softwares/cuda/include/cusparse.h"
typedef enum {
    CUSPARSE_MATRIX_TYPE_GENERAL=0,
    CUSPARSE_MATRIX_TYPE_SYMMETRIC=1,
    CUSPARSE_MATRIX_TYPE_HERMITIAN=2,
    CUSPARSE_MATRIX_TYPE_TRIANGULAR=3
} cusparseMatrixType_t;

typedef enum {
    CUSPARSE_FILL_MODE_LOWER=0,
    CUSPARSE_FILL_MODE_UPPER=1
} cusparseFillMode_t;

typedef enum {
    CUSPARSE_DIAG_TYPE_NON_UNIT=0,
    CUSPARSE_DIAG_TYPE_UNIT=1
} cusparseDiagType_t;

typedef enum {
    CUSPARSE_INDEX_BASE_ZERO=0,
    CUSPARSE_INDEX_BASE_ONE=1
} cusparseIndexBase_t;

typedef enum {
    CUSPARSE_OPERATION_NON_TRANSPOSE=0,
    CUSPARSE_OPERATION_TRANSPOSE=1,
    CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE=2
} cusparseOperation_t;

typedef enum {
    CUSPARSE_DIRECTION_ROW=0,
    CUSPARSE_DIRECTION_COLUMN=1
} cusparseDirection_t;


struct cusparseMatDescr;
typedef struct cusparseMatDescr *cusparseMatDescr_t;



cusparseStatus_t cusparseCreateMatDescr(cusparseMatDescr_t *descrA);
cusparseStatus_t cusparseDestroyMatDescr (cusparseMatDescr_t descrA);

cusparseStatus_t cusparseSetMatType(cusparseMatDescr_t descrA, cusparseMatrixType_t type);
cusparseMatrixType_t cusparseGetMatType(const cusparseMatDescr_t descrA);


cusparseStatus_t cusparseSetMatFillMode(cusparseMatDescr_t descrA, cusparseFillMode_t fillMode);
cusparseFillMode_t cusparseGetMatFillMode(const cusparseMatDescr_t descrA);


cusparseStatus_t cusparseSetMatDiagType(cusparseMatDescr_t descrA, cusparseDiagType_t diagType);
cusparseDiagType_t cusparseGetMatDiagType(const cusparseMatDescr_t descrA);

cusparseStatus_t cusparseSetMatIndexBase(cusparseMatDescr_t descrA, cusparseIndexBase_t base);
cusparseIndexBase_t cusparseGetMatIndexBase(const cusparseMatDescr_t descrA);



struct cusparseSolveAnalysisInfo;
typedef struct cusparseSolveAnalysisInfo *cusparseSolveAnalysisInfo_t;

cusparseStatus_t cusparseCreateSolveAnalysisInfo(cusparseSolveAnalysisInfo_t *info);

cusparseStatus_t cusparseDestroySolveAnalysisInfo(cusparseSolveAnalysisInfo_t info);






cusparseStatus_t cusparseSaxpyi(cusparseHandle_t handle,
                                            int nnz,
                                            float alpha,
                                            const float *xVal,
                                            const int *xInd,
                                            float *y,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDaxpyi(cusparseHandle_t handle,
                                            int nnz,
                                            double alpha,
                                            const double *xVal,
                                            const int *xInd,
                                            double *y,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseCaxpyi(cusparseHandle_t handle,
                                            int nnz,
                                            cuComplex alpha,
                                            const cuComplex *xVal,
                                            const int *xInd,
                                            cuComplex *y,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZaxpyi(cusparseHandle_t handle,
                                            int nnz,
                                            cuDoubleComplex alpha,
                                            const cuDoubleComplex *xVal,
                                            const int *xInd,
                                            cuDoubleComplex *y,
                                            cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseSdoti(cusparseHandle_t handle,
                                           int nnz,
                                           const float *xVal,
                                           const int *xInd,
                                           const float *y,
                                           float *resultHostPtr,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDdoti(cusparseHandle_t handle,
                                           int nnz,
                                           const double *xVal,
                                           const int *xInd,
                                           const double *y,
                                           double *resultHostPtr,
                                           cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseCdoti(cusparseHandle_t handle,
                                            int nnz,
                                            const cuComplex *xVal,
                                            const int *xInd,
                                            const cuComplex *y,
                                            cuComplex *resultHostPtr,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZdoti(cusparseHandle_t handle,
                                            int nnz,
                                            const cuDoubleComplex *xVal,
                                            const int *xInd,
                                            const cuDoubleComplex *y,
                                            cuDoubleComplex *resultHostPtr,
                                            cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseCdotci(cusparseHandle_t handle,
                                            int nnz,
                                            const cuComplex *xVal,
                                            const int *xInd,
                                            const cuComplex *y,
                                            cuComplex *resultHostPtr,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZdotci(cusparseHandle_t handle,
                                            int nnz,
                                            const cuDoubleComplex *xVal,
                                            const int *xInd,
                                            const cuDoubleComplex *y,
                                            cuDoubleComplex *resultHostPtr,
                                            cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseSgthr(cusparseHandle_t handle,
                                           int nnz,
                                           const float *y,
                                           float *xVal,
                                           const int *xInd,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDgthr(cusparseHandle_t handle,
                                           int nnz,
                                           const double *y,
                                           double *xVal,
                                           const int *xInd,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseCgthr(cusparseHandle_t handle,
                                           int nnz,
                                           const cuComplex *y,
                                           cuComplex *xVal,
                                           const int *xInd,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZgthr(cusparseHandle_t handle,
                                           int nnz,
                                           const cuDoubleComplex *y,
                                           cuDoubleComplex *xVal,
                                           const int *xInd,
                                           cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseSgthrz(cusparseHandle_t handle,
                                            int nnz,
                                            float *y,
                                            float *xVal,
                                            const int *xInd,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDgthrz(cusparseHandle_t handle,
                                            int nnz,
                                            double *y,
                                            double *xVal,
                                            const int *xInd,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseCgthrz(cusparseHandle_t handle,
                                            int nnz,
                                            cuComplex *y,
                                            cuComplex *xVal,
                                            const int *xInd,
                                            cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZgthrz(cusparseHandle_t handle,
                                            int nnz,
                                            cuDoubleComplex *y,
                                            cuDoubleComplex *xVal,
                                            const int *xInd,
                                            cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseSsctr(cusparseHandle_t handle,
                                           int nnz,
                                           const float *xVal,
                                           const int *xInd,
                                           float *y,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDsctr(cusparseHandle_t handle,
                                           int nnz,
                                           const double *xVal,
                                           const int *xInd,
                                           double *y,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseCsctr(cusparseHandle_t handle,
                                           int nnz,
                                           const cuComplex *xVal,
                                           const int *xInd,
                                           cuComplex *y,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZsctr(cusparseHandle_t handle,
                                           int nnz,
                                           const cuDoubleComplex *xVal,
                                           const int *xInd,
                                           cuDoubleComplex *y,
                                           cusparseIndexBase_t idxBase);



cusparseStatus_t cusparseSroti(cusparseHandle_t handle,
                                           int nnz,
                                           float *xVal,
                                           const int *xInd,
                                           float *y,
                                           float c,
                                           float s,
                                           cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDroti(cusparseHandle_t handle,
                                           int nnz,
                                           double *xVal,
                                           const int *xInd,
                                           double *y,
                                           double c,
                                           double s,
                                           cusparseIndexBase_t idxBase);






cusparseStatus_t cusparseScsrmv(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            float alpha,
                                            const cusparseMatDescr_t descrA,
                                            const float *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const float *x,
                                            float beta,
                                            float *y);

cusparseStatus_t cusparseDcsrmv(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            double alpha,
                                            const cusparseMatDescr_t descrA,
                                            const double *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const double *x,
                                            double beta,
                                            double *y);

cusparseStatus_t cusparseCcsrmv(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            cuComplex alpha,
                                            const cusparseMatDescr_t descrA,
                                            const cuComplex *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const cuComplex *x,
                                            cuComplex beta,
                                            cuComplex *y);

cusparseStatus_t cusparseZcsrmv(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            cuDoubleComplex alpha,
                                            const cusparseMatDescr_t descrA,
                                            const cuDoubleComplex *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const cuDoubleComplex *x,
                                            cuDoubleComplex beta,
                                            cuDoubleComplex *y);




cusparseStatus_t cusparseScsrsv_analysis(cusparseHandle_t handle,
                                                     cusparseOperation_t transA,
                                                     int m,
                                                     const cusparseMatDescr_t descrA,
                                                     const float *csrValA,
                                                     const int *csrRowPtrA,
                                                     const int *csrColIndA,
                                                     cusparseSolveAnalysisInfo_t info);

cusparseStatus_t cusparseDcsrsv_analysis(cusparseHandle_t handle,
                                                     cusparseOperation_t transA,
                                                     int m,
                                                     const cusparseMatDescr_t descrA,
                                                     const double *csrValA,
                                                     const int *csrRowPtrA,
                                                     const int *csrColIndA,
                                                     cusparseSolveAnalysisInfo_t info);

cusparseStatus_t cusparseCcsrsv_analysis(cusparseHandle_t handle,
                                                     cusparseOperation_t transA,
                                                     int m,
                                                     const cusparseMatDescr_t descrA,
                                                     const cuComplex *csrValA,
                                                     const int *csrRowPtrA,
                                                     const int *csrColIndA,
                                                     cusparseSolveAnalysisInfo_t info);

cusparseStatus_t cusparseZcsrsv_analysis(cusparseHandle_t handle,
                                                     cusparseOperation_t transA,
                                                     int m,
                                                     const cusparseMatDescr_t descrA,
                                                     const cuDoubleComplex *csrValA,
                                                     const int *csrRowPtrA,
                                                     const int *csrColIndA,
                                                     cusparseSolveAnalysisInfo_t info);


cusparseStatus_t cusparseScsrsv_solve(cusparseHandle_t handle,
                                                  cusparseOperation_t transA,
                                                  int m,
                                                  float alpha,
                                                  const cusparseMatDescr_t descrA,
                                                  const float *csrValA,
                                                  const int *csrRowPtrA,
                                                  const int *csrColIndA,
                                                  cusparseSolveAnalysisInfo_t info,
                                                  const float *x,
                                                  float *y);

cusparseStatus_t cusparseDcsrsv_solve(cusparseHandle_t handle,
                                                  cusparseOperation_t transA,
                                                  int m,
                                                  double alpha,
                                                  const cusparseMatDescr_t descrA,
                                                  const double *csrValA,
                                                  const int *csrRowPtrA,
                                                  const int *csrColIndA,
                                                  cusparseSolveAnalysisInfo_t info,
                                                  const double *x,
                                                  double *y);

cusparseStatus_t cusparseCcsrsv_solve(cusparseHandle_t handle,
                                                  cusparseOperation_t transA,
                                                  int m,
                                                  cuComplex alpha,
                                                  const cusparseMatDescr_t descrA,
                                                  const cuComplex *csrValA,
                                                  const int *csrRowPtrA,
                                                  const int *csrColIndA,
                                                  cusparseSolveAnalysisInfo_t info,
                                                  const cuComplex *x,
                                                  cuComplex *y);

cusparseStatus_t cusparseZcsrsv_solve(cusparseHandle_t handle,
                                                  cusparseOperation_t transA,
                                                  int m,
                                                  cuDoubleComplex alpha,
                                                  const cusparseMatDescr_t descrA,
                                                  const cuDoubleComplex *csrValA,
                                                  const int *csrRowPtrA,
                                                  const int *csrColIndA,
                                                  cusparseSolveAnalysisInfo_t info,
                                                  const cuDoubleComplex *x,
                                                  cuDoubleComplex *y);





cusparseStatus_t cusparseScsrmm(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            int k,
                                            float alpha,
                                            const cusparseMatDescr_t descrA,
                                            const float *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const float *B,
                                            int ldb,
                                            float beta,
                                            float *C,
                                            int ldc);

cusparseStatus_t cusparseDcsrmm(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            int k,
                                            double alpha,
                                            const cusparseMatDescr_t descrA,
                                            const double *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const double *B,
                                            int ldb,
                                            double beta,
                                            double *C,
                                            int ldc);

cusparseStatus_t cusparseCcsrmm(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            int k,
                                            cuComplex alpha,
                                            const cusparseMatDescr_t descrA,
                                            const cuComplex *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const cuComplex *B,
                                            int ldb,
                                            cuComplex beta,
                                            cuComplex *C,
                                            int ldc);

cusparseStatus_t cusparseZcsrmm(cusparseHandle_t handle,
                                            cusparseOperation_t transA,
                                            int m,
                                            int n,
                                            int k,
                                            cuDoubleComplex alpha,
                                            const cusparseMatDescr_t descrA,
                                            const cuDoubleComplex *csrValA,
                                            const int *csrRowPtrA,
                                            const int *csrColIndA,
                                            const cuDoubleComplex *B,
                                            int ldb,
                                            cuDoubleComplex beta,
                                            cuDoubleComplex *C,
                                            int ldc);


cusparseStatus_t cusparseSnnz(cusparseHandle_t handle,
                                          cusparseDirection_t dirA,
                                          int m,
                                          int n,
                                          const cusparseMatDescr_t descrA,
                                          const float *A,
                                          int lda,
                                          int *nnzPerRowCol,
                                          int *nnzHostPtr);

cusparseStatus_t cusparseDnnz(cusparseHandle_t handle,
                                          cusparseDirection_t dirA,
                                          int m,
                                          int n,
                                          const cusparseMatDescr_t descrA,
                                          const double *A,
                                          int lda,
                                          int *nnzPerRowCol,
                                          int *nnzHostPtr);

cusparseStatus_t cusparseCnnz(cusparseHandle_t handle,
                                          cusparseDirection_t dirA,
                                          int m,
                                          int n,
                                          const cusparseMatDescr_t descrA,
                                          const cuComplex *A,
                                          int lda,
                                          int *nnzPerRowCol,
                                          int *nnzHostPtr);

cusparseStatus_t cusparseZnnz(cusparseHandle_t handle,
                                          cusparseDirection_t dirA,
                                          int m,
                                          int n,
                                          const cusparseMatDescr_t descrA,
                                          const cuDoubleComplex *A,
                                          int lda,
                                          int *nnzPerRowCol,
                                          int *nnzHostPtr);


cusparseStatus_t cusparseSdense2csr(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const float *A,
                                                int lda,
                                                const int *nnzPerRow,
                                                float *csrValA,
                                                int *csrRowPtrA,
                                                int *csrColIndA);

cusparseStatus_t cusparseDdense2csr(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const double *A,
                                                int lda,
                                                const int *nnzPerRow,
                                                double *csrValA,
                                                int *csrRowPtrA,
                                                int *csrColIndA);

cusparseStatus_t cusparseCdense2csr(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuComplex *A,
                                                int lda,
                                                const int *nnzPerRow,
                                                cuComplex *csrValA,
                                                int *csrRowPtrA,
                                                int *csrColIndA);

cusparseStatus_t cusparseZdense2csr(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuDoubleComplex *A,
                                                int lda,
                                                const int *nnzPerRow,
                                                cuDoubleComplex *csrValA,
                                                int *csrRowPtrA,
                                                int *csrColIndA);


cusparseStatus_t cusparseScsr2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const float *csrValA,
                                                const int *csrRowPtrA,
                                                const int *csrColIndA,
                                                float *A,
                                                int lda);

cusparseStatus_t cusparseDcsr2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const double *csrValA,
                                                const int *csrRowPtrA,
                                                const int *csrColIndA,
                                                double *A,
                                                int lda);

cusparseStatus_t cusparseCcsr2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuComplex *csrValA,
                                                const int *csrRowPtrA,
                                                const int *csrColIndA,
                                                cuComplex *A,
                                                int lda);

cusparseStatus_t cusparseZcsr2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuDoubleComplex *csrValA,
                                                const int *csrRowPtrA,
                                                const int *csrColIndA,
                                                cuDoubleComplex *A,
                                                int lda);


cusparseStatus_t cusparseSdense2csc(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const float *A,
                                                int lda,
                                                const int *nnzPerCol,
                                                float *cscValA,
                                                int *cscRowIndA,
                                                int *cscColPtrA);

cusparseStatus_t cusparseDdense2csc(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const double *A,
                                                int lda,
                                                const int *nnzPerCol,
                                                double *cscValA,
                                                int *cscRowIndA,
                                                int *cscColPtrA);

cusparseStatus_t cusparseCdense2csc(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuComplex *A,
                                                int lda,
                                                const int *nnzPerCol,
                                                cuComplex *cscValA,
                                                int *cscRowIndA,
                                                int *cscColPtrA);

cusparseStatus_t cusparseZdense2csc(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const
                                                cusparseMatDescr_t descrA,
                                                const cuDoubleComplex *A,
                                                int lda,
                                                const int *nnzPerCol,
                                                cuDoubleComplex *cscValA,
                                                int *cscRowIndA,
                                                int *cscColPtrA);


cusparseStatus_t cusparseScsc2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const float *cscValA,
                                                const int *cscRowIndA,
                                                const int *cscColPtrA,
                                                float *A,
                                                int lda);

cusparseStatus_t cusparseDcsc2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const double *cscValA,
                                                const int *cscRowIndA,
                                                const int *cscColPtrA,
                                                double *A,
                                                int lda);

cusparseStatus_t cusparseCcsc2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuComplex *cscValA,
                                                const int *cscRowIndA,
                                                const int *cscColPtrA,
                                                cuComplex *A,
                                                int lda);

cusparseStatus_t cusparseZcsc2dense(cusparseHandle_t handle,
                                                int m,
                                                int n,
                                                const cusparseMatDescr_t descrA,
                                                const cuDoubleComplex *cscValA,
                                                const int *cscRowIndA,
                                                const int *cscColPtrA,
                                                cuDoubleComplex *A,
                                                int lda);


cusparseStatus_t cusparseXcoo2csr(cusparseHandle_t handle,
                                              const int *cooRowInd,
                                              int nnz,
                                              int m,
                                              int *csrRowPtr,
                                              cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseXcsr2coo(cusparseHandle_t handle,
                                              const int *csrRowPtr,
                                              int nnz,
                                              int m,
                                              int *cooRowInd,
                                              cusparseIndexBase_t idxBase);


cusparseStatus_t cusparseScsr2csc(cusparseHandle_t handle,
                                              int m,
                                              int n,
                                              const float *csrVal,
                                              const int *csrRowPtr,
                                              const int *csrColInd,
                                              float *cscVal,
                                              int *cscRowInd,
                                              int *cscColPtr,
                                              int copyValues,
                                              cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseDcsr2csc(cusparseHandle_t handle,
                                              int m,
                                              int n,
                                              const double *csrVal,
                                              const int *csrRowPtr,
                                              const int *csrColInd,
                                              double *cscVal,
                                              int *cscRowInd,
                                              int *cscColPtr,
                                              int copyValues,
                                              cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseCcsr2csc(cusparseHandle_t handle,
                                              int m,
                                              int n,
                                              const cuComplex *csrVal,
                                              const int *csrRowPtr,
                                              const int *csrColInd,
                                              cuComplex *cscVal,
                                              int *cscRowInd,
                                              int *cscColPtr,
                                              int copyValues,
                                              cusparseIndexBase_t idxBase);

cusparseStatus_t cusparseZcsr2csc(cusparseHandle_t handle,
                                              int m,
                                              int n,
                                              const cuDoubleComplex *csrVal,
                                              const int *csrRowPtr,
                                              const int *csrColInd,
                                              cuDoubleComplex *cscVal,
                                              int *cscRowInd,
                                              int *cscColPtr,
                                              int copyValues,
                                              cusparseIndexBase_t idxBase);
# 8 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h" 2


# 1 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h" 1







# 1 "/home/dpnkarthik/petsc-rnet/include/private/vecimpl.h" 1
# 9 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/dvecimpl.h" 1
# 14 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/dvecimpl.h"
typedef struct {
  PetscScalar *array; PetscScalar *array_allocated; PetscScalar *unplacedarray;
} Vec_Seq;

extern PetscErrorCode VecMDot_Seq(Vec,PetscInt,const Vec[],PetscScalar *);
extern PetscErrorCode VecMTDot_Seq(Vec,PetscInt,const Vec[],PetscScalar *);
extern PetscErrorCode VecMin_Seq(Vec,PetscInt*,PetscReal *);
extern PetscErrorCode VecSet_Seq(Vec,PetscScalar);
extern PetscErrorCode VecMAXPY_Seq(Vec,PetscInt,const PetscScalar *,Vec *);
extern PetscErrorCode VecAYPX_Seq(Vec,PetscScalar,Vec);
extern PetscErrorCode VecWAXPY_Seq(Vec,PetscScalar,Vec,Vec);
extern PetscErrorCode VecAXPBYPCZ_Seq(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
extern PetscErrorCode VecMaxPointwiseDivide_Seq(Vec,Vec,PetscReal*);
extern PetscErrorCode VecPlaceArray_Seq(Vec,const PetscScalar *);
extern PetscErrorCode VecResetArray_Seq(Vec);
extern PetscErrorCode VecReplaceArray_Seq(Vec,const PetscScalar *);
extern PetscErrorCode VecDot_Seq(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecTDot_Seq(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecScale_Seq(Vec,PetscScalar);
extern PetscErrorCode VecAXPY_Seq(Vec,PetscScalar,Vec);
extern PetscErrorCode VecAXPBY_Seq(Vec,PetscScalar,PetscScalar,Vec);
extern PetscErrorCode VecMax_Seq(Vec,PetscInt*,PetscReal *);
extern PetscErrorCode VecNorm_Seq(Vec,NormType,PetscReal*);
extern PetscErrorCode VecDestroy_Seq(Vec);
extern PetscErrorCode VecDuplicate_Seq(Vec,Vec*);
extern PetscErrorCode VecSetOption_Seq(Vec,VecOption,PetscBool);
extern PetscErrorCode VecGetValues_Seq(Vec,PetscInt,const PetscInt*,PetscScalar*);
extern PetscErrorCode VecSetValues_Seq(Vec,PetscInt,const PetscInt*,const PetscScalar*,InsertMode);
extern PetscErrorCode VecSetValuesBlocked_Seq(Vec,PetscInt,const PetscInt*,const PetscScalar*,InsertMode);
extern PetscErrorCode VecView_Seq(Vec,PetscViewer);
extern PetscErrorCode VecGetSize_Seq(Vec,PetscInt*);
extern PetscErrorCode VecCopy_Seq(Vec,Vec);
extern PetscErrorCode VecSwap_Seq(Vec,Vec);
extern PetscErrorCode VecConjugate_Seq(Vec);
extern PetscErrorCode VecDestroy_Seq(Vec);
extern PetscErrorCode VecSetRandom_Seq(Vec,PetscRandom);
extern PetscErrorCode VecPointwiseMult_Seq(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseMax_Seq(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseMaxAbs_Seq(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseMin_Seq(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseDivide_Seq(Vec,Vec,Vec);


extern PetscErrorCode VecCreate_Seq(Vec);

extern PetscErrorCode VecCreate_Seq_Private(Vec,const PetscScalar[]);
# 10 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h" 2


# 46 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"
 int mtNexts;
 uint mtNexti;
 uint s_seeds[624];
# 62 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"
static PetscErrorCode VecView_Seq_ASCII(Vec ,PetscViewer);
static PetscErrorCode PinnedMalloc(PetscScalar** x,PetscInt n);
static PetscErrorCode PinnedFree(PetscScalar* x);
extern PetscErrorCode VecDotNorm2_SeqGPU(Vec,Vec,PetscScalar *, PetscScalar *);
extern PetscErrorCode VecPointwiseDivide_SeqGPU(Vec,Vec,Vec);
extern PetscErrorCode VecMaxPointwiseDivide_SeqGPU(Vec,Vec,PetscReal*);
extern PetscErrorCode VecWAXPY_SeqGPU(Vec,PetscScalar,Vec,Vec);
extern PetscErrorCode VecMDot_SeqGPU(Vec,PetscInt,const Vec[],PetscScalar *);
extern PetscErrorCode VecSet_SeqGPU(Vec,PetscScalar);
extern PetscErrorCode VecMAXPY_SeqGPU(Vec,PetscInt,const PetscScalar *,Vec *);
extern PetscErrorCode VecAXPBYPCZ_SeqGPU(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
extern PetscErrorCode VecPointwiseMult_SeqGPU(Vec,Vec,Vec);
extern PetscErrorCode VecPlaceArray_SeqGPU(Vec,const PetscScalar *);
extern PetscErrorCode VecResetArray_SeqGPU(Vec);
extern PetscErrorCode VecReplaceArray_SeqGPU(Vec,const PetscScalar *);
extern PetscErrorCode VecDot_SeqGPU(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecTDot_SeqGPU(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecScale_SeqGPU(Vec,PetscScalar);
extern PetscErrorCode VecCopy_SeqGPU(Vec,Vec);
extern PetscErrorCode VecSwap_SeqGPU(Vec,Vec);
extern PetscErrorCode VecAXPY_SeqGPU(Vec,PetscScalar,Vec);
extern PetscErrorCode VecAXPBY_SeqGPU(Vec,PetscScalar,PetscScalar,Vec);
extern PetscErrorCode VecDuplicate_SeqGPU(Vec,Vec *);
extern PetscErrorCode VecNorm_SeqGPU(Vec,NormType,PetscReal*);
extern PetscErrorCode VecCreate_SeqGPU(Vec);
extern PetscErrorCode VecView_SeqGPU(Vec,PetscViewer);
extern PetscErrorCode VecDestroy_SeqGPU(Vec);
extern PetscErrorCode VecDestroyVecs_SeqGPU(PetscInt,Vec*);

extern PetscErrorCode VecSetRandom_SeqGPU(Vec,PetscRandom);
extern PetscErrorCode VecSetValues_SeqGPU(Vec,PetscInt,const PetscInt*,const PetscScalar *y, InsertMode);
extern PetscErrorCode VecCopyOverH2D(Vec,PetscScalar*);
extern PetscErrorCode VecCopyOverD2H(Vec,PetscScalar*);
extern PetscErrorCode VecCopyBlockD2H(Vec,PetscScalar*,PetscInt,PetscInt);
extern PetscErrorCode VecCopyBlockH2D(Vec,PetscScalar*,PetscInt,PetscInt);
extern PetscErrorCode VecCopyOverDevice(Vec,Vec);
extern PetscErrorCode VecCopyBlockDevice(Vec, Vec, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode VecView_SeqGPU(Vec,PetscViewer);
extern PetscErrorCode VecGPUAllocateCheck_Public(Vec);
extern PetscErrorCode VecCopyToGPUSome_Public(Vec,PetscInt);
extern PetscErrorCode VecGetArray_SeqGPU(Vec,PetscScalar**);
extern PetscErrorCode VecRestoreArray_SeqGPU(Vec,PetscScalar**);
extern PetscErrorCode VecCheckCUDAError(const char *);
extern PetscErrorCode VecCheckCUDAStatus(cudaError_t ,const char *);
extern PetscErrorCode VecGetLocalSize_SeqGPU(Vec , PetscInt *);
extern PetscErrorCode VecGetSize_SeqGPU(Vec , PetscInt *);
extern PetscErrorCode VecCompare_SeqGPU(Vec,Vec, PetscBool*, PetscInt, PetscInt);
extern PetscErrorCode VecCheck_SeqGPU(Vec);
extern void kernAXPBYPCZ(double*, double*, double*,int*);
extern void kernCheck(double*, int*);
extern void kernCopyLen(int*, int*);
extern void kernCODevice(double*,double*, int*);
extern void kernCompare(double*, double*,int*, int*, int*);
extern void kernSet(double*, int*);
extern void kernScale(double*, int*);
extern void kernCopy(double*, double*, int*, int*);
extern void kernDot(double*,double*,int*,int*,int*,double*);
extern void kernRedDot(int*,double*,double*);
extern void kernAXPY(double*, double*, int*,double*);
extern void kernWAXPY(double*, double*, int*, double*);
extern void kernWXPY(double*, double*, int*, double*);
extern void kernWXMY(double*, double*, int*, double*);
extern void kernXPY(double*, double*, int*);
extern void kernNorm2(double*,int*,int*,int*,double*);
extern void kernRedNorm(int*,double*,double*);
extern void kernPDIV(double*,double*,int*,double* );
extern void kernPMULT(double*,double*,int*,double* );
extern void kernMAXPDIV(double*,double*,int*,int*,int*,double*);
extern void kernMAX(int*,double*,double*);
extern void kernRand(double*, int*);
extern void kernRandS(uint *);

 void mt19937si(uint);
 void mt19937sai(uint*,uint);
 uint mt19937s(void);
 uint mt19937sl(void);
# 148 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"
extern PetscBool synchronizeGPU;
# 157 "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"
typedef enum {VEC_SINGLE,VEC_PERSIST,VEC_DEALLOC,VEC_COLLECT} VecUsageGPUFlag;
typedef enum {VEC_UNALLOC,VEC_ALLOC,VEC_GPU, VEC_CPU,VEC_SYNCHED} VecGPUFlag;
typedef struct{
  PetscScalar *array; PetscScalar *array_allocated; PetscScalar *unplacedarray;

  VecGPUFlag syncState;

  PetscInt* length;
  PetscScalar* scalar;
  PetscScalar* cpuptr;
  PetscScalar* devptr;
  PetscScalar* zval;
  cudaStream_t stream;
  PetscInt* offset;
  PetscInt* segment;

 }Vec_SeqGPU;






static inline PetscErrorCode VecGPUAllocateCheck(Vec v)
{



  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 185; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGPUAllocateCheck") && strcmp("VecGPUAllocateCheck","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",185,"VecGPUAllocateCheck","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecGPUAllocateCheck_SeqGPU\n");
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecCopyToGPU(Vec v)
{



  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 201; petscstack->currentsize++; } do { if (strcmp(__func__,"VecCopyToGPU") && strcmp("VecCopyToGPU","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",201,"VecCopyToGPU","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecCopyToGPU\n");





  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecCopyToGPUSome(Vec v,PetscInt x)
{



  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 222; petscstack->currentsize++; } do { if (strcmp(__func__,"VecCopyToGPUSome") && strcmp("VecCopyToGPUSome","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",222,"VecCopyToGPUSome","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecCopyToGPUSome\n");





  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}





static inline PetscErrorCode VecGPUGetArrayReadWrite(Vec v, PetscScalar *va)
{


  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 240; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGetArrayReadWrite") && strcmp("VecGetArrayReadWrite","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",240,"VecGetArrayReadWrite","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecGetArrayReadWrite\n");


  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecGPURestoreArrayReadWrite(Vec v, PetscScalar *va)
{


  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 257; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGPURestoreArrayReadWrite") && strcmp("VecGPURestoreArrayReadWrite","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",257,"VecGPURestoreArrayReadWrite","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecGetRestoreReadWrite\n");


  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecGPUGetArrayRead(Vec v, PetscScalar *va)
{


  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 274; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGPUGetArrayRead") && strcmp("VecGPUGetArrayRead","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",274,"VecGPUGetArrayRead","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecGetArrayRead\n");


  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecGPURestoreArrayRead(Vec v, PetscScalar *va)
{
  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 289; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGPURestoreArrayRead") && strcmp("VecGPURestoreArrayRead","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",289,"VecGPURestoreArrayRead","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecGetRestoreRead\n");


  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecGPUGetArrayWrite(Vec v,PetscScalar *va)
{


  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 306; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGPUGetArrayWrite") && strcmp("VecGPUGetArrayWrite","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",306,"VecGPUGetArrayWrite","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecGetArrayWrite\n");


  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







static inline PetscErrorCode VecGPURestoreArrayWrite(Vec v,PetscScalar *va)
{


  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 323; petscstack->currentsize++; } do { if (strcmp(__func__,"VecGPURestoreArrayWrite") && strcmp("VecGPURestoreArrayWrite","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","/home/dpnkarthik/petsc-rnet/include/../src/vec/vec/impls/seq/seqgpu/gpuvecimpl.h",323,"VecGPURestoreArrayWrite","__func__",__func__); } } while (0); } while (0);
  printf("Call to VecRestoreArrayWrite\n");


  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}


# 11 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h" 2

# 56 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h"
typedef struct {
  MatScalar *bdiag,*ibdiag,*ssor_work;
  PetscInt bdiagsize;
  PetscBool ibdiagvalid;

  PetscBool use;
  PetscInt node_count;
  PetscInt *size;
  PetscInt limit;
  PetscInt max_limit;
  PetscBool checked;
} Mat_SeqAIJ_Inode;

extern PetscErrorCode MatView_SeqAIJ_Inode(Mat,PetscViewer);
extern PetscErrorCode MatAssemblyEnd_SeqAIJ_Inode(Mat,MatAssemblyType);
extern PetscErrorCode MatDestroy_SeqAIJ_Inode(Mat);
extern PetscErrorCode MatCreate_SeqAIJ_Inode(Mat);
extern PetscErrorCode MatSetOption_SeqAIJ_Inode(Mat,MatOption,PetscBool );
extern PetscErrorCode MatDuplicate_SeqAIJ_Inode(Mat,MatDuplicateOption,Mat*);
extern PetscErrorCode MatDuplicateNoCreate_SeqAIJ(Mat,Mat,MatDuplicateOption,PetscBool );
extern PetscErrorCode MatLUFactorNumeric_SeqAIJ_Inode_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqAIJ_Inode(Mat,Mat,const MatFactorInfo*);


typedef struct{
  cusparseHandle_t handle;
  cusparseMatDescr_t descrip;
  PetscInt *dev_csrRowOffsets;
  PetscInt *dev_csrIndices;
  PetscScalar *dev_dataA;
}CudaSparseVars;



typedef struct {
  PetscBool roworiented; PetscInt nonew; PetscInt nounused; PetscBool singlemalloc; PetscInt maxnz; PetscInt *imax; PetscInt *ilen; PetscBool free_imax_ilen; PetscInt reallocs; PetscInt rmax; PetscBool keepnonzeropattern; PetscBool ignorezeroentries; PetscInt *xtoy,*xtoyB; Mat XtoY; PetscBool free_ij; PetscBool free_a; Mat_CompressedRow compressedrow; PetscInt nz; PetscInt *i; PetscInt *j; PetscInt *diag; PetscBool free_diag; MatScalar *a; PetscScalar *solve_work; IS row, col, icol; PetscBool pivotinblocks; Mat parent;
  Mat_SeqAIJ_Inode inode;
  MatScalar *saved_values;
  PetscScalar *idiag,*mdiag,*ssor_work;
  PetscBool idiagvalid;
  PetscScalar *ibdiag;
  PetscBool ibdiagvalid;
  PetscScalar fshift,omega;
  ISColoring coloring;
  CudaSparseVars cudav;
} Mat_SeqAIJ;
# 112 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h"
static inline PetscErrorCode MatSeqXAIJFreeAIJ(Mat AA,MatScalar **a,PetscInt **j,PetscInt **i)
{
  PetscErrorCode ierr;
  Mat_SeqAIJ *A = (Mat_SeqAIJ*) AA->data;
  if (A->singlemalloc) {
    ierr = ((*i)=0,(*j)=0,((*a) && ((*PetscTrFree)((void*)(*a),117,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/") || ((*a) = 0,0))));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),117,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    if (A->free_a) {ierr = ((*a) && ((*PetscTrFree)((void*)(*a),119,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/") || ((*a) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),119,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);}
    if (A->free_ij) {ierr = ((*j) && ((*PetscTrFree)((void*)(*j),120,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/") || ((*j) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),120,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);}
    if (A->free_ij) {ierr = ((*i) && ((*PetscTrFree)((void*)(*i),121,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/") || ((*i) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),121,__func__,"/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);}
  }
  return 0;
}
# 163 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/aij/seq/aij.h"

extern PetscErrorCode MatSeqAIJSetPreallocation_SeqAIJ(Mat,PetscInt,const PetscInt*);



extern PetscErrorCode MatMult_SeqAIJ(Mat A,Vec,Vec);


extern PetscErrorCode MatILUFactorSymbolic_SeqAIJ_inplace(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactorSymbolic_SeqAIJ(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactorSymbolic_SeqAIJ_ilu0(Mat,Mat,IS,IS,const MatFactorInfo*);

extern PetscErrorCode MatICCFactorSymbolic_SeqAIJ_inplace(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatICCFactorSymbolic_SeqAIJ(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorSymbolic_SeqAIJ_inplace(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorSymbolic_SeqAIJ(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorNumeric_SeqAIJ_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorNumeric_SeqAIJ(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatDuplicate_SeqAIJ(Mat,MatDuplicateOption,Mat*);
extern PetscErrorCode MatCopy_SeqAIJ(Mat,Mat,MatStructure);
extern PetscErrorCode MatMissingDiagonal_SeqAIJ(Mat,PetscBool *,PetscInt*);
extern PetscErrorCode MatMarkDiagonal_SeqAIJ(Mat);


extern PetscErrorCode MatMultAdd_SeqAIJ(Mat A,Vec,Vec,Vec);
extern PetscErrorCode MatMultTranspose_SeqAIJ(Mat A,Vec,Vec);
extern PetscErrorCode MatMultTransposeAdd_SeqAIJ(Mat A,Vec,Vec,Vec);
extern PetscErrorCode MatSOR_SeqAIJ(Mat,Vec,PetscReal,MatSORType,PetscReal,PetscInt,PetscInt,Vec);

extern PetscErrorCode MatSetColoring_SeqAIJ(Mat,ISColoring);
extern PetscErrorCode MatSetValuesAdic_SeqAIJ(Mat,void*);
extern PetscErrorCode MatSetValuesAdifor_SeqAIJ(Mat,PetscInt,void*);

extern PetscErrorCode MatGetSymbolicTranspose_SeqAIJ(Mat,PetscInt *[],PetscInt *[]);
extern PetscErrorCode MatGetSymbolicTransposeReduced_SeqAIJ(Mat,PetscInt,PetscInt,PetscInt *[],PetscInt *[]);
extern PetscErrorCode MatRestoreSymbolicTranspose_SeqAIJ(Mat,PetscInt *[],PetscInt *[]);
extern PetscErrorCode MatToSymmetricIJ_SeqAIJ(PetscInt,PetscInt*,PetscInt*,PetscInt,PetscInt,PetscInt**,PetscInt**);
extern PetscErrorCode MatLUFactorSymbolic_SeqAIJ_inplace(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorSymbolic_SeqAIJ(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqAIJ_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqAIJ(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqAIJ_InplaceWithPerm(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactor_SeqAIJ(Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatSolve_SeqAIJ_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqAIJ(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqAIJ_Inode_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqAIJ_Inode(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqAIJ_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqAIJ_NaturalOrdering(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqAIJ_InplaceWithPerm(Mat,Vec,Vec);
extern PetscErrorCode MatSolveAdd_SeqAIJ_inplace(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolveAdd_SeqAIJ(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqAIJ_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqAIJ(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTransposeAdd_SeqAIJ_inplace(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolveTransposeAdd_SeqAIJ(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMatSolve_SeqAIJ_inplace(Mat,Mat,Mat);
extern PetscErrorCode MatMatSolve_SeqAIJ(Mat,Mat,Mat);
extern PetscErrorCode MatEqual_SeqAIJ(Mat A,Mat B,PetscBool * flg);
extern PetscErrorCode MatFDColoringCreate_SeqAIJ(Mat,ISColoring,MatFDColoring);
extern PetscErrorCode MatLoad_SeqAIJ(Mat,PetscViewer);
extern PetscErrorCode RegisterApplyPtAPRoutines_Private(Mat);
extern PetscErrorCode MatMatMultSymbolic_SeqAIJ_SeqAIJ(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatMultNumeric_SeqAIJ_SeqAIJ(Mat,Mat,Mat);
extern PetscErrorCode MatPtAPSymbolic_SeqAIJ(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatPtAPNumeric_SeqAIJ(Mat,Mat,Mat);
extern PetscErrorCode MatPtAPSymbolic_SeqAIJ_SeqAIJ(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatPtAPNumeric_SeqAIJ_SeqAIJ(Mat,Mat,Mat);
extern PetscErrorCode MatMatMultTranspose_SeqAIJ_SeqAIJ(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatMatMultTransposeSymbolic_SeqAIJ_SeqAIJ(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatMultTransposeNumeric_SeqAIJ_SeqAIJ(Mat,Mat,Mat);
extern PetscErrorCode MatSetValues_SeqAIJ(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatGetRow_SeqAIJ(Mat,PetscInt,PetscInt*,PetscInt**,PetscScalar**);
extern PetscErrorCode MatRestoreRow_SeqAIJ(Mat,PetscInt,PetscInt*,PetscInt**,PetscScalar**);
extern PetscErrorCode MatAXPY_SeqAIJ(Mat,PetscScalar,Mat,MatStructure);
extern PetscErrorCode MatGetRowIJ_SeqAIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatRestoreRowIJ_SeqAIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatGetColumnIJ_SeqAIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatRestoreColumnIJ_SeqAIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,PetscInt *[],PetscInt *[],PetscBool *);
extern PetscErrorCode MatDestroy_SeqAIJ(Mat);
extern PetscErrorCode MatView_SeqAIJ(Mat,PetscViewer);

extern PetscErrorCode Mat_CheckInode(Mat,PetscBool );
extern PetscErrorCode Mat_CheckInode_FactorLU(Mat,PetscBool );

extern PetscErrorCode MatAXPYGetPreallocation_SeqAIJ(Mat,Mat,PetscInt*);


extern PetscErrorCode MatConvert_SeqAIJ_SeqSBAIJ(Mat,const char*,MatReuse,Mat*);
extern PetscErrorCode MatConvert_SeqAIJ_SeqBAIJ(Mat,const char*,MatReuse,Mat*);
extern PetscErrorCode MatConvert_SeqAIJ_SeqAIJPERM(Mat,const char*,MatReuse,Mat*);
extern PetscErrorCode MatReorderForNonzeroDiagonal_SeqAIJ(Mat,PetscReal,IS,IS);
extern PetscErrorCode MatMatMult_SeqAIJ_SeqAIJ(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatCreate_SeqAIJ(Mat);

extern PetscErrorCode MatAssemblyEnd_SeqAIJ(Mat,MatAssemblyType);
extern PetscErrorCode MatDestroy_SeqAIJ(Mat);
# 6 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/baij.h" 2
# 1 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/ftn-kernels/fsolvebaij.h" 1
# 7 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/baij.h" 2
# 28 "/home/dpnkarthik/petsc-rnet/include/../src/mat/impls/baij/seq/baij.h"
typedef struct {
  PetscBool roworiented; PetscInt nonew; PetscInt nounused; PetscBool singlemalloc; PetscInt maxnz; PetscInt *imax; PetscInt *ilen; PetscBool free_imax_ilen; PetscInt reallocs; PetscInt rmax; PetscBool keepnonzeropattern; PetscBool ignorezeroentries; PetscInt *xtoy,*xtoyB; Mat XtoY; PetscBool free_ij; PetscBool free_a; Mat_CompressedRow compressedrow; PetscInt nz; PetscInt *i; PetscInt *j; PetscInt *diag; PetscBool free_diag; MatScalar *a; PetscScalar *solve_work; IS row, col, icol; PetscBool pivotinblocks; Mat parent;
  PetscInt bs2; PetscInt mbs,nbs; PetscScalar *mult_work; PetscScalar *sor_work; MatScalar *saved_values; Mat sbaijMat; MatScalar *idiag; PetscBool idiagvalid;
} Mat_SeqBAIJ;


extern PetscErrorCode MatSeqBAIJSetPreallocation_SeqBAIJ(Mat,PetscInt,PetscInt,PetscInt*);

extern PetscErrorCode MatILUFactorSymbolic_SeqBAIJ_inplace(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactorSymbolic_SeqBAIJ(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatICCFactorSymbolic_SeqBAIJ(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorSymbolic_SeqBAIJ(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorNumeric_SeqBAIJ_N(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorNumeric_SeqBAIJ_N_NaturalOrdering(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatDuplicate_SeqBAIJ(Mat,MatDuplicateOption,Mat*);
extern PetscErrorCode MatMissingDiagonal_SeqBAIJ(Mat,PetscBool *,PetscInt*);
extern PetscErrorCode MatMarkDiagonal_SeqBAIJ(Mat);
extern PetscErrorCode MatILUDTFactor_SeqBAIJ(Mat,IS,IS,const MatFactorInfo*,Mat*);

extern PetscErrorCode MatLUFactorSymbolic_SeqBAIJ_inplace(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorSymbolic_SeqBAIJ(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactor_SeqBAIJ(Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatIncreaseOverlap_SeqBAIJ(Mat,PetscInt,IS*,PetscInt);
extern PetscErrorCode MatGetSubMatrix_SeqBAIJ(Mat,IS,IS,MatReuse,Mat*);
extern PetscErrorCode MatGetSubMatrices_SeqBAIJ(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat*[]);
extern PetscErrorCode MatMultTranspose_SeqBAIJ(Mat,Vec,Vec);
extern PetscErrorCode MatMultHermitianTranspose_SeqBAIJ(Mat,Vec,Vec);
extern PetscErrorCode MatMultTransposeAdd_SeqBAIJ(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultHermitianTransposeAdd_SeqBAIJ(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatScale_SeqBAIJ(Mat,PetscScalar);
extern PetscErrorCode MatNorm_SeqBAIJ(Mat,NormType,PetscReal *);
extern PetscErrorCode MatEqual_SeqBAIJ(Mat,Mat,PetscBool *);
extern PetscErrorCode MatGetDiagonal_SeqBAIJ(Mat,Vec);
extern PetscErrorCode MatDiagonalScale_SeqBAIJ(Mat,Vec,Vec);
extern PetscErrorCode MatGetInfo_SeqBAIJ(Mat,MatInfoType,MatInfo *);
extern PetscErrorCode MatZeroEntries_SeqBAIJ(Mat);
extern PetscErrorCode MatDestroy_SeqBAIJ(Mat);
extern PetscErrorCode MatAssemblyEnd_SeqBAIJ(Mat,MatAssemblyType);

extern PetscErrorCode MatSeqBAIJ_UpdateFactorNumeric_NaturalOrdering(Mat);

extern PetscErrorCode MatSolve_SeqBAIJ_1_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_1(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_1_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_1_NaturalOrdering(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_2_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_2(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_2_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_2_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolve_SeqBAIJ_3_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_3(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_3_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_3_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolve_SeqBAIJ_4_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_4(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_4_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_4_NaturalOrdering(Mat,Vec,Vec);





extern PetscErrorCode MatSolve_SeqBAIJ_5_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_5(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_5_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_5_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolve_SeqBAIJ_6_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_6(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_6_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_6_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolve_SeqBAIJ_7_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_7(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_7_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_7_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolve_SeqBAIJ_15_NaturalOrdering_ver1(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_15_NaturalOrdering_ver2(Mat,Vec,Vec);

extern PetscErrorCode MatSolve_SeqBAIJ_N_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_N(Mat,Vec,Vec);
extern PetscErrorCode MatSolve_SeqBAIJ_N_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_1_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_1(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_1_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_1_NaturalOrdering(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_2_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_2(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_2_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_2_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_3_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_3(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_3_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_3_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_4_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_4(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_4_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_4_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_5_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_5(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_5_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_5_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_6_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_6(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_6_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_6_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_7_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_7(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_7_NaturalOrdering_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_7_NaturalOrdering(Mat,Vec,Vec);

extern PetscErrorCode MatSolveTranspose_SeqBAIJ_N_inplace(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTranspose_SeqBAIJ_N(Mat,Vec,Vec);

extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_N(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_1_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_1(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_2_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_2(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_2_NaturalOrdering_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_2_NaturalOrdering(Mat,Mat,const MatFactorInfo*);

extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_3_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_3(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_3_NaturalOrdering_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_3_NaturalOrdering(Mat,Mat,const MatFactorInfo*);

extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_4_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_4(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_4_NaturalOrdering_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_4_NaturalOrdering(Mat,Mat,const MatFactorInfo*);





extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_5_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_5(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_5_NaturalOrdering_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_5_NaturalOrdering(Mat,Mat,const MatFactorInfo*);

extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_6_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_6(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_6_NaturalOrdering_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_6_NaturalOrdering(Mat,Mat,const MatFactorInfo*);

extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_7_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_7(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_7_NaturalOrdering_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_7_NaturalOrdering(Mat,Mat,const MatFactorInfo*);

extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_15_NaturalOrdering(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_N_inplace(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric_SeqBAIJ_N(Mat,Mat,const MatFactorInfo*);

extern PetscErrorCode MatMult_SeqBAIJ_1(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_2(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_3(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_4(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_5(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_6(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_7(Mat,Vec,Vec);

extern PetscErrorCode MatMult_SeqBAIJ_15_ver1(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_15_ver2(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_15_ver3(Mat,Vec,Vec);
extern PetscErrorCode MatMult_SeqBAIJ_15_ver4(Mat,Vec,Vec);

extern PetscErrorCode MatMult_SeqBAIJ_N(Mat,Vec,Vec);

extern PetscErrorCode MatMultAdd_SeqBAIJ_1(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_2(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_3(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_4(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_5(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_6(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_7(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultAdd_SeqBAIJ_N(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatLoad_SeqBAIJ(Mat,PetscViewer);
extern PetscErrorCode MatSeqBAIJSetNumericFactorization_inplace(Mat,PetscBool );
extern PetscErrorCode MatSeqBAIJSetNumericFactorization(Mat,PetscBool );
# 3 "baij2.c" 2
# 1 "/home/dpnkarthik/petsc-rnet/include/../src/mat/blockinvert.h" 1
# 15 "/home/dpnkarthik/petsc-rnet/include/../src/mat/blockinvert.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscblaslapack.h" 1
# 17 "/home/dpnkarthik/petsc-rnet/include/petscblaslapack.h"
# 1 "/home/dpnkarthik/petsc-rnet/include/petscblaslapack_uscore.h" 1
# 18 "/home/dpnkarthik/petsc-rnet/include/petscblaslapack.h" 2







extern void dgetrf_(PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
extern void dorgqr_(PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dgeqrf_(PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);

extern PetscReal ddot_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);
extern PetscReal dnrm2_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);
extern void dscal_(const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*);
extern void dcopy_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
extern void dswap_(const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
extern void daxpy_(const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
extern PetscReal dasum_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);
extern void dpttrf_(const PetscBLASInt*,PetscReal*,PetscScalar*,const PetscBLASInt*);
extern void dstein_(const PetscBLASInt*,PetscReal*,PetscReal*,const PetscBLASInt*,PetscReal*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*);
extern void dgesv_(const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);

extern void dpotrf_(const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dpotrs_(const char*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dgemv_(const char*,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,const PetscScalar *,const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*);
extern void dgetrs_(const char*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dtrmv_(const char*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
extern void dgemm_(const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
# 61 "/home/dpnkarthik/petsc-rnet/include/petscblaslapack.h"
extern void dgelss_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,const PetscReal*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
extern void dsyev_(const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dsyevx_(const char*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
extern void dsygv_(PetscBLASInt*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dsygvx_(PetscBLASInt*,const char*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
extern void dpttrs_(PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dstebz_(const char*,const char*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
 extern void dgerfs_(const char*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscInt*,PetscBLASInt*);
extern void dtrsen_(const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
extern void dtgsen_(PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
extern void dgges_(const char*,const char*,const char*,PetscBLASInt(*)(),PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
extern void dhseqr_(const char*,const char*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
# 86 "/home/dpnkarthik/petsc-rnet/include/petscblaslapack.h"
extern void dgeev_(const char*,const char*,PetscBLASInt *,PetscScalar *,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
extern void dgesvd_(const char*,const char*,PetscBLASInt *,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);


# 16 "/home/dpnkarthik/petsc-rnet/include/../src/mat/blockinvert.h" 2






extern PetscErrorCode LINPACKdgefa(MatScalar*,PetscInt,PetscInt*);
extern PetscErrorCode LINPACKdgedi(MatScalar*,PetscInt,PetscInt*,MatScalar*);
extern PetscErrorCode Kernel_A_gets_inverse_A_2(MatScalar*,PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_3(MatScalar*,PetscReal);
# 101 "/home/dpnkarthik/petsc-rnet/include/../src/mat/blockinvert.h"
extern PetscErrorCode Kernel_A_gets_inverse_A_4(MatScalar *,PetscReal);



extern PetscErrorCode Kernel_A_gets_inverse_A_5(MatScalar *,PetscInt*,MatScalar*,PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_6(MatScalar *,PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_7(MatScalar *,PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_9(MatScalar *,PetscReal);
extern PetscErrorCode Kernel_A_gets_inverse_A_15(MatScalar *,PetscInt*,MatScalar*,PetscReal);
# 4 "baij2.c" 2
# 1 "/home/dpnkarthik/petsc-rnet/include/petscbt.h" 1




# 30 "/home/dpnkarthik/petsc-rnet/include/petscbt.h"
typedef char* PetscBT;

extern char _BT_mask;
extern char _BT_c;
extern PetscInt _BT_idx;
# 78 "/home/dpnkarthik/petsc-rnet/include/petscbt.h"

# 5 "baij2.c" 2




PetscErrorCode MatIncreaseOverlap_SeqBAIJ(Mat A,PetscInt is_max,IS is[],PetscInt ov)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode ierr;
  PetscInt row,i,j,k,l,m,n,*nidx,isz,val,ival;
  const PetscInt *idx;
  PetscInt start,end,*ai,*aj,bs,*nidx2;
  PetscBT table;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 18; petscstack->currentsize++; } do { if (strcmp(__func__,"MatIncreaseOverlap_SeqBAIJ") && strcmp("MatIncreaseOverlap_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",18,"MatIncreaseOverlap_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  m = a->mbs;
  ai = a->i;
  aj = a->j;
  bs = A->rmap->bs;

  if (ov < 0) return PetscError(((MPI_Comm)0x44000001),24,__func__,"baij2.c","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"Negative overlap specified");

  ierr = (((((m)/8 +1)*sizeof(char) != 0) ? (*PetscTrMalloc)((((m)/8 +1)*sizeof(char)),26,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&(table))) : (*(&(table)) = 0,0) ) || PetscMemzero(table,sizeof(char)*((m)/8 +1)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),26,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (((m+1)*sizeof(PetscInt) != 0) ? (*PetscTrMalloc)(((m+1)*sizeof(PetscInt)),27,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&nidx)) : (*(&nidx) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),27,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (((A->rmap->N+1)*sizeof(PetscInt) != 0) ? (*PetscTrMalloc)(((A->rmap->N+1)*sizeof(PetscInt)),28,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&nidx2)) : (*(&nidx2) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),28,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  for (i=0; i<is_max; i++) {

    isz = 0;
    ierr = PetscMemzero(table,sizeof(char)*((m)/8 +1));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),33,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);


    ierr = ISGetIndices(is[i],&idx);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),36,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = ISGetLocalSize(is[i],&n);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),37,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);


    for (j=0; j<n ; ++j){
      ival = idx[j]/bs;
      if (ival>=m) return PetscError(((MPI_Comm)0x44000001),42,__func__,"baij2.c","src/mat/impls/baij/seq/",63,PETSC_ERROR_INITIAL,"index greater than mat-dim");
      if(!(_BT_idx = (ival)/8, _BT_c = table[_BT_idx], _BT_mask = (char)1 << ((ival)%8), table[_BT_idx] = _BT_c | _BT_mask, _BT_c & _BT_mask)) { nidx[isz++] = ival;}
    }
    ierr = ISRestoreIndices(is[i],&idx);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),45,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = ISDestroy(&is[i]);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),46,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

    k = 0;
    for (j=0; j<ov; j++){
      n = isz;
      for (; k<n ; k++){
        row = nidx[k];
        start = ai[row];
        end = ai[row+1];
        for (l = start; l<end ; l++){
          val = aj[l];
          if (!(_BT_idx = (val)/8, _BT_c = table[_BT_idx], _BT_mask = (char)1 << ((val)%8), table[_BT_idx] = _BT_c | _BT_mask, _BT_c & _BT_mask)) {nidx[isz++] = val;}
        }
      }
    }

    for (j=0; j<isz; j++) {
      for (k=0; k<bs; k++)
        nidx2[j*bs+k] = nidx[j]*bs+k;
    }
    ierr = ISCreateGeneral(((MPI_Comm)0x44000001),isz*bs,nidx2,PETSC_COPY_VALUES,is+i);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),66,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = ((table) && ((*PetscTrFree)((void*)(table),68,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((table) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),68,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ((nidx) && ((*PetscTrFree)((void*)(nidx),69,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((nidx) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),69,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ((nidx2) && ((*PetscTrFree)((void*)(nidx2),70,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((nidx2) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),70,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatGetSubMatrix_SeqBAIJ_Private(Mat A,IS isrow,IS iscol,MatReuse scall,Mat *B)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data,*c;
  PetscErrorCode ierr;
  PetscInt *smap,i,k,kstart,kend,oldcols = a->nbs,*lens;
  PetscInt row,mat_i,*mat_j,tcol,*mat_ilen;
  const PetscInt *irow,*icol;
  PetscInt nrows,ncols,*ssmap,bs=A->rmap->bs,bs2=a->bs2;
  PetscInt *aj = a->j,*ai = a->i;
  MatScalar *mat_a;
  Mat C;
  PetscBool flag,sorted;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 89; petscstack->currentsize++; } do { if (strcmp(__func__,"MatGetSubMatrix_SeqBAIJ_Private") && strcmp("MatGetSubMatrix_SeqBAIJ_Private","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",89,"MatGetSubMatrix_SeqBAIJ_Private","__func__",__func__); } } while (0); } while (0);
  ierr = ISSorted(iscol,&sorted);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),90,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (!sorted) return PetscError(((MPI_Comm)0x44000001),91,__func__,"baij2.c","src/mat/impls/baij/seq/",73,PETSC_ERROR_INITIAL,"IS is not sorted");

  ierr = ISGetIndices(isrow,&irow);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),93,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISGetIndices(iscol,&icol);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),94,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISGetLocalSize(isrow,&nrows);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),95,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISGetLocalSize(iscol,&ncols);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),96,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  ierr = (((1+oldcols)*sizeof(PetscInt) != 0) ? (*PetscTrMalloc)(((1+oldcols)*sizeof(PetscInt)),98,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&smap)) : (*(&smap) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),98,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ssmap = smap;
  ierr = (((1+nrows)*sizeof(PetscInt) != 0) ? (*PetscTrMalloc)(((1+nrows)*sizeof(PetscInt)),100,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&lens)) : (*(&lens) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),100,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = PetscMemzero(smap,oldcols*sizeof(PetscInt));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),101,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  for (i=0; i<ncols; i++) smap[icol[i]] = i+1;

  for (i=0; i<nrows; i++) {
    kstart = ai[irow[i]];
    kend = kstart + a->ilen[irow[i]];
    lens[i] = 0;
      for (k=kstart; k<kend; k++) {
        if (ssmap[aj[k]]) {
          lens[i]++;
        }
      }
    }

  if (scall == MAT_REUSE_MATRIX) {
    c = (Mat_SeqBAIJ *)((*B)->data);

    if (c->mbs!=nrows || c->nbs!=ncols || (*B)->rmap->bs!=bs) return PetscError(((MPI_Comm)0x44000001),118,__func__,"baij2.c","src/mat/impls/baij/seq/",60,PETSC_ERROR_INITIAL,"Submatrix wrong size");
    ierr = PetscMemcmp(c->ilen,lens,c->mbs *sizeof(PetscInt),&flag);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),119,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    if (!flag) return PetscError(((MPI_Comm)0x44000001),120,__func__,"baij2.c","src/mat/impls/baij/seq/",60,PETSC_ERROR_INITIAL,"Cannot reuse matrix. wrong no of nonzeros");
    ierr = PetscMemzero(c->ilen,c->mbs*sizeof(PetscInt));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),121,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    C = *B;
  } else {
    ierr = MatCreate(((PetscObject)A)->comm,&C);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),124,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = MatSetSizes(C,nrows*bs,ncols*bs,-1,-1);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),125,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = MatSetType(C,((PetscObject)A)->type_name);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),126,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = MatSeqBAIJSetPreallocation_SeqBAIJ(C,bs,0,lens);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),127,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  c = (Mat_SeqBAIJ *)(C->data);
  for (i=0; i<nrows; i++) {
    row = irow[i];
    kstart = ai[row];
    kend = kstart + a->ilen[row];
    mat_i = c->i[i];
    mat_j = c->j + mat_i;
    mat_a = c->a + mat_i*bs2;
    mat_ilen = c->ilen + i;
    for (k=kstart; k<kend; k++) {
      if ((tcol=ssmap[a->j[k]])) {
        *mat_j++ = tcol - 1;
        ierr = PetscMemcpy(mat_a,a->a+k*bs2,bs2*sizeof(MatScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),141,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
        mat_a += bs2;
        (*mat_ilen)++;
      }
    }
  }


  ierr = ISRestoreIndices(iscol,&icol);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),149,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ((smap) && ((*PetscTrFree)((void*)(smap),150,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((smap) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),150,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ((lens) && ((*PetscTrFree)((void*)(lens),151,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((lens) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),151,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),152,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),153,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  ierr = ISRestoreIndices(isrow,&irow);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),155,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  *B = C;
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatGetSubMatrix_SeqBAIJ(Mat A,IS isrow,IS iscol,MatReuse scall,Mat *B)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  IS is1,is2;
  PetscErrorCode ierr;
  PetscInt *vary,*iary,nrows,ncols,i,bs=A->rmap->bs,count;
  const PetscInt *irow,*icol;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 170; petscstack->currentsize++; } do { if (strcmp(__func__,"MatGetSubMatrix_SeqBAIJ") && strcmp("MatGetSubMatrix_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",170,"MatGetSubMatrix_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  ierr = ISGetIndices(isrow,&irow);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),171,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISGetIndices(iscol,&icol);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),172,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISGetLocalSize(isrow,&nrows);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),173,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISGetLocalSize(iscol,&ncols);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),174,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);



  ierr = ((*(&iary) = 0,(((a->mbs)*sizeof(PetscInt)+(a->mbs)*sizeof(PetscInt)+(16 -1) != 0) ? (*PetscTrMalloc)(((a->mbs)*sizeof(PetscInt)+(a->mbs)*sizeof(PetscInt)+(16 -1)),178,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&vary)) : (*(&vary) = 0,0) )) || (*(&iary) = (PetscInt*)(void*)((((uintptr_t)(*(&vary)+a->mbs))+(16 -1)) & ~(16 -1)),0));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),178,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = PetscMemzero(vary,a->mbs*sizeof(PetscInt));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),179,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  for (i=0; i<nrows; i++) vary[irow[i]/bs]++;
  count = 0;
  for (i=0; i<a->mbs; i++) {
    if (vary[i]!=0 && vary[i]!=bs) return PetscError(((MPI_Comm)0x44000001),183,__func__,"baij2.c","src/mat/impls/baij/seq/",60,PETSC_ERROR_INITIAL,"Index set does not match blocks");
    if (vary[i]==bs) iary[count++] = i;
  }
  ierr = ISCreateGeneral(((MPI_Comm)0x44000001),count,iary,PETSC_COPY_VALUES,&is1);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),186,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  ierr = PetscMemzero(vary,(a->mbs)*sizeof(PetscInt));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),188,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  for (i=0; i<ncols; i++) vary[icol[i]/bs]++;
  count = 0;
  for (i=0; i<a->mbs; i++) {
    if (vary[i]!=0 && vary[i]!=bs) return PetscError(((MPI_Comm)0x44000001),192,__func__,"baij2.c","src/mat/impls/baij/seq/",77,PETSC_ERROR_INITIAL,"Internal error in PETSc");
    if (vary[i]==bs) iary[count++] = i;
  }
  ierr = ISCreateGeneral(((MPI_Comm)0x44000001),count,iary,PETSC_COPY_VALUES,&is2);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),195,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISRestoreIndices(isrow,&irow);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),196,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISRestoreIndices(iscol,&icol);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),197,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ((iary)=0, ((vary) && ((*PetscTrFree)((void*)(vary),198,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((vary) = 0,0))));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),198,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  ierr = MatGetSubMatrix_SeqBAIJ_Private(A,is1,is2,scall,B);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),200,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISDestroy(&is1);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),201,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = ISDestroy(&is2);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),202,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatGetSubMatrices_SeqBAIJ(Mat A,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *B[])
{
  PetscErrorCode ierr;
  PetscInt i;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 213; petscstack->currentsize++; } do { if (strcmp(__func__,"MatGetSubMatrices_SeqBAIJ") && strcmp("MatGetSubMatrices_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",213,"MatGetSubMatrices_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  if (scall == MAT_INITIAL_MATRIX) {
    ierr = (((n+1)*sizeof(Mat) != 0) ? (*PetscTrMalloc)(((n+1)*sizeof(Mat)),215,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(B)) : (*(B) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),215,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }

  for (i=0; i<n; i++) {
    ierr = MatGetSubMatrix_SeqBAIJ(A,irow[i],icol[i],scall,&(*B)[i]);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),219,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
# 231 "baij2.c"
PetscErrorCode MatMult_SeqBAIJ_1(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z,sum;
  const PetscScalar *x;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,n,nonzerorow=0;
  const PetscInt *idx,*ii,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 242; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_1") && strcmp("MatMult_SeqBAIJ_1","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",242,"MatMult_SeqBAIJ_1","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),243,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),244,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    ierr = PetscMemzero(z,mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),250,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    mbs = a->mbs;
    ii = a->i;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0];
    v = a->a + ii[0];
    idx = a->j + ii[0];
    ii++;
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+1*n),*_end = (const char*)((v+1*n)+(1*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    sum = 0.0;
    { PetscInt __i;for(__i=0;__i<n;__i++) sum += v[__i] * x[idx[__i]];};
    if (usecprow){
      z[ridx[i]] = sum;
    } else {
      nonzerorow += (n>0);
      z[i] = sum;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),272,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),273,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)2.0*a->nz - nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),274,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMult_SeqBAIJ_2(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,*zarray;
  const PetscScalar *x,*xb;
  PetscScalar x1,x2;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 291; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_2") && strcmp("MatMult_SeqBAIJ_2","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",291,"MatMult_SeqBAIJ_2","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),292,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),293,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    sum1 = 0.0; sum2 = 0.0;
    nonzerorow += (n>0);
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+4*n),*_end = (const char*)((v+4*n)+(4*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 2*(*idx++); x1 = xb[0]; x2 = xb[1];
      sum1 += v[0]*x1 + v[2]*x2;
      sum2 += v[1]*x1 + v[3]*x2;
      v += 4;
    }
    if (usecprow) z = zarray + 2*ridx[i];
    z[0] = sum1; z[1] = sum2;
    if (!usecprow) z += 2;
  }
  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),323,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),324,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)8.0*a->nz - 2.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),325,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMult_SeqBAIJ_3(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,x1,x2,x3,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;






  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 346; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_3") && strcmp("MatMult_SeqBAIJ_3","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",346,"MatMult_SeqBAIJ_3","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),347,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),348,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0;
    nonzerorow += (n>0);
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+9*n),*_end = (const char*)((v+9*n)+(9*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 3*(*idx++); x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      sum1 += v[0]*x1 + v[3]*x2 + v[6]*x3;
      sum2 += v[1]*x1 + v[4]*x2 + v[7]*x3;
      sum3 += v[2]*x1 + v[5]*x2 + v[8]*x3;
      v += 9;
    }
    if (usecprow) z = zarray + 3*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3;
    if (!usecprow) z += 3;
  }
  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),379,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),380,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)18.0*a->nz - 3.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),381,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMult_SeqBAIJ_4(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,x1,x2,x3,x4,*zarray;
  const PetscScalar *x,*xb;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 397; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_4") && strcmp("MatMult_SeqBAIJ_4","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",397,"MatMult_SeqBAIJ_4","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),398,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),399,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0;
    nonzerorow += (n>0);
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+16*n),*_end = (const char*)((v+16*n)+(16*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 4*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      sum1 += v[0]*x1 + v[4]*x2 + v[8]*x3 + v[12]*x4;
      sum2 += v[1]*x1 + v[5]*x2 + v[9]*x3 + v[13]*x4;
      sum3 += v[2]*x1 + v[6]*x2 + v[10]*x3 + v[14]*x4;
      sum4 += v[3]*x1 + v[7]*x2 + v[11]*x3 + v[15]*x4;
      v += 16;
    }
    if (usecprow) z = zarray + 4*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4;
    if (!usecprow) z += 4;
  }
  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),432,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),433,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)32.0*a->nz - 4.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),434,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMult_SeqBAIJ_5(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar sum1,sum2,sum3,sum4,sum5,x1,x2,x3,x4,x5,*z = 0,*zarray;
  const PetscScalar *xb,*x;
  const MatScalar *v;
  PetscErrorCode ierr;
  const PetscInt *idx,*ii,*ridx=0;
  PetscInt mbs,i,j,n,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 451; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_5") && strcmp("MatMult_SeqBAIJ_5","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",451,"MatMult_SeqBAIJ_5","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),452,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),453,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0;
    nonzerorow += (n>0);
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+25*n),*_end = (const char*)((v+25*n)+(25*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 5*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4];
      sum1 += v[0]*x1 + v[5]*x2 + v[10]*x3 + v[15]*x4 + v[20]*x5;
      sum2 += v[1]*x1 + v[6]*x2 + v[11]*x3 + v[16]*x4 + v[21]*x5;
      sum3 += v[2]*x1 + v[7]*x2 + v[12]*x3 + v[17]*x4 + v[22]*x5;
      sum4 += v[3]*x1 + v[8]*x2 + v[13]*x3 + v[18]*x4 + v[23]*x5;
      sum5 += v[4]*x1 + v[9]*x2 + v[14]*x3 + v[19]*x4 + v[24]*x5;
      v += 25;
    }
    if (usecprow) z = zarray + 5*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5;
    if (!usecprow) z += 5;
  }
  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),487,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),488,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)50.0*a->nz - 5.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),489,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}




PetscErrorCode MatMult_SeqBAIJ_6(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,sum5,sum6;
  const PetscScalar *x,*xb;
  PetscScalar x1,x2,x3,x4,x5,x6,*zarray;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 507; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_6") && strcmp("MatMult_SeqBAIJ_6","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",507,"MatMult_SeqBAIJ_6","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),508,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),509,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0;
    nonzerorow += (n>0);
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+36*n),*_end = (const char*)((v+36*n)+(36*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 6*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5];
      sum1 += v[0]*x1 + v[6]*x2 + v[12]*x3 + v[18]*x4 + v[24]*x5 + v[30]*x6;
      sum2 += v[1]*x1 + v[7]*x2 + v[13]*x3 + v[19]*x4 + v[25]*x5 + v[31]*x6;
      sum3 += v[2]*x1 + v[8]*x2 + v[14]*x3 + v[20]*x4 + v[26]*x5 + v[32]*x6;
      sum4 += v[3]*x1 + v[9]*x2 + v[15]*x3 + v[21]*x4 + v[27]*x5 + v[33]*x6;
      sum5 += v[4]*x1 + v[10]*x2 + v[16]*x3 + v[22]*x4 + v[28]*x5 + v[34]*x6;
      sum6 += v[5]*x1 + v[11]*x2 + v[17]*x3 + v[23]*x4 + v[29]*x5 + v[35]*x6;
      v += 36;
    }
    if (usecprow) z = zarray + 6*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6;
    if (!usecprow) z += 6;
  }

  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),545,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),546,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)72.0*a->nz - 6.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),547,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMult_SeqBAIJ_7(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7;
  const PetscScalar *x,*xb;
  PetscScalar x1,x2,x3,x4,x5,x6,x7,*zarray;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 564; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_7") && strcmp("MatMult_SeqBAIJ_7","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",564,"MatMult_SeqBAIJ_7","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),565,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),566,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    nonzerorow += (n>0);
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+49*n),*_end = (const char*)((v+49*n)+(49*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 7*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      sum1 += v[0]*x1 + v[7]*x2 + v[14]*x3 + v[21]*x4 + v[28]*x5 + v[35]*x6 + v[42]*x7;
      sum2 += v[1]*x1 + v[8]*x2 + v[15]*x3 + v[22]*x4 + v[29]*x5 + v[36]*x6 + v[43]*x7;
      sum3 += v[2]*x1 + v[9]*x2 + v[16]*x3 + v[23]*x4 + v[30]*x5 + v[37]*x6 + v[44]*x7;
      sum4 += v[3]*x1 + v[10]*x2 + v[17]*x3 + v[24]*x4 + v[31]*x5 + v[38]*x6 + v[45]*x7;
      sum5 += v[4]*x1 + v[11]*x2 + v[18]*x3 + v[25]*x4 + v[32]*x5 + v[39]*x6 + v[46]*x7;
      sum6 += v[5]*x1 + v[12]*x2 + v[19]*x3 + v[26]*x4 + v[33]*x5 + v[40]*x6 + v[47]*x7;
      sum7 += v[6]*x1 + v[13]*x2 + v[20]*x3 + v[27]*x4 + v[34]*x5 + v[41]*x6 + v[48]*x7;
      v += 49;
    }
    if (usecprow) z = zarray + 7*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    if (!usecprow) z += 7;
  }

  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),603,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),604,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)98.0*a->nz - 7.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),605,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}






PetscErrorCode MatMult_SeqBAIJ_15_ver1(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar *zarray,xv;
  const MatScalar *v;
  PetscErrorCode ierr;
  const PetscInt *ii,*ij=a->j,*idx;
  PetscInt mbs,i,j,k,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 626; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_15_ver1") && strcmp("MatMult_SeqBAIJ_15_ver1","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",626,"MatMult_SeqBAIJ_15_ver1","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),627,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),628,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[i+1] - ii[i];
    idx = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    nonzerorow += (n>0);
    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);

      for(k=0;k<15;k++){
        xv = xb[k];
 sum1 += v[0]*xv;
        sum2 += v[1]*xv;
 sum3 += v[2]*xv;
 sum4 += v[3]*xv;
 sum5 += v[4]*xv;
        sum6 += v[5]*xv;
 sum7 += v[6]*xv;
 sum8 += v[7]*xv;
        sum9 += v[8]*xv;
        sum10 += v[9]*xv;
 sum11 += v[10]*xv;
 sum12 += v[11]*xv;
 sum13 += v[12]*xv;
        sum14 += v[13]*xv;
 sum15 += v[14]*xv;
 v += 15;
      }
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),678,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),679,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)450.0*a->nz - 15.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),680,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}




PetscErrorCode MatMult_SeqBAIJ_15_ver2(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar x1,x2,x3,x4,*zarray;
  const MatScalar *v;
  PetscErrorCode ierr;
  const PetscInt *ii,*ij=a->j,*idx;
  PetscInt mbs,i,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 699; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_15_ver2") && strcmp("MatMult_SeqBAIJ_15_ver2","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",699,"MatMult_SeqBAIJ_15_ver2","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),700,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),701,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[i+1] - ii[i];
    idx = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    nonzerorow += (n>0);
    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];

      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3 + v[45]*x4;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3 + v[46]*x4;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3 + v[47]*x4;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3 + v[48]*x4;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3 + v[49]*x4;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3 + v[50]*x4;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3 + v[51]*x4;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3 + v[52]*x4;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3 + v[53]*x4;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3 + v[54]*x4;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4;

      v += 60;

      x1 = xb[4]; x2 = xb[5]; x3 = xb[6]; x4 = xb[7];

      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3 + v[45]*x4;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3 + v[46]*x4;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3 + v[47]*x4;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3 + v[48]*x4;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3 + v[49]*x4;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3 + v[50]*x4;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3 + v[51]*x4;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3 + v[52]*x4;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3 + v[53]*x4;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3 + v[54]*x4;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4;
      v += 60;

      x1 = xb[8]; x2 = xb[9]; x3 = xb[10]; x4 = xb[11];
      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3 + v[45]*x4;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3 + v[46]*x4;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3 + v[47]*x4;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3 + v[48]*x4;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3 + v[49]*x4;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3 + v[50]*x4;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3 + v[51]*x4;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3 + v[52]*x4;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3 + v[53]*x4;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3 + v[54]*x4;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4;
      v += 60;

      x1 = xb[12]; x2 = xb[13]; x3 = xb[14];
      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3;
      v += 45;
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),805,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),806,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)450.0*a->nz - 15.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),807,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}




PetscErrorCode MatMult_SeqBAIJ_15_ver3(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar x1,x2,x3,x4,x5,x6,x7,x8,*zarray;
  const MatScalar *v;
  PetscErrorCode ierr;
  const PetscInt *ii,*ij=a->j,*idx;
  PetscInt mbs,i,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 826; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_15_ver3") && strcmp("MatMult_SeqBAIJ_15_ver3","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",826,"MatMult_SeqBAIJ_15_ver3","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),827,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),828,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[i+1] - ii[i];
    idx = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    nonzerorow += (n>0);
    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      x8 = xb[7];

      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3 + v[45]*x4 + v[60]*x5 + v[75]*x6 + v[90]*x7 + v[105]*x8;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3 + v[46]*x4 + v[61]*x5 + v[76]*x6 + v[91]*x7 + v[106]*x8;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3 + v[47]*x4 + v[62]*x5 + v[77]*x6 + v[92]*x7 + v[107]*x8;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3 + v[48]*x4 + v[63]*x5 + v[78]*x6 + v[93]*x7 + v[108]*x8;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3 + v[49]*x4 + v[64]*x5 + v[79]*x6 + v[94]*x7 + v[109]*x8;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3 + v[50]*x4 + v[65]*x5 + v[80]*x6 + v[95]*x7 + v[110]*x8;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3 + v[51]*x4 + v[66]*x5 + v[81]*x6 + v[96]*x7 + v[111]*x8;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3 + v[52]*x4 + v[67]*x5 + v[82]*x6 + v[97]*x7 + v[112]*x8;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3 + v[53]*x4 + v[68]*x5 + v[83]*x6 + v[98]*x7 + v[113]*x8;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3 + v[54]*x4 + v[69]*x5 + v[84]*x6 + v[99]*x7 + v[114]*x8;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4 + v[70]*x5 + v[85]*x6 + v[100]*x7 + v[115]*x8;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4 + v[71]*x5 + v[86]*x6 + v[101]*x7 + v[116]*x8;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4 + v[72]*x5 + v[87]*x6 + v[102]*x7 + v[117]*x8;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4 + v[73]*x5 + v[88]*x6 + v[103]*x7 + v[118]*x8;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4 + v[74]*x5 + v[89]*x6 + v[104]*x7 + v[119]*x8;
      v += 120;

      x1 = xb[8]; x2 = xb[9]; x3 = xb[10]; x4 = xb[11]; x5 = xb[12]; x6 = xb[13]; x7 = xb[14];

      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3 + v[45]*x4 + v[60]*x5 + v[75]*x6 + v[90]*x7;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3 + v[46]*x4 + v[61]*x5 + v[76]*x6 + v[91]*x7;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3 + v[47]*x4 + v[62]*x5 + v[77]*x6 + v[92]*x7;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3 + v[48]*x4 + v[63]*x5 + v[78]*x6 + v[93]*x7;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3 + v[49]*x4 + v[64]*x5 + v[79]*x6 + v[94]*x7;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3 + v[50]*x4 + v[65]*x5 + v[80]*x6 + v[95]*x7;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3 + v[51]*x4 + v[66]*x5 + v[81]*x6 + v[96]*x7;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3 + v[52]*x4 + v[67]*x5 + v[82]*x6 + v[97]*x7;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3 + v[53]*x4 + v[68]*x5 + v[83]*x6 + v[98]*x7;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3 + v[54]*x4 + v[69]*x5 + v[84]*x6 + v[99]*x7;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4 + v[70]*x5 + v[85]*x6 + v[100]*x7;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4 + v[71]*x5 + v[86]*x6 + v[101]*x7;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4 + v[72]*x5 + v[87]*x6 + v[102]*x7;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4 + v[73]*x5 + v[88]*x6 + v[103]*x7;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4 + v[74]*x5 + v[89]*x6 + v[104]*x7;
      v += 105;
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),896,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),897,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)450.0*a->nz - 15.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),898,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}





PetscErrorCode MatMult_SeqBAIJ_15_ver4(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *z = 0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sum13,sum14,sum15;
  const PetscScalar *x,*xb;
  PetscScalar x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,*zarray;
  const MatScalar *v;
  PetscErrorCode ierr;
  const PetscInt *ii,*ij=a->j,*idx;
  PetscInt mbs,i,j,n,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 918; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_15_ver4") && strcmp("MatMult_SeqBAIJ_15_ver4","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",918,"MatMult_SeqBAIJ_15_ver4","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),919,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),920,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[i+1] - ii[i];
    idx = ij + ii[i];
    sum1 = 0.0; sum2 = 0.0; sum3 = 0.0; sum4 = 0.0; sum5 = 0.0; sum6 = 0.0; sum7 = 0.0;
    sum8 = 0.0; sum9 = 0.0; sum10 = 0.0; sum11 = 0.0; sum12 = 0.0; sum13 = 0.0; sum14 = 0.0;sum15 = 0.0;

    nonzerorow += (n>0);
    for (j=0; j<n; j++) {
      xb = x + 15*(idx[j]);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      x8 = xb[7]; x9 = xb[8]; x10 = xb[9]; x11 = xb[10]; x12 = xb[11]; x13 = xb[12]; x14 = xb[13];x15 = xb[14];

      sum1 += v[0]*x1 + v[15]*x2 + v[30]*x3 + v[45]*x4 + v[60]*x5 + v[75]*x6 + v[90]*x7 + v[105]*x8 + v[120]*x9 + v[135]*x10 + v[150]*x11 + v[165]*x12 + v[180]*x13 + v[195]*x14 + v[210]*x15;
      sum2 += v[1]*x1 + v[16]*x2 + v[31]*x3 + v[46]*x4 + v[61]*x5 + v[76]*x6 + v[91]*x7 + v[106]*x8 + v[121]*x9 + v[136]*x10 + v[151]*x11 + v[166]*x12 + v[181]*x13 + v[196]*x14 + v[211]*x15;
      sum3 += v[2]*x1 + v[17]*x2 + v[32]*x3 + v[47]*x4 + v[62]*x5 + v[77]*x6 + v[92]*x7 + v[107]*x8 + v[122]*x9 + v[137]*x10 + v[152]*x11 + v[167]*x12 + v[182]*x13 + v[197]*x14 + v[212]*x15;
      sum4 += v[3]*x1 + v[18]*x2 + v[33]*x3 + v[48]*x4 + v[63]*x5 + v[78]*x6 + v[93]*x7 + v[108]*x8 + v[123]*x9 + v[138]*x10 + v[153]*x11 + v[168]*x12 + v[183]*x13 + v[198]*x14 + v[213]*x15;
      sum5 += v[4]*x1 + v[19]*x2 + v[34]*x3 + v[49]*x4 + v[64]*x5 + v[79]*x6 + v[94]*x7 + v[109]*x8 + v[124]*x9 + v[139]*x10 + v[154]*x11 + v[169]*x12 + v[184]*x13 + v[199]*x14 + v[214]*x15;
      sum6 += v[5]*x1 + v[20]*x2 + v[35]*x3 + v[50]*x4 + v[65]*x5 + v[80]*x6 + v[95]*x7 + v[110]*x8 + v[125]*x9 + v[140]*x10 + v[155]*x11 + v[170]*x12 + v[185]*x13 + v[200]*x14 + v[215]*x15;
      sum7 += v[6]*x1 + v[21]*x2 + v[36]*x3 + v[51]*x4 + v[66]*x5 + v[81]*x6 + v[96]*x7 + v[111]*x8 + v[126]*x9 + v[141]*x10 + v[156]*x11 + v[171]*x12 + v[186]*x13 + v[201]*x14 + v[216]*x15;
      sum8 += v[7]*x1 + v[22]*x2 + v[37]*x3 + v[52]*x4 + v[67]*x5 + v[82]*x6 + v[97]*x7 + v[112]*x8 + v[127]*x9 + v[142]*x10 + v[157]*x11 + v[172]*x12 + v[187]*x13 + v[202]*x14 + v[217]*x15;
      sum9 += v[8]*x1 + v[23]*x2 + v[38]*x3 + v[53]*x4 + v[68]*x5 + v[83]*x6 + v[98]*x7 + v[113]*x8 + v[128]*x9 + v[143]*x10 + v[158]*x11 + v[173]*x12 + v[188]*x13 + v[203]*x14 + v[218]*x15;
      sum10 += v[9]*x1 + v[24]*x2 + v[39]*x3 + v[54]*x4 + v[69]*x5 + v[84]*x6 + v[99]*x7 + v[114]*x8 + v[129]*x9 + v[144]*x10 + v[159]*x11 + v[174]*x12 + v[189]*x13 + v[204]*x14 + v[219]*x15;
      sum11 += v[10]*x1 + v[25]*x2 + v[40]*x3 + v[55]*x4 + v[70]*x5 + v[85]*x6 + v[100]*x7 + v[115]*x8 + v[130]*x9 + v[145]*x10 + v[160]*x11 + v[175]*x12 + v[190]*x13 + v[205]*x14 + v[220]*x15;
      sum12 += v[11]*x1 + v[26]*x2 + v[41]*x3 + v[56]*x4 + v[71]*x5 + v[86]*x6 + v[101]*x7 + v[116]*x8 + v[131]*x9 + v[146]*x10 + v[161]*x11 + v[176]*x12 + v[191]*x13 + v[206]*x14 + v[221]*x15;
      sum13 += v[12]*x1 + v[27]*x2 + v[42]*x3 + v[57]*x4 + v[72]*x5 + v[87]*x6 + v[102]*x7 + v[117]*x8 + v[132]*x9 + v[147]*x10 + v[162]*x11 + v[177]*x12 + v[192]*x13 + v[207]*x14 + v[222]*x15;
      sum14 += v[13]*x1 + v[28]*x2 + v[43]*x3 + v[58]*x4 + v[73]*x5 + v[88]*x6 + v[103]*x7 + v[118]*x8 + v[133]*x9 + v[148]*x10 + v[163]*x11 + v[178]*x12 + v[193]*x13 + v[208]*x14 + v[223]*x15;
      sum15 += v[14]*x1 + v[29]*x2 + v[44]*x3 + v[59]*x4 + v[74]*x5 + v[89]*x6 + v[104]*x7 + v[119]*x8 + v[134]*x9 + v[149]*x10 + v[164]*x11 + v[179]*x12 + v[194]*x13 + v[209]*x14 + v[224]*x15;
      v += 225;
    }
    if (usecprow) z = zarray + 15*ridx[i];
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    z[7] = sum8; z[8] = sum9; z[9] = sum10; z[10] = sum11; z[11] = sum12; z[12] = sum13; z[13] = sum14;z[14] = sum15;

    if (!usecprow) z += 15;
  }

  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),969,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),970,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)450.0*a->nz - 15.0*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),971,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}







PetscErrorCode MatMult_SeqBAIJ_N(Mat A,Vec xx,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*z = 0,*xb,*work,*workt,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,bs=A->rmap->bs,j,n,bs2=a->bs2;
  PetscInt ncols,k,*ridx=0,nonzerorow=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 991; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMult_SeqBAIJ_N") && strcmp("MatMult_SeqBAIJ_N","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",991,"MatMult_SeqBAIJ_N","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),992,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),993,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  if (!a->mult_work) {
    k = (((A->rmap->n)<(A->cmap->n)) ? (A->cmap->n) : (A->rmap->n));
    ierr = (((k+1)*sizeof(PetscScalar) != 0) ? (*PetscTrMalloc)(((k+1)*sizeof(PetscScalar)),1009,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&a->mult_work)) : (*(&a->mult_work) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1009,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  work = a->mult_work;
  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    ncols = n*bs;
    workt = work;
    nonzerorow += (n>0);
    for (j=0; j<n; j++) {
      xb = x + bs*(*idx++);
      for (k=0; k<bs; k++) workt[k] = xb[k];
      workt += bs;
    }
    if (usecprow) z = zarray + bs*ridx[i];
    { PetscScalar _one = 1.0,_zero = 0.0; PetscBLASInt _ione = 1,_bbs,_bncols; _bbs = bs; _bncols = ncols; dgemv_("N",&(_bbs),&_bncols,&_one,v,&(_bbs),work,&_ione,&_zero,z,&_ione); };

    v += n*bs2;
    if (!usecprow) z += bs;
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1028,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1029,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)2.0*a->nz*bs2 - bs*nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1030,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_1(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  const PetscScalar *x;
  PetscScalar *y,*z,sum;
  const MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,n,*ridx=0,nonzerorow=0;
  const PetscInt *idx,*ii;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1047; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_1") && strcmp("MatMultAdd_SeqBAIJ_1","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1047,"MatMultAdd_SeqBAIJ_1","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1048,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&y);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1049,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1051,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    z = y;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(z,y,mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1060,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0];
    ii++;
    if (!usecprow){
      nonzerorow += (n>0);
      sum = y[i];
    } else {
      sum = y[ridx[i]];
    }
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+n),*_end = (const char*)((v+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    { PetscInt __i;for(__i=0;__i<n;__i++) sum += v[__i] * x[idx[__i]];};
    v += n;
    idx += n;
    if (usecprow){
      z[ridx[i]] = sum;
    } else {
      z[i] = sum;
    }
  }
  ierr = VecRestoreArrayRead(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1089,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&y);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1090,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1092,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)2.0*a->nz - nonzerorow),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1094,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_2(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*y = 0,*z = 0,*xb,sum1,sum2;
  PetscScalar x1,x2,*yarray,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1110; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_2") && strcmp("MatMultAdd_SeqBAIJ_2","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1110,"MatMultAdd_SeqBAIJ_2","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1111,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1112,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1114,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    zarray = yarray;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,2*mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1123,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,a->mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1129,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
  } else {
    ii = a->i;
    y = yarray;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow){
      z = zarray + 2*ridx[i];
      y = yarray + 2*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1];
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+4*n),*_end = (const char*)((v+4*n)+(4*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 2*(*idx++); x1 = xb[0]; x2 = xb[1];
      sum1 += v[0]*x1 + v[2]*x2;
      sum2 += v[1]*x1 + v[3]*x2;
      v += 4;
    }
    z[0] = sum1; z[1] = sum2;
    if (!usecprow){
      z += 2; y += 2;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1157,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1158,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1160,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)4.0*a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1162,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_3(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*y = 0,*z = 0,*xb,sum1,sum2,sum3,x1,x2,x3,*yarray,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1177; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_3") && strcmp("MatMultAdd_SeqBAIJ_3","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1177,"MatMultAdd_SeqBAIJ_3","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1178,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1179,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1181,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    zarray = yarray;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,3*mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1190,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y = yarray;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow){
      z = zarray + 3*ridx[i];
      y = yarray + 3*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2];
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+9*n),*_end = (const char*)((v+9*n)+(9*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 3*(*idx++); x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      sum1 += v[0]*x1 + v[3]*x2 + v[6]*x3;
      sum2 += v[1]*x1 + v[4]*x2 + v[7]*x3;
      sum3 += v[2]*x1 + v[5]*x2 + v[8]*x3;
      v += 9;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3;
    if (!usecprow){
      z += 3; y += 3;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1222,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1223,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1225,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)18.0*a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1227,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_4(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*y = 0,*z = 0,*xb,sum1,sum2,sum3,sum4,x1,x2,x3,x4,*yarray,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1242; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_4") && strcmp("MatMultAdd_SeqBAIJ_4","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1242,"MatMultAdd_SeqBAIJ_4","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1243,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1244,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1246,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    zarray = yarray;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,4*mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1255,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y = yarray;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow){
      z = zarray + 4*ridx[i];
      y = yarray + 4*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3];
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+16*n),*_end = (const char*)((v+16*n)+(16*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 4*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      sum1 += v[0]*x1 + v[4]*x2 + v[8]*x3 + v[12]*x4;
      sum2 += v[1]*x1 + v[5]*x2 + v[9]*x3 + v[13]*x4;
      sum3 += v[2]*x1 + v[6]*x2 + v[10]*x3 + v[14]*x4;
      sum4 += v[3]*x1 + v[7]*x2 + v[11]*x3 + v[15]*x4;
      v += 16;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4;
    if (!usecprow){
      z += 4; y += 4;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1289,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1290,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1292,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)32.0*a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1294,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_5(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*y = 0,*z = 0,*xb,sum1,sum2,sum3,sum4,sum5,x1,x2,x3,x4,x5;
  PetscScalar *yarray,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1310; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_5") && strcmp("MatMultAdd_SeqBAIJ_5","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1310,"MatMultAdd_SeqBAIJ_5","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1311,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1312,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1314,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    zarray = yarray;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,5*mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1323,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y = yarray;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow){
      z = zarray + 5*ridx[i];
      y = yarray + 5*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4];
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+25*n),*_end = (const char*)((v+25*n)+(25*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 5*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4];
      sum1 += v[0]*x1 + v[5]*x2 + v[10]*x3 + v[15]*x4 + v[20]*x5;
      sum2 += v[1]*x1 + v[6]*x2 + v[11]*x3 + v[16]*x4 + v[21]*x5;
      sum3 += v[2]*x1 + v[7]*x2 + v[12]*x3 + v[17]*x4 + v[22]*x5;
      sum4 += v[3]*x1 + v[8]*x2 + v[13]*x3 + v[18]*x4 + v[23]*x5;
      sum5 += v[4]*x1 + v[9]*x2 + v[14]*x3 + v[19]*x4 + v[24]*x5;
      v += 25;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5;
    if (!usecprow){
      z += 5; y += 5;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1358,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1359,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1361,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)50.0*a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1363,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}


PetscErrorCode MatMultAdd_SeqBAIJ_6(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*y = 0,*z = 0,*xb,sum1,sum2,sum3,sum4,sum5,sum6;
  PetscScalar x1,x2,x3,x4,x5,x6,*yarray,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1378; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_6") && strcmp("MatMultAdd_SeqBAIJ_6","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1378,"MatMultAdd_SeqBAIJ_6","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1379,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1380,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1382,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    zarray = yarray;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,6*mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1391,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y = yarray;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow){
      z = zarray + 6*ridx[i];
      y = yarray + 6*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4]; sum6 = y[5];
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+36*n),*_end = (const char*)((v+36*n)+(36*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 6*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5];
      sum1 += v[0]*x1 + v[6]*x2 + v[12]*x3 + v[18]*x4 + v[24]*x5 + v[30]*x6;
      sum2 += v[1]*x1 + v[7]*x2 + v[13]*x3 + v[19]*x4 + v[25]*x5 + v[31]*x6;
      sum3 += v[2]*x1 + v[8]*x2 + v[14]*x3 + v[20]*x4 + v[26]*x5 + v[32]*x6;
      sum4 += v[3]*x1 + v[9]*x2 + v[15]*x3 + v[21]*x4 + v[27]*x5 + v[33]*x6;
      sum5 += v[4]*x1 + v[10]*x2 + v[16]*x3 + v[22]*x4 + v[28]*x5 + v[34]*x6;
      sum6 += v[5]*x1 + v[11]*x2 + v[17]*x3 + v[23]*x4 + v[29]*x5 + v[35]*x6;
      v += 36;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6;
    if (!usecprow){
      z += 6; y += 6;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1427,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1428,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1430,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)72.0*a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1432,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_7(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*y = 0,*z = 0,*xb,sum1,sum2,sum3,sum4,sum5,sum6,sum7;
  PetscScalar x1,x2,x3,x4,x5,x6,x7,*yarray,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs=a->mbs,i,*idx,*ii,j,n,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1448; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_7") && strcmp("MatMultAdd_SeqBAIJ_7","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1448,"MatMultAdd_SeqBAIJ_7","__func__",__func__); } } while (0); } while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1449,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1450,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1452,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    zarray = yarray;
  }

  idx = a->j;
  v = a->a;
  if (usecprow){
    if (zz != yy){
      ierr = PetscMemcpy(zarray,yarray,7*mbs*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1461,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    }
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    ii = a->i;
    y = yarray;
    z = zarray;
  }

  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    if (usecprow){
      z = zarray + 7*ridx[i];
      y = yarray + 7*ridx[i];
    }
    sum1 = y[0]; sum2 = y[1]; sum3 = y[2]; sum4 = y[3]; sum5 = y[4]; sum6 = y[5]; sum7 = y[6];
    do { const char *_p = (const char*)(idx+n),*_end = (const char*)((idx+n)+(n)); for ( ; _p < _end; _p += 64) ; } while (0);
    do { const char *_p = (const char*)(v+49*n),*_end = (const char*)((v+49*n)+(49*n)); for ( ; _p < _end; _p += 64) ; } while (0);
    for (j=0; j<n; j++) {
      xb = x + 7*(*idx++);
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3]; x5 = xb[4]; x6 = xb[5]; x7 = xb[6];
      sum1 += v[0]*x1 + v[7]*x2 + v[14]*x3 + v[21]*x4 + v[28]*x5 + v[35]*x6 + v[42]*x7;
      sum2 += v[1]*x1 + v[8]*x2 + v[15]*x3 + v[22]*x4 + v[29]*x5 + v[36]*x6 + v[43]*x7;
      sum3 += v[2]*x1 + v[9]*x2 + v[16]*x3 + v[23]*x4 + v[30]*x5 + v[37]*x6 + v[44]*x7;
      sum4 += v[3]*x1 + v[10]*x2 + v[17]*x3 + v[24]*x4 + v[31]*x5 + v[38]*x6 + v[45]*x7;
      sum5 += v[4]*x1 + v[11]*x2 + v[18]*x3 + v[25]*x4 + v[32]*x5 + v[39]*x6 + v[46]*x7;
      sum6 += v[5]*x1 + v[12]*x2 + v[19]*x3 + v[26]*x4 + v[33]*x5 + v[40]*x6 + v[47]*x7;
      sum7 += v[6]*x1 + v[13]*x2 + v[20]*x3 + v[27]*x4 + v[34]*x5 + v[41]*x6 + v[48]*x7;
      v += 49;
    }
    z[0] = sum1; z[1] = sum2; z[2] = sum3; z[3] = sum4; z[4] = sum5; z[5] = sum6; z[6] = sum7;
    if (!usecprow){
      z += 7; y += 7;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1498,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(yy,&yarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1499,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (zz != yy) {
    ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1501,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)98.0*a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1503,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultAdd_SeqBAIJ_N(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *x,*z = 0,*xb,*work,*workt,*zarray;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,bs=A->rmap->bs,j,n,bs2=a->bs2;
  PetscInt ncols,k,*ridx=0;
  PetscBool usecprow=a->compressedrow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1519; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultAdd_SeqBAIJ_N") && strcmp("MatMultAdd_SeqBAIJ_N","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1519,"MatMultAdd_SeqBAIJ_N","__func__",__func__); } } while (0); } while (0);
  ierr = VecCopy(yy,zz);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1520,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1521,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1522,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = a->compressedrow.nrows;
    ii = a->compressedrow.i;
    ridx = a->compressedrow.rindex;
  } else {
    mbs = a->mbs;
    ii = a->i;
    z = zarray;
  }

  if (!a->mult_work) {
    k = (((A->rmap->n)<(A->cmap->n)) ? (A->cmap->n) : (A->rmap->n));
    ierr = (((k+1)*sizeof(PetscScalar) != 0) ? (*PetscTrMalloc)(((k+1)*sizeof(PetscScalar)),1538,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&a->mult_work)) : (*(&a->mult_work) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1538,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  work = a->mult_work;
  for (i=0; i<mbs; i++) {
    n = ii[1] - ii[0]; ii++;
    ncols = n*bs;
    workt = work;
    for (j=0; j<n; j++) {
      xb = x + bs*(*idx++);
      for (k=0; k<bs; k++) workt[k] = xb[k];
      workt += bs;
    }
    if (usecprow) z = zarray + bs*ridx[i];
    { PetscScalar _one = 1.0; PetscBLASInt _ione = 1,_bbs,_bncols; _bbs = bs; _bncols = ncols; dgemv_("N",&(_bbs),&(_bncols),&_one,v,&(_bbs),work,&_ione,&_one,z,&_ione); };

    v += n*bs2;
    if (!usecprow){
      z += bs;
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1558,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&zarray);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1559,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)2.0*a->nz*bs2),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1560,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultHermitianTranspose_SeqBAIJ(Mat A,Vec xx,Vec zz)
{
  PetscScalar zero = 0.0;
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1571; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultHermitianTranspose_SeqBAIJ") && strcmp("MatMultHermitianTranspose_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1571,"MatMultHermitianTranspose_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  ierr = VecSet(zz,zero);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1572,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = MatMultHermitianTransposeAdd_SeqBAIJ(A,xx,zz,zz);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1573,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultTranspose_SeqBAIJ(Mat A,Vec xx,Vec zz)
{
  PetscScalar zero = 0.0;
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1584; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultTranspose_SeqBAIJ") && strcmp("MatMultTranspose_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1584,"MatMultTranspose_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  ierr = VecSet(zz,zero);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1585,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = MatMultTransposeAdd_SeqBAIJ(A,xx,zz,zz);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1586,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultHermitianTransposeAdd_SeqBAIJ(Mat A,Vec xx,Vec yy,Vec zz)

{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *zb,*x,*z,*xb = 0,x1,x2,x3,x4,x5;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,rval,bs=A->rmap->bs,j,n,bs2=a->bs2,*ib,*ridx=0;
  Mat_CompressedRow cprow = a->compressedrow;
  PetscBool usecprow=cprow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1603; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultHermitianTransposeAdd_SeqBAIJ") && strcmp("MatMultHermitianTransposeAdd_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1603,"MatMultHermitianTransposeAdd_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  if (yy != zz) { ierr = VecCopy(yy,zz);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1604,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0); }
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1605,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1606,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = cprow.nrows;
    ii = cprow.i;
    ridx = cprow.rindex;
  } else {
    mbs=a->mbs;
    ii = a->i;
    xb = x;
  }

  switch (bs) {
  case 1:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + ridx[i];
      x1 = xb[0];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j];
        z[rval] += (*v) * x1;
        v++;
      }
      if (!usecprow) xb++;
    }
    break;
  case 2:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 2*ridx[i];
      x1 = xb[0]; x2 = xb[1];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*2;
        z[rval++] += (v[0])*x1 + (v[1])*x2;
        z[rval++] += (v[2])*x1 + (v[3])*x2;
        v += 4;
      }
      if (!usecprow) xb += 2;
    }
    break;
  case 3:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 3*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*3;
        z[rval++] += (v[0])*x1 + (v[1])*x2 + (v[2])*x3;
        z[rval++] += (v[3])*x1 + (v[4])*x2 + (v[5])*x3;
        z[rval++] += (v[6])*x1 + (v[7])*x2 + (v[8])*x3;
        v += 9;
      }
      if (!usecprow) xb += 3;
    }
    break;
  case 4:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 4*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*4;
        z[rval++] += (v[0])*x1 + (v[1])*x2 + (v[2])*x3 + (v[3])*x4;
        z[rval++] += (v[4])*x1 + (v[5])*x2 + (v[6])*x3 + (v[7])*x4;
        z[rval++] += (v[8])*x1 + (v[9])*x2 + (v[10])*x3 + (v[11])*x4;
        z[rval++] += (v[12])*x1 + (v[13])*x2 + (v[14])*x3 + (v[15])*x4;
        v += 16;
      }
      if (!usecprow) xb += 4;
    }
    break;
  case 5:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 5*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      x4 = xb[3]; x5 = xb[4];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*5;
        z[rval++] += (v[0])*x1 + (v[1])*x2 + (v[2])*x3 + (v[3])*x4 + (v[4])*x5;
        z[rval++] += (v[5])*x1 + (v[6])*x2 + (v[7])*x3 + (v[8])*x4 + (v[9])*x5;
        z[rval++] += (v[10])*x1 + (v[11])*x2 + (v[12])*x3 + (v[13])*x4 + (v[14])*x5;
        z[rval++] += (v[15])*x1 + (v[16])*x2 + (v[17])*x3 + (v[18])*x4 + (v[19])*x5;
        z[rval++] += (v[20])*x1 + (v[21])*x2 + (v[22])*x3 + (v[23])*x4 + (v[24])*x5;
        v += 25;
      }
      if (!usecprow) xb += 5;
    }
    break;
  default: {
      PetscInt ncols,k;
      PetscScalar *work,*workt,*xtmp;

      return PetscError(((MPI_Comm)0x44000001),1706,__func__,"baij2.c","src/mat/impls/baij/seq/",56,PETSC_ERROR_INITIAL,"block size larger than 5 is not supported yet");
      if (!a->mult_work) {
        k = (((A->rmap->n)<(A->cmap->n)) ? (A->cmap->n) : (A->rmap->n));
        ierr = (((k+1)*sizeof(PetscScalar) != 0) ? (*PetscTrMalloc)(((k+1)*sizeof(PetscScalar)),1709,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&a->mult_work)) : (*(&a->mult_work) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1709,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
      }
      work = a->mult_work;
      xtmp = x;
      for (i=0; i<mbs; i++) {
        n = ii[1] - ii[0]; ii++;
        ncols = n*bs;
        ierr = PetscMemzero(work,ncols*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1716,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
        if (usecprow) {
          xtmp = x + bs*ridx[i];
        }
        { PetscScalar _one = 1.0; PetscBLASInt _ione = 1,_bbs,_bncols; _bbs = bs; _bncols = ncols; dgemv_("T",&_bbs,&_bncols,&_one,v,&_bbs,xtmp,&_ione,&_one,work,&_ione); };

        v += n*bs2;
        if (!usecprow) xtmp += bs;
        workt = work;
        for (j=0; j<n; j++) {
          zb = z + bs*(*idx++);
          for (k=0; k<bs; k++) zb[k] += workt[k] ;
          workt += bs;
        }
      }
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1733,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1734,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)2.0*a->nz*a->bs2),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1735,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatMultTransposeAdd_SeqBAIJ(Mat A,Vec xx,Vec yy,Vec zz)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscScalar *zb,*x,*z,*xb = 0,x1,x2,x3,x4,x5;
  MatScalar *v;
  PetscErrorCode ierr;
  PetscInt mbs,i,*idx,*ii,rval,bs=A->rmap->bs,j,n,bs2=a->bs2,*ib,*ridx=0;
  Mat_CompressedRow cprow = a->compressedrow;
  PetscBool usecprow=cprow.use;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1751; petscstack->currentsize++; } do { if (strcmp(__func__,"MatMultTransposeAdd_SeqBAIJ") && strcmp("MatMultTransposeAdd_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1751,"MatMultTransposeAdd_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  if (yy != zz) { ierr = VecCopy(yy,zz);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1752,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0); }
  ierr = VecGetArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1753,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1754,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  idx = a->j;
  v = a->a;
  if (usecprow){
    mbs = cprow.nrows;
    ii = cprow.i;
    ridx = cprow.rindex;
  } else {
    mbs=a->mbs;
    ii = a->i;
    xb = x;
  }

  switch (bs) {
  case 1:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + ridx[i];
      x1 = xb[0];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j];
        z[rval] += *v * x1;
        v++;
      }
      if (!usecprow) xb++;
    }
    break;
  case 2:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 2*ridx[i];
      x1 = xb[0]; x2 = xb[1];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*2;
        z[rval++] += v[0]*x1 + v[1]*x2;
        z[rval++] += v[2]*x1 + v[3]*x2;
        v += 4;
      }
      if (!usecprow) xb += 2;
    }
    break;
  case 3:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 3*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*3;
        z[rval++] += v[0]*x1 + v[1]*x2 + v[2]*x3;
        z[rval++] += v[3]*x1 + v[4]*x2 + v[5]*x3;
        z[rval++] += v[6]*x1 + v[7]*x2 + v[8]*x3;
        v += 9;
      }
      if (!usecprow) xb += 3;
    }
    break;
  case 4:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 4*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2]; x4 = xb[3];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*4;
        z[rval++] += v[0]*x1 + v[1]*x2 + v[2]*x3 + v[3]*x4;
        z[rval++] += v[4]*x1 + v[5]*x2 + v[6]*x3 + v[7]*x4;
        z[rval++] += v[8]*x1 + v[9]*x2 + v[10]*x3 + v[11]*x4;
        z[rval++] += v[12]*x1 + v[13]*x2 + v[14]*x3 + v[15]*x4;
        v += 16;
      }
      if (!usecprow) xb += 4;
    }
    break;
  case 5:
    for (i=0; i<mbs; i++) {
      if (usecprow) xb = x + 5*ridx[i];
      x1 = xb[0]; x2 = xb[1]; x3 = xb[2];
      x4 = xb[3]; x5 = xb[4];
      ib = idx + ii[0];
      n = ii[1] - ii[0]; ii++;
      for (j=0; j<n; j++) {
        rval = ib[j]*5;
        z[rval++] += v[0]*x1 + v[1]*x2 + v[2]*x3 + v[3]*x4 + v[4]*x5;
        z[rval++] += v[5]*x1 + v[6]*x2 + v[7]*x3 + v[8]*x4 + v[9]*x5;
        z[rval++] += v[10]*x1 + v[11]*x2 + v[12]*x3 + v[13]*x4 + v[14]*x5;
        z[rval++] += v[15]*x1 + v[16]*x2 + v[17]*x3 + v[18]*x4 + v[19]*x5;
        z[rval++] += v[20]*x1 + v[21]*x2 + v[22]*x3 + v[23]*x4 + v[24]*x5;
        v += 25;
      }
      if (!usecprow) xb += 5;
    }
    break;
  default: {
      PetscInt ncols,k;
      PetscScalar *work,*workt,*xtmp;

      if (!a->mult_work) {
        k = (((A->rmap->n)<(A->cmap->n)) ? (A->cmap->n) : (A->rmap->n));
        ierr = (((k+1)*sizeof(PetscScalar) != 0) ? (*PetscTrMalloc)(((k+1)*sizeof(PetscScalar)),1856,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&a->mult_work)) : (*(&a->mult_work) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1856,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
      }
      work = a->mult_work;
      xtmp = x;
      for (i=0; i<mbs; i++) {
        n = ii[1] - ii[0]; ii++;
        ncols = n*bs;
        ierr = PetscMemzero(work,ncols*sizeof(PetscScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1863,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
        if (usecprow) {
          xtmp = x + bs*ridx[i];
        }
        { PetscScalar _one = 1.0; PetscBLASInt _ione = 1,_bbs,_bncols; _bbs = bs; _bncols = ncols; dgemv_("T",&_bbs,&_bncols,&_one,v,&_bbs,xtmp,&_ione,&_one,work,&_ione); };

        v += n*bs2;
        if (!usecprow) xtmp += bs;
        workt = work;
        for (j=0; j<n; j++) {
          zb = z + bs*(*idx++);
          for (k=0; k<bs; k++) zb[k] += workt[k] ;
          workt += bs;
        }
      }
    }
  }
  ierr = VecRestoreArray(xx,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1880,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecRestoreArray(zz,&z);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1881,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)2.0*a->nz*a->bs2),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1882,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatScale_SeqBAIJ(Mat inA,PetscScalar alpha)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)inA->data;
  PetscInt totalnz = a->bs2*a->nz;
  PetscScalar oalpha = alpha;
  PetscErrorCode ierr;
  PetscBLASInt one = 1,tnz = totalnz;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1896; petscstack->currentsize++; } do { if (strcmp(__func__,"MatScale_SeqBAIJ") && strcmp("MatScale_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1896,"MatScale_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  dscal_(&tnz,&oalpha,a->a,&one);
  ierr = (_TotalFlops += 1.0*((PetscLogDouble)totalnz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1898,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatNorm_SeqBAIJ(Mat A,NormType type,PetscReal *norm)
{
  PetscErrorCode ierr;
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  MatScalar *v = a->a;
  PetscReal sum = 0.0;
  PetscInt i,j,k,bs=A->rmap->bs,nz=a->nz,bs2=a->bs2,k1;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1912; petscstack->currentsize++; } do { if (strcmp(__func__,"MatNorm_SeqBAIJ") && strcmp("MatNorm_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1912,"MatNorm_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  if (type == NORM_FROBENIUS) {
    for (i=0; i< bs2*nz; i++) {



      sum += (*v)*(*v); v++;

    }
    *norm = sqrt(sum);
  } else if (type == NORM_1) {
    PetscReal *tmp;
    PetscInt *bcol = a->j;
    ierr = (((A->cmap->n+1)*sizeof(PetscReal) != 0) ? (*PetscTrMalloc)(((A->cmap->n+1)*sizeof(PetscReal)),1925,__func__,"baij2.c","src/mat/impls/baij/seq/",(void**)(&tmp)) : (*(&tmp) = 0,0) );do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1925,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = PetscMemzero(tmp,A->cmap->n*sizeof(PetscReal));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1926,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    for (i=0; i<nz; i++){
      for (j=0; j<bs; j++){
        k1 = bs*(*bcol) + j;
        for (k=0; k<bs; k++){
          tmp[k1] += PetscAbsScalar(*v); v++;
        }
      }
      bcol++;
    }
    *norm = 0.0;
    for (j=0; j<A->cmap->n; j++) {
      if (tmp[j] > *norm) *norm = tmp[j];
    }
    ierr = ((tmp) && ((*PetscTrFree)((void*)(tmp),1940,__func__,"baij2.c","src/mat/impls/baij/seq/") || ((tmp) = 0,0)));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1940,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else if (type == NORM_INFINITY) {
    *norm = 0.0;
    for (k=0; k<bs; k++) {
      for (j=0; j<a->mbs; j++) {
        v = a->a + bs2*a->i[j] + k;
        sum = 0.0;
        for (i=0; i<a->i[j+1]-a->i[j]; i++) {
          for (k1=0; k1<bs; k1++){
            sum += PetscAbsScalar(*v);
            v += bs;
          }
        }
        if (sum > *norm) *norm = sum;
      }
    }
  } else return PetscError(((MPI_Comm)0x44000001),1956,__func__,"baij2.c","src/mat/impls/baij/seq/",56,PETSC_ERROR_INITIAL,"No support for this norm yet");
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}




PetscErrorCode MatEqual_SeqBAIJ(Mat A,Mat B,PetscBool * flg)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ *)A->data,*b = (Mat_SeqBAIJ *)B->data;
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 1968; petscstack->currentsize++; } do { if (strcmp(__func__,"MatEqual_SeqBAIJ") && strcmp("MatEqual_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",1968,"MatEqual_SeqBAIJ","__func__",__func__); } } while (0); } while (0);

  if ((A->rmap->N != B->rmap->N) || (A->cmap->n != B->cmap->n) || (A->rmap->bs != B->rmap->bs)|| (a->nz != b->nz)) {
    *flg = PETSC_FALSE;
    do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
  }


  ierr = PetscMemcmp(a->i,b->i,(a->mbs+1)*sizeof(PetscInt),flg);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1976,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (!*flg) {
    do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
  }


  ierr = PetscMemcmp(a->j,b->j,(a->nz)*sizeof(PetscInt),flg);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1982,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (!*flg) {
    do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
  }

  ierr = PetscMemcmp(a->a,b->a,(a->nz)*(A->rmap->bs)*(B->rmap->bs)*sizeof(PetscScalar),flg);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),1987,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);

}



PetscErrorCode MatGetDiagonal_SeqBAIJ(Mat A,Vec v)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode ierr;
  PetscInt i,j,k,n,row,bs,*ai,*aj,ambs,bs2;
  PetscScalar *x,zero = 0.0;
  MatScalar *aa,*aa_j;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 2002; petscstack->currentsize++; } do { if (strcmp(__func__,"MatGetDiagonal_SeqBAIJ") && strcmp("MatGetDiagonal_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",2002,"MatGetDiagonal_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  if (A->factortype) return PetscError(((MPI_Comm)0x44000001),2003,__func__,"baij2.c","src/mat/impls/baij/seq/",73,PETSC_ERROR_INITIAL,"Not for factored matrix");
  bs = A->rmap->bs;
  aa = a->a;
  ai = a->i;
  aj = a->j;
  ambs = a->mbs;
  bs2 = a->bs2;

  ierr = VecSet(v,zero);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2011,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetArray(v,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2012,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = VecGetLocalSize(v,&n);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2013,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (n != A->rmap->N) return PetscError(((MPI_Comm)0x44000001),2014,__func__,"baij2.c","src/mat/impls/baij/seq/",60,PETSC_ERROR_INITIAL,"Nonconforming matrix and vector");
  for (i=0; i<ambs; i++) {
    for (j=ai[i]; j<ai[i+1]; j++) {
      if (aj[j] == i) {
        row = i*bs;
        aa_j = aa+j*bs2;
        for (k=0; k<bs2; k+=(bs+1),row++) x[row] = aa_j[k];
        break;
      }
    }
  }
  ierr = VecRestoreArray(v,&x);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2025,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}



PetscErrorCode MatDiagonalScale_SeqBAIJ(Mat A,Vec ll,Vec rr)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  const PetscScalar *l,*r,*li,*ri;
  PetscScalar x;
  MatScalar *aa, *v;
  PetscErrorCode ierr;
  PetscInt i,j,k,lm,rn,M,m,n,mbs,tmp,bs,bs2,iai;
  const PetscInt *ai,*aj;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 2041; petscstack->currentsize++; } do { if (strcmp(__func__,"MatDiagonalScale_SeqBAIJ") && strcmp("MatDiagonalScale_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",2041,"MatDiagonalScale_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  ai = a->i;
  aj = a->j;
  aa = a->a;
  m = A->rmap->n;
  n = A->cmap->n;
  bs = A->rmap->bs;
  mbs = a->mbs;
  bs2 = a->bs2;
  if (ll) {
    ierr = VecGetArrayRead(ll,&l);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2051,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = VecGetLocalSize(ll,&lm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2052,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    if (lm != m) return PetscError(((MPI_Comm)0x44000001),2053,__func__,"baij2.c","src/mat/impls/baij/seq/",60,PETSC_ERROR_INITIAL,"Left scaling vector wrong length");
    for (i=0; i<mbs; i++) {
      M = ai[i+1] - ai[i];
      li = l + i*bs;
      v = aa + bs2*ai[i];
      for (j=0; j<M; j++) {
        for (k=0; k<bs2; k++) {
          (*v++) *= li[k%bs];
        }
      }
    }
    ierr = VecRestoreArrayRead(ll,&l);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2064,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = (_TotalFlops += 1.0*((PetscLogDouble)a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2065,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }

  if (rr) {
    ierr = VecGetArrayRead(rr,&r);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2069,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = VecGetLocalSize(rr,&rn);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2070,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    if (rn != n) return PetscError(((MPI_Comm)0x44000001),2071,__func__,"baij2.c","src/mat/impls/baij/seq/",60,PETSC_ERROR_INITIAL,"Right scaling vector wrong length");
    for (i=0; i<mbs; i++) {
      iai = ai[i];
      M = ai[i+1] - iai;
      v = aa + bs2*iai;
      for (j=0; j<M; j++) {
        ri = r + bs*aj[iai+j];
        for (k=0; k<bs; k++) {
          x = ri[k];
          for (tmp=0; tmp<bs; tmp++) v[tmp] *= x;
          v += bs;
        }
      }
    }
    ierr = VecRestoreArrayRead(rr,&r);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2085,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    ierr = (_TotalFlops += 1.0*((PetscLogDouble)a->nz),0);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2086,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}




PetscErrorCode MatGetInfo_SeqBAIJ(Mat A,MatInfoType flag,MatInfo *info)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 2098; petscstack->currentsize++; } do { if (strcmp(__func__,"MatGetInfo_SeqBAIJ") && strcmp("MatGetInfo_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",2098,"MatGetInfo_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  info->block_size = a->bs2;
  info->nz_allocated = a->bs2*a->maxnz;
  info->nz_used = a->bs2*a->nz;
  info->nz_unneeded = (double)(info->nz_allocated - info->nz_used);
  info->assemblies = A->num_ass;
  info->mallocs = A->info.mallocs;
  info->memory = ((PetscObject)A)->mem;
  if (A->factortype) {
    info->fill_ratio_given = A->info.fill_ratio_given;
    info->fill_ratio_needed = A->info.fill_ratio_needed;
    info->factor_mallocs = A->info.factor_mallocs;
  } else {
    info->fill_ratio_given = 0;
    info->fill_ratio_needed = 0;
    info->factor_mallocs = 0;
  }
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}




PetscErrorCode MatZeroEntries_SeqBAIJ(Mat A)
{
  Mat_SeqBAIJ *a = (Mat_SeqBAIJ*)A->data;
  PetscErrorCode ierr;

  do { if (petscstack && (petscstack->currentsize < 64)) { petscstack->function[petscstack->currentsize] = __func__; petscstack->file[petscstack->currentsize] = "baij2.c"; petscstack->directory[petscstack->currentsize] = "src/mat/impls/baij/seq/"; petscstack->line[petscstack->currentsize] = 2126; petscstack->currentsize++; } do { if (strcmp(__func__,"MatZeroEntries_SeqBAIJ") && strcmp("MatZeroEntries_SeqBAIJ","User provided function")) { (*PetscErrorPrintf)("%s%s:%d: __FUNCT__=\"%s\" does not agree with %s=\"%s\"\n","src/mat/impls/baij/seq/","baij2.c",2126,"MatZeroEntries_SeqBAIJ","__func__",__func__); } } while (0); } while (0);
  ierr = PetscMemzero(a->a,a->bs2*a->i[a->mbs]*sizeof(MatScalar));do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2127,__func__,"baij2.c","src/mat/impls/baij/seq/",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  do { if (petscstack && petscstack->currentsize > 0) { petscstack->currentsize--; petscstack->function[petscstack->currentsize] = 0; petscstack->file[petscstack->currentsize] = 0; petscstack->directory[petscstack->currentsize] = 0; petscstack->line[petscstack->currentsize] = 0; } return(0);} while (0);
}
