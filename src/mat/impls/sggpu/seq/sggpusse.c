/// SGGPU Matrix Type

#define PETSCMAT_DLL

#include "petsc-private/matimpl.h"

#ifdef __AVX__ //Use 256 AVX intrinsics
#include <immintrin.h>
#define _VEC4
#elif defined(__SSE2__) //Use 128 bit SSE intrinsics
#include <emmintrin.h>
#define _VEC2
#else
#define _VEC1
#endif

#define BETAMAX      512

// DEBUGGING
#define _TIME 0
#define _PRINT_UPDATE_RESULT 0

#define LOAD_TRI_MATRIX4(hd,offset) \
  mc0 = _mm_loadu_pd( hd + offset );		\
  mc1 = _mm_loadu_pd( hd + offset + 2 );	\
  mc2 = _mm_loadu_pd( hd + offset + 4 );	\
  mc3 = _mm_loadu_pd( hd + offset + 6 );	\
  mc4 = _mm_loadu_pd( hd + offset + 8 );	\
  mc5 = _mm_loadu_pd( hd + offset + 10 );	\
  mc6 = _mm_loadu_pd( hd + offset + 12 );	\
  mc7 = _mm_loadu_pd( hd + offset + 14 )

#define LOAD_RHS4(rhs,offset) \
  msum0 = _mm_loadu_pd( rhs + (offset) ); \
  msum1 = _mm_loadu_pd( rhs + (offset) + 2 )

#define LOAD_CLONE_KNOWN4(x,b) \
  xk0 = _mm_load1_pd ( &(x[(b)]) ); \
  xk1 = _mm_load1_pd ( &(x[(b) + 1]) ); \
  xk2 = _mm_load1_pd ( &(x[(b) + 2]) ); \
  xk3 = _mm_load1_pd ( &(x[(b) + 3]) )

#define MUL_SUB_BLOCK						\
	  msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc0, xk0));	\
	  msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc1, xk0));	\
	  msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc2, xk1));	\
	  msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc3, xk1));	\
	  msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc4, xk2));	\
	  msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc5, xk2));	\
	  msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc6, xk3));	\
	  msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc7, xk3))

#define LOAD_LOWER_MAIN_DIAG_MATRIX4(hd,offset)		\
	  mc0 = _mm_loadu_pd( hd + offset );		\
	  mc1 = _mm_loadu_pd( hd + offset + 2 );	\
	  mc2 = _mm_loadu_pd( hd + offset + 4 );	\
	  mc3 = _mm_loadu_pd( hd + offset + 6 );	\
	  mc5 = _mm_loadu_pd( hd + offset + 10 );	\
	  mc7 = _mm_loadu_pd( hd + offset + 14 )

#define LOAD_UPPER_MAIN_DIAG_MATRIX4(hd,offset)		\
	  mc0 = _mm_loadu_pd( hd + offset );		\
	  mc2 = _mm_loadu_pd( hd + offset + 4 );	\
	  mc4 = _mm_loadu_pd( hd + offset + 8 );	\
	  mc5 = _mm_loadu_pd( hd + offset + 10 );	\
	  mc6 = _mm_loadu_pd( hd + offset + 12 );	\
	  mc7 = _mm_loadu_pd( hd + offset + 14 )

//------------------------------------------------------
// general timer function using unix system call
// dlowell ANL-MCS
//------------------------------------------------------
double getclock() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU_sse4"
//-------------------------------------------------
//
// sse4 Version -- for dof==4 use vector intrinsics
//
//-------------------------------------------------
PetscErrorCode MatSolve_SeqSGGPU_sse4(int m, int n, int p, int dof, int dim, PetscScalar *hostData,
				      int numDiags, int *diagOffsets, PetscScalar *x, const PetscScalar *b)
{
  int          j, k, nz;
  int          gridSize = m * n * p;
  int          numElements = gridSize * dof;
  int          chunkSize, chunk, numChunks;
  int          rowCoord, block, colCoord;
  int          sbIndex, di;
  int          offset0, offset1, offset2;
  PetscScalar  beta[BETAMAX];
  double       dbg128d[4];

  // vector register values
  __m128d xk0, xk1, xk2, xk3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, msum0, msum1, msum2, msum3;

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif

  // set chunkSize -- this determines how many rows of the solve to do as a unit
  chunkSize = m;

#if _PRINT_UPDATE_RESULT
  FILE *dbgfile = fopen("/home/jeisenlohr/RF/petsc-rnet/sse4-debug.txt","a");
  fprintf(dbgfile,"Input:\n");
  for ( j = 0; j < numElements; j++ )
    fprintf(dbgfile,"b[%d] = %f\n", j, b[j]);
  fprintf(dbgfile,"\n\nLOWER SOLVE\n");
#endif

  // LOWER TRIANGULAR SOLVE
  // determine the diagonals that contribute already solved values to the RHS for this thread
  int center_diag = numDiags / 2;
  int max_solved_diag = center_diag;

  // initialize the RHS with the components of b
  // PetscMemcpy( iluSolveVec, b, numElements*sizeof(PetscScalar) );    
  PetscMemcpy( x, b, numElements*sizeof(PetscScalar) );    

  // PetscScalar sum[MAXBLOCKSIZE];

  numChunks = ((numElements/dof) + (chunkSize-1)) / chunkSize;
#if _TRACE
  printf("matSize = %d, chunkSize = %d, numChunks = %d, dof = %d\n", numElements, chunkSize, numChunks, dof);
#endif

  // proceed a chunk of rows at a time
  // result of LT solve is stored in iluSolveVec
  while ( (0 <= max_solved_diag) && ( diagOffsets[max_solved_diag] > -chunkSize ) )
    --max_solved_diag;

  for ( chunk = 0; chunk < numChunks; chunk++ ) {
    // offset into portion of iluSolveVec being solved for this chunk
    offset2 = (chunk*chunkSize)*dof;

    //---------------------------------------------
    // Update the RHS with elements of the solution
    // vector solved in previous chunks
    //---------------------------------------------
    for ( di = 0; di < center_diag; di++ ) {  //for ( int di = 0; di <= max_solved_diag; di++ ) {
      int d = diagOffsets[di];
      int startBlock = chunk*chunkSize + d >= 0 ? 0 : -(chunk*chunkSize + d);
      int endBlock = chunkSize > -d ? -d : chunkSize;

      // offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
      offset1 = (chunk*chunkSize + d)*dof;
      offset0 = di*dof*numElements + chunk*chunkSize*dof;
	
      for ( block = startBlock; block < endBlock; block++ ) {

	// load and clone the four components of the known block into 4 registers
	xk0 = _mm_load1_pd ( &(x[offset1 + block*dof]) );
	xk1 = _mm_load1_pd ( &(x[offset1 + block*dof + 1]) );
	xk2 = _mm_load1_pd ( &(x[offset1 + block*dof + 2]) );
	xk3 = _mm_load1_pd ( &(x[offset1 + block*dof + 3]) );

	// load the 4x4 block into 8 128 bit vec registers
	mc0 = _mm_loadu_pd( &(hostData[offset0 + block*dof]) );
	mc1 = _mm_loadu_pd( &(hostData[offset0 + block*dof + 2]) );
	mc2 = _mm_loadu_pd( &(hostData[offset0 + block*dof + numElements]) );
	mc3 = _mm_loadu_pd( &(hostData[offset0 + block*dof + numElements + 2]) );
	mc4 = _mm_loadu_pd( &(hostData[offset0 + block*dof + 2*numElements]) );
	mc5 = _mm_loadu_pd( &(hostData[offset0 + block*dof + 2*numElements + 2]) );
	mc6 = _mm_loadu_pd( &(hostData[offset0 + block*dof + 3*numElements]) );
	mc7 = _mm_loadu_pd( &(hostData[offset0 + block*dof + 3*numElements + 2]) );

	// load the 4 components of the RHS vector into registers
	msum0 = _mm_loadu_pd ( x + offset2 + block*dof );
	msum1 = _mm_loadu_pd ( x + offset2 + block*dof + 2 );
	
	// multiply the matrix by the vector
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc0, xk0));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc1, xk0));
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc2, xk1));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc3, xk1));			      				       
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc4, xk2));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc5, xk2));
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc6, xk3));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc7, xk3));

	_mm_storeu_pd(x + offset2 + block*dof     , msum0);
	_mm_storeu_pd(x + offset2 + block*dof + 2 , msum1);
      }
    }


    for ( block = 0; block < chunkSize; block++ ) {
      offset0 = center_diag*dof*numElements + (chunk*chunkSize + block)*dof;
      offset2 = (chunk*chunkSize + block)*dof;
      // solve diagBlock * x = beta

      // load and clone the four components of the RHS into 4 registers
      xk0 = _mm_load1_pd ( &(x[offset2]) );
      xk1 = _mm_load1_pd ( &(x[offset2 + 1]) );
      xk2 = _mm_load1_pd ( &(x[offset2 + 2]) );
      xk3 = _mm_load1_pd ( &(x[offset2 + 3]) );

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"block %d, modified RHS before solve\n",chunk*chunkSize + block);
    fprintf(dbgfile,"modRHS[%d] = (%f,%f,%f,%f)\n",
	    block,x[offset2], x[offset2 + 1], x[offset2 + 2], x[offset2 + 3] );
#endif

      // load the 4x4 block into 8 128 bit vec registers
      mc0 = _mm_loadu_pd( &(hostData[offset0]) );
      mc1 = _mm_loadu_pd( &(hostData[offset0 + 2]) );

      mc2 = _mm_loadu_pd( &(hostData[offset0 + numElements]) );
      mc3 = _mm_loadu_pd( &(hostData[offset0 + numElements + 2]) );
      //mc4 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements]) );
      mc5 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements + 2]) );
      //mc6 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements]) );
      mc7 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements + 2]) );

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"block %d, main diag matrix block RHS\n",block);
    _mm_store_pd(dbg128d, mc0);
    _mm_store_pd(dbg128d+2, mc1);
    fprintf(dbgfile,"col1[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
    _mm_store_pd(dbg128d, mc2);
    _mm_store_pd(dbg128d+2, mc3);
    fprintf(dbgfile,"col2[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
    _mm_store_pd(dbg128d, mc5);
    _mm_store_pd(dbg128d+2, mc7);
    fprintf(dbgfile,"col3+4[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
#endif

      // multiply the matrix by the vector
      msum0 = _mm_setzero_pd();
      msum1 = _mm_setzero_pd();
	
      msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc0, xk0));
      msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc1, xk0));
      msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc2, xk1));
      msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc3, xk1));			      				       
      //msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc4, xk2));
      msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc5, xk2));
      //msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc6, xk3));
      msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc7, xk3));

      // get and clone the four components of the block for which we just solved into 4 registers
      // do the update (below) before storing back to x+offset2[0..3]
      xk0 = _mm_unpacklo_pd( msum0, msum0 );
      xk1 = _mm_unpackhi_pd( msum0, msum0 );
      xk2 = _mm_unpacklo_pd( msum1, msum1 );
      xk3 = _mm_unpackhi_pd( msum1, msum1 );

      // update blocks below this diag block
      for ( di = max_solved_diag+1; di < center_diag; di++ ) {
	int d = diagOffsets[di];
	// is this thread in a diagonal block?
	if ( block - d < chunkSize ) {
	  offset0 = di*dof*numElements + (chunk*chunkSize + block - d)*dof;
	  offset1 = (chunk*chunkSize + block - d)*dof;

	  // load the appropriate 4x4 block of this diagonal into 8 128 bit vec registers
	  mc0 = _mm_loadu_pd( &(hostData[offset0]) );
	  mc1 = _mm_loadu_pd( &(hostData[offset0 + 2]) );
	  mc2 = _mm_loadu_pd( &(hostData[offset0 + numElements]) );
	  mc3 = _mm_loadu_pd( &(hostData[offset0 + numElements + 2]) );
	  mc4 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements]) );
	  mc5 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements + 2]) );
	  mc6 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements]) );
	  mc7 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements + 2]) );

	  // load the 4 components of the RHS vector into registers
	  msum2 = _mm_loadu_pd ( x + offset1 );
	  msum3 = _mm_loadu_pd ( x + offset1 + 2 );
	
	  // multiply the matrix by the vector
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc0, xk0));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc1, xk0));
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc2, xk1));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc3, xk1));			      				       
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc4, xk2));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc5, xk2));
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc6, xk3));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc7, xk3));

	  // store the updated RHS components
	  _mm_storeu_pd(x + offset1     , msum2);
	  _mm_storeu_pd(x + offset1 + 2 , msum3);
	}
      }

      // now store the values for which we just solved
      _mm_storeu_pd(x + offset2     , msum0);
      _mm_storeu_pd(x + offset2 + 2 , msum1);

#if _PRINT_UPDATE_RESULT
	for ( colCoord = 0; colCoord < dof; colCoord++ )
	  fprintf(dbgfile,"%f\n", x[offset2 + colCoord]);
#endif

    }
  }

#if _PRINT_UPDATE_RESULT
  fprintf(dbgfile,"\n  ***************  UPPER SOLVE  ***************\n\n");
#endif
  //-----------------------
  // UPPER TRIANGULAR SOLVE
  //-----------------------
  // PetscMemcpy( x, iluSolveVec, numElements*sizeof(PetscScalar) );    
  
  // proceed a chunk of rows at a time
  int min_solved_diag = center_diag+1;
  while ( (min_solved_diag < numDiags) && (diagOffsets[min_solved_diag] < chunkSize ) )
    ++min_solved_diag;

  for ( chunk = numChunks-1; chunk >= 0; chunk-- ) {

    //---------------------------------------------
    // Update the RHS with elements of the solution
    // vector solved in previous chunks
    //---------------------------------------------
    for ( di = center_diag+1; di < numDiags; di++ ) {
      int d = diagOffsets[di];
      int startBlock = d < chunkSize ? chunkSize - d : 0;
      int endBlock = (chunk+1)*chunkSize -1 + d < gridSize ? chunkSize - 1 : gridSize - chunk*chunkSize - 1 - d;
	
      for ( block = endBlock; block >= startBlock; block-- ) {

	// offset into matrix coefficients to top of stripe in chunk
	offset0 = ((di+1)*numElements + chunk*chunkSize + block) * dof;
	// offset into already solved portion of iluSolveVec needed for this chunk/diagonal combination 
	offset1 = (chunk*chunkSize + block + d) * dof;
	// offset into portion of x being solved for this chunk
	offset2 = (chunk*chunkSize + block)*dof;

	//for ( colCoord = 0; colCoord < dof; colCoord++ ) {
	//  x[offset2 + colCoord] -= hostData[ offset0 + colCoord ] * x[offset1 + rowCoord];
	//}

	// load and clone the four components of the known block into 4 registers
	xk0 = _mm_load1_pd ( &(x[offset1]) );
	xk1 = _mm_load1_pd ( &(x[offset1 + 1]) );
	xk2 = _mm_load1_pd ( &(x[offset1 + 2]) );
	xk3 = _mm_load1_pd ( &(x[offset1 + 3]) );

	// load the 4x4 block into 8 128 bit vec registers
	mc0 = _mm_loadu_pd( &(hostData[offset0]) );
	mc1 = _mm_loadu_pd( &(hostData[offset0 + 2]) );
	mc2 = _mm_loadu_pd( &(hostData[offset0 + numElements]) );
	mc3 = _mm_loadu_pd( &(hostData[offset0 + numElements + 2]) );
	mc4 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements]) );
	mc5 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements + 2]) );
	mc6 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements]) );
	mc7 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements + 2]) );

	// load the 4 components of the RHS vector into registers
	msum0 = _mm_loadu_pd ( x + offset2 );
	msum1 = _mm_loadu_pd ( x + offset2 + 2 );
	
	// multiply the matrix by the vector
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc0, xk0));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc1, xk0));
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc2, xk1));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc3, xk1));			      				       
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc4, xk2));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc5, xk2));
	msum0 = _mm_sub_pd(msum0 , _mm_mul_pd(mc6, xk3));
	msum1 = _mm_sub_pd(msum1 , _mm_mul_pd(mc7, xk3));

	_mm_storeu_pd(x + offset2     , msum0);
	_mm_storeu_pd(x + offset2 + 2 , msum1);

#if _PRINT_UPDATE_RESULT
	for ( colCoord = 0; colCoord < dof; colCoord++ )
	  fprintf(dbgfile,"%f\n", x[offset2 + colCoord]);
#endif

      }

    }

#if 0
    for ( block = chunkSize-1; block >= 0; block-- ) {
      offset0 = (center_diag+1)*dof*numElements + (chunk*chunkSize + block)*dof;
      offset2 = (chunk*chunkSize + block)*dof;
      // solve diagBlock * x = beta
      PetscMemcpy( beta, &(x[offset2]), dof*sizeof(PetscScalar) );
      for ( j = dof-1; j >= 0; j-- ) {
	x[offset2 + j] *= hostData[offset0 + j*numElements + j];
	for ( k = j+1; k < dof; k++ )
	  x[offset2 + j] += ( hostData[offset0 + k*numElements + j] * beta[k] );
      }

      // update blocks above this diag block
      for ( di = min_solved_diag-1; di > center_diag; di-- ) {
	int d = diagOffsets[di];
	// is this thread in a diagonal block?
	if ( block - d >= 0 ) {
	  offset0 = (di+1)*dof*numElements + (chunk*chunkSize + block - d)*dof;
	  offset1 = (chunk*chunkSize + block)*dof;
	  offset2 = (chunk*chunkSize + block - d)*dof;
	  for ( j = 0; j < dof; j++ ) {
	    for ( k = 0; k < dof; k++ ) {
	      x[offset2 + j] -= hostData[ offset0 + k*numElements + j ] * x[offset1 + k];
	    }	
	  }
	}
      }
    }
#endif

    for ( block = chunkSize-1; block >= 0; block-- ) {
      offset0 = ( (center_diag+1)*numElements + chunk*chunkSize + block ) * dof;
      offset2 = (chunk*chunkSize + block)*dof;

      // solve diagBlock * x = beta

      // load and clone the four components of the RHS into 4 registers
      xk0 = _mm_load1_pd ( &(x[offset2]) );
      xk1 = _mm_load1_pd ( &(x[offset2 + 1]) );
      xk2 = _mm_load1_pd ( &(x[offset2 + 2]) );
      xk3 = _mm_load1_pd ( &(x[offset2 + 3]) );

      // load the 4x4 block into 8 128 bit vec registers
      mc0 = _mm_loadu_pd( &(hostData[offset0]) );
      //mc1 = _mm_loadu_pd( &(hostData[offset0 + 2]) );
      mc2 = _mm_loadu_pd( &(hostData[offset0 + numElements]) );
      //mc3 = _mm_loadu_pd( &(hostData[offset0 + numElements + 2]) );
      mc4 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements]) );
      mc5 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements + 2]) );
      mc6 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements]) );
      mc7 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements + 2]) );

      // multiply the matrix by the vector
      msum0 = _mm_setzero_pd();
      msum1 = _mm_setzero_pd();
	
      msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc0, xk0));
      //msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc1, xk0));
      msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc2, xk1));
      //msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc3, xk1));			      				       
      msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc4, xk2));
      msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc5, xk2));
      msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc6, xk3));
      msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc7, xk3));

      // get and clone the four components of the block for which we just solved into 4 registers
      // then store back to x+offset2[0..3]
      xk0 = _mm_unpacklo_pd( msum0, msum0 );
      xk1 = _mm_unpackhi_pd( msum0, msum0 );
      xk2 = _mm_unpacklo_pd( msum1, msum1 );
      xk3 = _mm_unpackhi_pd( msum1, msum1 );
      _mm_storeu_pd(x + offset2     , msum0);
      _mm_storeu_pd(x + offset2 + 2 , msum1);

      // update blocks above this diag block
      for ( di = min_solved_diag-1; di > center_diag; di-- ) {
	int d = diagOffsets[di];
	// is this thread in a diagonal block?
	if ( block - d >= 0 ) {
	  offset0 = ( (di+1)*numElements + chunk*chunkSize + block - d) * dof;
	  offset1 = (chunk*chunkSize + block - d)*dof;
	  // load the appropriate 4x4 block of this diagonal into 8 128 bit vec registers
	  mc0 = _mm_loadu_pd( &(hostData[offset0]) );
	  mc1 = _mm_loadu_pd( &(hostData[offset0 + 2]) );
	  mc2 = _mm_loadu_pd( &(hostData[offset0 + numElements]) );
	  mc3 = _mm_loadu_pd( &(hostData[offset0 + numElements + 2]) );
	  mc4 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements]) );
	  mc5 = _mm_loadu_pd( &(hostData[offset0 + 2*numElements + 2]) );
	  mc6 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements]) );
	  mc7 = _mm_loadu_pd( &(hostData[offset0 + 3*numElements + 2]) );

	  // load the 4 components of the RHS vector into registers
	  msum2 = _mm_loadu_pd ( x + offset1 );
	  msum3 = _mm_loadu_pd ( x + offset1 + 2 );
	
	  // multiply the matrix by the vector
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc0, xk0));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc1, xk0));
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc2, xk1));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc3, xk1));			      				       
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc4, xk2));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc5, xk2));
	  msum2 = _mm_sub_pd(msum2 , _mm_mul_pd(mc6, xk3));
	  msum3 = _mm_sub_pd(msum3 , _mm_mul_pd(mc7, xk3));

	  // store the updated RHS components
	  _mm_storeu_pd(x + offset1     , msum2);
	  _mm_storeu_pd(x + offset1 + 2 , msum3);
	}
      }      
    }

  }

  //ierr = PetscFree(testXcpu); CHKERRQ(ierr);

#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("sse4 solve time %lf\n",elapsed);
#endif
#if _PRINT_UPDATE_RESULT
  fclose(dbgfile);
#endif
  //ierr = PetscLogFlops(2*nz - A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__  
#define __FUNCT__ "MatSolve_SeqSGGPU_sse4_nochunk"
//-------------------------------------------------
//
// sse4 Version -- for dof==4 use vector intrinsics
//
//-------------------------------------------------
PetscErrorCode MatSolve_SeqSGGPU_sse4_nochunk(int m, int n, int p, int dof, int dim, PetscScalar *hostData,
					      int numDiags, int *diagOffsets, PetscScalar *x, const PetscScalar *b)
{
  int          j, k, nz;
  int          gridSize = m * n * p;
  int          numElements = gridSize * dof;
  int          diagSize = numElements * dof;
  int          block, di;
  PetscInt     startDiag,endDiag;
  int          offset0, offset1,halfRowSize,dofdof;
  PetscScalar  sum, beta[BETAMAX];
  PetscScalar *uTri;
  double       dbg128d[4];

  // vector register values
  __m128d xk0, xk1, xk2, xk3, mc0, mc1, mc2, mc3, mc4, mc5, mc6, mc7, msum0, msum1, msum2, msum3;

#if _TIME
  double t_start, t_end, elapsed;
  t_start = getclock();
#endif

#if _PRINT_UPDATE_RESULT
  int colCoord;
  FILE *dbgfile = fopen("/home/jeisenlohr/RF/petsc-rnet/sse4-nochunk-debug.txt","a");
  fprintf(dbgfile,"Input:\n");
  for ( j = 0; j < numElements; j++ )
    fprintf(dbgfile,"b[%d] = %f\n", j, b[j]);
  fprintf(dbgfile,"\n\nLOWER SOLVE\n");
#endif

  // LOWER TRIANGULAR SOLVE
  startDiag = 0;
  endDiag = numDiags / 2;
  dofdof = 16;
  halfRowSize = dofdof * ((numDiags+1) / 2);

  for ( block = 0; block < gridSize; block++ ) {

    // load the 4 components of the RHS vector into registers
    //msum0 = _mm_loadu_pd ( b + block*dof );
    //msum1 = _mm_loadu_pd ( b + block*dof + 2 );
    LOAD_RHS4(b,4*block);

    for ( di = startDiag; di < endDiag; di++ ) {
      // offset to next non-zero block of matrix in this row
      offset0 = block * halfRowSize + dofdof * di;
      // offset to corresponding known block of solution vector
      offset1 = dof*(block + diagOffsets[di]);
      if ( offset1 > -1 ) {
	// for dof=1, this is:
	// sum -= x[offset1] * hostData[offset0];

	// load and clone the four components of the known block into 4 registers
	LOAD_CLONE_KNOWN4(x,offset1);

	// load the 4x4 block into 8 128 bit vec registers
	LOAD_TRI_MATRIX4(hostData,offset0);
	
	// multiply the matrix by the vector and subtract product
	MUL_SUB_BLOCK;
      }
    }

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"block %d, modified RHS before solve\n",block);
    _mm_store_pd(dbg128d, msum0);
    _mm_store_pd(dbg128d+2, msum1);
    fprintf(dbgfile,"modRHS[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
#endif

    // solve the system given by the block on the main diagonal
    // and the updated RHS stored in (msum0,msum1)

    // load and clone (msum0,msum1) into (xk0,xk1,xk2,xk3)
    xk0 = _mm_shuffle_pd(msum0,msum0,_MM_SHUFFLE2(0,0));
    xk1 = _mm_shuffle_pd(msum0,msum0,_MM_SHUFFLE2(1,1));
    xk2 = _mm_shuffle_pd(msum1,msum1,_MM_SHUFFLE2(0,0));
    xk3 = _mm_shuffle_pd(msum1,msum1,_MM_SHUFFLE2(1,1));

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"loaded and cloned RHS:\n");
    _mm_store_pd(dbg128d, xk0);
    fprintf(dbgfile,"xk0[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
    _mm_store_pd(dbg128d, xk1);
    fprintf(dbgfile,"xk1[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
    _mm_store_pd(dbg128d, xk2);
    fprintf(dbgfile,"xk2[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
    _mm_store_pd(dbg128d, xk3);
    fprintf(dbgfile,"xk3[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
#endif

    // offset to the main diagonal block for this row
    offset0 = block * halfRowSize + dofdof * endDiag;

    // load the 4x4 block into 8 128 bit vec registers
    LOAD_LOWER_MAIN_DIAG_MATRIX4(hostData,offset0);

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"block %d, main diag matrix block RHS\n",block);
    _mm_store_pd(dbg128d, mc0);
    _mm_store_pd(dbg128d+2, mc1);
    fprintf(dbgfile,"main diag block offset %d\n",offset0);
    fprintf(dbgfile,"col1[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
    _mm_store_pd(dbg128d, mc2);
    _mm_store_pd(dbg128d+2, mc3);
    fprintf(dbgfile,"col2[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
    _mm_store_pd(dbg128d, mc5);
    _mm_store_pd(dbg128d+2, mc7);
    fprintf(dbgfile,"col3+4[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
#endif

    // multiply the matrix by the vector
    msum0 = _mm_setzero_pd();
    msum1 = _mm_setzero_pd();
	
    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc0, xk0));
    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc1, xk0));
    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc2, xk1));
    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc3, xk1));			      				       
    //msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc4, xk2));
    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc5, xk2));
    //msum1 = _mm_add_pd(msum0 , _mm_mul_pd(mc6, xk3));
    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc7, xk3));

    _mm_storeu_pd(x + block*dof     , msum0);
    _mm_storeu_pd(x + block*dof + 2 , msum1);

#if _PRINT_UPDATE_RESULT
	for ( colCoord = 0; colCoord < dof; colCoord++ )
	  fprintf(dbgfile,"%f\n", x[block*dof + colCoord]);
#endif

  }
 
#if _PRINT_UPDATE_RESULT
  fprintf(dbgfile,"\n  ***************  UPPER SOLVE  ***************\n\n");
#endif
  // UPPER TRIANGULAR SOLVE
  startDiag = endDiag + 1;
  endDiag = numDiags;
  uTri = hostData + startDiag * diagSize;

  for ( block = gridSize-1; block >= 0; block-- ) {

    // load the 4 components of the RHS vector into registers
    //msum0 = _mm_loadu_pd ( x + block*dof );
    //msum1 = _mm_loadu_pd ( x + block*dof + 2 );
    LOAD_RHS4(x,4*block);

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"block %d, original RHS\n",block);
    _mm_store_pd(dbg128d, msum0);
    _mm_store_pd(dbg128d+2, msum1);
    fprintf(dbgfile,"msum1[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
#endif

    for ( di = startDiag; di < endDiag; di++ ) {
      // offset to next non-zero block of matrix in this row
      offset0 = block * halfRowSize + dofdof * (di-startDiag+1);
      // offset to corresponding known block of solution vector
      offset1 = block + diagOffsets[di];
      if ( offset1 < gridSize ) {
	offset1 *= dof;
	// for dof==1, this is:
	// sum -= x[offset1] * hostData[offset0];

	// load and clone the four components of the known block into 4 registers
	LOAD_CLONE_KNOWN4(x,offset1);

	// load the 4x4 block into 8 128 bit vec registers
	LOAD_TRI_MATRIX4(uTri,offset0);

	// multiply the matrix by the vector and subtract product
	MUL_SUB_BLOCK;
      }

#if _PRINT_UPDATE_RESULT
      fprintf(dbgfile,"block %d, diagonal %d\n",
	      block, diagOffsets[di]);
      fprintf(dbgfile,"offset0 %d, offset1 %d, modified RHS\n",
	      offset0, offset1);
      fprintf(dbgfile,"known[%d] = (%f,%f,%f,%f)\n",
	      block + diagOffsets[di],
	      x[offset1],x[offset1+1],x[offset1+2],x[offset1+3]);
      _mm_store_pd(dbg128d, msum0);
      _mm_store_pd(dbg128d+2, msum1);
      fprintf(dbgfile,"msum1[%d] = (%f,%f,%f,%f)\n",
	      block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
#endif
    }

    // solve the system given by the block on the main diagonal
    // and the updated RHS stored in (msum0,msum1)

    // load and clone (msum0,msum1) into (xk0,xk1,xk2,xk3)
    xk0 = _mm_shuffle_pd(msum0,msum0,_MM_SHUFFLE2(0,0));
    xk1 = _mm_shuffle_pd(msum0,msum0,_MM_SHUFFLE2(1,1));
    xk2 = _mm_shuffle_pd(msum1,msum1,_MM_SHUFFLE2(0,0));
    xk3 = _mm_shuffle_pd(msum1,msum1,_MM_SHUFFLE2(1,1));

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"loaded and cloned RHS:\n");
    _mm_store_pd(dbg128d, xk0);
    fprintf(dbgfile,"xk0[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
    _mm_store_pd(dbg128d, xk1);
    fprintf(dbgfile,"xk1[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
    _mm_store_pd(dbg128d, xk2);
    fprintf(dbgfile,"xk2[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
    _mm_store_pd(dbg128d, xk3);
    fprintf(dbgfile,"xk3[%d] = (%f,%f)\n", block,dbg128d[0],dbg128d[1]);
#endif

    // offset to the main diagonal block for this row
    offset0 = block * halfRowSize;

    // load the 4x4 block into 8 128 bit vec registers
    LOAD_UPPER_MAIN_DIAG_MATRIX4(uTri,offset0);

#if _PRINT_UPDATE_RESULT
    fprintf(dbgfile,"block %d, main diag matrix block RHS\n",block);
    _mm_store_pd(dbg128d, mc0);
    _mm_store_pd(dbg128d+2, mc2);
    fprintf(dbgfile,"main diag block offset %d\n",offset0);
    fprintf(dbgfile,"col1+2[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
    _mm_store_pd(dbg128d, mc4);
    _mm_store_pd(dbg128d+2, mc5);
    fprintf(dbgfile,"col3[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
    _mm_store_pd(dbg128d, mc6);
    _mm_store_pd(dbg128d+2, mc7);
    fprintf(dbgfile,"col4[%d] = (%f,%f,%f,%f)\n",
	    block,dbg128d[0],dbg128d[1],block,dbg128d[2],dbg128d[3]);
#endif

    // multiply the matrix by the vector
    msum0 = _mm_setzero_pd();
    msum1 = _mm_setzero_pd();
	
    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc0, xk0));
    //msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc1, xk0));
    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc2, xk1));
    //msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc3, xk1));			      				       
    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc4, xk2));
    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc5, xk2));
    msum0 = _mm_add_pd(msum0 , _mm_mul_pd(mc6, xk3));
    msum1 = _mm_add_pd(msum1 , _mm_mul_pd(mc7, xk3));

    _mm_storeu_pd(x + block*dof     , msum0);
    _mm_storeu_pd(x + block*dof + 2 , msum1);

#if _PRINT_UPDATE_RESULT
    for ( colCoord = 0; colCoord < dof; colCoord++ )
      fprintf(dbgfile,"%f\n", x[block*dof + colCoord]);
#endif

  }


#if _TIME
  t_end = getclock();
  elapsed = t_end - t_start;
  printf("sse4 solve time %lf\n",elapsed);
#endif
#if _PRINT_UPDATE_RESULT
  fclose(dbgfile);
#endif
  //ierr = PetscLogFlops(2*nz - A->cmap->n);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
