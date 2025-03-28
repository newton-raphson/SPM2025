
/**
  @file parUtils.C
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include "mpi.h"
#include "binUtils.h"
#include "dtypes.h"
#include "parUtils.h"

#ifdef __DEBUG__
#ifndef __DEBUG_PAR__
#define __DEBUG_PAR__
#endif
#endif

namespace par {

  int DBG_rProc, DBG_nProc;

  unsigned int splitCommBinary( MPI_Comm orig_comm, MPI_Comm *new_comm) {
    int npes, rank;

    MPI_Group  orig_group, new_group;

    MPI_Comm_size(orig_comm, &npes);
    MPI_Comm_rank(orig_comm, &rank);

    unsigned int splitterRank = binOp::getPrevHighestPowerOfTwo(npes);

    int *ranksAsc, *ranksDesc;
    //Determine sizes for the 2 groups 
    ranksAsc = new int[splitterRank];
    ranksDesc = new int[( npes - splitterRank)];

    int numAsc = 0;
    int numDesc = ( npes - splitterRank - 1);

    //This is the main mapping between old ranks and new ranks.
    for(int i=0; i<npes; i++) {
      if( static_cast<unsigned int>(i) < splitterRank) {
        ranksAsc[numAsc] = i;
        numAsc++;
      }else {
        ranksDesc[numDesc] = i;
        numDesc--;
      }
    }//end for i

    MPI_Comm_group(orig_comm, &orig_group);

    /* Divide tasks into two distinct groups based upon rank */
    if (static_cast<unsigned int>(rank) < splitterRank) {
      MPI_Group_incl(orig_group, splitterRank, ranksAsc, &new_group);
    }else {
      MPI_Group_incl(orig_group, (npes-splitterRank), ranksDesc, &new_group);
    }

    MPI_Comm_create(orig_comm, new_group, new_comm);

    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);

    delete [] ranksAsc;
    ranksAsc = NULL;
    
    delete [] ranksDesc;
    ranksDesc = NULL;

    return splitterRank;
  }//end function

  unsigned int splitCommBinaryNoFlip( MPI_Comm orig_comm, MPI_Comm *new_comm) {
    int npes, rank;

    MPI_Group  orig_group, new_group;

    MPI_Comm_size(orig_comm, &npes);
    MPI_Comm_rank(orig_comm, &rank);

    unsigned int splitterRank =  binOp::getPrevHighestPowerOfTwo(npes);

    int *ranksAsc, *ranksDesc;
    //Determine sizes for the 2 groups 
    ranksAsc = new int[splitterRank];
    ranksDesc = new int[( npes - splitterRank)];

    int numAsc = 0;
    int numDesc = 0; //( npes - splitterRank - 1);

    //This is the main mapping between old ranks and new ranks.
    for(int i = 0; i < npes; i++) {
      if(static_cast<unsigned int>(i) < splitterRank) {
        ranksAsc[numAsc] = i;
        numAsc++;
      }else {
        ranksDesc[numDesc] = i;
        numDesc++;
      }
    }//end for i

    MPI_Comm_group(orig_comm, &orig_group);

    /* Divide tasks into two distinct groups based upon rank */
    if (static_cast<unsigned int>(rank) < splitterRank) {
      MPI_Group_incl(orig_group, splitterRank, ranksAsc, &new_group);
    }else {
      MPI_Group_incl(orig_group, (npes-splitterRank), ranksDesc, &new_group);
    }

    MPI_Comm_create(orig_comm, new_group, new_comm);


    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);


    delete [] ranksAsc;
    ranksAsc = NULL;
    
    delete [] ranksDesc;
    ranksDesc = NULL;

    return splitterRank;
  }//end function

  //create Comm groups and remove empty processors...
  int splitComm2way(bool iAmEmpty, MPI_Comm * new_comm, MPI_Comm comm) {
#ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_SPLIT_COMM_2WAY_BEGIN

      MPI_Group  orig_group, new_group;
    int size;
    MPI_Comm_size(comm, &size);

    bool* isEmptyList = new bool[size];
    par::Mpi_Allgather<bool>(&iAmEmpty, isEmptyList, 1, comm);

    int numActive=0, numIdle=0;
    for(int i = 0; i < size; i++) {
      if(isEmptyList[i]) {
        numIdle++;
      }else {
        numActive++;
      }
    }//end for i

    int* ranksActive = new int[numActive];
    int* ranksIdle = new int[numIdle];

    numActive=0;
    numIdle=0;
    for(int i = 0; i < size; i++) {
      if(isEmptyList[i]) {
        ranksIdle[numIdle] = i;
        numIdle++;
      }else {
        ranksActive[numActive] = i;
        numActive++;
      }
    }//end for i

    delete [] isEmptyList;	
    isEmptyList = NULL;

    /* Extract the original group handle */
    MPI_Comm_group(comm, &orig_group);

    /* Divide tasks into two distinct groups based upon rank */
    if (!iAmEmpty) {
      MPI_Group_incl(orig_group, numActive, ranksActive, &new_group);
    }else {
      MPI_Group_incl(orig_group, numIdle, ranksIdle, &new_group);
    }

    /* Create new communicator */
    MPI_Comm_create(comm, new_group, new_comm);

    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);

    delete [] ranksActive;
    ranksActive = NULL;
    
    delete [] ranksIdle;
    ranksIdle = NULL;

    PROF_SPLIT_COMM_2WAY_END
  }//end function

  int splitCommUsingSplittingRank(int splittingRank, MPI_Comm* new_comm,
      MPI_Comm comm) {
#ifdef __PROFILE_WITH_BARRIER__
    MPI_Barrier(comm);
#endif
    PROF_SPLIT_COMM_BEGIN

      MPI_Group  orig_group, new_group;
    int size;
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int* ranksActive = new int[splittingRank];
    int* ranksIdle = new int[size - splittingRank];

    for(int i = 0; i < splittingRank; i++) {
      ranksActive[i] = i;
    }

    for(int i = splittingRank; i < size; i++) {
      ranksIdle[i - splittingRank] = i;
    }

    /* Extract the original group handle */
    MPI_Comm_group(comm, &orig_group);

    /* Divide tasks into two distinct groups based upon rank */
    if (rank < splittingRank) {
      MPI_Group_incl(orig_group, splittingRank, ranksActive, &new_group);
    }else {
      MPI_Group_incl(orig_group, (size - splittingRank), ranksIdle, &new_group);
    }

    /* Create new communicator */
    MPI_Comm_create(comm, new_group, new_comm);

    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);

    delete [] ranksActive;
    ranksActive = NULL;
    
    delete [] ranksIdle;
    ranksIdle = NULL;

    PROF_SPLIT_COMM_END
  }//end function

  //create Comm groups and remove empty processors...
  int splitComm2way(const bool* isEmptyList, MPI_Comm * new_comm, MPI_Comm comm) {
      
    MPI_Group  orig_group, new_group;
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int numActive=0, numIdle=0;
    for(int i = 0; i < size; i++) {
      if(isEmptyList[i]) {
        numIdle++;
      }else {
        numActive++;
      }
    }//end for i

    int* ranksActive = new int[numActive];
    int* ranksIdle = new int[numIdle];

    numActive=0;
    numIdle=0;
    for(int i = 0; i < size; i++) {
      if(isEmptyList[i]) {
        ranksIdle[numIdle] = i;
        numIdle++;
      }else {
        ranksActive[numActive] = i;
        numActive++;
      }
    }//end for i

    /* Extract the original group handle */
    MPI_Comm_group(comm, &orig_group);

    /* Divide tasks into two distinct groups based upon rank */
    if (!isEmptyList[rank]) {
      MPI_Group_incl(orig_group, numActive, ranksActive, &new_group);
    }else {
      MPI_Group_incl(orig_group, numIdle, ranksIdle, &new_group);
    }

    /* Create new communicator */
    MPI_Comm_create(comm, new_group, new_comm);

    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);


    delete [] ranksActive;
    ranksActive = NULL;
    
    delete [] ranksIdle;
    ranksIdle = NULL;

    return 0;
  }//end function



  //
  // mpi_next_if
  //
  std::pair<int, int> mpi_next_if(bool included, MPI_Comm comm)
  {
    // @author Masado Ishii
    // predecessor = exclusive_scan[max](included ? rank : -1).
    // successor = reversed_value(reversed(predecessor(reversed pattern))).

    int comm_size, comm_rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &comm_rank);
    const int max_rank = comm_size - 1;

    const auto mpi_reverse = [&](int value) {
      const int reverse_rank = max_rank - comm_rank;
      if (reverse_rank == comm_rank)
        return value;
      int rev_value;
      par::Mpi_Sendrecv(&value, 1, reverse_rank, int{},
                        &rev_value, 1, reverse_rank, int{},
                        comm, MPI_STATUS_IGNORE);
      return rev_value;
    };

    enum { SELF, REVERSE };
    int effective_rank[2];
    effective_rank[SELF] =    included              ? comm_rank : -1;
    effective_rank[REVERSE] = mpi_reverse(included) ? comm_rank : -1;

    int predecessor[2] = {-1, -1};
    par::Mpi_Exscan(effective_rank, predecessor, 2, MPI_MAX, comm);
    int successor = max_rank - mpi_reverse(predecessor[REVERSE]);

    if (predecessor[SELF] == -1)
      predecessor[SELF] = MPI_PROC_NULL;
    if (successor == comm_size)
      successor = MPI_PROC_NULL;

    return {predecessor[SELF], successor};
  }




  void waitForDebugger(MPI_Comm comm, int commRank)
  {
    if (getenv("USE_MPI_DEBUGGER") != NULL && commRank == 0)
    {
      volatile int goAhead = 0;
      fprintf(stderr, "pid %ld waiting for debugger\n", (long) getpid());
      while (goAhead == 0)
      {
        /* Wait for change in `goAhead' by the debugger. */
      }
    }
    MPI_Barrier(comm);
  }


}// end namespace

