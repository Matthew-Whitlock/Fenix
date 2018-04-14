/*
//@HEADER
// ************************************************************************
//
//
//            _|_|_|_|  _|_|_|_|  _|      _|  _|_|_|  _|      _|
//            _|        _|        _|_|    _|    _|      _|  _|
//            _|_|_|    _|_|_|    _|  _|  _|    _|        _|
//            _|        _|        _|    _|_|    _|      _|  _|
//            _|        _|_|_|_|  _|      _|  _|_|_|  _|      _|
//
//
//
//
// Copyright (C) 2016 Rutgers University and Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Author Marc Gamell, Eric Valenzuela, Keita Teranishi, Manish Parashar
//        and Michael Heroux
//
// Questions? Contact Keita Teranishi (knteran@sandia.gov) and
//                    Marc Gamell (mgamell@cac.rutgers.edu)
//
// ************************************************************************
//@HEADER
*/

#if 1

#include "fenix_data_recovery.h"
#include "fenix_opt.h"
#include "fenix_process_recovery.h"
#include "fenix_util.h"
#include "fenix_ext.h"
#include "fenix_metadata.h"

inline void __fenix_init_group_metadata ( fenix_group_entry_t *gentry, int groupid, MPI_Comm comm, int timestamp,
                                    int depth  )
{
   gentry->groupid = groupid;
   gentry->comm = comm;
   gentry->timestart = timestamp;
   gentry->timestamp = timestamp;
   gentry->depth = depth + 1;
   gentry->state = OCCUPIED;
}

inline void __fenix_reinit_group_metadata ( fenix_group_entry_t *gentry  )

{
  gentry->current_rank = __fenix_get_current_rank( gentry->comm );
  gentry->comm_size    = __fenix_get_world_size( gentry->comm );
  gentry->in_rank      = ( gentry->current_rank + gentry->comm_size - gentry->rank_separation ) % gentry->comm_size;
  gentry->out_rank     = ( gentry->current_rank + gentry->comm_size + gentry->rank_separation ) % gentry->comm_size;
}

inline void __fenix_data_member_init_metadata ( fenix_member_entry_t *mentry, int memberid, void *data, int count, MPI_Datatype datatype )

{
    mentry->memberid = memberid;
    mentry->state = OCCUPIED;
    mentry->user_data = data;
    mentry->current_count = count;
    mentry->current_datatype = datatype;
    int dsize;
    MPI_Type_size(datatype, &dsize);

    mentry->datatype_size = mentry->current_size = dsize;
}


inline void __fenix_data_member_init_store_packet ( fenix_member_store_packet_t *lentry_packet, fenix_buffer_entry_t *lentry, int flag )
{
  if ( flag == 0 ) {
    lentry_packet->rank = lentry->origin_rank;
    lentry_packet->datatype = lentry->datatype;
    lentry_packet->entry_count = lentry->count;
    lentry_packet->entry_size  = lentry->datatype_size;
    lentry_packet->entry_real_count  = lentry->count;
    lentry_packet->num_blocks  = 0;
  } else if ( flag == 1 ) {
    lentry_packet->rank = lentry->origin_rank;
    lentry_packet->datatype = lentry->datatype;
    lentry_packet->entry_count = lentry->count;
    lentry_packet->entry_size  = lentry->datatype_size;
    lentry_packet->entry_real_count  = 0;
    lentry_packet->num_blocks  = 0;
  } else if (flag == 2 ) { /* Subset */


  }
}
#endif
