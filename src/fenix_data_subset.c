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

#include "mpi.h"
#include "fenix-config.h"
#include "fenix_ext.h"
#include "fenix_data_subset.h"


int __fenix_data_subset_init(int num_blocks, Fenix_Data_subset* subset){
   int retval = -1;
   if(num_blocks <= 0){
      debug_print("ERROR __fenix_data_subset_init: num_regions <%d> must be positive\n",
                num_blocks);
   } else {
      subset->start_offsets = (int*) s_malloc(sizeof(int) * num_blocks);
      subset->end_offsets = (int*) s_malloc(sizeof(int) * num_blocks);
      subset->num_repeats = (int*) s_malloc(sizeof(int) * num_blocks);
      subset->num_blocks = num_blocks;
   }
   return retval;
}


/**
 * @brief
 * @param num_blocks
 * @param start_offset
 * @param end_offset
 * @param stride
 * @param subset_specifier
 *
 * This routine creates 
 */
int __fenix_data_subset_create(int num_blocks, int start_offset, int end_offset, int stride,
                       Fenix_Data_subset *subset_specifier) {
  int retval = -1;
  if (num_blocks <= 0) {
    debug_print("ERROR Fenix_Data_subset_create: num_blocks <%d> must be positive\n",
                num_blocks);
    retval = FENIX_ERROR_SUBSET_NUM_BLOCKS;
  } else if (start_offset < 0) {
    debug_print("ERROR Fenix_Data_subset_create: start_offset <%d> must be positive\n",
                start_offset);
    retval = FENIX_ERROR_SUBSET_START_OFFSET;
  } else if (end_offset <= 0) {
    debug_print("ERROR Fenix_Data_subset_create: end_offset <%d> must be positive\n",
                end_offset);
    retval = FENIX_ERROR_SUBSET_END_OFFSET;
  } else if (stride <= 0) {
    debug_print("ERROR Fenix_Data_subset_create: stride <%d> must be positive\n", stride);
    retval = FENIX_ERROR_SUBSET_STRIDE;
  } else {
    //This is a simple subset with a single region descriptor that simply
    //repeats num_blocks times.
    __fenix_data_subset_init(1 /*Only 1 block, repeated*/, subset_specifier);

    subset_specifier->start_offsets[0] = start_offset;
    subset_specifier->end_offsets[0] = end_offset;
    subset_specifier->num_repeats[0] = num_blocks-1;
    subset_specifier->stride = stride;
    subset_specifier->specifier = __FENIX_SUBSET_CREATE;
    retval = FENIX_SUCCESS;
  }
  return retval;
}

/**
 * @brief
 * @param num_blocks
 * @param array_start_offsets
 * @param array_end_offsets
 * @param subset_specifier
 */
int __fenix_data_subset_createv(int num_blocks, int *array_start_offsets, int *array_end_offsets,
                        Fenix_Data_subset *subset_specifier) {

  int retval = -1;
  if (num_blocks <= 0) {
    debug_print("ERROR Fenix_Data_subset_createv: num_blocks <%d> must be positive\n",
                num_blocks);
    retval = FENIX_ERROR_SUBSET_NUM_BLOCKS;
  } else if (array_start_offsets == NULL) {
    debug_print( "ERROR Fenix_Data_subset_createv: array_start_offsets %s must be at least of size 1\n", "");
    retval = FENIX_ERROR_SUBSET_START_OFFSET;
  } else if (array_end_offsets == NULL) {
    debug_print( "ERROR Fenix_Data_subset_createv: array_end_offsets %s must at least of size 1\n", "");
    retval = FENIX_ERROR_SUBSET_END_OFFSET;
  } else {

    // first check that the start offsets and end offsets are valid
    int index;
    int invalid_index = -1;
    int found_invalid_index = 0;
    for (index = 0; found_invalid_index != 1 && (index < num_blocks); index++) {
      if (array_start_offsets[index] > array_end_offsets[index]) {
        invalid_index = index;
        found_invalid_index = 1;
      }
    }

    if (found_invalid_index != 1) { // if not true (!= 1)
      __fenix_data_subset_init(num_blocks, subset_specifier);

      //Createv type, so each region is never repeated, just occurs once.
      for(int block_index = 0; block_index < num_blocks; block_index++){
         subset_specifier->num_repeats[block_index] = 0;
      }

      memcpy(subset_specifier->start_offsets, array_start_offsets, ( num_blocks * sizeof(int))); // deep copy
      memcpy(subset_specifier->end_offsets, array_end_offsets, ( num_blocks * sizeof(int))); // deep copy
      
      subset_specifier->specifier = __FENIX_SUBSET_CREATEV;
      retval = FENIX_SUCCESS;
    } else {
      debug_print(
              "ERROR Fenix_Data_subset_createv: array_end_offsets[%d] must be less than array_start_offsets[%d]\n",
              invalid_index, invalid_index);
      retval = FENIX_ERROR_SUBSET_END_OFFSET;
    }
  }
  return retval;
}

//This should only be used to copy to a currently non-inited subset
// If the destination already has memory allocated in the num_blocks/offsets regions
// then this can lead to double-mallocs or memory leaks.
void __fenix_data_subset_deep_copy(Fenix_Data_subset* from, Fenix_Data_subset* to){
   if(from->specifier == __FENIX_SUBSET_FULL || from->specifier == __FENIX_SUBSET_EMPTY){
      to->specifier = from->specifier;
   } else {
      __fenix_data_subset_init(from->num_blocks, to);
      memcpy(to->num_repeats, from->num_repeats, to->num_blocks);
      memcpy(to->start_offsets, from->start_offsets, to->num_blocks);
      memcpy(to->end_offsets, from->end_offsets, to->num_blocks);
      to->specifier = from->specifier;
      to->stride = from->stride;
   }
}

//This function checks for any overlapping regions and removes them.
void __fenix_data_subset_simplify_regions(Fenix_Data_subset* ss){
   int space_allocated = ss->num_blocks;
   
   if(ss->specifier == __FENIX_SUBSET_CREATE){
      //We will handle this by viewing the data as regions of size stride.
      //Each block will be broken into a value dictating which regions it is
      //within, and what data within each region it is within.
      //
      //If two blocks do not overlap within regions, there is no overlap.
      //If they overlap within regions, but the regions they touch do not overlap,
      //there is no overlap. etc.
      
      for(int i = 0; i < ss->num_blocks-1; i++){
         int did_merge = 0;
         size_t start_region_i = ss->start_offsets[i] / ss->stride;
         size_t end_region_i = start_region_i + ss->num_repeats[i];
         size_t start_index_i = ss->start_offsets[i] % ss->stride;
         size_t end_index_i = ss->end_offsets[i] % ss->stride;
         
         for(int j = i+1; j < ss->num_blocks; j++){
            size_t start_region_j = ss->start_offsets[j] / ss->stride;
            size_t end_region_j = start_region_j + ss->num_repeats[j];
            size_t start_index_j = ss->start_offsets[j] % ss->stride;
            size_t end_index_j = ss->end_offsets[j] % ss->stride;
            
            //First, do these have the potential to overlap?
            if( !( (start_index_i <= start_index_j && end_index_i >= start_index_j)
                  || (start_index_j <= start_index_i && end_index_j >= start_index_i) ) ){
               //Even if they touch the same regions, there is no overlap.
               continue;
            }

            //Now, do they overlap on which regions they touch?
            if( !( (start_region_i <= start_region_j && end_region_i >= start_region_j)
                  || (start_region_j <= start_region_i && end_region_j >= start_region_i) ) ){
               //They do not touch the same regions ever.
               continue;
            } 
            
            
            //Now, in which regions do these overlap?
            //We will simplify the logic by switching from i and j referencing
            //to viewing the two blocks in the order that they touch regions.
            int first_block; 
            int first_block_start;
            int first_block_end;

            int second_block;
            int second_block_start;
            int second_block_end;

            if(start_region_i < start_region_j){
               first_block = i;
               first_block_start = start_region_i;
               first_block_end = end_region_i;
               second_block = j;
               second_block_start = start_region_j;
               second_block_end = end_region_j;
            } else {
               first_block = j;
               first_block_start = start_region_j;
               first_block_end = end_region_j;
               second_block = i;
               second_block_start = start_region_i;
               second_block_end = end_region_i;
            }

            int length_first_only = first_block_start - second_block_start;
            int length_both = first_block_end < second_block_end ?
                  second_block_start - first_block_end :
                  second_block_start - second_block_end;

            //Could be negative, but that's fixed later.
            int length_second_only = second_block_end - first_block_end;

            
            
            //Now we know what the overlap is, so we make the changes to the data subset.
            ss->num_repeats[first_block] = length_first_only - 1;
            
            ss->num_repeats[second_block] = length_second_only - 1;
            ss->start_offsets[second_block] = ss->start_offsets[second_block] % ss->stride +
                  ss->stride * (first_block_end+1);
            ss->end_offsets[second_block] = ss->end_offsets[second_block] % ss->stride +
                  ss->stride * (first_block_end+1);

            //Now add the overlapped region into the subset
            int merged_dest = -1;
            if(ss->num_repeats[first_block] < 0){
               //Place the overlapping area in the now-nonexistant first block's space.
               merged_dest = first_block;
            } else if(ss->num_repeats[second_block] < 0){
               //Place the overlapping area in the now-nonexistant second block's space.
               merged_dest = second_block;
            } else {
               //Place the overlapping area in a new block.
               merged_dest = ss->num_blocks;

               //Make sure there is space, and increment num_blocks.
               ss->num_blocks++;
               if(ss->num_blocks > space_allocated){
                  
                  ss->end_offsets = (int*) s_realloc(ss->end_offsets,
                                (space_allocated * 2) * sizeof(int));
                  ss->start_offsets = (int*) s_realloc(ss->start_offsets,
                                (space_allocated * 2) * sizeof(int));
                  ss->num_repeats = (int*) s_realloc(ss->num_repeats,
                                (space_allocated * 2) * sizeof(int));
                  space_allocated *= 2;
               }
            }

            ss->num_repeats[merged_dest] = length_both-1;
            ss->start_offsets[merged_dest] = ss->stride*second_block_start +
                  (start_index_i<start_index_j ? start_index_i : start_index_j);
            ss->end_offsets[merged_dest] = ss->stride*second_block_start +
                  (end_index_i>end_index_j ? end_index_i : end_index_j);
            

            //Check if num_repeats[second_block] < 0, if so remove it.
            //This could occur if both blocks can be perfectly minimized to a single block.
            if(ss->num_repeats[second_block] < 0){
               if(second_block == ss->num_blocks-1){
                  //Don't need to move anything.
                  ss->num_blocks--;
               } else {
                  //We need to move everything over by one.
                  memmove(ss->num_repeats + second_block, ss->num_repeats + second_block + 1, 
                        ss->num_blocks - second_block - 1);
                  memmove(ss->start_offsets + second_block, ss->start_offsets + second_block + 1, 
                        ss->num_blocks - second_block - 1);
                  memmove(ss->end_offsets + second_block, ss->end_offsets + second_block + 1, 
                        ss->num_blocks - second_block - 1);
                  ss->num_blocks--;
               }
            } 
            
            //Update this for the new block i. 
            start_region_i = ss->start_offsets[i] / ss->stride;
            end_region_i = start_region_i + ss->num_repeats[i];
            start_index_i = ss->start_offsets[i] % ss->stride;
            end_index_i = ss->end_offsets[i] % ss->stride;
            
            did_merge = 1;
         }

         //If we merged w/ anything, recheck w/ new merged block.
         if(did_merge) i--;
      }
   } else if(ss->specifier == __FENIX_SUBSET_CREATEV){
      //This is much simpler than with CREATE type, since we don't have to
      //worry about repetition.
      for(int i = 0; i < ss->num_blocks-1; i++){
         int did_merge = 0;
         
         for(int j = i+1; j < ss->num_blocks; j++){
            if(   ( ss->start_offsets[i] < ss->start_offsets[j]  &&
                     ss->end_offsets[i] > ss->start_offsets[j] ) 
                  || 
                  ( ss->start_offsets[j] < ss->start_offsets[i] &&
                     ss->end_offsets[j] > ss->start_offsets[i] )){
               did_merge = 1;

               ss->start_offsets[i] = (ss->start_offsets[i] < ss->start_offsets[j]) ?
                     ss->start_offsets[i] :
                     ss->start_offsets[j];

               ss->end_offsets[i] = (ss->end_offsets[i] > ss->end_offsets[j]) ?
                     ss->end_offsets[i] :
                     ss->end_offsets[j];
               
               //Move everything over to remove j
               memmove(ss->num_repeats + j, ss->num_repeats + j + 1, 
                     ss->num_blocks - j - 1);
               memmove(ss->start_offsets + j, ss->start_offsets + j + 1, 
                     ss->num_blocks - j - 1);
               memmove(ss->end_offsets + j, ss->end_offsets + j + 1, 
                     ss->num_blocks - j - 1);
               ss->num_blocks--;
            }
         }

         if(did_merge) i--;
      }
   }
}

//This should only be used to copy to a currently non-inited subset
// If the destination already has memory allocated in the num_blocks/offsets regions
// then this can lead to double-mallocs or memory leaks.
void __fenix_data_subset_merge(Fenix_Data_subset* first_subset, Fenix_Data_subset* second_subset,
      Fenix_Data_subset* output){
      
   //Simple cases first
   if(first_subset->specifier == __FENIX_SUBSET_FULL || 
         second_subset->specifier == __FENIX_SUBSET_FULL){
      //We don't need to populate anything else.
      output->specifier = __FENIX_SUBSET_FULL;
   
   } else if(first_subset->specifier == __FENIX_SUBSET_EMPTY){
      __fenix_data_subset_deep_copy(second_subset, output);
   
   } else if(second_subset->specifier == __FENIX_SUBSET_EMPTY){
      __fenix_data_subset_deep_copy(first_subset, output);

   } else if(first_subset->specifier == __FENIX_SUBSET_CREATE &&
         second_subset->specifier == __FENIX_SUBSET_CREATE &&
         first_subset->stride == second_subset->stride){
      //Output is just a CREATE type with combined descriptors. 
      //Start by making a list of all descriptors, then merge any with overlaps.
      output->stride = first_subset->stride;
      output->num_blocks = first_subset->num_blocks 
         + second_subset->num_blocks;
      __fenix_data_subset_init(output->num_blocks, output);
      output->specifier = __FENIX_SUBSET_CREATE; 
      
      memcpy(output->num_repeats, first_subset->num_repeats, first_subset->num_blocks);
      memcpy(output->num_repeats+first_subset->num_blocks, second_subset->num_repeats, 
            second_subset->num_blocks);

      memcpy(output->start_offsets, first_subset->start_offsets, first_subset->num_blocks);
      memcpy(output->start_offsets+first_subset->num_blocks, second_subset->start_offsets, 
            second_subset->num_blocks);
      
      memcpy(output->end_offsets, first_subset->end_offsets, first_subset->num_blocks);
      memcpy(output->end_offsets+first_subset->num_blocks, second_subset->end_offsets, 
            second_subset->num_blocks);
   
      //Now we have all of the regions, so we just need to simplify them.
      __fenix_data_subset_simplify_regions(output); 
   } else {
      output->specifier = __FENIX_SUBSET_CREATEV;
      
      output->num_blocks = first_subset->num_blocks + second_subset->num_blocks;
      if(first_subset->specifier == __FENIX_SUBSET_CREATE){
         for(int i = 0; i < first_subset->num_blocks; i++){
            output->num_blocks += first_subset->num_repeats[i];
         }
      }
      if(second_subset->specifier == __FENIX_SUBSET_CREATE){
         for(int i = 0; i < second_subset->num_blocks; i++){
            output->num_blocks += second_subset->num_repeats[i];
         }
      }

      //TODO: Interrupted here by forced PC update.

     
   }

}


int __fenix_data_subset_free( Fenix_Data_subset *subset_specifier ) {
  int  retval = FENIX_SUCCESS;
  free( subset_specifier->num_repeats );
  free( subset_specifier->start_offsets );
  free( subset_specifier->end_offsets );
  subset_specifier->specifier = __FENIX_SUBSET_UNDEFINED;
  return retval;
}

/**
 * @brief
 * @param subset_specifier
 */
int __fenix_data_subset_delete( Fenix_Data_subset *subset_specifier ) {
  __fenix_data_subset_free(subset_specifier);
  free(subset_specifier);
  return FENIX_SUCCESS;
}
