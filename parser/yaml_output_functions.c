/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL Flexible Modeling System (FMS).
 *
 * FMS is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FMS is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
 **********************************************************************/
#ifdef use_yaml

#include <stdio.h>
#include <yaml.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

 // Should always match values of string_len_parameter and lvl2_key_parameter in fms_yaml_output_mod
#define LVL2KEY_NUM 8
#define KEY_STR_LEN 255
#define LVL2KEY_SIZE LVL2KEY_NUM*KEY_STR_LEN

#define DEBUG 0

struct fmsyamloutkeys {
	char key1 [KEY_STR_LEN];
  char key2 [KEY_STR_LEN];
  char key3 [KEY_STR_LEN];
  char key4 [KEY_STR_LEN];
  char key5 [KEY_STR_LEN];
  char key6 [KEY_STR_LEN];
  char key7 [KEY_STR_LEN];
  char key8 [KEY_STR_LEN];
  char key9 [KEY_STR_LEN];
  char key10 [KEY_STR_LEN];
  char key11 [KEY_STR_LEN];
  char key12 [KEY_STR_LEN];
  char key13 [KEY_STR_LEN];
  char key14 [KEY_STR_LEN];
  char key15 [KEY_STR_LEN];
  char key16 [KEY_STR_LEN];
  int level2key_offset;
	char level2key [LVL2KEY_SIZE];
};
struct fmsyamloutvalues {
  char val1 [KEY_STR_LEN];
  char val2 [KEY_STR_LEN];
  char val3 [KEY_STR_LEN];
  char val4 [KEY_STR_LEN];
  char val5 [KEY_STR_LEN];
  char val6 [KEY_STR_LEN];
  char val7 [KEY_STR_LEN];
  char val8 [KEY_STR_LEN];
  char val9 [KEY_STR_LEN];
  char val10 [KEY_STR_LEN];
  char val11 [KEY_STR_LEN];
  char val12 [KEY_STR_LEN];
  char val13 [KEY_STR_LEN];
  char val14 [KEY_STR_LEN];
  char val15 [KEY_STR_LEN];
  char val16 [KEY_STR_LEN];
};

/* \breif Prints a warning message that the yaml was not written correctly
 * \param emitter The libyaml emitter for this file
 * \param event The libyaml eent pointer
 * \param yamlname The name of the yaml file
 * \param keys yamlout The file pointer for the yaml file
 */
void error(char * yamlname ,yaml_event_t * event, yaml_emitter_t * emitter, FILE * yamlout){
  /* Write a warning to stderr and srdout */
  fprintf(stderr, "WARNING: YAML_OUTPUT: No output %s written.  Failed to emit event %d: %s\n", yamlname, event->type, emitter->problem);
  fprintf(stdout, "WARNING: YAML_OUTPUT: No output %s written.  Failed to emit event %d: %s\n", yamlname, event->type, emitter->problem);
  yaml_emitter_delete(emitter);
  fclose(yamlout);
}
void keyerror(yaml_event_t * event, yaml_emitter_t * emitter){
  /* Write a warning to stderr and srdout */
  fprintf(stderr, "WARNING: YAML_OUTPUT: Failed to emit event %d: %s\n", event->type, emitter->problem);
  fprintf(stdout, "WARNING: YAML_OUTPUT: Failed to emit event %d: %s\n", event->type, emitter->problem);
}
/* \breif Writes the key/value pairs of the fmsyamloutkeys and fmsyamloutvalues structs
 * \param emitter The libyaml emitter for this file
 * \param event The libyaml eent pointer
 * \param aindex The index of keys and vals that are being written currently
 * \param keys The keys to be written
 * \param vals The values correcponding to keys
 */
void write_keys_vals_yaml (yaml_emitter_t * emitter, yaml_event_t * event , int aindex, struct fmsyamloutkeys *keys, struct fmsyamloutvalues *vals){
  /* Check if a key exists */
  if (keys[aindex].key1[0] !='\0')  {
    /* Write the key */
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
        (yaml_char_t *)keys[aindex].key1, strlen(keys[aindex].key1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    /* Check for errors */
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    /* Write the value */
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
        (yaml_char_t *)vals[aindex].val1, strlen(vals[aindex].val1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    /* check for errors */
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  /* Repeat for all possible keys */
  if (keys[aindex].key2[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
        (yaml_char_t *)keys[aindex].key2, strlen(keys[aindex].key2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
        (yaml_char_t *)vals[aindex].val2, strlen(vals[aindex].val2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key3[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key3, strlen(keys[aindex].key3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val3, strlen(vals[aindex].val3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key4[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key4, strlen(keys[aindex].key4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val4, strlen(vals[aindex].val4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key5[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key5, strlen(keys[aindex].key5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val5, strlen(vals[aindex].val5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key6[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key6, strlen(keys[aindex].key6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val6, strlen(vals[aindex].val6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key7[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key7, strlen(keys[aindex].key7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val7, strlen(vals[aindex].val7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key8[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key8, strlen(keys[aindex].key8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val8, strlen(vals[aindex].val8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key9[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key9, strlen(keys[aindex].key9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val9, strlen(vals[aindex].val9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key10[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key10, strlen(keys[aindex].key10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val10, strlen(vals[aindex].val10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key11[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key11, strlen(keys[aindex].key11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val11, strlen(vals[aindex].val11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key12[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key12, strlen(keys[aindex].key12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val12, strlen(vals[aindex].val12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key13[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key13, strlen(keys[aindex].key13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val13, strlen(vals[aindex].val13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key14[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key14, strlen(keys[aindex].key14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val14, strlen(vals[aindex].val14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key15[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key15, strlen(keys[aindex].key15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val15, strlen(vals[aindex].val15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }
  if (keys[aindex].key16[0] !='\0')  {
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key16, strlen(keys[aindex].key16), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
    yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val16, strlen(vals[aindex].val16), 1, 0, YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(emitter, event)){
      keyerror(event, emitter);
      return;
    }
  }

return ;
}
/* \description Writes a YAML.
 * In this function, the YAML can have one key per "level" that has an array of keys.
 * If a struct has a level2key, then that key is written as a sequence. If level2key is empty, no second level is written.
 * If a key is null, then the kay/value pair is not written.
 * Users must pass in the size of each struct array for level 2 or 3.
 * L3 arrays are essentially 2D arrays in 1D and the number of elements can change between each instance, so the user must pass
 * the array `neach` which has the number of elements per instance.
 * The top level size must be 1 (for now).
 * \param yamlname The name of the output yaml file
 * \param asize The size of the top level array (must be 1)
 * \param topkeys The Keys for the top level of the yaml
 * \param topvals The values corresponding to the topkeys
 * \param a2size The size of the second level arrays and n3each
 * \param l2keys The keys for the second level of the yaml
 * \param l2vals The values corresponding to l2keys
 * \param a3size The full size of the l3 arrays
 * \param n3each Array that has the number of elements for each l2 array's third level elements
 * \param l3keys The keys for the third level of the yaml, any lvl2 keys will not be printed
 * \param l3vals The values corresponding to l3keys
 * \param lvl2keyeach array to indicate how many structs to print per level2key for the top level keys, should be the size of LVL2KEY_NUM
 */
void write_yaml_from_struct_3 (char *yamlname, int asize, struct fmsyamloutkeys *topkeys, struct fmsyamloutvalues *topvals, int a2size, struct fmsyamloutkeys *l2keys,
                               struct fmsyamloutvalues *l2vals, int a3size, int * n3each, struct fmsyamloutkeys *l3keys, struct fmsyamloutvalues *l3vals,
                               int* lvl2keyeach){
  yaml_emitter_t emitter; /* libyaml emitter */
  yaml_event_t event; /* libyaml event for the output yaml */
  int s2count = 0; /* A counter to keep track of the number of level 2 arrays output */
  int s3count = 0; /* A counter to keep track of the number of level 3 arrays output */
  FILE * yamlout; /* The file for the YAML output. */
  int i_n3 = 0;/* index for the n3each argument*/

  // trim any trailing whitespace
  int ws_ind = strlen(yamlname)-1;
  while(isspace(*(yamlname+ws_ind))) ws_ind--;
  if( ws_ind != strlen(yamlname)-1) yamlname[ws_ind+1] = '\0';

  /* open the yaml output file. Only 1 core should do this */
  yamlout = fopen(yamlname,"w");
  /* Start the emmitter */
  yaml_emitter_initialize(&emitter);
  yaml_emitter_set_output_file(&emitter, yamlout);
  yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
  if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
  }

  yaml_document_start_event_initialize(&event, NULL, NULL, NULL, 0);
  if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
  }
  /* start the event (top level) */
  yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
                                        1, YAML_ANY_MAPPING_STYLE);
  if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
  }

  /* write the top level */
  write_keys_vals_yaml (&emitter, &event , 0, topkeys, topvals);
  char* curr_topkey = topkeys->level2key;

  /* loop through the top level 2 keys */
  int top_ind;
  for (top_ind=0; top_ind < topkeys->level2key_offset; top_ind++) {
    /* Start the secodn level event */
    yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      				 (yaml_char_t *)curr_topkey, strlen(curr_topkey), 1, 0,
      				 YAML_PLAIN_SCALAR_STYLE);
    if (!yaml_emitter_emit(&emitter, &event)){
 	    error(yamlname, &event, &emitter, yamlout);
	    return;
    }
    /* Start the sequencing */
    yaml_sequence_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_SEQ_TAG,
      					 1, YAML_ANY_SEQUENCE_STYLE);
    if (!yaml_emitter_emit(&emitter, &event)){
  	  error(yamlname, &event, &emitter, yamlout);
	    return;
    }
    /* loop through the structs for this key*/
    int s2;
    for (s2 = 0 ; s2 < lvl2keyeach[top_ind]; s2++){
      yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
        				  1, YAML_ANY_MAPPING_STYLE);
      if (!yaml_emitter_emit(&emitter, &event)){
 	      error(yamlname, &event, &emitter, yamlout);
	      return;
      }
      /* call the write function */
      write_keys_vals_yaml (&emitter, &event , s2count, l2keys, l2vals);

      /* Next level keys */
      char * curr_l2key = (&l2keys[s2count])->level2key;
      int l2_ind;
      for (l2_ind = 0; l2_ind < (&l2keys[s2count])->level2key_offset; l2_ind++) {
        /* Start the third level event */
     	  yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
     		                            (yaml_char_t *)curr_l2key, strlen(curr_l2key), 1, 0,
     		                             YAML_PLAIN_SCALAR_STYLE);
     	  if (!yaml_emitter_emit(&emitter, &event)){
 		      error(yamlname, &event, &emitter, yamlout);
		      return;
 	      }
     	  /* Start the sequencing */
     	  yaml_sequence_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_SEQ_TAG,
         		1, YAML_ANY_SEQUENCE_STYLE);
        if (!yaml_emitter_emit(&emitter, &event)){
 		      error(yamlname, &event, &emitter, yamlout);
		      return;
        }
    	  /* loop through the structs */
        int s3start = s3count;
        int s3end = s3start + n3each[i_n3];
        i_n3++;
        int s3;
        for (s3 = s3start ; s3 < s3end ; s3++){
          yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
                    1, YAML_ANY_MAPPING_STYLE);
          if (!yaml_emitter_emit(&emitter, &event)){
            error(yamlname, &event, &emitter, yamlout);
            return;
          }
          /* call the write function */
          write_keys_vals_yaml (&emitter, &event , s3, l3keys, l3vals);
          yaml_mapping_end_event_initialize(&event);
          if(!yaml_emitter_emit(&emitter, &event)){
            error(yamlname, &event, &emitter, yamlout);
            return;
          }
          s3count ++;
        } // for s3
        yaml_sequence_end_event_initialize(&event);
        if (!yaml_emitter_emit(&emitter, &event)){
 		      error(yamlname, &event, &emitter, yamlout);
		      return;
        }
        curr_l2key = (&l2keys[s2count])->level2key + ((l2_ind+1) * KEY_STR_LEN);
      } /* for l2_ind */
      yaml_mapping_end_event_initialize(&event);
      if (!yaml_emitter_emit(&emitter, &event)){
 	      error(yamlname, &event, &emitter, yamlout);
	      return;
      }

      s2count++;
    }/* for s2 loop */

    yaml_sequence_end_event_initialize(&event);
    if (!yaml_emitter_emit(&emitter, &event)){
 	    error(yamlname, &event, &emitter, yamlout);
      return;
    }
    curr_topkey = topkeys->level2key + ((top_ind+1) * KEY_STR_LEN);
  }/* for top_ind */

  /* end the emitter */
  yaml_mapping_end_event_initialize(&event);
  if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
  }

  yaml_document_end_event_initialize(&event, 0);
  if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
  }
  yaml_stream_end_event_initialize(&event);
  if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
  }
  yaml_emitter_delete(&emitter);
  fclose(yamlout);
  return;
}
/// @brief Adds a key to the level2key character array in the struct.
/// Uses an offset to store multiple keys, size(amount of strings) is set by LVL2KEY_NUM (avoids c to fortran 2d array issues)
/// @param key_name name of level 2 key to add
/// @param key_length string length of key_name
/// @param keys key struct to add level 2 key to
void add_level2key(int key_length, char* key_name, struct fmsyamloutkeys* keys){

  // local fixed length copy to try to mitigate any fortran to c string weirdness
  char kname_loc[key_length + 1];
  //memset(kname_loc, '\0', sizeof(kname_loc));
  strncpy(kname_loc, key_name, key_length);
  kname_loc[key_length] = '\0';
  // error checking
  if ( strlen(kname_loc) > KEY_STR_LEN){
    fprintf(stderr, "WARNING: YAML_OUTPUT: invalid level two key passed to add_level2key. Max string size is %d, passed in string: %s",  KEY_STR_LEN, key_name);
    fprintf(stdout, "WARNING: YAML_OUTPUT: invalid level two key passed to add_level2key. Max string size is %d, passed in string: %s",  KEY_STR_LEN, key_name);
  }
  if( keys->level2key_offset  >= LVL2KEY_NUM ){
    fprintf(stderr, "WARNING: YAML_OUTPUT: max amount of level 2 keys (%d) has been exceeded", LVL2KEY_NUM);
    fprintf(stdout, "WARNING: YAML_OUTPUT: max amount of level 2 keys (%d) has been exceeded", LVL2KEY_NUM);
  }
  // check if string is set to initialize offset count
  if ( keys->level2key[0] == '\0'){
    keys->level2key_offset = 0;
  }
  // calculate offset and copy into the level2key array
  int offset = keys->level2key_offset * KEY_STR_LEN;
  char* curr_key;
  curr_key = keys->level2key + offset;
  strcpy(curr_key, kname_loc);

  keys->level2key_offset++;

  if (DEBUG) {
    printf("key length: %d \n key_name:", key_length);
    printf(key_name);
    printf("\n");
    printf("offset: %d \n", offset);
    printf("kname_loc:");
    printf(kname_loc);
    printf("\n*l2key:");
    printf(keys->level2key);
    printf("\nl2key+offset:");
    printf(keys->level2key + offset);
    printf("\n");
  }
}

#endif

