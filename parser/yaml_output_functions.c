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

struct fmsyamloutkeys {
	char key1 [255];
        char key2 [255];
        char key3 [255];
        char key4 [255];
        char key5 [255];
        char key6 [255];
        char key7 [255];
        char key8 [255];
        char key9 [255];
        char key10 [255];
        char key11 [255];
        char key12 [255];
        char key13 [255];
        char key14 [255];
        char key15 [255];
	char level2key [255];
} fmsyamloutkeys;
struct fmsyamloutvalues {
        char val1 [255];
        char val2 [255];
        char val3 [255];
        char val4 [255];
        char val5 [255];
        char val6 [255];
        char val7 [255];
        char val8 [255];
        char val9 [255];
        char val10 [255];
        char val11 [255];
        char val12 [255];
        char val13 [255];
        char val14 [255];
        char val15 [255];
} fmsyamloutvalues;

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
 * \param l3keys The keys for the third level of the yaml
 * \param l3vals The values corresponding to l3keys
 */
void write_yaml_from_struct_3 (char *yamlname, int asize, struct fmsyamloutkeys *topkeys, struct fmsyamloutvalues *topvals, int a2size, struct fmsyamloutkeys *l2keys, struct fmsyamloutvalues *l2vals, int a3size, int * n3each, struct fmsyamloutkeys *l3keys, struct fmsyamloutvalues *l3vals){
  yaml_emitter_t emitter; /* libyaml emitter */
  yaml_event_t event; /* libyaml event for the output yaml */
  int s3count = 0; /* A counter to keep track of the number of level 3 arrays output */
  FILE * yamlout; /* The file for the YAML output. */
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
    /* Check for the next level key */
  if (topkeys->level2key[0] !='\0') {
  /* Start the secodn level event */
    yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      				 (yaml_char_t *)topkeys->level2key, strlen(topkeys->level2key), 1, 0,
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
    for (int s2 = 0 ; s2 < a2size ; s2++){
      yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
        				  1, YAML_ANY_MAPPING_STYLE);
      if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
      }
      /* call the write function */
      write_keys_vals_yaml (&emitter, &event , s2, l2keys, l2vals);
      /* Next level keys */
      if (l2keys->level2key[0] !='\0') {
        /* Start the third level event */
     	yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
     		(yaml_char_t *)l2keys->level2key, strlen(l2keys->level2key), 1, 0,
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
	int s3end = s3start + n3each[s2];
	for (int s3 = s3start ; s3 < s3end ; s3++){
      	  yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
          			1, YAML_ANY_MAPPING_STYLE);
      	  if (!yaml_emitter_emit(&emitter, &event)){
 		error(yamlname, &event, &emitter, yamlout);
		return;
	  }
      	  /* call the write function */
      	  write_keys_vals_yaml (&emitter, &event , s3, l3keys, l3vals);
	  yaml_mapping_end_event_initialize(&event);
      	  if (!yaml_emitter_emit(&emitter, &event)){
 		error(yamlname, &event, &emitter, yamlout);
		return;
	  }
	  s3count ++;
     	}
        yaml_sequence_end_event_initialize(&event);
        if (!yaml_emitter_emit(&emitter, &event)){
 		error(yamlname, &event, &emitter, yamlout);
		return;
        }
      } /* if (l2keys->level2key[0] !='\0') */
      yaml_mapping_end_event_initialize(&event);
      if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
      }
    }/* for s2 loop */
    yaml_sequence_end_event_initialize(&event);
    if (!yaml_emitter_emit(&emitter, &event)){
 	error(yamlname, &event, &emitter, yamlout);
	return;
    }

  }/* if (topkeys->level2key[0] !='\0') */

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

#endif

