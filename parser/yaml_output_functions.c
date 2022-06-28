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

void write_keys_vals_yaml (yaml_emitter_t * emitter, yaml_event_t * event , int aindex, struct fmsyamloutkeys *keys, struct fmsyamloutvalues *vals){
    if (keys[aindex].key1[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key1, strlen(keys[aindex].key1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val1, strlen(vals[aindex].val1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key2[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key2, strlen(keys[aindex].key2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val2, strlen(vals[aindex].val2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key3[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key3, strlen(keys[aindex].key3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val3, strlen(vals[aindex].val3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key4[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key4, strlen(keys[aindex].key4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val4, strlen(vals[aindex].val4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key5[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key5, strlen(keys[aindex].key5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val5, strlen(vals[aindex].val5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key6[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key6, strlen(keys[aindex].key6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val6, strlen(vals[aindex].val6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key7[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key7, strlen(keys[aindex].key7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val7, strlen(vals[aindex].val7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key8[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key8, strlen(keys[aindex].key8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val8, strlen(vals[aindex].val8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key9[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key9, strlen(keys[aindex].key9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val9, strlen(vals[aindex].val9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key10[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key10, strlen(keys[aindex].key10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val10, strlen(vals[aindex].val10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key11[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key11, strlen(keys[aindex].key11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val11, strlen(vals[aindex].val11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key12[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key12, strlen(keys[aindex].key12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val12, strlen(vals[aindex].val12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key13[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key13, strlen(keys[aindex].key13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val13, strlen(vals[aindex].val13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key14[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key14, strlen(keys[aindex].key14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val14, strlen(vals[aindex].val14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
    }
    if (keys[aindex].key15[0] !='\0')  {
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys[aindex].key15, strlen(keys[aindex].key15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      yaml_scalar_event_initialize(event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals[aindex].val15, strlen(vals[aindex].val15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(emitter, event)) goto error;
      }

return ; 
error:
    fprintf(stderr, "YAML_OUTPUT: No output YAML written.  Failed to emit event %d: %s\n", event->type, emitter->problem);
}
/* \brief Writes a YAML.  Writes out each key/value of the struct.
 * If a struct has a level2key, then that key is written as a sequence.   
 * If a key is null, then the kay/value pair is not written.
 * Users must pass in the size of each struct array for level 2. 3, or 4.  The top level size should be 1.
 */
void write_yaml_from_struct_3 (int asize, struct fmsyamloutkeys *topkeys, struct fmsyamloutvalues *topvals, int a2size, struct fmsyamloutkeys *l2keys, struct fmsyamloutvalues *l2vals, int a3size, int * n3each, struct fmsyamloutkeys *l3keys, struct fmsyamloutvalues *l3vals){
  yaml_emitter_t emitter;
  yaml_event_t event;
  char buffer[64];
  int s3count = 0;
  /* Start the emmitter */  
    yaml_emitter_initialize(&emitter);
    yaml_emitter_set_output_file(&emitter, stdout);

    yaml_stream_start_event_initialize(&event, YAML_UTF8_ENCODING);
if (!yaml_emitter_emit(&emitter, &event)) goto error;

    yaml_document_start_event_initialize(&event, NULL, NULL, NULL, 0);
if (!yaml_emitter_emit(&emitter, &event)) goto error;
/* start the event (top level) */
    yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
        1, YAML_ANY_MAPPING_STYLE);
    if (!yaml_emitter_emit(&emitter, &event)) goto error;

/* write the top level */
 write_keys_vals_yaml (&emitter, &event , 0, topkeys, topvals);

 /* Check for the next level key */
 if (topkeys->level2key[0] !='\0') {
  /* Start the secodn level event */
  yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
  (yaml_char_t *)topkeys->level2key, strlen(topkeys->level2key), 1, 0, 
  YAML_PLAIN_SCALAR_STYLE);
  if (!yaml_emitter_emit(&emitter, &event)) goto error;
  /* Start the sequencing */
  yaml_sequence_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_SEQ_TAG,
      1, YAML_ANY_SEQUENCE_STYLE);
  if (!yaml_emitter_emit(&emitter, &event)) goto error;
  /* loop through the structs */
  for (int s2 = 0 ; s2 < a2size ; s2++){
    	yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
        	1, YAML_ANY_MAPPING_STYLE);
    	if (!yaml_emitter_emit(&emitter, &event)) goto error;
    	/* call the write function */
    	write_keys_vals_yaml (&emitter, &event , s2, l2keys, l2vals);
    	/* Next level keys */
    	if (l2keys->level2key[0] !='\0') {
     		/* Start the third level event */
     		yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
     		(yaml_char_t *)l2keys->level2key, strlen(l2keys->level2key), 1, 0,
     		YAML_PLAIN_SCALAR_STYLE);
     		if (!yaml_emitter_emit(&emitter, &event)) goto error;
     		/* Start the sequencing */
     		yaml_sequence_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_SEQ_TAG,
         		1, YAML_ANY_SEQUENCE_STYLE);
     		if (!yaml_emitter_emit(&emitter, &event)) goto error;
    		/* loop through the structs */
		int s3start = s3count;
		int s3end = s3start + n3each[s2];
		for (int s3 = s3start ; s3 < s3end ; s3++){
      			yaml_mapping_start_event_initialize(&event, NULL, (yaml_char_t *)YAML_MAP_TAG,
          			1, YAML_ANY_MAPPING_STYLE);
      			if (!yaml_emitter_emit(&emitter, &event)) goto error;
      			/* call the write function */
      			write_keys_vals_yaml (&emitter, &event , s3, l3keys, l3vals);
      			yaml_mapping_end_event_initialize(&event);
      			if (!yaml_emitter_emit(&emitter, &event)) goto error;	
			s3count ++;
     		}
     		yaml_sequence_end_event_initialize(&event);
    	 	if (!yaml_emitter_emit(&emitter, &event)) goto error;
    	}
    	
	yaml_mapping_end_event_initialize(&event);
    	if (!yaml_emitter_emit(&emitter, &event)) goto error;
  }
  yaml_sequence_end_event_initialize(&event);
  if (!yaml_emitter_emit(&emitter, &event)) goto error;

 }

/* end the emitter */    
    yaml_mapping_end_event_initialize(&event);
if (!yaml_emitter_emit(&emitter, &event)) goto error;

    yaml_document_end_event_initialize(&event, 0);
if (!yaml_emitter_emit(&emitter, &event)) goto error;

    yaml_stream_end_event_initialize(&event);
if (!yaml_emitter_emit(&emitter, &event)) goto error;

    yaml_emitter_delete(&emitter);

    return;
error:
    fprintf(stderr, "YAML_OUTPUT: No output YAML written.  Failed to emit event %d: %s\n", event.type, emitter.problem);
    yaml_emitter_delete(&emitter);
  
}	

/*int main () {
	struct fmsyamloutkeys k;
	struct fmsyamloutvalues v;
        int fs = 2; /* files */
/*        struct fmsyamloutkeys k2[fs];
        struct fmsyamloutvalues v2[fs];
        int n31 = 2; /*n variables in file 1 */
/*	int n32 = 3; /*n variables in file 2 */
/*	int n3ar [fs]; /* array with the number of variables for each file */
/*	struct fmsyamloutkeys k3[n31+n32]; /*variables */
/*        struct fmsyamloutvalues v3[n31+n32];

	int i;

	i = 1;
	strcpy(k.key1,"English");
	strcpy(v.val1,"Hello");
	strcpy(k.key2,"Spanish");
        strcpy(v.val2,"Hola");
	strcpy(k.key3,"French");
        strcpy(v.val3,"Bonjour");
	strcpy(k.key15,"");
        strcpy(k.level2key,"Outro");

	strcpy(k2[0].key1,"English");
	strcpy(v2[0].val1,"Goodbye");
	strcpy(k2[0].key2,"Spanish");
        strcpy(v2[0].val2,"Adios");
        strcpy(k2[0].key7,"French");
        strcpy(v2[0].val7,"Au revoir");
	strcpy(k2[0].level2key,"nextLev");
	strcpy(k2[1].key1,"Eng");
        strcpy(v2[1].val1,"Goo");
        strcpy(k2[1].key2,"ish");
        strcpy(v2[1].val2,"dios");
        strcpy(k2[1].key13,"Frch");
        strcpy(v2[1].val13,"Auevoir");
        strcpy(k2[1].level2key,"nextLev");

	strcpy(k3[0].key1,"z");
        strcpy(v3[0].val1,"Z");
        strcpy(k3[0].key2,"y");
        strcpy(v3[0].val2,"y");
        strcpy(k3[0].key3,"X");
        strcpy(v3[0].val3,"x");
        strcpy(k3[0].level2key,"");
        strcpy(k3[1].key1,"v");
        strcpy(v3[1].val1,"V");
        strcpy(k3[1].key2,"w");
        strcpy(v3[1].val2,"W");
        strcpy(k3[1].key10,"u");
        strcpy(v3[1].val10,"U");
        strcpy(k3[1].level2key,"");


	strcpy(k3[2].key1,"a");
        strcpy(v3[2].val1,"1");
        strcpy(k3[2].key2,"b");
        strcpy(v3[2].val2,"2");
        strcpy(k3[2].key3,"c");
        strcpy(v3[2].val3,"3");
        strcpy(k3[2].level2key,"");
        strcpy(k3[4].key1,"d");
        strcpy(v3[4].val1,"4");
        strcpy(k3[4].key2,"e");
        strcpy(v3[4].val2,"5");
        strcpy(k3[4].key3,"f");
        strcpy(v3[4].val3,"f");
        strcpy(k3[4].level2key,"");
        strcpy(k3[3].key1,"g");
        strcpy(v3[3].val1,"G");
        strcpy(k3[3].key2,"h");
        strcpy(v3[3].val2,"H");
        strcpy(k3[3].key3,"I");
        strcpy(v3[3].val3,"i");
        strcpy(k3[3].level2key,"");
	n3ar[0] = n31 ;
	n3ar[1] = n32 ;

        write_yaml (i, &k, &v, fs, k2, v2, (n31+n32), n3ar, k3, v3);
 return 0;
}*/
/* ifdef use_yaml */
