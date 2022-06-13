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

void write_keys_vals_yaml (yaml_emitter_t * emitter, yaml_event_t * event , int asize, struct fmsyamloutkeys *keys, struct fmsyamloutvalues *vals){
    if (keys->key1[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key1, strlen(keys->key1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val1, strlen(vals->val1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key2[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key2, strlen(keys->key2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val2, strlen(vals->val2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key3[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key3, strlen(keys->key3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val3, strlen(vals->val3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key4[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key4, strlen(keys->key4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val4, strlen(vals->val4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key5[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key5, strlen(keys->key5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val5, strlen(vals->val5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key6[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key6, strlen(keys->key6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val6, strlen(vals->val6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key7[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key7, strlen(keys->key7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val7, strlen(vals->val7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key8[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key8, strlen(keys->key8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val8, strlen(vals->val8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key9[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key9, strlen(keys->key9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val9, strlen(vals->val9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key10[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key10, strlen(keys->key10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val10, strlen(vals->val10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key11[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key11, strlen(keys->key11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val11, strlen(vals->val11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key12[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key12, strlen(keys->key12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val12, strlen(vals->val12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key13[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key13, strlen(keys->key13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val13, strlen(vals->val13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key14[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key14, strlen(keys->key14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val14, strlen(vals->val14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (keys->key15[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)keys->key15, strlen(keys->key15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)vals->val15, strlen(vals->val15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      }
error:
    fprintf(stderr, "YAML_OUTPUT: No output YAML written.  Failed to emit event %d: %s\n", event.type, emitter.problem);
}
/* \brief Writes a 1 level YAML using the toplevel struct.  Writes out each key/value of the struct.
 * Ignores the level2key if one is given.  If a key is null, then the kay/value pair is not written.
 */
void write_yaml (int asize, struct fmsyamloutkeys *topkeys, struct fmsyamloutvalues *topvals, int a2size, struct fmsyamloutkeys *l2keys, struct fmsyamloutvalues *l2vals, int a3size, struct fmsyamloutkeys *l3keys, struct fmsyamloutvalues *l3vals,int a4size, struct fmsyamloutkeys *l4keys, struct fmsyamloutvalues *l4vals){
  yaml_emitter_t emitter;
  yaml_event_t event;
  char buffer[64];
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

/* start top level parsing */
 write_keys_vals_yaml (&emitter, &event , asize, topkeys, topvals);


/*
    if (topkeys->key1[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key1, strlen(topkeys->key1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val1, strlen(topvals->val1), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key2[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key2, strlen(topkeys->key2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val2, strlen(topvals->val2), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key3[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key3, strlen(topkeys->key3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val3, strlen(topvals->val3), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key4[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key4, strlen(topkeys->key4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val4, strlen(topvals->val4), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key5[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key5, strlen(topkeys->key5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val5, strlen(topvals->val5), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key6[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key6, strlen(topkeys->key6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val6, strlen(topvals->val6), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key7[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key7, strlen(topkeys->key7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val7, strlen(topvals->val7), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key8[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key8, strlen(topkeys->key8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val8, strlen(topvals->val8), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key9[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key9, strlen(topkeys->key9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val9, strlen(topvals->val9), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key10[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key10, strlen(topkeys->key10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val10, strlen(topvals->val10), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key11[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key11, strlen(topkeys->key11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val11, strlen(topvals->val11), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key12[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key12, strlen(topkeys->key12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val12, strlen(topvals->val12), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key13[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key13, strlen(topkeys->key13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val13, strlen(topvals->val13), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key14[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key14, strlen(topkeys->key14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val14, strlen(topvals->val14), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
    }
    if (topkeys->key15[0] !='\0')  {
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topkeys->key15, strlen(topkeys->key15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
       if (!yaml_emitter_emit(&emitter, &event)) goto error;
      yaml_scalar_event_initialize(&event, NULL, (yaml_char_t *)YAML_STR_TAG,
      (yaml_char_t *)topvals->val15, strlen(topvals->val15), 1, 0, YAML_PLAIN_SCALAR_STYLE);
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

int main () {
	struct fmsyamloutkeys k;
	struct fmsyamloutvalues v;
	int i;

	i = 1;
	strcpy(k.key1,"English");
	strcpy(v.val1,"Hello");
	strcpy(k.key2,"Spanish");
        strcpy(v.val2,"Hola");
	strcpy(k.key3,"French");
        strcpy(v.val3,"Bonjour");
	strcpy(k.key15,"");

        write_yaml (i, &k, &v, 0, NULL, NULL, 0, NULL, NULL, 0, NULL, NULL);
 return 0;
}
/* ifdef use_yaml */
