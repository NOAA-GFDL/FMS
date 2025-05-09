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
#include <stdbool.h>

/* Type to store info about key */
typedef struct {
   int key_number;        /* Id of this key */
   char key[255];         /* Name of the key */
   char value[255];       /* Value of the key */
   char parent_name[255]; /* Name of the block the key belongs to */
   int parent_key;        /* Id of the block the key belongs to */
}key_value_pairs;

/* Type to store all of the keys */
typedef struct  {
   int nkeys;
   key_value_pairs *keys;
}yaml_file;

/* Type to store all the yaml files that are opened */
typedef struct {
   yaml_file *files;
}file_type;

file_type my_files; /* Array of opened yaml files */
int nfiles = 0;     /* Number of files in the yaml file */

/* @brief  Private c function that gets the number of key-value pairs in a block
   @return Number of key-value pairs in this block */
int get_nkeys_binding(int *file_id, int *block_id)
{
  int nkeys = 0;    /* Number of key-value pairs */
  int i;            /* For loops */
  int j = *file_id; /* To minimize the typing :) */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(my_files.files[j].keys[i].parent_key == *block_id && !strcmp(my_files.files[j].keys[i].parent_name, "") ) nkeys = nkeys + 1;
  }

  return nkeys;

}

/* @brief Private c function that gets the ids of the key-value pairs in a block */
void get_key_ids_binding(int *file_id, int *block_id, int *key_ids)
{
  int i;              /* For loops */
  int key_count = -1; /* Number of key-value pairs */
  int j = *file_id;   /* To minimize the typing :) */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(my_files.files[j].keys[i].parent_key == *block_id && !strcmp(my_files.files[j].keys[i].parent_name, "") ){
        key_count = key_count + 1;
        key_ids[key_count] = i;
     }
  }

  return;
}

/* @brief Private c function that get the key from a key_id in a yaml file
   @return Name of the key obtained */
char *get_key(int *file_id, int *key_id)
{
  int j = *file_id;   /* To minimize the typing :) */
  return my_files.files[j].keys[*key_id].key;
}

/* @brief Private c function that get the value from a key_id in a yaml file
   @return String containing the value obtained */
char *get_value(int *file_id, int *key_id)
{
  int j = *file_id;   /* To minimize the typing :) */
  return my_files.files[j].keys[*key_id].value;
}

/* @brief Private c functions get gets the block name from a block id
   @return String containing the value obtained */
char *get_block(int *file_id, int *block_id)
{
  int j = *file_id;   /* To minimize the typing :) */
  return my_files.files[j].keys[*block_id].parent_name;
}

/* @brief Private c function that determines they value of a key in yaml_file
   @return c pointer with the value obtained */
char *get_value_from_key_wrap(int *file_id, int *block_id, char *key_name, int *sucess) /*, char *key_name) */
{
  int i;              /* For loops */
  int j = *file_id;   /* To minimize the typing :) */

  *sucess = 0;          /* Flag indicating if the search was sucessful */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if (my_files.files[j].keys[i].parent_key == *block_id)
     {
        if( strcmp(my_files.files[j].keys[i].key, key_name) == 0)
        {
           *sucess = 1;
           break;
        }
     }
  }
  if (*sucess == 1) {return my_files.files[j].keys[i].value;} else {return "";}
}

/* @brief Private c function that determines the number of blocks with block_name in the yaml file
   @return Number of blocks with block_name */
int get_num_blocks_all(int *file_id, char *block_name)
{
  int nblocks = 0;    /* Number of blocks */
  int i;              /* For loops */
  int j = *file_id;   /* To minimize the typing :) */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0) nblocks = nblocks + 1;
  }

  return nblocks;
}

/* @brief Private c function that determines the number of unique blocks (i.e diag_files, varlist, etc)
   @return The number of unique blocks */
int get_num_unique_blocks_bind(int *file_id, int *parent_block_id)
{
  int nblocks = 0;    /* Number of blocks */
  int i;              /* For loops */
  int j = *file_id;   /* To minimize the typing :) */
  char block_names[my_files.files[j].nkeys][255]; /* Array that stores the names of the unique blocks*/
  bool found;         /* True if the block name was already found (i.e it not unqiue)*/
  int k;              /* For loops */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
    if (my_files.files[j].keys[i].parent_key == *parent_block_id )
    {
      if (strcmp(my_files.files[j].keys[i].parent_name, "") == 0){
        continue;
      }
      found = false;
      for (k = 1; k <= nblocks; k++)
      {
        if (strcmp(block_names[k], my_files.files[j].keys[i].parent_name) == 0)
        {
          found = true;
          break;
        }
      }

      if (found) continue;

      nblocks = nblocks + 1;
      strcpy(block_names[nblocks], my_files.files[j].keys[i].parent_name);
      // printf("Block names: %s \n", block_names[nblocks]);
    }
  }
  return nblocks;
}

/* @brief Private c function that determines the ids of the unique blocks (i.e diag_files, varlist, etc)
   @return The ids of the unique blocks */
void get_unique_block_ids_bind(int *file_id, int *block_ids, int *parent_block_id)
{
  int nblocks = 0;    /* Number of blocks */
  int i;              /* For loops */
  int j = *file_id;   /* To minimize the typing :) */
  char block_names[my_files.files[j].nkeys][255]; /* Array that stores the names of the unique blocks*/
  bool found;         /* True if the block name was already found (i.e it not unqiue)*/
  int k;              /* For loops */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
    if (my_files.files[j].keys[i].parent_key == *parent_block_id )
    {
      if (strcmp(my_files.files[j].keys[i].parent_name, "") == 0){
        continue;
      }
      found = false;
      for (k = 1; k <= nblocks; k++)
      {
        if (strcmp(block_names[k], my_files.files[j].keys[i].parent_name) == 0)
        {
          found = true;
          break;
        }
      }

      if (found) continue;

      nblocks = nblocks + 1;
      block_ids[nblocks - 1] = my_files.files[j].keys[i].key_number;
      strcpy(block_names[nblocks], my_files.files[j].keys[i].parent_name);
      //printf("Block names: %s \n", block_names[nblocks]);
    }
  }
  return;
}
/* @brief Private c function that determines the number of blocks with block_name that belong to
   a parent block with parent_block_id in the yaml file
   @return Number of blocks with block_name */
int get_num_blocks_child(int *file_id, char *block_name, int *parent_block_id)
{
  int nblocks = 0;    /* Number of blocks */
  int i;              /* For loops */
  int j = *file_id;   /* To minimize the typing :) */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0 && my_files.files[j].keys[i].parent_key == *parent_block_id) nblocks = nblocks + 1;
  }

  return nblocks;
}


/* @brief Private c function that gets the the ids of the blocks with block_name in the yaml file */
void get_block_ids_all(int *file_id, char *block_name, int *block_ids)
{
  int i;              /* For loops */
  int nblocks = -1;   /* Number of blocks */
  int j = *file_id;   /* To minimize the typing :) */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0) {
        nblocks = nblocks + 1;
        block_ids[nblocks] = my_files.files[j].keys[i].key_number;
     }
  }
  return;
}

/* @brief Private c function that gets the the ids of the blocks with block_name and that
   belong to a parent block id in the yaml file */
void get_block_ids_child(int *file_id, char *block_name, int *block_ids, int *parent_key_id )
{
  int i;              /* For loops */
  int nblocks = -1;   /* Number of blocks */
  int j = *file_id;   /* To minimize the typing :) */

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0 && my_files.files[j].keys[i].parent_key == *parent_key_id) {
        nblocks = nblocks + 1;
        block_ids[nblocks] = my_files.files[j].keys[i].key_number;
     }
  }
  return;
}

/* @brief Private c function to determine if a block_id is valid */
bool is_valid_block_id(int *file_id, int *block_id)
{
   /* If the block id it not in the allowed range is not a valid block id */
   if (*block_id <= -1 || *block_id > my_files.files[*file_id].nkeys) {return false;}

   /* If the block id has an empty parent name then it is not a valid block id */
   if (*block_id != 0 && strcmp(my_files.files[*file_id].keys[*block_id].parent_name, "") == 0) {return false;}
   return true;
}

/* @brief Private c function to determine if a key_id is valid */
bool is_valid_key_id(int *file_id, int *key_id)
{
   if (*key_id > -1 && *key_id <= my_files.files[*file_id].nkeys) {return true;}
   else { return false;}
}

/* @brief Private c function to determine if a file_id is valid */
bool is_valid_file_id(int *file_id)
{
   if (*file_id > -1 && *file_id < nfiles) {return true;}
   else { return false;}
}

/* @brief Private c function that opens and parses a yaml file and saves it in a struct
   @return Flag indicating if the read was sucessful */
int open_and_parse_file_wrap(char *filename, int *file_id)
{
  yaml_parser_t parser;
  yaml_token_t  token;
  FILE *file;

  bool is_key = false;         /* Flag indicating if the current token in a key */
  char key_value[255];         /* Value of a key */
  int layer = 0;               /* Current layer (block level) */
  int key_count=0;             /* Current number of keys */
  int parent[10];              /* Ids of blocks */
  int current_parent;          /* Id of the current block */
  char layer_name[10][255];    /* Array of block names */
  char current_layername[255]; /* Name of the current block */
  int i;                       /* To minimize the typing :) */
  int j;                       /* To minimize the typing :) */

  if (nfiles == 0 )
  {
     my_files.files = (yaml_file*)calloc(1, sizeof(yaml_file));
  } else
  {
     my_files.files = realloc(my_files.files, (nfiles+1)*sizeof(yaml_file));
  }

  j = nfiles;
  *file_id =j;

/*  printf("Opening file: %s.\nThere are %i files opened.\n", filename, j); */
  file = fopen(filename, "r");
  if (file == NULL) return -1;

  if(!yaml_parser_initialize(&parser)) return -2;

  my_files.files[j].keys = (key_value_pairs*)calloc(1, sizeof(key_value_pairs));

  parent[0]=0;
  strcpy(layer_name[0], "TOP");
  /* Set input file */
  yaml_parser_set_input_file(&parser, file);
  do {
    if (!yaml_parser_scan(&parser, &token)) {
      return -3;
    }
    switch(token.type)
    {
    case YAML_KEY_TOKEN:
       {
        is_key = true;
        break;
       }
    case YAML_VALUE_TOKEN:
       {
        is_key = false;
        break;
       }
    case YAML_BLOCK_ENTRY_TOKEN:
       {
        layer = layer + 1;

        if (strcmp(key_value, ""))
        {
           strcpy(layer_name[layer], key_value);
        }
        key_count = key_count + 1;
        i = key_count;
        my_files.files[j].keys = realloc(my_files.files[j].keys, (i+1)*sizeof(key_value_pairs));
        my_files.files[j].keys[i].key_number=i;
        my_files.files[j].keys[i].parent_key = parent[layer-1];
        strcpy(my_files.files[j].keys[i].parent_name, layer_name[layer]);
        strcpy(my_files.files[j].keys[i].key, "");
        strcpy(my_files.files[j].keys[i].value, "");
        parent[layer]=key_count;
        /*printf("KEY:%i LAYER:%i NAME:%s for %s=%i\n", key_count, layer, layer_name[layer], layer_name[layer-1], parent[layer-1]); */

        break;
       }
    case YAML_BLOCK_END_TOKEN:
       {
        layer = layer - 1;
        break;
       }
    case YAML_SCALAR_TOKEN:
       {
        if ( ! is_key)
        {
         current_parent = parent[layer];
         strcpy(current_layername, "");
         key_count = key_count + 1;
         i = key_count;
         my_files.files[j].keys = realloc(my_files.files[j].keys, (i+1)*sizeof(key_value_pairs));
         my_files.files[j].keys[i].key_number=i;
         my_files.files[j].keys[i].parent_key = current_parent;
         strcpy(my_files.files[j].keys[i].parent_name, current_layername);
         strcpy(my_files.files[j].keys[i].key, key_value);
         strcpy(my_files.files[j].keys[i].value, token.data.scalar.value);
         my_files.files[j].nkeys = key_count;
         /* printf("----> LAYER:%i LAYER_NAME=%s PARENT:%i, KEYCOUNT:%i KEY: %s VALUE: %s \n", layer, current_layername, current_parent, key_count, key_value, token.data.scalar.value); */
         strcpy(key_value,"");
        }
        else
         {strcpy(key_value,token.data.scalar.value);}
        }
     break;
     }
    if(token.type != YAML_STREAM_END_TOKEN)
      yaml_token_delete(&token);
  } while(token.type != YAML_STREAM_END_TOKEN);
  yaml_token_delete(&token);
  yaml_parser_delete(&parser);

  /*
  for ( i = 1; i <= my_files.files[j].nkeys; i++ ) {
       printf("Key_number:%i Parent_key:%i Parent_name:%s Key:%s Value:%s \n", my_files.files[j].keys[i].key_number, my_files.files[j].keys[i].parent_key, my_files.files[j].keys[i].parent_name, my_files.files[j].keys[i].key, my_files.files[j].keys[i].value);
  }
  printf("/\n");
   */

  nfiles = nfiles + 1;
/*  printf("closing file: %s\n", filename); */
  fclose(file);

  return 1;
}

#endif
