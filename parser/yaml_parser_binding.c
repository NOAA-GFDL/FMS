#ifdef use_yaml

#include <stdio.h>
#include <yaml.h>
#include <stdbool.h>

typedef struct {
   int key_number;
   char key[255];
   char value[255];
   char parent_name[255];
   int parent_key;
}key_value_pairs;

typedef struct  {
   int nkeys;
   key_value_pairs *keys;
}yaml_file;

typedef struct {
   yaml_file *files;
}file_type;

file_type my_files;
int nfiles = 0;

int get_nkeys(int *file_id, int *block_id)
{
  int nkeys = 0;
  int i;
  int j = *file_id;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(my_files.files[j].keys[i].parent_key == *block_id && !strcmp(my_files.files[j].keys[i].parent_name, "") ) nkeys = nkeys + 1;
  }

  return nkeys;

}

void get_key_ids(int *file_id, int *block_id, int key_ids[*])
{
  int i;
  int nkeys = -1;
  int j = *file_id;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(my_files.files[j].keys[i].parent_key == *block_id && !strcmp(my_files.files[j].keys[i].parent_name, "") ){
        nkeys = nkeys + 1;
        key_ids[nkeys] = i;
     }
  }

  return;
}

char *get_key(int *file_id, int *key_id)
{
   char *key_name;
   int j = *file_id;

  key_name = malloc(sizeof(char) * (strlen(my_files.files[j].keys[*key_id].key) + 1));
  strcpy(key_name, my_files.files[j].keys[*key_id].key);

  return key_name;
}

char *get_value(int *file_id, int *key_id)
{
  char *key_value;
  int j = *file_id;

  key_value = malloc(sizeof(char) * (strlen(my_files.files[j].keys[*key_id].value) + 1));
  strcpy(key_value, my_files.files[j].keys[*key_id].value);

  return key_value;
}

char *get_value_from_key_wrap(int *file_id, int *block_id, char *key_name, bool *sucess) /*, char *key_name) */
{
  int i;
  int j = *file_id;
 
 char *key_value=NULL;
 *sucess = false;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if (my_files.files[j].keys[i].parent_key == *block_id)
     {
        if( strcmp(my_files.files[j].keys[i].key, key_name) == 0)
        {
           key_value = malloc(sizeof(char) * (strlen(my_files.files[j].keys[i].value) + 1));
           strcpy(key_value, my_files.files[j].keys[i].value);
           *sucess = true;
           break;
        }
     }
  }
  return key_value;
}

int get_num_blocks_all(int *file_id, char *block_name)
{
  int nblocks = 0;
  int i;
  int j = *file_id;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0) nblocks = nblocks + 1;
  }

  return nblocks;
}

int get_num_blocks_child(int *file_id, char *block_name, int *parent_key_id)
{
  int nblocks = 0;
  int i;
  int j = *file_id;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0 && my_files.files[j].keys[i].parent_key == *parent_key_id) nblocks = nblocks + 1;
  }

  return nblocks;
}

void get_block_ids_all(int *file_id, char *block_name, int block_ids[*])
{
  int i;
  int nblocks = -1;
  int j = *file_id;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0) {
        nblocks = nblocks + 1;
        block_ids[nblocks] = my_files.files[j].keys[i].key_number;
     }
  }
  return;
}

void get_block_ids_child(int *file_id, char *block_name, int block_ids[*], int *parent_key_id )
{
  int i;
  int nblocks = -1;
  int j = *file_id;

  for ( i = 1; i <= my_files.files[j].nkeys; i++ )
  {
     if(strcmp(my_files.files[j].keys[i].parent_name, block_name) == 0 && my_files.files[j].keys[i].parent_key == *parent_key_id) {
        nblocks = nblocks + 1;
        block_ids[nblocks] = my_files.files[j].keys[i].key_number;
     }
  }
  return;
}

bool open_and_parse_file_wrap(char *filename, int *file_id)
{
  yaml_parser_t parser;
  yaml_token_t  token;
  FILE *file;

  bool is_key = false;
  char key_value[255];
  int layer = 0;
  int key_count=0;
  int parent[10];
  int current_parent;
  char layer_name[10][255];
  char current_layername[255];
  int i;
  int j;

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
  if (file == NULL) return false;

  if(!yaml_parser_initialize(&parser)) return false;

  my_files.files[j].keys = (key_value_pairs*)calloc(1, sizeof(key_value_pairs));

  parent[0]=0;
  strcpy(layer_name[0], "TOP");
  /* Set input file */
  yaml_parser_set_input_file(&parser, file);
  do {
    yaml_parser_scan(&parser, &token);
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

  for ( i = 1; i <= my_files.files[j].nkeys; i++ ) {
       printf("Key_number:%i Parent_key:%i Parent_name:%s Key:%s Value:%s \n", my_files.files[j].keys[i].key_number, my_files.files[j].keys[i].parent_key, my_files.files[j].keys[i].parent_name, my_files.files[j].keys[i].key, my_files.files[j].keys[i].value);
  }
  printf("/\n");

  nfiles = nfiles + 1;
/*  printf("closing file: %s\n", filename); */
  fclose(file);

  return true;
}

#endif
