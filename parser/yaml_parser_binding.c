#include <stdio.h>
#include <yaml.h>
#include <stdbool.h>

struct key_value_pairs {
   int key_number;
   char key[50];
   char value[50];
   char parent_name[50];
   int parent_key;
};

struct yaml_file {
   int nkeys;
   struct key_value_pairs keys[50];
};

struct yaml_file my_file;

int get_num_blocks(char *block_name)
{
  int nblocks = 0;
  int i;

  for ( i = 1; i <= my_file.nkeys; i++ )
  {
     if(strcmp(my_file.keys[i].parent_name, block_name) == 0) nblocks = nblocks + 1;
  }

  return nblocks;
}

bool get_block_keys(char *block_name, int nfiles, int p[*])
{
  int i;

  for ( i = 1; i <= my_file.nkeys; i++ )
  {
     if(strcmp(my_file.keys[i].parent_name, block_name) == 0) {
        p[i] = my_file.keys[i].key_number;
     }
  }
  return true;
}

int main(void)
{
  yaml_parser_t parser;
  yaml_token_t  token;
  FILE *file;

  bool is_key = false;
  char key_value[50];
  int layer = 0;
  int key_count=0;
  int parent[10];
  int current_parent;
  char layer_name[10][50];
  char current_layername[50];
  int i;
  printf("opening file: %s\n", "diag_table.yaml");
  file = fopen("diag_table.yaml", "r");

  if(!yaml_parser_initialize(&parser))
    fputs("Failed to initialize parser!\n", stderr);

  bool is_new=false;
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

        is_new=true;
        if (strcmp(key_value, "")) strcpy(layer_name[layer], key_value);
        /* printf("LAYER:%i NAME:%s for %s=%i\n", layer, layer_name[layer], layer_name[layer-1], parent[layer-1]); */
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
         if (is_new) {
            parent[layer]=key_count;
            current_parent = parent[layer-1];
            strcpy(current_layername, layer_name[layer]);
            is_new = false;
         }
         i = key_count;
         my_file.keys[i].key_number=i;
         my_file.keys[i].parent_key = current_parent;
         strcpy(my_file.keys[i].parent_name, current_layername);
         strcpy(my_file.keys[i].key, key_value);
         strcpy(my_file.keys[i].value, token.data.scalar.value);
         my_file.nkeys = key_count;
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

  for ( i = 1; i <= my_file.nkeys; i++ ) {
       printf("Key_number:%i Parent_key:%i Parent_name:%s Key:%s Value:%s \n", my_file.keys[i].key_number, my_file.keys[i].parent_key, my_file.keys[i].parent_name, my_file.keys[i].key, my_file.keys[i].value);

  }
  printf("closing file: %s\n", "diag_table.yaml");
  fclose(file);

  int nfiles= get_num_blocks("diag_files");
  printf("diag_files = %i\n", nfiles);

  int block_keys[nfiles];
  if(get_block_keys("diag_files", nfiles, block_keys)){
     for ( i = 0; i <= nfiles-1; i++ ) {
        printf("nfile:%i block_key=%i \n", i, block_keys[i]);
     }
  }

  return true;
}
