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

// TODO: FMS_MAX_FILE_LEN is repeated from fms_platform.h.
// Consider consolidating constant definitions or using a shared header.
#ifndef FMS_MAX_FILE_LEN
#define FMS_MAX_FILE_LEN 255
#endif

#define FMS_FILE_LEN FMS_MAX_FILE_LEN

/**
 * @brief Represents a key-value pair within a YAML file.
 *
 * Includes metadata for tracking the key's structure and origin.
 */
typedef struct {
    int key_id;                     ///< ID of this key
    char key[FMS_FILE_LEN];         ///< Name of the key
    char value[FMS_FILE_LEN];       ///< Value associated with the key
    char parent_name[FMS_FILE_LEN]; ///< Name of the parent block
    int parent_key;                 ///< ID of the parent block
} KeyValuePairs;

/**
 * @brief Represents the contents of a YAML file as an array of key-value pairs.
 */
typedef struct {
    int nkeys;                ///< Number of keys defined in the YAML file
    KeyValuePairs *keys;      ///< Array of key-value pairs
} YamlFile;

/**
 * @brief Represents all YAML files that have been opened and parsed.
 */
typedef struct {
    YamlFile *files;          ///< Array of parsed YAML files
} FileType;

FileType my_files;    // Struct holding parsed YAML files
int nfiles = 0;       // Number of files opened and parsed so far

/**
 * @brief Private C function that gets the number of key-value pairs in a block.
 *
 * @param file_id Pointer to the file index.
 * @param block_id Pointer to the block identifier.
 * @return Number of key-value pairs in this block.
 */
int get_nkeys_binding(const int *file_id, const int *block_id)
{
  int nkeys = 0;
  int j = *file_id;

  for (int i = 1; i <= my_files.files[j].nkeys; i++)
  {
     if (my_files.files[j].keys[i].parent_key == *block_id &&
         my_files.files[j].keys[i].parent_name[0] == '\0')
         nkeys++;
  }

  return nkeys;
}

/**
 * @brief Private C function that gets the ids of the key-value pairs in a block.
 *
 * @param file_id Pointer to the file index.
 * @param block_id Pointer to the block identifier.
 * @param key_ids Output array to store the key ids.
*/
void get_key_ids_binding(const int *file_id, const int *block_id, int *key_ids)
{
    int j = *file_id;
    int key_count = 0;

    for (int i = 1; i <= my_files.files[j].nkeys; i++)
    {
        if (my_files.files[j].keys[i].parent_key == *block_id &&
            my_files.files[j].keys[i].parent_name[0] == '\0')
        {
            key_ids[key_count++] = i;
        }
    }
}

/**
 * @brief Gets the key name corresponding to a key ID in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param key_id Pointer to the index of the key within that file (read-only).
 * @return Pointer to the key name.
 *
 * @note Assumes file_id and key_id are valid and checked beforehand.
 */
char *get_key(const int *file_id, const int *key_id)
{
    return my_files.files[*file_id].keys[*key_id].key;
}

/**
 * @brief Gets the key value corresponding to a key ID in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param key_id Pointer to the index of the key within that file (read-only).
 * @return Pointer to the key value.
 *
 * @note Assumes file_id and key_id are valid and checked beforehand.
 */
char *get_value(const int *file_id, const int *key_id)
{
  return my_files.files[*file_id].keys[*key_id].value;
}

/**
 * @brief Gets the block name corresponding to a block id in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_id Pointer to the index of the block within that file (read-only).
 * @return Pointer to the block name.
 *
 * @note Assumes file_id and block_id are valid and checked beforehand.
 */
char *get_block(const int *file_id, const int *block_id)
{
  return my_files.files[*file_id].keys[*block_id].parent_name;
}

/**
 * @brief Searches for the value of a given key name within a specified block in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_id Pointer to the ID of the block to search within (read-only).
 * @param key_name The name of the key to search for (read-only).
 * @param success Pointer to an int flag that will be set to 1 if found, 0 otherwise.
 * @return Pointer to the value string if found, or NULL if not found.
 *
 * @note Assumes file_id and block_id are valid and checked beforehand.
 */
char *get_value_from_key_wrap(const int *file_id, const int *block_id, const char *key_name, int *success)
{
    int i;
    int j = *file_id;

    *success = 0;  // Initialize flag to failure

    for (i = 1; i <= my_files.files[j].nkeys; i++)
    {
        if (my_files.files[j].keys[i].parent_key == *block_id)
        {
            if (strcmp(my_files.files[j].keys[i].key, key_name) == 0)
            {
                *success = 1;
                return my_files.files[j].keys[i].value;
            }
        }
    }

    return "";  // Key not found
}

/**
 * @brief Counts the number of blocks with a given block name in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_name Name of the block to count (read-only).
 * @return Number of blocks matching the given block name.
 *
 * @note Assumes valid file_id and non-null block_name.
 *       The function skips the key at index 0 (if indexing starts at 1).
 */
int get_num_blocks_all(const int *file_id, const char *block_name)
{
    int nblocks = 0;
    int i;
    int fid = *file_id;

    /* Loop through keys, assuming keys are 1-based indexed */
    for (i = 1; i <= my_files.files[fid].nkeys; i++)
    {
        if (strcmp(my_files.files[fid].keys[i].parent_name, block_name) == 0)
        {
            nblocks++;
        }
    }

    return nblocks;
}

/**
 * @brief Counts the number of unique blocks with a given parent block ID.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param parent_block_id Pointer to the ID of the parent block (read-only).
 * @return Number of unique blocks found.
 */
int get_num_unique_blocks_bind(const int *file_id, const int *parent_block_id)
{
    int nblocks = 0;
    int i, k;
    int fid = *file_id;
    char block_names[my_files.files[fid].nkeys][FMS_FILE_LEN];  // Assuming 255 or defined length
    bool found;

    for (i = 1; i <= my_files.files[fid].nkeys; i++)
    {
        if (my_files.files[fid].keys[i].parent_key == *parent_block_id)
        {
            if (strcmp(my_files.files[fid].keys[i].parent_name, "") == 0)
                continue;

            found = false;
            for (k = 0; k < nblocks; k++)
            {
                if (strcmp(block_names[k], my_files.files[fid].keys[i].parent_name) == 0)
                {
                    found = true;
                    break;
                }
            }

            if (found)
                continue;

            strcpy(block_names[nblocks], my_files.files[fid].keys[i].parent_name);
            ++nblocks;
        }
    }

    return nblocks;
}

/**
 * @brief Gets the IDs of the unique blocks with a given parent block ID.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_ids Array to store the unique block IDs (output).
 * @param parent_block_id Pointer to the ID of the parent block (read-only).
 */
void get_unique_block_ids_bind(const int *file_id, int *block_ids, const int *parent_block_id)
{
    int nblocks = 0;
    int i, k;
    int fid = *file_id;
    char block_names[my_files.files[fid].nkeys][FMS_FILE_LEN];
    bool found;

    for (i = 1; i <= my_files.files[fid].nkeys; i++)
    {
        if (my_files.files[fid].keys[i].parent_key == *parent_block_id)
        {
            if (strcmp(my_files.files[fid].keys[i].parent_name, "") == 0)
                continue;

            found = false;
            for (k = 0; k < nblocks; k++)
            {
                if (strcmp(block_names[k], my_files.files[fid].keys[i].parent_name) == 0)
                {
                    found = true;
                    break;
                }
            }

            if (found)
                continue;

            block_ids[nblocks] = my_files.files[fid].keys[i].key_id;
            strcpy(block_names[nblocks], my_files.files[fid].keys[i].parent_name);
            ++nblocks;
        }
    }
}

/**
 * @brief Counts the number of blocks with a given name that belong to a specified parent block.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_name Name of the blocks to count.
 * @param parent_block_id Pointer to the ID of the parent block (read-only).
 * @return Number of blocks with the specified block_name and parent_block_id.
 */
int get_num_blocks_child(const int *file_id, const char *block_name, const int *parent_block_id)
{
    int nblocks = 0;          // Counter for matching blocks
    int i;
    int fid = *file_id;       // Local copy for clarity

    for (i = 1; i <= my_files.files[fid].nkeys; i++)
    {
        if (strcmp(my_files.files[fid].keys[i].parent_name, block_name) == 0 &&
            my_files.files[fid].keys[i].parent_key == *parent_block_id)
        {
            nblocks++;
        }
    }

    return nblocks;
}

/**
 * @brief Retrieves the IDs of all blocks with the specified block_name in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_name Name of the block to search for.
 * @param block_ids Output array to store the found block IDs.
 *
 * @note Assumes block_ids points to enough memory to hold all matching block IDs.
 */
void get_block_ids_all(const int *file_id, const char *block_name, int *block_ids)
{
    int i;
    int nblocks = 0;         // Number of matching blocks found
    int fid = *file_id;      // Local copy of file ID for convenience

    // Loop over keys; assuming keys index starts at 1
    for (i = 1; i <= my_files.files[fid].nkeys; i++)
    {
        if (strcmp(my_files.files[fid].keys[i].parent_name, block_name) == 0)
        {
            block_ids[nblocks] = my_files.files[fid].keys[i].key_id;
            nblocks++;
        }
    }
}

/**
 * @brief Finds the IDs of child blocks with a given name and parent block ID in a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_name Name of the block to search for (read-only).
 * @param block_ids Output array to store matching block IDs.
 * @param parent_key_id Pointer to the parent block ID (read-only).
 *
 * @note Assumes block_ids points to an array large enough to hold all matching IDs.
 *       Assumes valid file_id and non-null pointers.
 */
void get_block_ids_child(const int *file_id, const char *block_name, int *block_ids, const int *parent_key_id)
{
  int nblocks = 0;  // Tracks number of matches found
  int fid = *file_id;

  // no need to check i == 0
  for ( int i = 1; i <= my_files.files[fid].nkeys; i++ )
  {
     if(strcmp(my_files.files[fid].keys[i].parent_name, block_name) == 0 &&
        my_files.files[fid].keys[i].parent_key == *parent_key_id)
      {
        block_ids[nblocks++] = my_files.files[fid].keys[i].key_id;
      }
  }
}

/**
 * @brief Checks whether a given block ID is valid within a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param block_id Pointer to the block ID to validate (read-only).
 * @return true if the block ID is valid, false otherwise.
 *
 * @note A block ID is invalid if:
 *       - It is less than zero or greater than the number of keys in the file.
 *       - Its associated parent name is empty (except for block ID 0).
 */
bool is_valid_block_id(const int *file_id, const int *block_id)
{
    int fid = *file_id;
    int bid = *block_id;

    if (bid <= -1 || bid > my_files.files[fid].nkeys) {
        return false;
    }

    if (bid != 0 && strcmp(my_files.files[fid].keys[bid].parent_name, "") == 0) {
        return false;
    }

    return true;
}

/**
 * @brief Checks whether a given key ID is valid within a YAML file.
 *
 * @param file_id Pointer to the index of the YAML file (read-only).
 * @param key_id Pointer to the key ID to validate (read-only).
 * @return true if the key ID is valid, false otherwise.
 *
 * @note A key ID is valid if it is between 0 and the number of keys in the file (inclusive).
 */
bool is_valid_key_id(const int *file_id, const int *key_id)
{
    int fid = *file_id;
    int kid = *key_id;

    return (kid > -1 && kid <= my_files.files[fid].nkeys);
}

/**
 * @brief Checks whether a given file ID is valid.
 *
 * @param file_id Pointer to the file ID to validate (read-only).
 * @return true if the file ID is valid, false otherwise.
 *
 * @note A file ID is valid if it is between 0 (inclusive) and nfiles (exclusive).
 */
bool is_valid_file_id(const int *file_id)
{
    int fid = *file_id;

    return (fid > -1 && fid < nfiles);
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
     my_files.files = (YamlFile*)calloc(1, sizeof(YamlFile));
  } else
  {
     my_files.files = realloc(my_files.files, (nfiles+1)*sizeof(YamlFile));
  }

  j = nfiles;
  *file_id =j;

/*  printf("Opening file: %s.\nThere are %i files opened.\n", filename, j); */
  file = fopen(filename, "r");
  if (file == NULL) return -1;

  if(!yaml_parser_initialize(&parser)) return -2;

  my_files.files[j].keys = (KeyValuePairs*)calloc(1, sizeof(KeyValuePairs));

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
        my_files.files[j].keys = realloc(my_files.files[j].keys, (i+1)*sizeof(KeyValuePairs));
        my_files.files[j].keys[i].key_id=i;
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
         my_files.files[j].keys = realloc(my_files.files[j].keys, (i+1)*sizeof(KeyValuePairs));
         my_files.files[j].keys[i].key_id=i;
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
       printf("Key_number:%i Parent_key:%i Parent_name:%s Key:%s Value:%s \n", my_files.files[j].keys[i].key_id, my_files.files[j].keys[i].parent_key, my_files.files[j].keys[i].parent_name, my_files.files[j].keys[i].key, my_files.files[j].keys[i].value);
  }
  printf("/\n");
   */

  nfiles = nfiles + 1;
/*  printf("closing file: %s\n", filename); */
  fclose(file);

  return 1;
}

#endif
