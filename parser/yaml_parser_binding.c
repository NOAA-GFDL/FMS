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

// TODO: These constants are also defined in yaml_parser.F90
// Consider consolidating constant definitions or using a shared header.
const int MISSING_FILE = -1;
const int PARSER_INIT_ERROR = -2;
const int INVALID_YAML = -3;
const int INVALID_ALIAS = -4;
const int MAX_LEVELS_REACH = -5;
const int SUCCESSFUL = 1;

#ifndef MAX_LEVELS
#define MAX_LEVELS 10
#endif

// DEBUG_FMS_YAML_PARSER is a hidden macro that may be useful when debugging parser issues
#ifdef DEBUG_FMS_YAML_PARSER
  #define DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
  #define DEBUG_PRINT(...) // nothing
#endif

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
 * @brief Represents the contents of a YAML anchor as an array of key-value pairs.
 */
typedef struct {
    int nkeys;                          ///< Number of keys defined for the anchor
    KeyValuePairs *keys;                ///< Array of key-value pairs
    char anchor_name[FMS_FILE_LEN];     ///< Name of the anchor
    int nlevels;
    int pid[MAX_LEVELS];
    char parent_names[MAX_LEVELS][FMS_FILE_LEN];
} AnchorsType;

/**
 * @brief Represents the contents of a YAML file as an array of key-value pairs.
 */
typedef struct {
    int nkeys;                ///< Number of keys defined in the YAML file
    KeyValuePairs *keys;      ///< Array of key-value pairs
    AnchorsType *Anchors;     ///< Array of anchors in the yaml file
    int nanchors;              ///< Number of anchors that have been defined in the file
} YamlFile;

/**
 * @brief Represents all YAML files that have been opened and parsed.
 */
typedef struct {
    YamlFile *files;          ///< Array of parsed YAML files
} FileType;

// Local Variables:
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

/**
 * @brief Increments the number of levels, enforcing a maximum limit.
 *
 * Increments the value pointed to by `nlevels`. If the new value exceeds
 * the maximum allowed (`MAX_LEVELS`), the function returns an error code.
 *
 * @param nlevels   Pointer to the current number of levels to be incremented.
 *
 * @return SUCCESSFUL (typically 0) if increment is valid;
 *         MAX_LEVELS_REACH if the maximum level count is exceeded.
 */
int increment_nlevels(int *nlevels) {
    (*nlevels) ++;
    if (*nlevels > MAX_LEVELS){
        return MAX_LEVELS_REACH;
    }
    return SUCCESSFUL;
}

/**
 * @brief Initializes an AnchorsType instance.
 *
 * @param anchor Pointer to the AnchorsType to initialize.
 * @param name Name of the anchor (null-terminated string).
 */
void init_anchor(AnchorsType *anchor, const char *name, const char *parent_name)
{
    if (anchor == NULL || name == NULL) return;

    anchor->nkeys = 0;
    anchor->keys = (KeyValuePairs*)calloc(1, sizeof(KeyValuePairs));
    strcpy(anchor->anchor_name, name);
    anchor->nlevels = 0;
    anchor->pid[0] = 0;
    strcpy(anchor->parent_names[0], parent_name);
}

/**
 * @brief Populates a KeyValuePairs structure with provided key data.
 *
 * @param my_key        Pointer to a KeyValuePairs structure to populate.
 * @param key_id        Integer identifier for the key.
 * @param parent_key    Integer identifier for the parent key.
 * @param key           String representing the key name (can be NULL).
 * @param value         String representing the key value (can be NULL).
 * @param parent_name   String representing the parent key's name (can be NULL).
 */
void add_key(KeyValuePairs *my_key, const int key_id, const int parent_key,
             const char *key, const char *value, const char *parent_name) {

    my_key->key_id = key_id;
    my_key->parent_key = parent_key;

    if (key) {
        strncpy(my_key->key, key, sizeof(my_key->key) - 1);
        my_key->key[sizeof(my_key->key) - 1] = '\0';
    } else {
        my_key->key[0] = '\0';
    }

    if (value) {
        strncpy(my_key->value, value, sizeof(my_key->value) - 1);
        my_key->value[sizeof(my_key->value) - 1] = '\0';
    } else {
        my_key->value[0] = '\0';
    }

    if (parent_name) {
        strncpy(my_key->parent_name, parent_name, sizeof(my_key->parent_name) - 1);
        my_key->parent_name[sizeof(my_key->parent_name) - 1] = '\0';
    } else {
        my_key->parent_name[0] = '\0';
    }
}

/**
 * @brief Populates a KeyValuePairs structure in anchor with a key/value
 *
 * @param anchor Pointer to the AnchorsType to populate
 * @param key Key of the yaml key/value pair
 * @param value Value of the yaml key/value pair
 */
void add_anchor_key(AnchorsType *anchor, const char *key, const char *value)
{
    anchor->nkeys++;
    anchor->keys = realloc(anchor->keys, (anchor->nkeys+1)*sizeof(KeyValuePairs));

    KeyValuePairs *my_key = &anchor->keys[anchor->nkeys];
    add_key(my_key, anchor->nkeys, anchor->pid[anchor->nlevels],
      key, value, "");

    DEBUG_PRINT("ANCHOR :: Key_number: %i, parent_key: %i, %s:%s \n ", my_key->key_id, my_key->parent_key, my_key->key, my_key->value);
}

/**
 * @brief Populates a KeyValuePairs structure in anchor with a new block
 *
 * @param anchor Pointer to the AnchorsType to populate
 * @param key Name of the block
 * @return 1 if successful otherwise error code
 */
int add_anchor_parent(AnchorsType *anchor, const char *key)
{
    anchor->nkeys++;

    anchor->keys = realloc(anchor->keys, (anchor->nkeys+1)*sizeof(KeyValuePairs));
    int err_code = increment_nlevels(&anchor->nlevels);
    if (err_code =! SUCCESSFUL) return err_code;

    anchor->pid[anchor->nlevels] = anchor->nkeys;

    if (strcmp(key, "")) {
        strcpy(anchor->parent_names[anchor->nlevels],key );
    }
    KeyValuePairs *my_key = &anchor->keys[anchor->nkeys];
    add_key(my_key, anchor->nkeys, anchor->pid[anchor->nlevels -1],
      "", "", anchor->parent_names[anchor->nlevels]);
    DEBUG_PRINT("ANCHOR :: Key_number: %i, parent_key: %i, parent_name: %s \n ", my_key->key_id, my_key->parent_key, my_key->parent_name);

    return SUCCESSFUL;
}
/**
 * @brief Retrieves the index of an anchor by its alias name.
 *
 * Searches the YamlFile's list of anchors for a matching alias name.
 * If a match is found, returns the corresponding index (starting from 1).
 *
 * @param this         Pointer to the YamlFile structure containing anchors.
 * @param alias_name   The alias name to search for.
 *
 * @return Index of the matching anchor if found; otherwise, returns INVALID_ALIAS.
 */
int get_anchor_id(YamlFile *this, const char *alias_name)
{
    for (int i = 1; i < this->nanchors + 1; i++) {
        AnchorsType *my_anchor = &this->Anchors[i];
        if (strcmp(my_anchor->anchor_name, alias_name) == 0) {
            return i;
        }
    }
    return INVALID_ALIAS;
}

/**
 * @brief Opens and parses a YAML file.
 *
 * @param filename Pointer to the name of the YAML file (read-only).
 * @param file_id Pointer to an integer where the assigned file ID will be stored (output).
 * @return 1 if the file was read successfully, or < 0 if there was an error.
 */
int open_and_parse_file_wrap(const char *filename, int *file_id)
{
  yaml_parser_t parser;
  yaml_token_t  token;

  int fid;                       /* To minimize the typing :) */

  // Allocate space to store all the yaml file's info
  if (nfiles == 0 )
  {
     my_files.files = (YamlFile*)calloc(1, sizeof(YamlFile));
  } else
  {
     my_files.files = realloc(my_files.files, (nfiles+1)*sizeof(YamlFile));
  }

  // Assign the file id
  fid = nfiles;
  *file_id =fid;

  DEBUG_PRINT("Opening file: %s.\n There are %i files opened.\n", filename, nfiles);

  FILE *file;
  file = fopen(filename, "r");
  if (file == NULL) return MISSING_FILE;

  if(!yaml_parser_initialize(&parser)) return PARSER_INIT_ERROR;

  YamlFile *my_file = &my_files.files[fid];
  my_file->keys = (KeyValuePairs*)calloc(1, sizeof(KeyValuePairs));

  int nlevels = 0;
  int nkeys = 0;
  int pid[MAX_LEVELS]; // parent_ids
  char parent_names[MAX_LEVELS][FMS_FILE_LEN]; // the name of the parent at each level
  char key_value[FMS_FILE_LEN];
  bool defining_value;
  bool defining_anchor;

  // Initialize some variables
  defining_value = false;
  defining_anchor = false;
  pid[0]=0;
  strcpy(parent_names[0], "TOP");

  yaml_parser_set_input_file(&parser, file);
  do {
    if (!yaml_parser_scan(&parser, &token)) {
      return INVALID_YAML;
    }
    switch(token.type)
    {
    case YAML_KEY_TOKEN:
        DEBUG_PRINT("YAML_KEY_TOKEN \n");
        defining_value = false;
        break;

    case YAML_VALUE_TOKEN:
        DEBUG_PRINT("YAML_VALUE_TOKEN \n");
        defining_value = true;
        break;

    case YAML_SCALAR_TOKEN:
        DEBUG_PRINT("YAML_SCALAR_TOKEN \n");
        if (!defining_value) {
         strcpy(key_value, token.data.scalar.value);
         break;
        }
        if (defining_anchor) {
           AnchorsType *my_anchor = &my_file->Anchors[my_file->nanchors];
           add_anchor_key(my_anchor, key_value, token.data.scalar.value);
        } else {
           nkeys ++;
           my_file->keys = realloc(my_file->keys, (nkeys+1)*sizeof(KeyValuePairs));
           KeyValuePairs *my_key = &my_file->keys[nkeys];
           add_key(my_key, nkeys, pid[nlevels],
             key_value, token.data.scalar.value, "");
           DEBUG_PRINT("Key_number: %i, parent_key: %i, ----- %s:%s \n ",
             my_file->keys[nkeys].key_id, my_file->keys[nkeys].parent_key,
             my_file->keys[nkeys].key, my_file->keys[nkeys].value);
        }
        defining_value = false;
        strcpy(key_value, "" );
        break;

    case YAML_BLOCK_ENTRY_TOKEN:
        DEBUG_PRINT("YAML_BLOCK_ENTRY_TOKEN \n");
        if (defining_anchor) {
           AnchorsType *my_anchor = &my_file->Anchors[my_file->nanchors];

           int err_code = add_anchor_parent(my_anchor, key_value);
           if (err_code != SUCCESSFUL) return err_code;

        } else {
           int err_code = increment_nlevels(&nlevels);
           if (err_code == MAX_LEVELS_REACH) return MAX_LEVELS_REACH;

           nkeys ++;
           pid[nlevels] = nkeys;
           if (strcmp(key_value, "")) {
             strcpy(parent_names[nlevels], key_value);
           }
           my_file->keys = realloc(my_file->keys, (nkeys+1)*sizeof(KeyValuePairs));
           KeyValuePairs *my_key = &my_file->keys[nkeys];
           add_key(my_key, nkeys, pid[nlevels-1],
             "", "", parent_names[nlevels]);
           DEBUG_PRINT("Key_number: %i, parent_key: %i, parent_name: %s \n ", my_file->keys[nkeys].key_id, my_file->keys[nkeys].parent_key, my_file->keys[nkeys].parent_name);
        }
        defining_value = false;
        strcpy(key_value, "" );
        break;
    case YAML_BLOCK_END_TOKEN:
        DEBUG_PRINT("YAML_BLOCK_END_TOKEN \n");
        if (defining_anchor) {
            AnchorsType *my_anchor = &my_file->Anchors[my_file->nanchors];
            my_anchor->nlevels--;
            if (my_anchor->nlevels == -1) {
                defining_anchor = false;
                DEBUG_PRINT("FINISHED WITH ANCHOR :: ----------------------------- \n");
            }
        } else {
            nlevels --;
        }
        defining_value = false;
        strcpy(key_value, "" );
        break;
    case YAML_ANCHOR_TOKEN: {
        DEBUG_PRINT("YAML_ANCHOR_TOKEN \n");
        my_file->nanchors ++;
        if (my_file->nanchors == 1) {
           my_file->Anchors = (AnchorsType*)calloc(my_file->nanchors + 1, sizeof(AnchorsType));
        } else {
           my_file->Anchors = realloc(my_file->Anchors, (my_file->nanchors + 1)*sizeof(AnchorsType));
        }
        defining_anchor = true;
        defining_value = false;
        AnchorsType *my_anchor = &my_file->Anchors[my_file->nanchors];
        init_anchor(my_anchor, token.data.anchor.value, key_value);
        break;
    }
    case YAML_ALIAS_TOKEN: {
        DEBUG_PRINT("YAML_ALIAS_TOKEN \n");
        int top_key = nkeys;

        int aid = get_anchor_id(my_file, token.data.alias.value);
        if (aid == INVALID_ALIAS) return INVALID_ALIAS;

        AnchorsType *my_anchor = &my_file->Anchors[aid];
        for (int i = 2; i < my_anchor->nkeys + 1; i++) {
            nkeys ++;
            my_file->keys = realloc(my_file->keys, (nkeys+1)*sizeof(KeyValuePairs));
            KeyValuePairs *my_key = &my_file->keys[nkeys];
            int parent_key_id;
            char parent_name[FMS_FILE_LEN];
            if (my_anchor->keys[i].parent_key == 0){
                strcpy(parent_name, parent_names[nlevels]);
                parent_key_id = pid[nlevels-1];
            } else {
                strcpy(parent_name, my_anchor->keys[i].parent_name);
                parent_key_id = top_key + my_anchor->keys[i].parent_key - 1;
            }
            add_key(my_key, nkeys, parent_key_id,
              my_anchor->keys[i].key, my_anchor->keys[i].value, parent_name);
            if (strcmp(my_file->keys[nkeys].key, "")) {
                DEBUG_PRINT("**:: Key_number: %i, parent_key: %i, ----- %s:%s \n ",
                    my_file->keys[nkeys].key_id, my_file->keys[nkeys].parent_key,
                    my_file->keys[nkeys].key, my_file->keys[nkeys].value);
            }else {
                DEBUG_PRINT("**:: Key_number: %i, parent_key: %i, parent_name: %s \n ",
                    my_file->keys[nkeys].key_id, my_file->keys[nkeys].parent_key,
                    my_file->keys[nkeys].parent_name);
            }
        }
        nlevels --;

        break;
    }
    } //end of switch
    if(token.type != YAML_STREAM_END_TOKEN)
      yaml_token_delete(&token);
  } while(token.type != YAML_STREAM_END_TOKEN);
  yaml_token_delete(&token);
  yaml_parser_delete(&parser);

  my_file->nkeys = nkeys;
  nfiles = nfiles + 1;
  DEBUG_PRINT("closing file: %s\n", filename);
  fclose(file);

  return SUCCESSFUL;
}

#endif
