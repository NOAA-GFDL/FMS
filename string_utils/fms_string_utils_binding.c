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
#include <stdlib.h>
#include <string.h>
#include <math.h>

// struct to store a string and id associated with that string
typedef struct{
  char arr_name[255];
  int id;
}my_type;

// Compares two my_type types by the arr_name
static int arr_name_sorter(const void* p1, const void* p2)
{
  const my_type *the_type1 = p1;
  const my_type *the_type2 = p2;

  return strcmp(the_type1->arr_name, the_type2->arr_name);
}

// Sorts an array of strings in alphabetical order
// Implements a binary search to search for a string in an array of strings
// arr -> pointer of character array
// n -> length of the array
// id - > indices of the character array
void fms_sort_this(char **arr, int* n, int* id)
{
  int i; // For do loops
  my_type *the_type;

  // Save the array and the id into a struct
  the_type = (my_type*)calloc(*n, sizeof(my_type));
    for(i=0; i<*n; i++){
      the_type[i].id = id[i];
      strcpy(the_type[i].arr_name, arr[i]);
    }

  qsort(the_type, *n, sizeof(my_type), arr_name_sorter);

  // Copy the sorted array and the sorted ids
  for(i=0; i<*n; i++){
    id[i] = the_type[i].id;
    strcpy(arr[i], the_type[i].arr_name);
  }
}

// Implements a binary search to search for a string in an array of strings
// arr -> pointer of character array
// n -> length of the array
// find me -> string to find
// np -> the number of times the string was found
// returns a string with the indices ;)
char* fms_find_my_string_binding(char** arr, int *n, char *find_me, int *np)
{
  int L= 0;     // Left bound
  int R = *n;   // Right bound
  int m;        // Middle of the bound
  int mm;       // Index currently looking at
  int is_found; // Result from strcmp: 0 if string was found  <0 if the string is "less" >0 if the string is "greater"
  int *p;       // Array to store the indices
  int i;        // For do loops

  *np = 0;
  is_found = -1;
  while(L != R){
    // Start looking in the midle of the array
	  m = ceil((L + R) / 2);
    //printf("L is set to %i from L=%i and R=%i \n", m, L, R);

    //printf("Checking %i:%s \n", m, (arr[m]));
	  is_found = strcmp(find_me,(arr[m]));
	  if (is_found == 0)
	  {
      *np = 1;
      p = malloc(sizeof(int) * *np);
      p[*np-1] = m + 1; //Because fortran indices start at 1 ;)
      //printf("Array found at %i %i %i \n", *np, m, p[*np-1]);

      // The string can be found in multiple indices of the array, so look to the left of the index where the string
      // was initially found
      mm = m;
      while (is_found == 0) {
        if (mm != 0) { // Only look to the left if m is not the begining of the array
			    mm= mm -1;
          is_found = strcmp(find_me,(arr[mm]));
			    if (is_found == 0 ) {
            *np = *np + 1;
            p = realloc(p, sizeof(int) * *np);
            p[*np-1] = mm + 1;
            //printf("Array found at %i %i %i\n", *np, mm, p[*np-1]);
          }
        } else {is_found = -999;} //Done looking
		  }
      // The string can be found in multiple indices of the array, so look to the right of the index where the string was
      // initially found
		  mm = m;
		  is_found =0;
      while (is_found == 0) {
        if (mm != *n-1) { // Only look to the right if m is not the end of the array
          mm = mm + 1;
          is_found = strcmp(find_me,(arr[mm]));
          if (is_found == 0 ) {
            *np = *np + 1;
            p = realloc(p, sizeof(int) * *np);
		        p[*np-1] = mm + 1;
            //printf("Array found at %i %i %i\n", *np, mm, p[*np-1]);
         }
        } else {is_found = -999;} //Done looking}
      }
		  L = R;
      // If find_me is greater than arr[m] (i.e find_me="potato" is greater than arr[m]="banana")
	  } else if (is_found > 0) {
      // Set the lower bound to start in m (ignore the first half)
  	  L = m + 1;
      //printf("L is set to %i \n", L);
    } else
    // If find_me is less than arr[m] (i.e find_me="potato" is less than arr[m] = "soccer")
    {
      // Set the upper bound to start in m (ignore the lower half)
      R = m;
      //printf("R is set to %i \n", R);
	}

}

  // This is the magical part:
  // Save the array of indices where the string was found into a string
  // The fortran side is going to allocate the array to the correct size and read the string into the array
  // The alternative (normal) way is to have a seperate function that gets the number of times the string is found
  // The fortran side will allocate the array to the correct size and send that into another function that
  // fill in that array. That will require you to search through the array twice ...
  char string[255];
  char *string_p;

  strcpy(string, "");

  for(i=0; i<*np; i++){
    if (i == *np-1) {sprintf( &string[ strlen(string) ],  "%d ", p[i] );}
    else {sprintf( &string[ strlen(string) ],  "%d ,", p[i] );}
  }

  string_p = (char*) malloc((strlen(string)+1)*sizeof(char));
  strcpy(string_p, string);
  return string_p;
}

// Finds the number of unique strings in an array
int fms_find_unique(char** arr, int *n)
{
  int i;
  int nfind;

  nfind=1;
  //printf("n is %i", *n);
  for(i=1; i<*n; i++){
    //printf("Comparing %s and %s \n",arr[i], arr[i-1]);
    if (strcmp(arr[i], arr[i-1]) != 0){ nfind = nfind + 1;}
  }

  //printf("nfind=%i",nfind);
  return nfind;
}
