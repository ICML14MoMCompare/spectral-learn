/**************************************************************************/
/*   treba - probabilistic FSM and HMM training and decoding              */
/*   Copyright © 2013 Mans Hulden                                         */

/*   This file is part of treba.                                          */

/*   Treba is free software: you can redistribute it and/or modify        */
/*   it under the terms of the GNU General Public License version 2 as    */
/*   published by the Free Software Foundation.                           */

/*   Treba is distributed in the hope that it will be useful,             */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/*   GNU General Public License for more details.                         */

/*   You should have received a copy of the GNU General Public License    */
/*   along with treba.  If not, see <http://www.gnu.org/licenses/>.       */
/**************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "treba.h"

extern int g_input_format;
extern int g_output_format;

PROB output_convert(PROB x) {
    /* Internal format is log2 */
    switch (g_output_format) {
        case FORMAT_LOG2:   return(x);
        case FORMAT_LOG10:  return(0.30102999566398119521 * x);
        case FORMAT_LN:     return(0.69314718055994530942 * x);
        case FORMAT_REAL:   return(x <= SMRZERO_LOG ? 0 : EXP(x));
        case FORMAT_NLOG2:  return(x == 0 ? 0 : -x);
        case FORMAT_NLOG10: return(x == 0 ? 0 : 0.30102999566398119521 * -x);
        case FORMAT_NLN:    return(x == 0 ? 0 : 0.69314718055994530942 * -x);
    }
    return(x);
}

PROB input_convert(PROB x) {
    /* Internal format is log2 */
    switch (g_input_format) {
        case FORMAT_LOG2:   return(x);
        case FORMAT_LOG10:  return(3.32192809488736234787 * x);
        case FORMAT_LN:     return(1.44269504088896340736 * x);
        case FORMAT_REAL:   return(x == 0 ? SMRZERO_LOG : log2(x));
        case FORMAT_NLOG2:  return(-x);
        case FORMAT_NLOG10: return(x == 0 ? 0 : 3.32192809488736234787 * -x);
        case FORMAT_NLN:    return(x == 0 ? 0 : 1.44269504088896340736 * -x);
    }
    return(x);
}

char *file_to_mem(char *name) {
    FILE    *infile;
    size_t    numbytes;
    char *buffer;
    infile = fopen(name, "rb");
    if(infile == NULL) {
        fprintf(stderr,"Error opening file '%s'\n",name);
        return NULL;
    }
    fseek(infile, 0L, SEEK_END);
    numbytes = ftell(infile);
    fseek(infile, 0L, SEEK_SET);
    buffer = (char*)malloc((numbytes+1) * sizeof(char));
    if(buffer == NULL) {
        perror("Error reading file");
        return NULL;
    }
    if (fread(buffer, sizeof(char), numbytes, infile) != numbytes) {
        perror("Error reading file");
        return NULL;
    }
    fclose(infile);
    *(buffer+numbytes)='\0';
    return(buffer);
}

void hmm_print(struct hmm *hmm) {
    PROB thisprob;
    int i,j;
    for (i = 0; i < hmm->num_states; i++) {
	for (j = 0; j < hmm->num_states; j++) {
	    thisprob = *HMM_TRANSITION_PROB(hmm,i,j);
	    if (thisprob > SMRZERO_LOG) {
		if (g_output_format != FORMAT_REAL || output_convert(thisprob) > 0) {
		    printf("%i > %i %.17g\n", i,j, output_convert(thisprob));
		}
	    }
	}
    }
    for (i = 0; i < hmm->num_states; i++) {
	for (j = 0; j < hmm->alphabet_size; j++) {
	    thisprob = *HMM_EMISSION_PROB(hmm,i,j);
	    if (thisprob > SMRZERO_LOG) {
		if (g_output_format != FORMAT_REAL || output_convert(thisprob) > 0) {
		    printf("%i %i %.17g\n", i,j, output_convert(thisprob));
		}
	    }
	}
    }
}

void wfsa_print(struct wfsa *fsm) {
    int i,j,k;
    PROB thisprob;
    for (i = 0; i < fsm->num_states; i++) {
	for (j = 0; j < fsm->alphabet_size; j++) {
	    for (k = 0; k < fsm->num_states; k++) {
		thisprob = *(TRANSITION(fsm,i,j,k));
		if (thisprob > SMRZERO_LOG) {
		    if (g_output_format != FORMAT_REAL || output_convert(thisprob) > 0) {
			printf("%i %i %i %.17g\n", i,k,j, output_convert(thisprob));
		    }
		}
	    }
	}
    }
    for (i=0; i < fsm->num_states; i++) {
	thisprob = *(fsm->final_table + i);
	if (thisprob > SMRZERO_LOG) {
	    if (g_output_format != FORMAT_REAL || output_convert(thisprob) > 0) {
		printf("%i %.17g\n",i,output_convert(thisprob));
	    }
	}
    }
}

int char_in_array(char c, char *array) {
    int i;
    for (i = 0; *(array+i) != '\0'; i++) {
	if (c == *(array+i)) {
	    return 1;
	}
    }
    return 0;
}

int line_count_elements(char **ptr) {
    int i, elements;
    char seps[] = {'0','1','2','3','4','5','6','7','8','9','.','-','e','E','>','\0'};
    for (i = 0, elements = 0; *(*ptr+i) != '\n' && *(*ptr+i) != '\0'; i++) {
	if (!char_in_array(*(*ptr+i), seps)) {
	    continue;
	}
	if (!char_in_array(*(*ptr+i+1), seps)) {
	    elements++;
	}
    }
    if (**ptr == '#') {
	elements = 0;
    }
    if (*(*ptr+i) == '\0') {
	*ptr = *ptr+i;
	elements = -1;
    } else {
	*ptr = *ptr+i+1;
    }
    return(elements);
}

char *line_to_int_array(char *ptr, int **line, int *size) {
    /* Reads (destructively) a line of integers (separated by non-integers) and returns a malloced array  */
    /* of numbers (ints) with the line in it + a size count, and also a pointer to the next line.         */
    char *nextline, *startptr;
    int i, j, elements, lastnumber;
    if (*ptr == '\n') {
	*line = NULL;
	*size = 0;
	return ptr+1;
    }
    for (i = 0, elements = 0; ptr[i] != '\n' && ptr[i] != '\0'; i++) {
	/* If ptr+i is a digit and ptr+i+1 is not, we've seen a number */
	if (!isdigit(ptr[i])) {
	    continue;
	}
	if (!isdigit(ptr[i+1])) {
	    elements++;
	}
    }
    if (ptr[i] == '\0') {
	nextline = ptr+i;
    } else {
	nextline = ptr+i+1;
    }
    *size = elements;
    if (!elements) {
	return nextline;
    }
    *line = malloc(sizeof(int) * elements);
    for (i = 0, j = 0, lastnumber = 0; ; i++, j++) {
	/* Find next digit */
	startptr = ptr + i;
	while (!isdigit(ptr[i]) && ptr[i] != '\n' && ptr[i] != '\0') {
	    i++;
	    startptr = ptr+i;
	}
	if (ptr[i] == '\n' || ptr[i] == '\0') {
	    break;
	}
	/* Find number end */
	while (isdigit(ptr[i])) {
	    i++;
	}
	if (ptr[i] == '\n' || ptr[i] == '\0') {
	    lastnumber = 1;
	}
	ptr[i] = '\0';
	*(*line+j) = atoi(startptr);
	if (lastnumber)
	    break;
    }
    return(nextline);
}

