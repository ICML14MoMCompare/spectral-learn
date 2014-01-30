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

#include <stdlib.h>

#include "treba.h"

/* Data structures for storing observation sequences */

int obssortcmp(struct observations **a, struct observations **b) {
    int *one, *two, i;
    one = (*a)->data; two = (*b)->data;
    for (i = 0; i < (*a)->size && i < (*b)->size; i++) { 
	if (*(one+i) < *(two+i)) {
	    return -1;
	}
	if (*(one+i) > *(two+i)) {
	    return 1;
	}
    }
    if ((*a)->size > (*b)->size) {
	return -1;
    }
    if ((*a)->size < (*b)->size) {
	return 1;
    }
    return 0;
}

int observations_alphabet_size(struct observations *ohead) {
    int i, maxsigma;
    for (maxsigma = -1; ohead != NULL ; ohead = ohead->next) {
	for (i = 0; i < ohead->size; i++) {
	    maxsigma = *(ohead->data+i) > maxsigma ? *(ohead->data+i) : maxsigma;		
	}
    }
    return(maxsigma+1);
}

struct observations **observations_to_array(struct observations *ohead, int *numobs) {
    int i;
    struct observations **obsarray, *o;
    for (o = ohead, i = 0; o != NULL; o = o->next, i++) { }
    *numobs = i;
    obsarray = malloc(sizeof(struct observations *) * i) ;
    for (o = ohead, i = 0; o != NULL; o = o->next, i++) {
	*(obsarray+i) = o;
    }
    return obsarray;
}

struct observations *observations_uniq(struct observations *ohead) {
    struct observations *o, *onext;
    for (o = ohead; o != NULL && o->next != NULL; ) {
	onext = o->next;
	if (onext != NULL) {
	    if (obssortcmp(&o,&onext) == 0) {
		o->next = onext->next;
		free(onext->data);
		free(onext);
		o->occurrences = o->occurrences + 1;
		continue;
	    }
	}
	o = o->next;
    }
    return (ohead);
}

struct observations *observations_sort(struct observations *ohead) {
    struct observations *o, **sptr, *curro, *lasto;
    int i, numobs;
    int (*sorter)() = obssortcmp;
    for (o = ohead, numobs = 0; o != NULL; o = o->next) {
	numobs++;
    }
    if (numobs < 1) {
	return ohead;
    }
    sptr = malloc(sizeof(struct observations *) * numobs);
    for (o = ohead, i = 0; o != NULL; o = o->next, i++) {
	*(sptr+i) = o;
    }
    qsort(sptr, numobs, sizeof(struct observations *), sorter);
    for (i = 1, ohead = lasto = *sptr, lasto->next = NULL; i < numobs; i++) {	
	curro = *(sptr+i);
	curro->next = NULL;
	lasto->next = curro;
	lasto = curro;
    }
    return(ohead);
}


void observations_destroy(struct observations *ohead) {
    struct observations *o, *olast;
    for (o = ohead; o != NULL; ) {
	olast = o;
	o = o->next;
	free(olast->data);
	free(olast);
    }
}

struct observations *observations_read(char *filename) {
    char *obs_char_data, *optr;
    struct observations *ohead, *o, *olast;
    int *line, size, comment;
    if ((obs_char_data = file_to_mem(filename)) == NULL) {
	return(NULL);
    }
    ohead = olast = NULL;
    for (optr = obs_char_data, o = ohead ; ;) {
	if (*optr == '\0') {
	    break;
	}
	comment = 0;
	if (*optr == '#')
	    comment = 1;
	optr = line_to_int_array(optr, &line, &size);
	if (optr == NULL) {
	    break;
	}
	if (comment)
	    continue;
	if (olast == NULL) {
	    ohead = olast = malloc(sizeof(struct observations));
	    o = olast;
	} else {
	    o = malloc(sizeof(struct observations));
	    olast->next = o;
	    olast = o;
	}
	o->size = size;
	o->data = line;	
	o->occurrences = 1;	
	o->next = NULL;	
    }
    free(obs_char_data);
    return(ohead);
}
