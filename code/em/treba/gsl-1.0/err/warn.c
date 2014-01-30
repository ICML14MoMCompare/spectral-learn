/* err/warn.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_message.h>

int gsl_warnings_off = 0;

void
gsl_warning (const char * reason, const char * file, int line, int gsl_errno)
{
  if (!gsl_warnings_off)
    {
      const char * error_string = gsl_strerror(gsl_errno) ;
      gsl_errno = 0;		/* stop complaints about unused variables */
      gsl_stream_printf ("WARNING", file, line, reason);
      gsl_stream_printf ("(ERRNO)", file, line, error_string);
    }
}


