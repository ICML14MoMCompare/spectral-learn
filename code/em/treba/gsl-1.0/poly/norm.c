/* poly/norm.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

static double norm_companion (const double *m, size_t nc);

static double
norm_companion (const double *m, size_t nc)
{
  size_t i;
  double A, s2 = 0;

  for (i = 1; i < nc - 1; i++)
    {
      const double t = MAT (m, i + 1, i, nc);
      s2 += t * t;
    }

  for (i = 0; i < nc; i++)
    {
      const double t = MAT (m, i, nc, nc);
      s2 += t * t;
    }

  A = sqrt (s2);

  return A;
}

