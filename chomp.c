/*
  (C) 2010 Iñaki García Etxebarria

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <string.h>
#include "chomp.h"

#define wrong_input(buf) do {\
  fprintf(stderr, "[%s:%d] Wrong input!!\n", __FILE__, __LINE__);\
  fprintf(stderr, "<<%s>>\n", buf);				 \
  abort();} while(0)

/* Returns nonzero iff buf starts with str. Both strings should be
   zero terminated. */
static int startswith(const char *buf, const char *str)
{
  while (1) {
    if (*str == 0)
      return 1;

    if (*buf == 0 || *buf != *str)
      return 0;

    buf++;
    str++;
  }
}

static int process_output(FILE *fd)
{
  char *line = NULL;
  size_t nline = 0;
  int result=0;

  if (getline(&line, &nline, fd) < 0)
    wrong_input(line);

  if (!startswith(line, "HOMCHAIN"))
    wrong_input(line);

  while (1) {
    char *p;
    if (getline(&line, &nline, fd) < 0)
      break;

    if (!startswith(line, "H_1 = "))
      continue;

    p = line + strlen("H_1 = ");

    if (startswith(p, "Z^")) {
      char *p2;
      
      p += strlen("Z^");

      result = strtol(p, &p2, 10);

      if (p == p2)
	wrong_input(line);
      
      break;
    } else if (!strcmp(p, "Z\n")) {
      result = 1;
      break;
    } else if (startswith(p, "0")) {
      result = 0;
      break;
    } else
      wrong_input(line);
  }

  free(line);

  /* If we reach this we didn't find H^1 in the output, meaning that
     it is 0. */
  return result;
}

int run_chomp(const char *fname)
{
  FILE *homchain;
  char *cmd;
  int result;

  assert(asprintf(&cmd, "homchain -d %s", fname) > 0);

  homchain = popen(cmd, "r");

  result = process_output(homchain);

  pclose(homchain);

  free(cmd);

  return result;
}
