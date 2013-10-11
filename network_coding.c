/*
Serval DNA network coding functions
Copyright (C) 2013 Paul Gardner-Stephen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


#include "serval.h"

// You need a separate nc structure for TX and RX
struct nc {
  // Define the size parameters of the network coding data structures
  uint32_t window_size; // limited to the number of bits we can fit in a uint32_t
  uint32_t datagram_size; // number of bytes in each fixed sized unit
  uint32_t queue_length; // maximum queued packets

  // For RX only
  uint32_t *linear_combinations; // bitmap of which packets are in each received random linear combination of packets

  // For both RX & TX

  uint32_t window_start;   // ID of first datagram in window
  uint32_t window_used;    // # of datagrams used in the window

  // Buffers for holding sent or received datagrams
  uint8_t **datagram_buffers; // queue_length long list of datagram_size byte buffers
  // Vector of which buffer is used for each datagram.
  // This indirection is used for RX where as we do the Gauss-Jordan reduction the order of
  // received datagrams may get rearranged, and we don't want to waste time exchanging buffer contents.
  uint32_t *buffer_numbers; // queue_length long list
};

void nc_free(struct nc *n)
{
  int i;

  if (!n) return;

  if (n->datagram_buffers) {
    for(i=0;i<n->queue_length;i++) free(n->datagram_buffers[i]);
    free(n->datagram_buffers);
    n->datagram_buffers=NULL;
  }
  if (n->buffer_numbers) {
    free(n->buffer_numbers);
    n->buffer_numbers=NULL;
  }
  free(n);
  return;
}

struct nc *nc_new(uint32_t window_size, uint32_t datagram_size, uint32_t queue_length)
{
  // Sanity check inputs
  if (window_size<0||window_size>32) return NULL;
  if (datagram_size<0||datagram_size>65536) return NULL;
  if (queue_length<(window_size+1)||queue_length>256) return NULL;

  // Allocate structure pre-zeroed
  struct nc *n=calloc(sizeof(struct nc),1);
  if (!n) return NULL;
  
  n->window_size=window_size;
  n->datagram_size=datagram_size;
  n->queue_length=queue_length;

  int i;

  // Allocate vector of pointers to buffers
  n->datagram_buffers=malloc(sizeof(uint32_t)*n->queue_length);
  if (!n->datagram_buffers) { nc_free(n); return NULL; }
  // Buffer numbers vector also
  n->buffer_numbers=malloc(sizeof(uint32_t)*n->queue_length);
  if (!n->buffer_numbers) { nc_free(n); return NULL; }
  // Initialise buffer numbers in cardinal order
  for(i=0;i<n->queue_length;i++) n->buffer_numbers[i]=i;

  // Allocate the datagram buffers themselves
  for(i=0;i<n->queue_length;i++)
    if ((n->datagram_buffers[i]=malloc(n->datagram_size))==NULL) 
      { nc_free(n); return NULL; }

  // Structure is now fully initialised, so return it.
  return n;
}

#define ROWS 40

int count=0;
int print_set(char *msg, unsigned int set[ROWS])
{
  int i;
  printf(">> %04d %s:\n",count++,msg);
  for(i=0;i<ROWS;i++) {
    printf("  %02d: 0x%08x ",
	   i,set[i]);
    int j;
    for(j=0;j<32;j++) printf("%0d",(set[i]>>(31-j))&1);
    printf("\n");
  }
  return 0;
}

// sort in reverse order
int comp_rows(const void *a,const void *b)
{
  const unsigned int *aa=a,*bb=b;
  if (*aa<*bb) return 1; 
  if (*aa>*bb) return -1;
  return 0;
}

int sort_set(unsigned int set[ROWS])
{
  qsort(set,ROWS,sizeof(unsigned int),comp_rows);
  return 0;
}

int reduce_set(unsigned int set[ROWS])
{
  int i,j;
  for(i=0;i<ROWS;i++)
    for(j=i+1;j<ROWS;j++)
      {
	if ((set[i]^set[j])<set[i]) { 
	  set[i]=set[i]^set[j];
	  // TODO: enable break and combine sort into process to make it 
	  // faster (probably). With separate sort, it is faster to XOR
	  // everything if it keeps making it smaller.
	  // break; 
	}
      }
  return 0;
}

int main(int argc,char **argv)
{
  unsigned int set[ROWS];
  int i;
  for(i=0;i<ROWS;i++) set[i]=random()|((random()&1)<<31);

  print_set("Original set",set);
  while(1) {
    sort_set(set);
    print_set("after sort",set);
    reduce_set(set);
    print_set("after reduce",set);
  }
}
