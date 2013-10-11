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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#define uint32_t unsigned int
#define uint8_t unsigned char

// You need a separate nc structure for TX and RX
struct nc {
  // Define the size parameters of the network coding data structures
  uint32_t window_size; // limited to the number of bits we can fit in a uint32_t
  uint32_t datagram_size; // number of bytes in each fixed sized unit
  uint32_t max_queue_size; // maximum queued packets

  // For RX only
  uint32_t *linear_combinations; // bitmap of which packets are in each received random linear combination of packets
  uint32_t queue_size; // number of remaining queued datagrams being decoded

  uint32_t max_recent_datagrams; // buffer of recently decoded datagrams
  uint32_t recent_datagrams_start;
  uint32_t recent_datagrams_count;
  uint8_t **recent_datagram_buffers; // max_recent_datagrams long list of datagram_size byte buffers
  uint32_t *recent_datagram_buffer_numbers; // max_recent_datagrams long list

  // For both RX & TX

  uint32_t window_start;   // ID of first datagram in window
  uint32_t window_used;    // # of datagrams used in the window

  // Buffers for holding sent or received datagrams
  uint8_t **datagram_buffers; // max_queue_size long list of datagram_size byte buffers
  // Vector of which buffer is used for each datagram.
  // This indirection is used for RX where as we do the Gauss-Jordan reduction the order of
  // received datagrams may get rearranged, and we don't want to waste time exchanging buffer contents.
  uint32_t *buffer_numbers; // max_queue_size long list
};

void nc_free(struct nc *n)
{
  int i;

  if (!n) return;

  if (n->datagram_buffers) {
    for(i=0;i<n->max_queue_size;i++) free(n->datagram_buffers[i]);
    free(n->datagram_buffers);
    n->datagram_buffers=NULL;
  }
  if (n->buffer_numbers) {
    free(n->buffer_numbers);
    n->buffer_numbers=NULL;
  }
  if (n->recent_datagram_buffers) {
    for(i=0;i<n->max_recent_datagrams;i++) free(n->recent_datagram_buffers[i]);
    free(n->recent_datagram_buffers);
    n->recent_datagram_buffers=NULL;
  }
  if (n->recent_datagram_buffer_numbers) {
    free(n->recent_datagram_buffer_numbers);
    n->recent_datagram_buffer_numbers=NULL;
  }
  free(n);
  return;
}

struct nc *nc_new(uint32_t window_size, uint32_t datagram_size, 
		  uint32_t max_queue_size, uint32_t recent_datagram_count)
{
  // Sanity check inputs
  if (window_size<0||window_size>32) return NULL;
  if (datagram_size<0||datagram_size>65536) return NULL;
  // max_queue_size MUST be at least window_size both for the RX algorithm to work,
  // and also to ensure that we allocate the right number of buffers.  Failure
  // would result in memory corruption.
  if (recent_datagram_count) {
    if (max_queue_size<(window_size+1)||max_queue_size>256) return NULL;
  } else {
    if (max_queue_size!=window_size||max_queue_size>256) return NULL;
  }
  if (recent_datagram_count<0||recent_datagram_count>window_size) return NULL;
  
  // Allocate structure pre-zeroed
  struct nc *n=calloc(sizeof(struct nc),1);
  if (!n) return NULL;
  
  n->window_size=window_size;
  n->datagram_size=datagram_size;
  n->max_queue_size=max_queue_size;

  int i;

  // Allocate vector of pointers to buffers
  n->datagram_buffers=malloc(sizeof(uint32_t*)*n->max_queue_size);
  if (!n->datagram_buffers) { nc_free(n); return NULL; }
  // Buffer numbers vector also
  n->buffer_numbers=malloc(sizeof(uint32_t*)*n->max_queue_size);
  if (!n->buffer_numbers) { nc_free(n); return NULL; }
  // Initialise buffer numbers in cardinal order
  for(i=0;i<n->max_queue_size;i++) n->buffer_numbers[i]=i;

  // Allocate the datagram buffers themselves
  for(i=0;i<n->max_queue_size;i++)
    if ((n->datagram_buffers[i]=malloc(n->datagram_size))==NULL) 
      { nc_free(n); return NULL; }

  if (recent_datagram_count) {
    // Allocate vector of pointers to buffers
    n->recent_datagram_buffers=malloc(sizeof(uint32_t*)*n->max_recent_datagrams);
    if (!n->recent_datagram_buffers) { nc_free(n); return NULL; }
    // Buffer numbers vector also
    n->recent_datagram_buffer_numbers
      =malloc(sizeof(uint32_t*)*n->max_recent_datagrams);
    if (!n->recent_datagram_buffer_numbers) { nc_free(n); return NULL; }
    // Initialise buffer numbers in cardinal order
    for(i=0;i<n->max_recent_datagrams;i++) n->recent_datagram_buffer_numbers[i]=i;
    
    // Allocate the datagram buffers themselves
    for(i=0;i<n->max_recent_datagrams;i++)
      if ((n->recent_datagram_buffers[i]=malloc(n->datagram_size))==NULL) 
	{ nc_free(n); return NULL; }
  }

  // Structure is now fully initialised, so return it.
  return n;
}

int nc_tx_enqueue_datagram(struct nc *n,unsigned char *d,int len)
{
  if (!n) return -1;
  if (!d) return -1;
  if (len!=n->datagram_size) return -1;
  // On the TX side the maximum number of queued packets
  // is the MINIMUM of the maximum queue size and the window
  // size.
  if (n->window_used>=n->max_queue_size) return -1;
  if (n->window_used>=n->window_size) return -1;
  
  // Add datagram to queue
  bcopy(d,n->datagram_buffers[n->buffer_numbers[n->window_used]],n->datagram_size);
  n->window_used++;
  return 0;
}

int nc_tx_ack_dof(struct nc *n,uint32_t latest_dof)
{
  if (!n) return -1;
  // Acknowledgement is for DOF that is yet to exist
  if (latest_dof>(n->window_start+n->window_size)) return -1;

  // Nothing to do if we have already acknowledged this degree of freedom
  if (latest_dof<n->window_start) return 0;

  // DOF is in window, so shift out all datagrams that preceed the datagram
  // indicated by the DOF
  int i;
  do {
    // Release buffer for the datagram in question, and shuffle the rest down one.
    // This could be done more efficiently by shuffling down many at once, but
    // the current code is clear and simple, and the computational cost is low.
    uint32_t freed_datagram_buffer=n->buffer_numbers[0];
    for(i=0;i<n->max_queue_size-1;i++) n->buffer_numbers[i]=n->buffer_numbers[i+1];
    n->buffer_numbers[n->max_queue_size-1]=freed_datagram_buffer;
    n->window_start++;
    n->window_used--;
  } while(n->window_start!=latest_dof);
  return 0;
}

int nc_tx_random_linear_combination(struct nc *n,uint8_t *combined_datagram,
				    uint32_t buffer_size,uint32_t *output_size)
{
  if (!n) return -1;
  if (!combined_datagram) return -1;
  // TODO: Don't waste more bytes than we need to on the bitmap and sequence number
  if (buffer_size<n->datagram_size+sizeof(uint32_t)+sizeof(uint32_t)) return -1;

  if (!n->window_used) {
    // Nothing to send, so just return
    *output_size=0;
    return 1;
  }

  // TODO: Check that combination is linearly independent of recently produced
  // combinations, i.e., that it contributes information.

  // get 32 bit random number.  random() only returns 31 bits, 
  // hence the double call and shift
  uint32_t combination=random()^(random()<<1);

  // restrict set bits to only those in the window
  // i.e., zero lower (32-n->window_used) bits
  combination=(combination>>(32-n->window_used))<<(32-n->window_used);

  // Never send all zeros, since that conveys no information.  
  if (!combination) {
    combination=0xffffffff;
    // restrict set bits to only those in the window
    // i.e., zero lower (32-n->window_used) bits
    combination=(combination>>(32-n->window_used))<<(32-n->window_used);
  }
  // should never be zero
  assert(combination!=0);

  // Now produce the output packet
  // Write out the current start of window
  combined_datagram[0]=(n->window_start>>24)&0xff;
  combined_datagram[1]=(n->window_start>>16)&0xff;
  combined_datagram[2]=(n->window_start>> 8)&0xff;
  combined_datagram[3]=(n->window_start>> 0)&0xff;
  // Write out bitmap of combinations involved
  combined_datagram[4]=(combination>>24)&0xff;
  combined_datagram[5]=(combination>>16)&0xff;
  combined_datagram[6]=(combination>> 8)&0xff;
  combined_datagram[7]=(combination>> 0)&0xff;
  // Produce linear combination
  bzero(&combined_datagram[sizeof(uint32_t)+sizeof(uint32_t)],n->datagram_size);
  int i,j;
  for(i=0;i<n->window_used;i++) {
    if ((combination<<i)&0x80000000) {
      // Combination includes this packet, so XOR with it
      for(j=0;j<n->datagram_size;j++)
	combined_datagram[sizeof(uint32_t)+sizeof(uint32_t)+j]
	  ^=n->datagram_buffers[n->buffer_numbers[i]][j];
    }
  }
  *output_size=n->datagram_size+sizeof(uint32_t)+sizeof(uint32_t);
  return 0;
}

int nc_release_recent_datagrams(struct nc *n,uint32_t combination_window_start)
{
  if (!n) return -1;

  if (combination_window_start>=n->recent_datagrams_start
      &&combination_window_start<=n->window_start)
    {
      // TODO: Shift down multiple slots at once when warranted to 
      // improve efficiency.
      while(combination_window_start>=n->recent_datagrams_start) {
	uint32_t freed_datagram_buffer=n->recent_datagram_buffer_numbers[0];
	int i;
	for(i=0;i<n->max_recent_datagrams-1;i++) 
	  n->recent_datagram_buffer_numbers[i]
	    =n->recent_datagram_buffer_numbers[i+1];
	n->recent_datagram_buffer_numbers[n->max_queue_size-1]
	  =freed_datagram_buffer;
	n->recent_datagrams_start++;
      }
    }
  return 0;
}

int nc_rx_linear_combination(struct nc *n,uint8_t *combination,int len)
{
  if (!n) return -1;
  if (!combination) return -1;
  if (len!=n->datagram_size+sizeof(uint32_t)+sizeof(uint32_t)) return -1;

  // Fail if there is no space
  if (n->queue_size>=n->max_queue_size) return -1;
  
  // Translate combination into our current window frame.
  // If it contains previously received datagrams, we should check if we
  // have them on hand, and if so we can reduce those from the combination.
  // If they are not on hand, then we cannot process the datagram.
  uint32_t combination_window_start
    =(combination[0]<<24)|(combination[1]<<16)|(combination[2]<<8)|combination[3];
  uint32_t combination_bitmap
    =(combination[4]<<24)|(combination[5]<<16)|(combination[6]<<8)|combination[7];

  nc_release_recent_datagrams(n,combination_window_start);

  while(combination_window_start<n->window_start) {
    if (combination_bitmap&0x80000000) {
      // TODO: Look through recently decoded datagrams, and reduce from this
      // combination. If we cannot, then we cannot decode this 
      // Until we can do that, we must simply fail datagrams that start before
      // the current window.
	if (combination_window_start==n->recent_datagrams_start
	    &&(combination_window_start<
	       (n->recent_datagrams_start+n->recent_datagrams_count)))
	  {
	    // We have the datagram - reduce it from the current combination
	    int i;
	    int buffer
	      =n->recent_datagram_buffer_numbers[combination_window_start
						 -n->recent_datagrams_start];
	    for(i=0;i<n->datagram_size;i++) 
	      combination[8+i]^=n->recent_datagram_buffers[buffer][i];
	    // Then advance the combination start and shift the bitmap accordingly
	    combination_bitmap=combination_bitmap<<1;
	    combination_window_start++;
	  }
	else
	  return -1;
    }    
    // If there is no information left in the combination, then stop processing it
    if (!combination_bitmap) return -1;
  }

  // Attempt to reduce datagram before inserting.
  // We know that previous calls leave the list sorted in descending order,
  // so we can perform the reduction efficiently.

  // Add to queue
  //  uint32_t buffer_number=n->buffer_numbers[n->buffer_numbers[n->

  return -1;
}

/* TODO: Tests that should be written.
   1. nc_new() works, and rejects bad input.
   2. nc_free() works, including on partially initialised structures.
   3. nc_tx_enqueue_datagram() works, including failing on bad input and when the
      queue is full.
   4. nc_tx_ack_dof() works, rejects bad input, and correctly releases buffers.
   5. nc_tx_random_linear_combination() works, rejects bad input, and produces valid
      linear combinations of the enqueued datagrams, and never produces all zeroes.
   6. nc_rx_linear_combination() works, rejects bad input
   7. nc_rx_linear_combination() rejects when RX queue full, when combination starts
      before current window.
*/

int nc_test_dump(char *name, unsigned char *addr, int len)
{
  int i,j;
  if (name)
    fprintf(stderr,"Dump of %s\n",name);
  for(i=0;i<len;i+=16){
    fprintf(stderr,"  %04x :",i);
    for(j=0;j<16&&(i+j)<len;j++) 
      fprintf(stderr," %02x",addr[i+j]);
    for(;j<16;j++) 
      fprintf(stderr,"   ");
    fprintf(stderr,"    ");
    for(j=0;j<16&&(i+j)<len;j++)
      fprintf(stderr,"%c",addr[i+j]>=' '&&addr[i+j]<0x7f?addr[i+j]:'.');
    fprintf(stderr,"\n");
  }
  return 0;
}

int nc_test_random_datagram(uint8_t *d,int len)
{
  int i;
  for(i=0;i<len;i++) d[i]=random()&0xff;
  return 0;
}

int nc_test()
{
  struct nc *rx,*tx;

  rx=nc_new(32,200,40,32);
  if (!rx) {
    fprintf(stderr,"FAIL: Failed to create validly defined nc struct for RX.\n");
    fprintf(stderr,"FATAL: Cannot continue tests.\n");
    return -1;
  } else fprintf(stderr,"PASS: Created validly defined nc struct for RX.\n");
  tx=nc_new(32,200,32,0);
  if (!tx) {
    fprintf(stderr,"FAIL: Failed to create validly defined nc struct for TX.\n");
    fprintf(stderr,"FATAL: Cannot continue tests.\n");
    return -1;
  } else fprintf(stderr,"PASS: Created validly defined nc struct for TX.\n");

  // Prepare some random datagrams for subsequent tests
  int i;
  uint8_t adatagram[200];
  uint8_t bdatagram[200];
  uint8_t cdatagram[200];
  uint8_t ddatagram[200];
  uint8_t edatagram[200];
  nc_test_random_datagram(adatagram,200);
  nc_test_random_datagram(bdatagram,200);
  nc_test_random_datagram(cdatagram,200);
  nc_test_random_datagram(ddatagram,200);
  nc_test_random_datagram(edatagram,200);

  // Test inserting datagrams into the queue

  if (nc_tx_enqueue_datagram(tx,adatagram,200)) {
    fprintf(stderr,"FAIL: Failed to enqueue datagram for TX\n");
    return -1;
  } else fprintf(stderr,"PASS: Enqueued datagram for TX\n");
  if (tx->window_used==1) 
    fprintf(stderr,"PASS: Enqueueing datagram increases window_used\n");
  else
    fprintf(stderr,"FAIL: Enqueueing datagram increases window_used\n");

  int j,fail=0;
  for(i=0;i<10;i++) {
    uint8_t outbuffer[4+4+200];
    int len=4+4+200;
    uint32_t written=0;
    fail=nc_tx_random_linear_combination(tx,outbuffer,len,&written);
    if (fail) { 
      fprintf(stderr,"FAIL: Produce random linear combination of single packet for TX\n");
      break;
    }
    if (!outbuffer[4]) {
      fprintf(stderr,"FAIL: Should not produce empty linear combination bitmap\n");
      nc_test_dump("output combination with headers",outbuffer,written);
      fail=1;
      break;
    }
    for(j=0;j<200;j++) if (outbuffer[8+j]!=adatagram[j]) break;
    if (j<200) {
      fprintf(stderr,"FAIL: Output identity datagram when only one in queue\n");
      nc_test_dump("output combination with headers",outbuffer,written);
      nc_test_dump("original datagram",adatagram,200);
      printf("Error was in %dth byte\n",j);
      fail=1;
      break;
    }
  }
  if (fail) return -1;
  fprintf(stderr,"PASS: Produce random linear combination of single packet for TX\n");
  fprintf(stderr,"PASS: Should not produce empty linear combination bitmap\n");
  fprintf(stderr,"PASS: Output identity datagram when only one in queue\n");

  return 0;
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
  return nc_test();

  /*
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
  */
}
