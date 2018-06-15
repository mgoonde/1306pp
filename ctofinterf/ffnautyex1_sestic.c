/* This program prints generators for the automorphism group of an
   n-vertex polygon, where n is a number supplied by the user.

   This version uses a fixed limit for MAXN.
*/

#define MAXN 1000    /* Define this before including nauty.h */
#include "nauty.h"   /* which includes <stdio.h> and other system files */
#include "naututil.h" /* include hashgraph */

void ffnautyex1_sestic(int *coco, int connect[*coco][*coco], int flab[*coco], int color[*coco], int color_true[*coco], int *hash_val1, int *hash_val2, int *hash_val3)
{
    graph g[MAXN*MAXM], canong[MAXN*MAXN];
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    int n,m,v,w;
    long zseed,colorhash;
/*    long hash_val; */

    n=*coco;

/*    printf("input value of n: %d (from C) \n",n); */

 /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
    Here we change those options that we want to be different from the
    defaults.  writeautoms=TRUE causes automorphisms to be written.     */

/*    options.writeautoms = TRUE; */
    options.writeautoms = FALSE;

    options.defaultptn= FALSE;
    
    options.getcanon= TRUE;

     /* The nauty parameter m is a value such that an array of
        m setwords is sufficient to hold n bits.  The type setword
        is defined in nauty.h.  The number of bits in a setword is
        WORDSIZE, which is 16, 32 or 64.  Here we calculate
        m = ceiling(n/WORDSIZE).                                  */

        m = SETWORDSNEEDED(n);

     /* The following optional call verifies that we are linking
        to compatible versions of the nauty routines.            */

        nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

     /* Now we create the cycle.  First we zero the graph, than for
        each v, we add the edge (v,v+1), where values are mod n. */
        
        EMPTYGRAPH(canong,m,n);
  
        EMPTYGRAPH(g,m,n);
        for (v = 0; v < n; ++v)  {
             for(w=0;w<n; ++w) {
		if(connect[w][v]==1) ADDONEEDGE(g,v,w,m);
              }
           }

/* filling ptn the color vetor*/
         for (v = 0; v < n; ++v)  {
      	    ptn[v]=color[v];
            lab[v]=flab[v];
/*            printf("%d -- %d \n",lab[v],ptn[v]);  */
           }
      
/*        printf("Generators for Aut(C[%d]) (from C):\n",n);  */

     /* we require canonical labelling, NULL pointer is replaced 
  	by the graph type canong*/

        densenauty(g,lab,ptn,orbits,&options,&stats,m,n,canong);

     /* The size of the group is returned in stats.grpsize1 and
        stats.grpsize2. */

/*        printf("Automorphism group size = ");
        writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
        printf("\n"); */

     /* print lab, which should be the canonical labelling of vertex in canongraph*/
/*        printf("should be canonical lab (from C) \n"); */
        for (v = 0; v < n; ++v) {
           flab[v] = lab[v];
/*           printf("%d,%d \n",v, flab[v]);  */
        }
/*        printf("canong (from C)\n");
        for (v=0; v<m*(size_t)n; ++v) {
           printf("%lu \n",canong[v]);
        }

        printf("colors in nauty format\n");
        for (v=0;v<n; ++v) {
           printf("%d \n", color[v]);
        }  */

/* Get the hash of color vector, call is: listhash(vector, nx, key), where
   nx is the element up to which to get the hash. Also, which key do we use?? */
/*        printf("colors in true format\n");
        for (v=0;v<n; ++v) {
           printf("%d \n", color_true[v]);
        } */

/*        printf("listhash\n"); */
        colorhash = listhash(color_true,n,489317L);
/*        printf("%lu \n",colorhash); */
/*        printf("listhash\n"); 
        printf("%lu \n",colorhash); */

/*        printf("hash1 \n");  */
/*       printf("%lu \n",zseed);  */


/* Call hashgraph 3 times, like in dreadnaut. Why?
   Each time it overwrites the "zseed" value anyway, so why... 
   Ans: They move around with 3 hash values just in case, I guess. No description
   is given in dreadnaut manual.
   In order to make final graph hash color dependent, add "colorhash" to the key.
   Might aswell just use the "colorhash" as the key itself, no? */

        zseed = hashgraph(canong,m,n,2922320L+colorhash); 
/*        zseed = hashgraph(canong,m,n,colorhash); */
/*        printf("%lu \n",zseed);  */
        *hash_val1 = zseed;  

        zseed = hashgraph(canong,m,n,19883109L+colorhash); 
/*        zseed = hashgraph(canong,m,n,colorhash); */
/*        printf(" %lu \n",zseed); */
        *hash_val2 = zseed;  

        zseed = hashgraph(canong,m,n,489317L+colorhash);  
/*        zseed = hashgraph(canong,m,n,colorhash); */
/*        printf(" %lu\n",zseed);  */
        *hash_val3 = zseed;  

 }


