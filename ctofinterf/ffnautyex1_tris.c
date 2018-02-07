/* This program prints generators for the automorphism group of an
   n-vertex polygon, where n is a number supplied by the user.

   This version uses a fixed limit for MAXN.
*/

#define MAXN 1000    /* Define this before including nauty.h */
#include "nauty.h"   /* which includes <stdio.h> and other system files */

void ffnautyex1_tris(int *coco, int connect[*coco][*coco], int flab[*coco], int color[*coco])
{
    graph g[MAXN*MAXM], canong[MAXN*MAXN];
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    int n,m,v,w;

    n=*coco;

    printf("input value of n: %d \n",n);

 /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
    Here we change those options that we want to be different from the
    defaults.  writeautoms=TRUE causes automorphisms to be written.     */

    options.writeautoms = TRUE;

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
            printf("%d -- %d \n",lab[v],ptn[v]);
           }
      
        printf("Generators for Aut(C[%d]):\n",n);

     /* we require canonical labelling, NULL pointer is replaced 
  	by the graph type canong*/

        densenauty(g,lab,ptn,orbits,&options,&stats,m,n,canong);

     /* The size of the group is returned in stats.grpsize1 and
        stats.grpsize2. */

        printf("Automorphism group size = ");
        writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
        printf("\n");

     /* print lab, which should be the canonical labelling of vertex in canongraph*/
        printf("cococo \n");
        for (v = 0; v < n; ++v) {
           printf("%d,%d \n",v, lab[v]);
        }
 }

