#include "mex.h"


#include <stdio.h>

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#define square(a) ((a)*(a))
#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define F0 stdout
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define PI 3.141592653


#define NS            5000
#define NBRANCH       (NS*2-2)
#define NCODE         64
#define NCATG         40

struct CommonInfo {
   char *z[2*NS-1], spname[NS][10], daafile[96], cleandata;
   int ns,ls,npatt,*fpatt,np,ntime,ncode,clock,rooted,model,icode;
   int seqtype, *pose, ncatG, npi0;
   double kappa, omega, alpha, pi[64],*rates, *lkl, daa[20*20], pi_sqrt[NCODE];
   double freqK[NCATG],rK[NCATG];
}  com;
struct TREEB {
   int nbranch, nnode, root, branches[NBRANCH][2];
}  tree;
struct TREEN {
   int father, nson, sons[NS], ibranch;
   double branch, divtime, omega, label, *lkl, daa[20*20];
}  *nodes;


int IncludeNodeLabel=0;


void BranchToNode (void)
{
/* tree.root need to be specified before calling this
*/
   int i, from, to;
   
   tree.nnode=tree.nbranch+1;

/*
   if (tree.root<0 || tree.root>com.ns*2-2) 
      { printf ("root at %d", tree.root+1); error2("tree root"); }
*/
   FOR (i,tree.nnode)
      { nodes[i].father=nodes[i].ibranch=-1;  nodes[i].nson=0; }
   FOR (i,tree.nbranch) {
      from=tree.branches[i][0];
      to  =tree.branches[i][1];
      nodes[from].sons[nodes[from].nson++]=to;
      nodes[to].father=from;
      nodes[to].ibranch=i;
   }
/*
   printf("\nNode\n%7s%7s%7s%7s%7s\n","father","node","branch","nson:","sons");
   FOR (i, tree.nnode) {
      printf ("\n%7d%7d%7d:%6d  ",
         nodes[i].father+1, i+1, nodes[i].ibranch+1, nodes[i].nson);
      FOR(j, nodes[i].nson) printf (" %2d", nodes[i].sons[j]+1);
   }
*/
}


int MakeTreeIb (int ns, int Ib[], int rooted)
{
/* construct tree from Ib[] using the algorithm of adding species
   Ib[k] marks the branch to which the (k+3)_th species (or the root) 
   is added.  Ib[k] should be in the range [0,k+3]
*/
   int i,j,k, is,it;

   tree.nbranch=3;
   for (i=0; i<3; i++)  { tree.branches[i][0]=3;  tree.branches[i][1]=i; }
   for (k=0; k<ns-3; k++) {
      is=k+3;       /* add species (is) to branch Ib[k] */

      for (i=0; i<tree.nbranch; i++)  for (j=0; j<2; j++)
         if (tree.branches[i][j]>=is) tree.branches[i][j]+=2;
      it=tree.branches[Ib[k]][1];
      tree.branches[Ib[k]][1]=is+1;
      tree.branches[tree.nbranch][0]=is+1;
      tree.branches[tree.nbranch++][1]=it;
      tree.branches[tree.nbranch][0]=is+1;
      tree.branches[tree.nbranch++][1]=is;
   }
   tree.root=tree.branches[0][0];
   BranchToNode ();
   
   if (rooted) {
      it=tree.branches[Ib[k]][0];
      tree.branches[Ib[k]][0]=tree.branches[tree.nbranch][0]=2*ns-2;
      tree.branches[tree.nbranch][1]=it;
      for (; it!=tree.root;  it=nodes[it].father) {
         tree.branches[nodes[it].ibranch][0]=it;
         tree.branches[nodes[it].ibranch][1]=nodes[it].father;
      }
      tree.root=2*ns-2;  tree.nbranch++;
      BranchToNode ();
   }
   return (0);
}


int OutSubTreeN (FILE *fout, int inode, int spnames, int branchlen)
{
   int i,ison;

   fputc ('(', fout);
   for(i=0; i<nodes[inode].nson; i++) {
      ison=nodes[inode].sons[i];
      if(nodes[ison].nson==0) {
         if(spnames) {
            if(IncludeNodeLabel) fprintf(fout, "%d_",ison+1);
/*            fprintf(fout,"%s",com.spname[ison]); */
         }
         else 
            fprintf(fout,"%d",ison+1);
      }
      else
         OutSubTreeN(fout,ison,spnames,branchlen);
      if(IncludeNodeLabel && nodes[ison].nson) fprintf(fout," %d ",ison+1);
      if(branchlen) { 

/*  Add branch labels to be read by Rod Page's TreeView. */
#ifdef CODEML
/*         if(com.verbose>1 && com.seqtype==1 && com.model && com.model!=7 && !com.NSsites) {
            fprintf(fout," #%.4f ", nodes[ison].omega);
*/         }
#endif

/*
if(IncludeNodeLabel && ison>=com.ns) fprintf(fout, " #%8.4f",nodes[ison].label);
*/
         fprintf(fout,": %.5f", nodes[ison].branch);
         /* fprintf(fout,": %.8e", nodes[ison].branch); */
      }
      if(i<nodes[inode].nson-1) fprintf(fout,", ");
   }
   fputc (')', fout);

   /* output node times for treeview */
/*
   if(noisy>=9 && com.ns>9 && com.clock && branchlen)
      fprintf(fout," #%.1f ", nodes[inode].divtime);

   if(noisy>=9 && com.verbose && com.ns>9)
      fprintf(fout," #%d ", inode+1);
*/
   return (0);
}



int OutaTreeN (FILE *fout, int spnames, int branchlen)
{
/* IncludeNodeLabel=1 will label the nodes
*/
   OutSubTreeN(fout,tree.root,spnames,branchlen);  
   if(IncludeNodeLabel) fprintf(fout," %d ", tree.root+1);
   if(branchlen && nodes[tree.root].branch>0) 
      fprintf(fout,": %.5f", nodes[tree.root].branch);

   fputc(';',fout);
   return(0);
}






int ListTrees (FILE* fout, int ns, int rooted)
{
/* list trees by adding species, works fine with large ns
*/
   int NTrees, NTreeRoot=3;
   int i, Ib[NS-2], ns1=ns+rooted, nM=ns1-3,finish;

   if(ns<=12) {
      printf ("%20s%20s%20s\n", "Taxa", "Unrooted trees", "Rooted trees");
      for (i=4,NTrees=1; i<=ns; i++)  
         printf ("%20d%20d%20d\n", i, (NTrees*=2*i-5), (NTreeRoot*=2*i-3));
      fprintf (fout, "%10d %10d\n", ns, (!rooted?NTrees:NTreeRoot));
   }
   for (i=0;i<nM;i++) Ib[i]=0;
   for (NTrees=0; ; ) {
      MakeTreeIb(ns, Ib, rooted);
      OutaTreeN(fout, 0, 0);
      fprintf(fout, " %7d\n", NTrees++);

      for (i=nM-1,Ib[nM-1]++,finish=0; i>=0; i--) {
         if (Ib[i]<2*i+3) break;
         if (i==0) { finish=1; break; }
         Ib[i]=0; Ib[i-1]++; 
      }
      if (finish) break;
   }
   FPN(fout);
   return (0);
}



/*
int main(int argc, char *argv[])
{
	FILE *fp;
	int rooted;
	int ns; int i;

	printf("Hello, world\n");

    fp=(FILE*)fopen("filename", "w");
   if(fp==NULL) {
      printf("\nerror when opening file %s\n", "filename");
      exit(-1);
   }
   rooted=1;
   ns=5;

   i=(ns*2-1)*sizeof(struct TREEN);
   if((nodes=(struct TREEN*)malloc(i))==NULL) exit(1);

    ListTrees(fp,ns,rooted);
    FPN(fp);
	return 0;
}
*/

void gotrees(double ns, double rooted)
{
FILE *fp;
int i;
    fp=(FILE*)fopen("outfile", "w");
   if(fp==NULL) {
      printf("\nerror when opening file %s\n", "outfile");
      exit(-1);
   }

   i=(ns*2-1)*sizeof(struct TREEN);
   if((nodes=(struct TREEN*)malloc(i))==NULL) exit(1);
    ListTrees(fp,ns,rooted);
     fclose(fp);
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double x;
  double y;
  char *input_buf;
  int   buflen,status;

  /* Create a 1-by-1 matrix for the return argument. */
  /* plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); */

  /* Get the scalar value of the input x. */
  /* Note: mxGetScalar returns a value, not a pointer. */
  if (nrhs !=2) {
    mexErrMsgTxt("Two input argument required.");
    mexErrMsgTxt("listtree(ns,rooted).");
  }

  x = mxGetScalar(prhs[0]);
  y = mxGetScalar(prhs[1]);
  if (x<3) {
    mexErrMsgTxt("ns<3, no need to do this?");
  }
  if (x>10) {
    mexErrMsgTxt("ns>10, too big output?");
  }

  /* Call the timestwo_alt subroutine. */
  gotrees(x,y);
}

