#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "lsd.h"

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

#define STACK_MAX 10000
/*----------------------------------------------------------------------------*/
/** Fatal error, print a message to standard-error output and exit.
 */

double *gSegments;

static void error(char * msg)
{
  fprintf(stderr,"%s\n",msg);
  exit(EXIT_FAILURE);
}


typedef struct {
  int *array;
  size_t used;
  size_t size;
} Array;

void initArray(Array *a, size_t initialSize) {
  a->array = (int *)malloc(initialSize * sizeof(int));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Array *a, int element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (int *)realloc(a->array, a->size * sizeof(int));
  }
  a->array[a->used++] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

struct Stack {
    int     data[STACK_MAX];
    int     size;
};

typedef struct Stack Stack;


void Stack_Init(Stack *S)
{
    S->size = 0;
}

int Stack_Top(Stack *S)
{
    if (S->size == 0) {
        printf("Error: stack empty\n");
        return -1;
    } 

    return S->data[S->size-1];
}

void Stack_Push(Stack *S, int d)
{
    if (S->size < STACK_MAX)
        S->data[S->size++] = d;
    else
        printf("Error: stack full\n");
}

void Stack_Pop(Stack *S)
{
    if (S->size == 0)
        printf("Error: stack empty\n");
    else
        S->size--;
}

/*----------------------------------------------------------------------------*/
/** Read a ASCII number from a PGM file.
 */
static int get_num(FILE * f)
{
  int num,c;

  while(isspace(c=getc(f)));
  if(!isdigit(c)) error("Error: corrupted PGM file.");
  num = c - '0';
  while( isdigit(c=getc(f)) ) num = 10 * num + c - '0';
  if( c != EOF && ungetc(c,f) == EOF )
    error("Error: unable to 'ungetc' while reading PGM file.");

  return num;
}

/*----------------------------------------------------------------------------*/
/** Skip white characters and comments in a PGM file.
 */
static void skip_whites_and_comments(FILE * f)
{
  int c;
  do
    {
      while(isspace(c=getc(f))); /* skip spaces */
      if(c=='#') /* skip comments */
        while( c!='\n' && c!='\r' && c!=EOF )
          c=getc(f);
    }
  while( c == '#' || isspace(c) );
  if( c != EOF && ungetc(c,f) == EOF )
    error("Error: unable to 'ungetc' while reading PGM file.");
}

static double * read_pgm_image_double(int * X, int * Y, char * name)
{
  FILE * f;
  int c,bin;
  int xsize,ysize,depth,x,y;
  double * image;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdin;
  else f = fopen(name,"rb");
  if( f == NULL ) error("Error: unable to open input image file.");

  /* read header */
  if( getc(f) != 'P' ) error("Error: not a PGM file!");
  if( (c=getc(f)) == '2' ) bin = FALSE;
  else if( c == '5' ) bin = TRUE;
  else error("Error: not a PGM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);            /* X size */
  if(xsize<=0) error("Error: X size <=0, invalid PGM file\n");
  skip_whites_and_comments(f);
  ysize = get_num(f);            /* Y size */
  if(ysize<=0) error("Error: Y size <=0, invalid PGM file\n");
  skip_whites_and_comments(f);
  depth = get_num(f);            /* depth */
  if(depth<=0) fprintf(stderr,"Warning: depth<=0, probably invalid PGM file\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("Error: corrupted PGM file.");

  /* get memory */
  image = (double *) calloc( (size_t) (xsize*ysize), sizeof(double) );
  if( image == NULL ) error("Error: not enough memory.");

  /* read data */
  for(y=0;y<ysize;y++)
    for(x=0;x<xsize;x++)
      image[ x + y * xsize ] = bin ? (double) getc(f)
                                   : (double) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
    error("Error: unable to close file while reading PGM file.");

  /* return image */
  *X = xsize;
  *Y = ysize;
  return image;
}

static int dist(double x1, double y1, double x2, double y2)
{
  //printf("%f %f %f %f \n", x1, y1, x2, y2);
  return (int)sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}

static void write_eps( double * segs, int n, int dim,
                       char * filename, int xsize, int ysize, double width )
{
  FILE * eps;
  int i;

  /* check input */
  if( segs == NULL || n < 0 || dim <= 0 )
    error("Error: invalid line segment list in write_eps.");
  if( xsize <= 0 || ysize <= 0 )
    error("Error: invalid image size in write_eps.");

  /* open file */
  eps = fopen(filename,"w");
  if( eps == NULL ) error("Error: unable to open EPS output file.");

  /* write EPS header */
  fprintf(eps,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(eps,"%%%%BoundingBox: 0 0 %d %d\n",xsize,ysize);
  fprintf(eps,"%%%%Creator: LSD, Line Segment Detector\n");
  fprintf(eps,"%%%%Title: (%s)\n",filename);
  fprintf(eps,"%%%%EndComments\n");

  printf("X size %d Y size %d n %d \n",xsize,ysize,n);

  /* write line segments */
  for(i=0;i<n;i++)
  {
    //if(i>=2) continue;
    int length = (int)segs[i*dim+6];
    if(length <= 0) continue;
    // printf( "newpath %f %f moveto %f %f lineto %f setlinewidth stroke size %f \n",
    //          segs[i*dim+0],
    //          (double) ysize - segs[i*dim+1],
    //          segs[i*dim+2],
    //          (double) ysize - segs[i*dim+3],
    //          width <= 0.0 ? segs[i*dim+4] : width, segs[i*dim+6]);
    fprintf( eps,"newpath %f %f moveto %f %f lineto %f setlinewidth stroke\n",
             segs[i*dim+0],
             (double) ysize - segs[i*dim+1],
             segs[i*dim+2],
             (double) ysize - segs[i*dim+3],
             width <= 0.0 ? segs[i*dim+4] : width);
  }


  /* close EPS file */
  fprintf(eps,"showpage\n");
  fprintf(eps,"%%%%EOF\n");
  if( eps != stdout && fclose(eps) == EOF )
    error("Error: unable to close file while writing EPS file.");
}

static void write_eps_BS( Array * val, double * segs, int n, int dim,
                       char * filename, int xsize, int ysize, double width , int belongsto)
{
  FILE * eps;
  int i;

  /* check input */
  if( segs == NULL || n < 0 || dim <= 0 )
    error("Error: invalid line segment list in write_eps.");
  if( xsize <= 0 || ysize <= 0 )
    error("Error: invalid image size in write_eps.");

  if(val->used == 1) return;
  /* open file */
  eps = fopen(filename,"w");
  if( eps == NULL ) error("Error: unable to open EPS output file.");

  /* write EPS header */
  fprintf(eps,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(eps,"%%%%BoundingBox: 0 0 %d %d\n",xsize,ysize);
  fprintf(eps,"%%%%Creator: LSD, Line Segment Detector\n");
  fprintf(eps,"%%%%Title: (%s)\n",filename);
  fprintf(eps,"%%%%EndComments\n");

  //printf("X size %d Y size %d n %d parallel segments found %zu \n",xsize,ysize,n, val->used);

  /* write line segments */
  for(i=0;i<val->used;i++)
  {
    int length = (int)segs[val->array[i]*dim+6];
    double localWidth = width;
    
    if(length <= 0) continue;//last line of defence
    int v = val->array[i];
    if(v==belongsto) localWidth = 6.0; //darken the index line
    //printf("%d \n",v);
    //if(v == 0 ||v == 1 || v== 4 || v==5 || v==154 || v== 751 || v== 888) 
    //if(v==0 || v == 10061 || v==8009 || v==6301 || v==4693 || v==2894 || v==751)
    //if(v==0 || v==2894)
    {
      // printf( "newpath %f %f moveto %f %f lineto %f setlinewidth stroke length %f line num %d \n",
      //          segs[val->array[i]*dim+0],
      //          (double) ysize - segs[val->array[i]*dim+1],
      //          segs[val->array[i]*dim+2],
      //          (double) ysize - segs[val->array[i]*dim+3],
      //          localWidth <= 0.0 ? segs[val->array[i]*dim+4] : localWidth, segs[val->array[i]*dim+6],i);
      fprintf( eps,"newpath %f %f moveto %f %f lineto %f setlinewidth stroke\n",
               segs[val->array[i]*dim+0],
               (double) ysize - segs[val->array[i]*dim+1],
               segs[val->array[i]*dim+2],
               (double) ysize - segs[val->array[i]*dim+3],
               localWidth <= 0.0 ? segs[val->array[i]*dim+4] : localWidth);
    }
  }


  /* close EPS file */
  fprintf(eps,"showpage\n");
  fprintf(eps,"%%%%EOF\n");
  if( eps != stdout && fclose(eps) == EOF )
    error("Error: unable to close file while writing EPS file.");
}

static void write_eps_BS1(double segs [][5], int n, int dim,
                       char * filename, int xsize, int ysize )
{
  FILE * eps;
  int i;

  /* check input */
  if( segs == NULL || n < 0 || dim <= 0 )
    error("Error: invalid line segment list in write_eps.");
  if( xsize <= 0 || ysize <= 0 )
    error("Error: invalid image size in write_eps.");

  /* open file */
  eps = fopen(filename,"w");
  if( eps == NULL ) error("Error: unable to open EPS output file.");

  /* write EPS header */
  fprintf(eps,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(eps,"%%%%BoundingBox: 0 0 %d %d\n",xsize,ysize);
  fprintf(eps,"%%%%Creator: LSD, Line Segment Detector\n");
  fprintf(eps,"%%%%Title: (%s)\n",filename);
  fprintf(eps,"%%%%EndComments\n");

  //printf("X size %d Y size %d n %d parallel segments found %zu \n",xsize,ysize,n, val->used);

  /* write line segments */
  for(i=0;i<n;i++)
  {
    // int length = (int)segs[val->array[i]*dim+6];
     double localWidth = 2.0;
    if(segs[i][4] == 1) continue;
    // if(length <= 0) continue;//last line of defence
    // int v = val->array[i];
    // if(v==belongsto) localWidth = 6.0; //darken the index line
    //printf("%d \n",v);
    //if(v == 0 ||v == 1 || v== 4 || v==5 || v==154 || v== 751 || v== 888) 
    //if(v==0 || v == 10061 || v==8009 || v==6301 || v==4693 || v==2894 || v==751)
    //if(v==0 || v==2894)
    {
      // printf( "newpath %f %f moveto %f %f lineto %f setlinewidth stroke length %f line num %d \n",
      //          segs[val->array[i]*dim+0],
      //          (double) ysize - segs[val->array[i]*dim+1],
      //          segs[val->array[i]*dim+2],
      //          (double) ysize - segs[val->array[i]*dim+3],
      //          localWidth <= 0.0 ? segs[val->array[i]*dim+4] : localWidth, segs[val->array[i]*dim+6],i);
      fprintf( eps,"newpath %f %f moveto %f %f lineto %f setlinewidth stroke\n",
               segs[i][0],
               (double) ysize - segs[i][1],
               segs[i][2],
               (double) ysize - segs[i][3],
               localWidth);
    }
  }


  /* close EPS file */
  fprintf(eps,"showpage\n");
  fprintf(eps,"%%%%EOF\n");
  if( eps != stdout && fclose(eps) == EOF )
    error("Error: unable to close file while writing EPS file.");
}

int cmpfunc (const void * a, const void * b) //descending order sort
{
  return (int)( *(double *)(b + sizeof(double) * (6)) - *(double*)(a + sizeof(double) * (6)) );
}

int cmpfunc1 (const void * a, const void * b) //ascending order sort
{
  return (int)( *(double *)(a)  - *(double*)(b) );
}

//To compare parallel to  x axis
int cmpfunc2 (const void * a, const void * b) //ascending order sort
{
  int line1 = *((int *)a);
  int line2 = *((int *)b);
  //We do this since segment coordinates may not be in order
  double yStart1 = gSegments[line1*7+0] >= gSegments[line1*7+2] ? gSegments[line1*7+2] : gSegments[line1*7+0];
  double yStart2 = gSegments[line2*7+0] >= gSegments[line2*7+2] ? gSegments[line2*7+2] : gSegments[line2*7+0];
  return  (int)(yStart1 - yStart2);
}

//To compare parallel to  y axis
int cmpfunc3 (const void * a, const void * b) //ascending order sort
{
  int line1 = *((int *)a);
  int line2 = *((int *)b);
  //We do this since segment coordinates may not be in order
  double yStart1 = gSegments[line1*7+1] >= gSegments[line1*7+3] ? gSegments[line1*7+3] : gSegments[line1*7+1];
  double yStart2 = gSegments[line2*7+1] >= gSegments[line2*7+3] ? gSegments[line2*7+3] : gSegments[line2*7+3];
  return  (int)(yStart1 - yStart2);
}

void getPerpendicular2D(double x1, double y1, double *x2, double *y2){
  *x2 = -y1;
  *y2 = x1;
}

//
//before calling this function filter based on length of lines greater than a threshold
double findDistance(double ls1[], double ls2[], double xsize, double ysize)
{
    //Begin and ending points of First line segment .
    double X1 = ls1[0];
    //double Y1 = (double) ysize - ls1[1];
    double Y1 = ls1[1];
    double X2 = ls1[2];
    //double Y2 = (double) ysize - ls1[3];
    double Y2 = ls1[3];
    //Ending point of line segment
    double X3 = ls2[2];
    double Y4 = ls2[3];
    /*double Y4 = (double) ysize -ls2[3];*/

    double vecx = (X2-X1);
    double vecy = (Y2-Y1);
    double lineLength = ls1[6];//this is the length of that line

    //printf("startx %f starty %f endx %f endy %f lengt [%f] \n",X1,Y1,X2,Y2,lineLength);
  
    double normx, normy;
    getPerpendicular2D(vecx/lineLength, vecy/lineLength, &normx, &normy);
    
    //printf("normx [%f] normy [%f] \n",normx, normy);
    double dotProduct =((X3-X1) *normx +(Y4-Y1)*normy);

    return fabs(dotProduct);
}

//computer distance array n * n for each line segment O(n^2)

void computeDistances(double **distArray, double *out, int n, double X, double Y)
{
  for(int idx=0; idx<n; idx++)
  {
    //printf("outer [%d] length is [%f] ",idx, out[idx*7+6]);
    for(int idx1=0; idx1<n; idx1++)
    {
      if(idx == idx1)
        distArray[idx][idx1] = 0;//distance between a point to same point is same
      else
      {
        if(distArray[idx][idx1] == 0) 
        {
          distArray[idx][idx1] = findDistance(out+(idx*7), out+(idx1*7), X, Y );
          distArray[idx1][idx] = distArray[idx][idx1];
        }
      }
    }
  }
}

double retAngle(double s_x1, double s_y1, double e_x2, double e_y2, double e_x3, double e_y3)// start is same for both segments
{
  double dx21 = e_x2 - s_x1;
  double dx31 = e_x3 - s_x1;
  double dy21 = e_y2 - s_y1;
  double dy31 = e_y3 - s_y1;
  double m12 = sqrt(dx21*dx21 + dy21*dy21);
  double m13 = sqrt(dx31*dx31 + dy31*dy31);
  double rad = acos((dx21*dx31 + dy21*dy31) / (m12*m13));
  double theta = rad * 180.0 / M_PI; // convert radians to angle
  return theta;
}

double calculateTheta(double * segs, double slope, int cur, int w, int h)
{
  double theta; 
  if(slope >= 0)//then we calculate parallel to X axis the normal line and take angle between them two.
    theta = retAngle(segs[cur*7+0], segs[cur*7+1], segs[cur*7+2], segs[cur*7+3], w - segs[cur*7+0], segs[cur*7+1] );
  else //then we calculate parallel to Y axis the normal line and take angle between them two.
    theta = retAngle(segs[cur*7+0], segs[cur*7+1], segs[cur*7+2], segs[cur*7+3], segs[cur*7+0] , 0);
  return theta;
}

//When we have overlap distance, eucledian distance and slopes calculated , we can classify using these criteria
//All segments that have overlap -1 means they are parallel. We are trying to reduce false positives by 
//checking slope and distance between lines to group them together.

void classifySegments(double ** dist, double *slopeOfLines, double **overlap, int n, Array ** parallelLines, Array ** singleLines, int T, int lengthT, double *out)
{
  // FILE * classified;
  // classified = fopen("classified_segments.txt","w");

  *parallelLines = (Array *)(malloc(sizeof( Array )* n));
  *singleLines = (Array *)(malloc(sizeof( Array )* n));

  Array * parallelLineEntries = *parallelLines;
  Array * singleLineEntries = *singleLines;
  for(int outer=0; outer<n; outer++) //points to frame - line to compare with
  {
    //printf("start %d \n",outer);
    Array * singleLineEntry = (Array *)malloc(sizeof(Array ) );
    initArray(singleLineEntry, 50);//first entry stores the line itself
    insertArray(singleLineEntry, outer);
    Array * parallelLineEntry = (Array *)malloc(sizeof(Array ) );
    initArray(parallelLineEntry, 50);//first entry stores the line itself
    insertArray(parallelLineEntry, outer);
    

    for(int inner=0; inner<n; inner++) // points to ref - line who is compared with frame
    {
      if(inner == outer) continue;
      //if(out[inner*7+6] < lengthT) continue;
      //if(outer !=0 || ) break;
      //im taking all int since it makes it easier to compare. anyways all double data is rounded usinf fabs so we wont lose any data here.
      int slopeF = slopeOfLines[outer]; // when slope is 0 it need not mean parallel to x. since i rounded it off.
      int slopeR = slopeOfLines[inner];
      int distFromFrameToRef = (int)dist[outer][inner];
      int overlapDistFromFrameToRef = (int)overlap[outer][inner];

      int t=5; // slopeThreshold - do we need this ?
      int t2=5; // distThreshold1
      int t3=10; // distThreshold2
      int t4=50;//overlap distance threshold
      int t5=5; //upto 5 degree variation
      //those whose overlap = -1 , can be considered they are not parallel to each other . so that condition can be used
      //rather than t1 threshold check.

      //fprintf(classified,"Classify segments pair(%d,%d) SLOPE DIFF= [%d], DIST= [%d], OVERLAP= [%d],  \n",outer,inner,abs(slopeR-slopeF),distFromFrameToRef,overlapDistFromFrameToRef);

      //if ((absolute difference in slope(LS1) and slope(LS2)) < t1 ) && (distance between them less than < t2) && (overlap between them == 0) 
      // then it can be classified as one single line segment
      if(overlap[outer][inner] != -1 && abs(slopeF-slopeR) <= t5  && distFromFrameToRef < t2 && overlapDistFromFrameToRef <t4 /*&& !(out[inner*7+6] < lengthT)*/)
      //if(overlap[outer][inner] != -1 && abs(slopeF-slopeR) <= t5  && distFromFrameToRef < t2 && overlapDistFromFrameToRef <t4 && !(out[inner*7+6] < lengthT))
        insertArray(singleLineEntry, inner);

      //if ((absolute difference in slope(LS1) and slope(LS2)) < t1 ) && (t2 <= distance between them less than > t3) && (overlap between them > 0)
      // then it can be classified as parallel line segment .
      //else
        //if(overlap[outer][inner] != -1 && abs(slopeF-slopeR) <= t5 && (distFromFrameToRef >t3 && distFromFrameToRef < T) /*&& overlapDistFromFrameToRef > 0*/)
      if(/*overlap[outer][inner] != -1 &&*/ abs(slopeF-slopeR) <= t5 && !(out[inner*7+6] < lengthT) /*&& distFromFrameToRef < t4 */)
        insertArray(parallelLineEntry, inner);
    }
    parallelLineEntries[outer] = *parallelLineEntry;
    singleLineEntries[outer] = *singleLineEntry;
  }
  //printf("done \n");
}

void findOverLapDistance(double ** dist, Array * innerLines, double * slopeOfLines, double * segs, int cur, double w, double h, int num_segs)
{

  //double mFrame = slopeOfLines[cur];
  double thetaFrame =calculateTheta(segs, slopeOfLines[cur], cur, w, h); 
  for(int l = 0; l < num_segs; l++)
    dist[cur][l] = -1;//-1 means since they were not parallel in either axis , hence we dint calculate overlap .

  for(int outer=0;/*outer<innerLines->used*/ outer<num_segs;outer++) //we are looping through all lines that are parallel to cur index line
  {
    //outer index is just to iterate through the parallel lines array, dont use it for any thing else .

    int lineNum = /*innerLines->array[outer]*/ outer; //lineNum is nothing but array index num that was stored to look up in segs n*5 d array
    
    double thetaRef =calculateTheta(segs, slopeOfLines[lineNum], lineNum, w, h);
    //assume both lines project on X axis as 4 points X1, X2, X3, X4
    //assume both lines project on Y axis as 4 points Y1, Y2, Y3, Y4

    //Lets use new logic as per prof

    //unfortunately all line segments given by LSD are not in the format of start and end so we can find the start segment by choosing the one that has least X

    double d= -1.0;
    if(outer != cur)
    {
    //points 2894 (5017.796397,4508.583200) (5005.967515,5005.967515)
      //For X pair(0,751) , 0, 1329, 555, 1273, d=807 
      //For X pair(0,2894) , 0, 1329, 0, 402, d=7054 

      if(slopeOfLines[cur] >=0 && slopeOfLines[lineNum] >=0)
      {
        double frame_s_X1 = segs[cur*7+0] >= segs[cur*7+2] ? segs[cur*7+2] : segs[cur*7+0];//s means start
        double frame_e_X2 = segs[cur*7+0] >= segs[cur*7+2] ? segs[cur*7+0] : segs[cur*7+2];//e means end
        double ref_s_X3 = segs[lineNum*7+0] >= segs[lineNum*7+2] ? segs[lineNum*7+2] : segs[lineNum*7+0];
        double ref_e_X4 = segs[lineNum*7+0] >= segs[lineNum*7+2] ? segs[lineNum*7+0] : segs[lineNum*7+2];
        // int frame_s_X1 = segs[cur*7+0];
        // int frame_e_X2 = segs[cur*7+2];
        // int ref_s_X3 = segs[lineNum*7+0];
        // int ref_e_X4 = segs[lineNum*7+2];
        if( frame_s_X1 >= ref_s_X3 && frame_e_X2 <= ref_e_X4) //it means frame is between and smaller than reference line seg
          d = fabs((frame_e_X2 - frame_s_X1)/cos(thetaFrame));
        else if( frame_s_X1 <= ref_s_X3 && frame_e_X2 >= ref_e_X4) //it means ref is between and smaller than frame line seg
          d = fabs((ref_e_X4 - ref_s_X3)/cos(thetaRef));
        //means they may some what overlap and don't at all- we need to consider 2 combinations of ref and frame arrangement
        else if(frame_s_X1 <= ref_e_X4)
          d = fabs((ref_e_X4 - frame_s_X1)/cos(thetaFrame));
        else if(frame_e_X2 >= ref_s_X3)
          d = fabs((ref_s_X3 - frame_e_X2)/cos(thetaRef));
        else
          d=0.0;
        //-1.0 means there is no overlap
        //printf("For X pair(%d,%d) , %d, %d, %d, %d, d=%d \n",cur,lineNum,frame_s_X1,frame_e_X2,ref_s_X3, ref_e_X4,d);
      }
      else if(slopeOfLines[cur] <0 && slopeOfLines[lineNum] <0)
      //these lines are oriented to Y axis , note we have not checked for slope infinity lines but thos lines may not appear
      {
        double frame_s_Y1 = segs[cur*7+1] >= segs[cur*7+3] ? segs[cur*7+3] : segs[cur*7+1];//s means start
        double frame_e_Y2 = segs[cur*7+1] >= segs[cur*7+3] ? segs[cur*7+1] : segs[cur*7+3];//e means end
        double ref_s_Y3 = segs[lineNum*7+1] >= segs[lineNum*7+3] ? segs[lineNum*7+3] : segs[lineNum*7+1];
        double ref_e_Y4 = segs[lineNum*7+1] >= segs[lineNum*7+3] ? segs[lineNum*7+1] : segs[lineNum*7+3];
        //need to take into account where segment start or end may be equal.
        if( frame_s_Y1 >= ref_s_Y3 && frame_e_Y2 <= ref_e_Y4) //it means frame is between and smaller than reference line seg
          d = fabs((frame_e_Y2 - frame_s_Y1)/sin(thetaFrame));
        else if( frame_s_Y1 <= ref_s_Y3 && frame_e_Y2 >= ref_e_Y4) //it means ref is between and smaller than frame line seg
          d = fabs((ref_e_Y4 - ref_s_Y3)/sin(thetaRef));
        //means they may some what overlap and don't at all - we need to consider 2 combinations of ref and frame arrangement
        else if(frame_s_Y1 <= ref_e_Y4)
          d = fabs((ref_e_Y4 - frame_s_Y1)/sin(thetaFrame));
        else if(frame_e_Y2 >= ref_s_Y3)
          d = fabs((ref_s_Y3 - frame_e_Y2)/sin(thetaRef));
        else
          d=0.0;
        //-1.0 means there is no overlap
        //printf("For Y pair(%d,%d) %d, %d, %d, %d, d=%d \n",cur,lineNum,frame_s_Y1,frame_e_Y2,ref_s_Y3, ref_e_Y4,d);
      }
    }
    dist[cur][lineNum] = d;
  }
}

//returns parallel lines and also calculates slope of lines
void findParallelLines(Array ** parallelLines , double * segs, int n, int X, int Y, int dim, double * slopeOfLines, int bins)
{
  printf("Enter findParallelLines bins [%d]\n",bins);
  *parallelLines = (Array *)(malloc(sizeof(Array )* n));
  
  Array * innerLines = *parallelLines;
  for(int outer = 0; outer <n ; outer++)
  {
    double X1 = segs[outer*dim+0];
    //double Y1 = (double) Y - segs[outer*dim+1];
    double Y1 = segs[outer*dim+1];
    double X2 = segs[outer*dim+2];
    //double Y2 = (double) Y - segs[outer*dim+3];
    double Y2 = segs[outer*dim+3];
    slopeOfLines[outer] = (Y2 -Y1) / (X2-X1); 

    double thetaLine = 180 * atan2(Y2-Y1,X2-X1)/M_PI + 90;
    int iTheta = thetaLine; 
    if(iTheta >= 180) iTheta = (iTheta + 180)%360;
    slopeOfLines[outer] = (iTheta/5.0 + 0.5) * bins;

    //Min slope is [-11650.748350] and max slope is [6320.331073]
    Array * parallels = (Array *)malloc(sizeof(Array));
    initArray(parallels, 50);
    insertArray(parallels, outer);
    for(int inner=0; inner <n ; inner++)
    {
      if(inner == outer) continue;

      double X3 = segs[inner*dim+0];
      double Y3 = (double)Y - segs[inner*dim+1];
      //double Y3 = segs[inner*dim+1];
      double X4 = segs[inner*dim+2];
      double Y4 = (double)Y - segs[inner*dim+3];
      //double Y4 = segs[inner*dim+3];
      
      double thetaLine1 = 180 * atan2(Y4-Y4,X3-X3)/M_PI + 90;
      int iTheta1 = thetaLine1; 
      if(iTheta1 >= 180) iTheta1 = (iTheta1 + 180)%360;
      double slope1 = (iTheta1/5.0 + 0.5) * bins;

      //if(fmax(X1,X2) < fmin(X3,X4)) continue; // there is no mutual abcissa so continue

      //int abcissa_1 = (int)(Y1-Y2)/(X1-X2);
      //int abcissa_2 = (int)(Y3-Y4)/(X3-X4);

      if(fabs(slopeOfLines[outer]-slope1) <= 5 ) //it means its parallel keeping 5 threshold
      {
        insertArray(parallels, inner);
      }
    }
    innerLines[outer] = *parallels;
  }

  //binning
  // qsort(slopeOfLines, n , sizeof(double), cmpfunc1);
  // double min = slopeOfLines[0];
  // double max = slopeOfLines[n-1];
  // int binSize = (int)(fabs(min) + fabs(max))/bins; // each bin can contain line segments within this slope range
  // printf("BIN SIZE [%d] [%f] [%f]\n",binSize,min,max);
  // for(int idx=0 ;idx <n; idx++)//Loop through all slopes O(n*binSize)
  // {
  //   //We need to bin the slopes to some bucket.
  //   double slope= slopeOfLines[idx];
    
  //   double bestBin ;
  //   int isFound = 0;
  //   double start=min;
  //   //Improvize using binary search O(n*log(binSize))
  //   for( ; start <= max; start += binSize)
  //   {
  //     if(slope > start){
  //       isFound = 1;
  //       bestBin = start;
  //     }
  //     else{
  //       break;
  //     }
  //   }
  //   if(isFound == 1)
  //     slopeOfLines[idx] = bestBin;
  //   else
  //     slopeOfLines[idx] = start;//it could be in beginning or end then worst case.
  // }

}

void merge_lines(double merged_lines[][5] , int cur, int n , double * segments, int *lines, double *slopeOfLines, int w, int h)
{
  int set = 0;
  int t4=50;//overlap distance threshold
  int t5=5; //upto 5 degree variation
  int t2=5; // distThreshold1

  for(int ni=0; ni < n; ni++)
  {
    int idx = lines[ni]; // get line index
    if(set == 1) //It means already have a line segment assigned
    {
      if(idx < cur)// it means we may have already merged idx line segment
      {
        if(merged_lines[idx][4] == 0) //check if its valid
        {
          //First lets get its perpendicular distance and overlap distance
           double thetaMerged =calculateTheta(segments, slopeOfLines[idx], idx, w, h);
           double thetaCur =calculateTheta(segments, slopeOfLines[cur], cur, w, h);
           double overlap = -1;
            if(slopeOfLines[cur] >=0 && slopeOfLines[idx] >=0)
            {
              double merged_s_X1 = segments[idx*7+0] >= segments[idx*7+2] ? segments[idx*7+2] : segments[idx*7+0];
              double merged_e_X2 = segments[idx*7+0] >= segments[idx*7+2] ? segments[idx*7+0] : segments[idx*7+2];
              double cur_s_X3 = segments[cur*7+0] >= segments[cur*7+2] ? segments[cur*7+2] : segments[cur*7+0];//s means start
              double cur_e_X4 = segments[cur*7+0] >= segments[cur*7+2] ? segments[cur*7+0] : segments[cur*7+2];//e means end

              if( merged_s_X1 >= cur_s_X3 && merged_e_X2 <= cur_e_X4) //it means frame is between and smaller than reference line seg
                overlap = fabs((merged_e_X2 - merged_s_X1)/cos(thetaMerged));
              else if( merged_s_X1 <= cur_s_X3 && merged_e_X2 >= cur_e_X4) //it means ref is between and smaller than frame line seg
                overlap = fabs((cur_e_X4 - cur_s_X3)/cos(thetaCur));
              //means they may some what overlap and don't at all- we need to consider 2 combinations of ref and frame arrangement
              else if(merged_s_X1 <= cur_e_X4)
                overlap = fabs((cur_e_X4 - merged_s_X1)/cos(thetaMerged));
              else if(merged_e_X2 >= cur_s_X3)
                overlap = fabs((cur_s_X3 - merged_e_X2)/cos(thetaCur));
              else
                overlap=0.0;
            }
            else if(slopeOfLines[cur] <0 && slopeOfLines[idx] <0)
            //these lines are oriented to Y axis , note we have not checked for slope infinity lines but thos lines may not appear
            {
              double merged_s_Y1 = segments[idx*7+1] >= segments[idx*7+3] ? segments[idx*7+3] : segments[idx*7+1];
              double merged_e_Y2 = segments[idx*7+1] >= segments[idx*7+3] ? segments[idx*7+1] : segments[idx*7+3];
              double cur_s_Y3 = segments[cur*7+1] >= segments[cur*7+3] ? segments[cur*7+3] : segments[cur*7+1];//s means start
              double cur_e_Y4 = segments[cur*7+1] >= segments[cur*7+3] ? segments[cur*7+1] : segments[cur*7+3];//e means end
              //need to take into account where segment start or end may be equal.
              if( merged_s_Y1 >= cur_s_Y3 && merged_e_Y2 <= cur_e_Y4) //it means frame is between and smaller than reference line seg
                overlap = fabs((merged_e_Y2 - merged_s_Y1)/sin(thetaMerged));
              else if( merged_s_Y1 <= cur_s_Y3 && merged_e_Y2 >= cur_e_Y4) //it means ref is between and smaller than frame line seg
                overlap = fabs((cur_e_Y4 - cur_s_Y3)/sin(thetaCur));
              //means they may some what overlap and don't at all - we need to consider 2 combinations of ref and frame arrangement
              else if(merged_s_Y1 <= cur_e_Y4)
                overlap = fabs((cur_e_Y4 - merged_s_Y1)/sin(thetaMerged));
              else if(merged_e_Y2 >= cur_s_Y3)
                overlap = fabs((cur_s_Y3 - merged_e_Y2)/sin(thetaCur));
              else
                overlap=0.0;
            }

            //Now lets find perpendicular distance
            double perpend_dist = findDistance(&segments[cur*7+0], &segments[idx*7+0] , w, h);
            if(overlap != -1 && fabs(slopeOfLines[cur]-slopeOfLines[idx]) <= t5  && perpend_dist < t2 && overlap <t4 )
            { // it means the other merged linesegment is overshadowed by this so make it invalid
              if((merged_lines[cur][1] < segments[idx*7+1]) && (merged_lines[cur][3] < segments[idx*7+3])  )
              {//copy all the coordinates
                //lets copy only the end coordinates
                merged_lines[cur][2] = segments[idx*7+2];
                merged_lines[cur][3] = segments[idx*7+3];   
              }

              merged_lines[idx][4]=1; // mark it invalid means not required
            }
        }
        set=1;
        continue;
      }
      else //its a new line so just take
      {
        //We come here to Check if its a valid candidate to merge
        if((merged_lines[cur][1] < segments[idx*7+1]) && (merged_lines[cur][3] < segments[idx*7+3])  )
        {
          //lets copy only the end coordinates
          merged_lines[cur][2] = segments[idx*7+2];
          merged_lines[cur][3] = segments[idx*7+3];   
        }
      }
      set=1;
    }
    else { //lets copy the line segment data
        merged_lines[cur][0] = segments[idx*7+0];
        merged_lines[cur][1] = segments[idx*7+1];
        merged_lines[cur][2] = segments[idx*7+2];
        merged_lines[cur][3] = segments[idx*7+3]; 
        merged_lines[cur][4] = 0;//make it a valid line
        set =1;
    }
  }
}

int main(int argc, char **argv)
{
  double * image;
  double * out;
  int x,y,i,j,n;
  int X = 128;  /* x image size */
  int Y = 128;  /* y image size */
  int dim=7;
  int bins = -1;
  int threshold = -1;
  int lengthThreshold=-1;

  /* create a simple image: left half black, right half gray */
  //image = (double *) malloc( X * Y * sizeof(double) );

  //if(argc != 3)
    //exit(EXIT_FAILURE);
  if(argc == 4) //get the number of bins 
    bins = atoi(argv[3]);
  if(argc == 5) //get the number of bins  and threshold T for distance 
  {
    bins = atoi(argv[3]);
    threshold = atoi(argv[4]);
  }
  if(argc == 6) //get the number of bins , threshold T for distance and threshold t6 for length
  {
    bins = atoi(argv[3]);
    threshold = atoi(argv[4]);
    lengthThreshold=atoi(argv[5]);
  }
  if(bins == -1) bins = 5; //lets default if there is no bins option
  if(lengthThreshold ==-1)lengthThreshold =15; //default pixel distance

  char *name = argv[1];

  printf("lenght thresld is %d \n",lengthThreshold);

  printf("File Name is %s \n",name);
  image = (double *)read_pgm_image_double(&X, &Y, name);
  if( image == NULL )
  {
    fprintf(stderr,"error: not enough memory\n");
    exit(EXIT_FAILURE);
  }

  printf(" File dimensions %d, %d", X,Y);

  /* LSD call */
  out = lsd(&n,image,X,Y);

  printf("op file name %s \n",argv[2]);

  FILE * outputFile;
  outputFile = fopen("raw-op.txt","w");

  // FILE * outputFile1;
  // outputFile1 = fopen("raw-op1.txt","w");

  /* print output */
  printf("%d line segments found:\n",n);
  //double * ptr = &out[10*dim+6];
  //Just store length in the last cell since im not using it
  for(i=0;i<n;i++)
  {
    /*fill length in last cell of the double array*/
    double length  = (double)dist(out[i*dim+0], (double) Y - out[i*dim+1], out[i*dim+2], (double) Y - out[i*dim+3]);
    out[i*dim+6] = length;
    out[i*dim+0] = fabs(out[i*dim+0]);
    out[i*dim+1] = fabs(out[i*dim+1]);
    out[i*dim+2] = fabs(out[i*dim+2]);
    out[i*dim+3] = fabs(out[i*dim+3]);
  }

  //printf("the dist of line 10 before sorting is [%f] \n",*ptr);

  // for(i=0;i<n;i++)
  // {
  //   for(j=0;j<dim;j++)
  //     fprintf(outputFile1,"%f ",out[i*dim+j]);
  //   //fprintf(outputFile,"%f ",slopeOfLines[i]);//print slope we calculated
  //   fprintf(outputFile1,"\n");
  // }
  // fclose(outputFile1);

  //sort based on line segment distance that we computed- which is needed to find parallel lines
  qsort(out , n, sizeof(double *) * dim, cmpfunc);

  //printf("the dist of line 10 after sorting is [%f] \n",*ptr);

  //find parallel lines and also calculate slope
  Array * parallelLines = NULL;
  double *slopeOfLines = calloc(n, sizeof(double));
  findParallelLines(&parallelLines, out, n, X, Y, dim, slopeOfLines, bins);

  //TEST
  //qsort(slopeOfLines, n , sizeof(double), cmpfunc1);

  //printf("Min slope is [%f] and max slope is [%f]",slopeOfLines[0],slopeOfLines[n-1]);

  //Print all linesegment details along with calculated slope
  for(i=0;i<n;i++)
  {
    for(j=0;j<dim;j++){
      fprintf(outputFile,"%f ,",out[i*dim+j]);
    }
    fprintf(outputFile,"%f ",slopeOfLines[i]);//print slope we calculated
    fprintf(outputFile,"\n");
  }
  fclose(outputFile);

  //Calculate eucledian distance between each pair of line segment
  double **distArray = (double **) calloc(n , sizeof(double *));
  for (int i=0; i<n; i++)
    distArray[i] = (double *)calloc(n , sizeof(double));
  computeDistances(distArray, out, n, X, Y);
  
  // for(int idx=0; idx <n ; idx++)
  // {
  //   if(idx !=0) continue;
  //   //printf("slope of this line is %f where length is %f ",slopeOfLines[idx], out[idx*dim+6]);
  //   Array *obj = &parallelLines[idx];
  //   //printf("parallel lines for %d is %d \n",idx,obj->used );
  //   write_eps_BS(obj,out, n, dim, "2.eps", X, Y, 2.0);
  // }

  //printf("distance betwen first and second is %f \n", distArray[6][6]);
  //write_eps(out, n, dim, "2.eps", X, Y, 2.0);

  
  //overlap array stores all overlap distance
  double **overlapArray = (double **) malloc(n * sizeof(double *));
  for (int i=0; i<n; i++)
    overlapArray[i] = (double *)malloc(n * sizeof(double));

  for(int idx=0; idx <n ; idx++)
  {
    findOverLapDistance(overlapArray, &parallelLines[idx], slopeOfLines, out, idx, X, Y,n);
  }


  //now we have overlap distance(overlapArray) , slope(slopeOfLines) and distance(distArray)
  // now we can classify and group accordingly.
  Array * allParallelLinesGrouped = NULL;
  Array * allSingleLinesGrouped = NULL;
  classifySegments(distArray, slopeOfLines, overlapArray, n, &allParallelLinesGrouped, &allSingleLinesGrouped, threshold, lengthThreshold, out);
  
  for(int ni =0; ni< n; ni++) {
    char str[50]={0};
    // sprintf(str,"parallels/%d-.eps",ni);
    // write_eps_BS(&allParallelLinesGrouped[ni],out, ni, dim, str, X, Y, 2.0, ni);  
  }
  for(int ni =0; ni< n; ni++) {
    char str[50]={0};
    sprintf(str,"singles/%d-.eps",ni);
    write_eps_BS(&allSingleLinesGrouped[ni],out, ni, dim, str, X, Y, 2.0, ni);  
  }

  /*Merge all lines*/
  gSegments = out; //exposed for comparator

  /*this array stores all the new lne segments that are merged*/
  double merged_lines[n][5]; //0 startx 1 starty 2 endx 3 endy 4 invalid
  for(int ni=0; ni < n; ni++)
    merged_lines[ni][4] = 1;//all are invalid by default
  int lines_of_tower [] = {7,15,16,19,22,24,108,207,222,229,230,271,295,299,337,345,478,477,528,535,671,1002,1235};
  for(int ni =0; ni< 23; ni++)
  {
    if(slopeOfLines[lines_of_tower[ni]] <0)
      continue;
    int * parallel_lines = allSingleLinesGrouped[lines_of_tower[ni]].array;
    qsort(parallel_lines, allSingleLinesGrouped[lines_of_tower[ni]].used , sizeof(int), cmpfunc3);//only y axis for now
    merged_lines[lines_of_tower[ni]][4] = 0;//Make it valid
    merge_lines(merged_lines, lines_of_tower[ni], allSingleLinesGrouped[lines_of_tower[ni]].used, out ,parallel_lines, slopeOfLines, X, Y); //force coersion here
  }

  write_eps_BS1(merged_lines, n, dim, "final.eps", X,  Y);
  /* free memory */
  free( (void *) image );
  free( (void *) out );

  return EXIT_SUCCESS;
}