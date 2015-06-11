/*********************************************************************
MPEG2MOV.C

This is the main file for the MPEG to Matlab movie decoder.

**********************************************************************/

/* MPEG Decoder (mpeg_play) copyright notice */

/*
 * Copyright (c) 1992 The Regents of the University of California.
 * All rights reserved.
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose, without fee, and without written agreement is
 * hereby granted, provided that the above copyright notice and the following
 * two paragraphs appear in all copies of this software.
 * 
 * IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
 * OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
 * CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
 * ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */


/************************************************************************
                              R C S Information
*************************************************************************/

/* $Log: mpgread.c,v $
 * Revision 1.7  1994/01/14  15:19:54  daf
 * Added RCS info header
 *
 * revision 1.6    locked by: daf;
 * date: 1994/01/14 15:13:04;  author: daf;  state: Exp;  lines: +16 -15
 * Fixed calls to MATLAB to get matlab file path
 * 
 * revision 1.5
 * date: 1994/01/12 09:05:59;  author: daf;  state: Exp;  lines: +2 -2
 * Changed includes so string.h is always included.
 *
 * revision 1.4
 * date: 1994/01/11 20:14:02;  author: daf;  state: Exp;  lines: +27 -3
 * Modified for PC compatiblity
 * 
 * revision 1.3
 * date: 1994/01/11 19:02:53;  author: daf;  state: Exp;  lines: +27 -1
 * Now adds .mpg extension if file is not found without it.
 * 
 * revision 1.2
 * date: 1994/01/11 14:48:41;  author: daf;  state: Exp;  lines: +120 -52
 * Added [R,G,B] = ... format
 * 
 * revision 1.1
 * date: 1994/01/07 16:52:26;  author: daf;  state: Exp;
 * Initial revision
*/


#include "video.h"
#include "proto.h"
#include <stdio.h>
#include <sys/types.h>
#include <signal.h>
#ifndef WIN32
#include <memory.h>
#endif
#include <string.h>


#ifdef WIN32
    /* need to reverse byte order */
    /* note -- we assume here that htonl is called on a variable, not a
     * constant; thus, this is not for general use.
     */
#define htonl(x)    \
    ((((unsigned char *)(&x))[0] << 24) |	\
     (((unsigned char *)(&x))[1] << 16) |	\
     (((unsigned char *)(&x))[2] << 8) |	\
     (((unsigned char *)(&x))[3]))
#define ntohl(x)    htonl(x)
#define htons(x)    \
    ((((unsigned char *)(&x))[0] << 8) |	\
     ((unsigned char *)(&x))[1])
#define ntohs(x)    htons(x)
#else
#ifndef MIPS
#include <netinet/in.h>
#else
#include <bsd/netinet/in.h>
#endif /* MIPS */
#endif  /* WIN32 */

#include "util.h"
#include "dither.h"

#include "mex.h"


/* Define buffer length. */

#define BUF_LENGTH 80000

/* External function declarations */
extern int mpegVidRsrc();
extern VidStream *NewVidStream();
extern void ConvertColor();

/* global variable to indicate start of new stream */
int start_decode;

/* Declaration of global variable to hold dither info. */

int ditherType;

/* Global file pointer to incoming data. */
FILE *input;

/* End of File flag. */
static int EOF_flag = 0;

/* Loop flag. */
int loopFlag = 0;

/* Shared memory flag. */
int shmemFlag = 0;

/* Quiet flag. */
int quietFlag = 0;

/* Display image on screen? */
int noDisplayFlag = 0;

typedef enum {ELITTLE_ENDIAN, EBIG_ENDIAN} ByteOrder;

ByteOrder MOVIE_GetClientByteOrder(void)
{
    const short one = 1;
    return((*((char *) &one) == 1) ? ELITTLE_ENDIAN : EBIG_ENDIAN );
}

#define FLIPBYTES(x)                                                     \
(                                                                         \
((x&0xff) << 24) | ((x & 0xff00) << 8) | ((x & 0xff0000) >> 8) | (x >> 24) \
)


/*
 *--------------------------------------------------------------
 *
 * get_more_data --
 *
 *	Called by correct_underflow in bit parsing utilities to
 *      read in more data.
 *
 * Results:
 *	Input buffer updated, buffer length updated.
 *      Returns 1 if data read, 0 if EOF, -1 if error.
 *
 * Side effects:
 *      None.
 *
 *--------------------------------------------------------------
 */

int 
get_more_data(unsigned int *buf_start,
	      int max_length,
	      int *length_ptr,
	      unsigned int **buf_ptr)
{
  
  int length, num_read, i, request;
  unsigned char *buffer, *mark;
  unsigned int *lmark;

  if (EOF_flag) return 0;

  length = *length_ptr;
  buffer = (unsigned char *) *buf_ptr;

  if (length > 0) {
    memcpy((unsigned char *) buf_start, buffer, (length*4));
    mark = ((unsigned char *) (buf_start + length));
  }
  else {
    mark = (unsigned char *) buf_start;
    length = 0;
  }

  request = (max_length-length)*4;
  
  num_read = fread( mark, 1, request, input);

  /* Paulo Villegas - 26/1/1993: Correction for 4-byte alignment */
  {
    int num_read_rounded;
    unsigned char *index;
 
    num_read_rounded = 4*(num_read/4);
 
    /* this can happen only if num_read<request; i.e. end of file reached */
    if( num_read_rounded < num_read )
      { 
 	num_read_rounded = 4*( num_read/4+1 );
 	/* fill in with zeros */
 	for( index=mark+num_read; index<mark+num_read_rounded; *(index++)=0 );
 	/* advance to the next 4-byte boundary */
 	num_read = num_read_rounded;
      }
  }
  
  if   (num_read < 0) {
    return -1;
  }
  else if (num_read == 0) {
    *buf_ptr = buf_start;
    
    /* Make 32 bits after end equal to 0 and 32
       bits after that equal to seq end code
       in order to prevent messy data from infinite
       recursion.
    */

    *(buf_start + length) = 0x0;
    *(buf_start + length+1) = SEQ_END_CODE;

    EOF_flag = 1;
    return 0;
  }

  lmark = (unsigned int *) mark;

  num_read = num_read/4;

  for (i=0; i<num_read; i++) {
    *lmark = htonl((*lmark));
    lmark++;
  }

  *buf_ptr = buf_start;
  *length_ptr = length + num_read;
 
  return 1;
}



/*
======================================================================
mexFunction

This is the Matlab interface function for the MPEG to Matlab decoder.
This routine initializes variables, creates the output matrices and
calls the MPEG encoder functions to decode the MPEG file.

======================================================================
*/
void
mexFunction(int nlhs,
	    mxArray *plhs[],
	    int nrhs,
	    const mxArray *prhs[])
{
  char filename[256];             /* input MPEG file name */
  char err_text[80];              /* string to hold error message */
  static VidStream *theStream;    /* pointer to Video Stream struct. */
  int h,i,j,k;                    /* loop counters */
  double *cm;                     /* pointer to colormap matrix values */
  double *movie;                  /* pointer to movie matrix values */
  double *red;                    /* pointer to R,G,B matrix values */
  double *green;
  double *blue;
  int movie_rows;                 /* number of rows in movie matrix */
  int num_pix_doubles;            /* number of doubles used by each frame */
  int n_colors;                   /* number of colors in colormap */
  unsigned char r,g,b;            /* hold RGB values converted from YUV */
  unsigned char *m_ptr;           /* pointer to bytes in movie matrix */
  double *r_ptr;                  /* pointers to doubles in R,G,B matrices */
  double *g_ptr;
  double *b_ptr;
  unsigned int lum, crom;         /* indices into luminance and crominance */
  int vid_status;                 /* holds return from main decoder call */
  int last_frame;                 /* holds last frame requested by user */
  int n_req_frames;               /* number of requested frames */
  double *req_frames;             /* points to array of frame numbers */
  int frames_stored;              /* counts # of frames stored in matrix */
  int no_colormap;                /* whether or not to return colormap */
  int rgb_format;                 /* RGB if true, otherwise colormapped */
  int NewOutput = 0;              /* 0: Pre-5.3 movie format, 1: 5.3+ movie format */
  int NewTrueColorOutput = 0;     /* 0: Indexed movie, 1: TrueColor movie */
  unsigned int temp_int;          /* used as temp storage for flipping bytes*/
  mxArray *fid[2], *filestr[2];    /* matrices used in searching Matlab path */


  /* check arguments */
  no_colormap = 1;
  rgb_format = 0;
  if (nlhs == 2)
    no_colormap = 0;
  else if (nlhs == 3)
    rgb_format = 1;

  else if (nlhs > 3)
    mexErrMsgTxt("Too many output arguments.");
  else if (nlhs < 1)
    mexErrMsgTxt("Too few output arguments.");

  if (nrhs < 1)
    mexErrMsgTxt("Too few input arguments.");
  if (nrhs > 3)
    mexErrMsgTxt("Too many input arguments.");

  if (!mxIsChar(prhs[0]))
    mexErrMsgTxt("First argument (file name) must be a string.");

  mxGetString(prhs[0], filename, 256);

  /* get list of frame numbers if passed*/
  if ((nrhs >= 2) && (!mxIsEmpty(prhs[1])))
  {
    if (mxGetN(prhs[1]) > 1)
    {
      n_req_frames = mxGetN(prhs[1]);
      if (mxGetM(prhs[1]) > 1)
        mexErrMsgTxt("Second argument must be a vector");
    }
    else
    {
      n_req_frames = mxGetM(prhs[1]);
    }
    req_frames = mxGetPr(prhs[1]);

    /* find the highest number frame requested */
    last_frame = 1;
    for (j = 0; j < n_req_frames; j++)
    {
      if ((int)req_frames[j] > last_frame)
        last_frame = (int)req_frames[j];
    }
  }
  else
  {
    last_frame = 0;
    n_req_frames = 0;
  }

  if (nrhs == 3) {
    char str[100];
    mxGetString(prhs[2], str, 100);
    if (strcmp(str, "truecolor") == 0) {
      NewTrueColorOutput = 1;
    }
    NewOutput = 1;
  }

  /* open the MPEG file */
  /* call Matlab to search the Matlab path */
  filestr[0] = mxCreateString(filename);
  mexCallMATLAB(1,fid,1,filestr,"fopen");
  if (*(mxGetPr(fid[0])) == -1)
  {
    if (strchr(filename, '.') == NULL)
    {
      strcat(filename, ".mpg");
      mxDestroyArray(filestr[0]);
      mxDestroyArray(fid[0]);
      filestr[0] = mxCreateString(filename);
      mexCallMATLAB(1,fid,1,filestr,"fopen");
      if (*mxGetPr(fid[0]) == -1)
        mexErrMsgTxt("Could not open file.");
    }
    else
      mexErrMsgTxt("Could not open file.");
  }
  mxDestroyArray(filestr[0]);
  mexCallMATLAB(1,filestr,1,fid,"fopen");
  mxGetString(filestr[0], filename, 255);
  mxDestroyArray(filestr[0]);
  mexCallMATLAB(1,filestr,1,fid,"fclose");
  mxDestroyArray(fid[0]);
  mxDestroyArray(filestr[0]);
#ifdef WIN32
  input = fopen(filename, "rb");
#else
  input = fopen(filename, "r");
#endif
  if (input == NULL) {
    mexErrMsgTxt("Could not open file.");
  }

  /* initialization */
  ditherType = ORDERED2_DITHER;
  LUM_RANGE = 8;
  CR_RANGE = CB_RANGE = 4;
  noDisplayFlag = 0;

#ifdef SH_MEM
  shmemFlag = 1;
#endif

  init_tables();

  InitColor();
  InitOrdered2Dither();

  
  EOF_flag = 0;
  curBits = 0;
  bitOffset = 0;
  bufLength = 0;
  bitBuffer = NULL;
  totNumFrames = 0;
  
  start_decode = 1;
  curVidStream = NULL;
  theStream = NewVidStream(BUF_LENGTH);
  if (theStream == NULL)
    mexErrMsgTxt("Out of memory.");

  /* parse the MPEG header */
  mpegVidRsrc(0, theStream);

  /* create the movie and colormap Matrices */
  n_colors = LUM_RANGE*CB_RANGE*CR_RANGE;

  /*realTimeStart = ReadSysClock();*/

  if (n_req_frames == 0)
  {
    /* if user did not specify frames to get, one pass is needed to determine
     * the number of frames before actually putting them into the matrix */
    vid_status = 0;
    while (vid_status != 1)
    {
      vid_status = mpegVidRsrc(0, theStream);
    }
    n_req_frames = last_frame = totNumFrames;
    req_frames = (double *)mxCalloc(n_req_frames, sizeof(double));
    for (i = 0; i < n_req_frames; i++)
      req_frames[i] = i+1;
    DestroyVidStream(theStream);
    rewind(input);
    EOF_flag = 0;
    curBits = 0;
    bitOffset = 0;
    bufLength = 0;
    bitBuffer = NULL;
    totNumFrames = 0;
    theStream = NewVidStream(BUF_LENGTH);
    mpegVidRsrc(0, theStream);
  }

  /* create the output matrices */
  if (NewOutput) 
    {
      char *fieldnames[] = {"cdata", "colormap"};
      plhs[0] = mxCreateStructMatrix(1, n_req_frames, 2, (const char **)fieldnames);
      if (NewTrueColorOutput) 
	{
	  int i;
	  int dims[3];
	  dims[0] = theStream->v_size; dims[1] = theStream->h_size; dims[2] = 3;
	  for (i = 0; i < n_req_frames; i++) 
	    {
	      mxArray *cdata = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
	      mxSetField(plhs[0], i, "cdata", cdata);
	    }
	}
      else
	{
	  int i;
	  int dims[2];
	  dims[0] = theStream->v_size; dims[1] = theStream->h_size;
	  for (i = 0; i < n_req_frames; i++) 
	    {
	      mxArray *cdata = mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL);
	      mxArray *colormap = mxCreateDoubleMatrix(n_colors, 3, mxREAL);
	      mxSetField(plhs[0], i, "cdata", cdata);
	      mxSetField(plhs[0], i, "colormap", colormap);
	      
	      {
		/* Fill in the colormap */
		double *cm = mxGetPr(colormap);
		int j;
		for (j = 0; j < n_colors; j++)
		  {
		    /* fprintf(stderr, "Color: %d\n", j); */
		    ConvertColor(lum_values[(j/(CR_RANGE*CB_RANGE))%LUM_RANGE],
				 cr_values[(j/CB_RANGE)%CR_RANGE],
				 cb_values[j%CB_RANGE], 
				 &r, &g, &b);
		    cm[j] = (double)r/255;
		    cm[n_colors + j] = (double)g/255;
		    cm[n_colors * 2 + j] = (double)b/255;
		  }
	      }
	    }
	}
    }
  else
    {
      if (rgb_format)
	{
	  plhs[0] = mxCreateDoubleMatrix(theStream->v_size, theStream->h_size*n_req_frames,
					 mxREAL);
	  if (plhs[0] == NULL)
	    mexErrMsgTxt("Out of memory.");

	  plhs[1] = mxCreateDoubleMatrix(theStream->v_size, theStream->h_size*n_req_frames,
					 mxREAL);
	  if (plhs[1] == NULL)
	    mexErrMsgTxt("Out of memory.");

	  plhs[2] = mxCreateDoubleMatrix(theStream->v_size, theStream->h_size*n_req_frames,
					 mxREAL);
	  if (plhs[2] == NULL)
	    mexErrMsgTxt("Out of memory.");

	  red = mxGetPr(plhs[0]);
	  green = mxGetPr(plhs[1]);
	  blue = mxGetPr(plhs[2]);
	  movie_rows = mxGetM(plhs[0]);
	}
      else
	{
	  num_pix_doubles = (int)((theStream->h_size*theStream->v_size + 7)/8);
  
	  if (n_req_frames > 0)
	    plhs[0] = mxCreateDoubleMatrix(388 + num_pix_doubles, n_req_frames, mxREAL);

	  if (plhs[0] == NULL)
	    mexErrMsgTxt("Out of memory.");
	  movie = mxGetPr(plhs[0]);
	  movie_rows = mxGetM(plhs[0]);
	}

      if (!no_colormap)
	{
	  plhs[1] = mxCreateDoubleMatrix(n_colors, 3, mxREAL);
	  if (plhs[1] == NULL)
	    mexErrMsgTxt("Out of memory.");
	  cm = mxGetPr(plhs[1]);
	  /* create a Matlab colormap */
	  /* fprintf(stderr, "About to create colormap\n");*/
	  for (j = 0; j < n_colors; j++)
	    {
	      /* fprintf(stderr, "Color: %d\n", j); */
	      ConvertColor(lum_values[(j/(CR_RANGE*CB_RANGE))%LUM_RANGE],
			   cr_values[(j/CB_RANGE)%CR_RANGE],
			   cb_values[j%CB_RANGE], 
			   &r, &g, &b);
	      cm[j] = (double)r/255;
	      cm[n_colors + j] = (double)g/255;
	      cm[n_colors * 2 + j] = (double)b/255;
	    }
	}
    }
  /*fprintf(stderr, "movie rows: %d\n", movie_rows);*/

  /* get requested frames and store them in output martrix (matrices).*/
  frames_stored = 0;
  for (j = 1; j <= last_frame; j++)
  {
    /*fprintf(stderr, "Generating frame #%d\n", j);*/
    while (totNumFrames < j)
    {
      vid_status = mpegVidRsrc(0, theStream);
      if (vid_status == 1)
      {
        sprintf(err_text, "Frame(s) requested beyond last frame (%d).",
                totNumFrames);
	mexErrMsgTxt(err_text);
      }
    } 

    /* store this frame wherever it has been requested */   
    for (i = 0; i < n_req_frames; i++)
    {
      if ((int)req_frames[i] == j)
      {
        frames_stored++;
        
	if (NewOutput) 
	  {
	    mxArray *cdata = mxGetField(plhs[0], i, "cdata");
	    unsigned char *pr = (unsigned char *)mxGetPr(cdata);
	    if (NewTrueColorOutput) 
	      {
		unsigned char *pg = pr + theStream->v_size * theStream->h_size;
		unsigned char *pb = pg + theStream->v_size * theStream->h_size;

		lum = 0;
		crom = 0;
		for (h = 0; h < theStream->h_size; h++) 
		  {
		    lum = h;
		    crom = h/2;
		    for (k = 0; k < theStream->v_size; k++) 
		      {
			ConvertColor(*(curVidStream->current->luminance + lum),
				     *(curVidStream->current->Cr + crom),
				     *(curVidStream->current->Cb + crom),
				     &r, &g, &b);
			*pr = (unsigned char)r;
			*pg = (unsigned char)g;
			*pb = (unsigned char)b;
			pr++; pg++; pb++;

			lum += theStream->h_size;
			if ((k % 2) == 1)
			  {
			    crom += theStream->h_size / 2;
			  }
		      }
		  }
	      }
	    else
	      {
		int index = 0;
		for (h = 0; h < theStream->h_size; h++) 
		  {
		    for (k = 0; k < theStream->v_size; k++) 
		      {
			index = k * theStream->h_size + h;
			*pr = *(curVidStream->current->display + index);
			pr++;
		      }
		  }
	      }
	  }
	else
	  {
	    if (rgb_format)
	      {
		lum = 0;
		crom = 0;
		for (h = 0; h < theStream->v_size; h++)
		  {
		    r_ptr = red + theStream->v_size*theStream->h_size * i + h;
		    g_ptr = green + theStream->v_size*theStream->h_size * i + h;
		    b_ptr = blue + theStream->v_size*theStream->h_size * i + h;
		    for (k = 0; k < theStream->h_size; k++)
		      {
			ConvertColor(*(curVidStream->current->luminance + lum),
				     *(curVidStream->current->Cr + crom),
				     *(curVidStream->current->Cb + crom),
				     &r, &g, &b);
			*r_ptr = (double)r/255;
			*g_ptr = (double)g/255;
			*b_ptr = (double)b/255;
              
			r_ptr += theStream->v_size;
			g_ptr += theStream->v_size;
			b_ptr += theStream->v_size;

			lum++;
			if (k%2 == 0)
			  crom++;
		      }
		    if (h%2 == 0)
		      crom -= theStream->h_size / 2;
		  }
	      }
	    else
	      {
		/* create Matlab movie frame header */
		m_ptr = (unsigned char *)(movie + (movie_rows * i));
		*((double *)m_ptr) = curVidStream->h_size;
		m_ptr += sizeof(double);
		*((double *)m_ptr) = curVidStream->v_size;
		m_ptr += sizeof(double);
		*((double *)m_ptr) = 0;
		m_ptr += sizeof(double);
		if (MOVIE_GetClientByteOrder() == ELITTLE_ENDIAN)
		  {
		    *((int *)m_ptr) = 0;
		    m_ptr += sizeof(int);
		    *((int *)m_ptr) = LUM_RANGE*CB_RANGE*CR_RANGE;
		    m_ptr += sizeof(int);
		  }
		else
		  {
		    *((int *)m_ptr) = LUM_RANGE*CB_RANGE*CR_RANGE;
		    m_ptr += sizeof(int);
		    *((int *)m_ptr) = 0;
		    m_ptr += sizeof(int);
		  }
		m_ptr += 256*3*sizeof(float);

          /* copy frame data into movie matrix */
		memcpy(m_ptr, curVidStream->current->display, 
		       (theStream->h_size*theStream->v_size));
		if (MOVIE_GetClientByteOrder() == ELITTLE_ENDIAN)
		  {
		    for (k = 0; k < (movie_rows - 388); k++)
		      {
			temp_int = FLIPBYTES(*((unsigned int *)m_ptr));
			*((unsigned int *)m_ptr) = 
			  FLIPBYTES(*((unsigned int *)m_ptr + 1));
			*((unsigned int *)m_ptr + 1) = temp_int;
			m_ptr += sizeof(double);
		      }
		  }
	      }
	  }
      }
    }
  }


  /* destroy the video stream */
  DestroyVidStream(theStream);
  fclose(input);
}


/*
 *--------------------------------------------------------------
 *
 * DoDitherImage --
 *
 *	Called when image needs to be dithered. Selects correct
 *      dither routine based on info in ditherType.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	None.
 *
 *--------------------------------------------------------------
 */

void
DoDitherImage(unsigned char *l,
	      unsigned char *Cr, 
	      unsigned char *Cb,
	      unsigned char *disp,
	      int h,
	      int w)
{
  Ordered2DitherImage(l, Cr, Cb, disp, h, w);
}



















