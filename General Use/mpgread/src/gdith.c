/* gdith.c
   General functions related to dithering
*/

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

/***********************************************************************
                          R C S Information
************************************************************************/

/* $Log: gdith.c,v $
 * Revision 1.3  1994/01/14  15:45:44  daf
 * Added RCS info header
 * /
 *
 * revision 1.2    locked by: daf;
 * date: 1994/01/07 16:45:41;  author: daf;  state: Exp;  lines: +21 -467
 * Modified for use as .mex file
 *   XWindows calls removed
 * 
 * revision 1.1
 * date: 1994/01/07 16:13:52;  author: daf;  state: Exp;
 * Initial revision
*/


#include <math.h>
#include "video.h"
#include "proto.h"
#include "dither.h"
#include "mex.h"

/* Range values for lum, cr, cb. */
int LUM_RANGE;
int CR_RANGE;
int CB_RANGE;

/* Array that remaps color numbers to actual pixel values used by X server. */

unsigned char pixel[256];

/* Arrays holding quantized value ranged for lum, cr, and cb. */

int *lum_values;
int *cr_values;
int *cb_values;

/* Declaration of global variable containing dither type. */

extern int ditherType;

/* Structures used by the X server. */
/*
Display *display;

static XImage *ximage = NULL;
static Colormap cmap;
static Window window;
static GC gc;
*/


/*
 *--------------------------------------------------------------
 *
 * InitColor --
 *
 *	Initialized lum, cr, and cb quantized range value arrays.
 *
 * Results: 
 *      None.
 *
 * Side effects:
 *      None.
 *
 *--------------------------------------------------------------
 */

void
InitColor(void)
{
  int i;


  lum_values = (int *) mxCalloc(LUM_RANGE,sizeof(int));
  cr_values = (int *) mxCalloc(CR_RANGE,sizeof(int));
  cb_values = (int *) mxCalloc(CB_RANGE,sizeof(int));

  for (i=0; i<LUM_RANGE; i++) {
    lum_values[i] = ((i * 256) / (LUM_RANGE)) + (256/(LUM_RANGE*2));
  }

  for (i=0; i<CR_RANGE; i++) {
    cr_values[i] = ((i * 256) / (CR_RANGE)) + (256/(CR_RANGE*2));
  }

  for (i=0; i<CB_RANGE; i++) {
    cb_values[i] = ((i * 256) / (CB_RANGE)) + (256/(CB_RANGE*2));
  }

  /* set up pixel color array used by dithering routine */
  /* This is required to replace the X stuff in InitDisplay. */
  for (i = 0; i < (LUM_RANGE*CR_RANGE*CB_RANGE); i++)
    pixel[i] = (unsigned char)i;

}


/*
 *--------------------------------------------------------------
 *
 * ConvertColor --
 *
 *	Given a l, cr, cb tuple, converts it to r,g,b.
 *
 * Results:
 *	r,g,b values returned in pointers passed as parameters.
 *
 * Side effects:
 *      None.
 *
 *--------------------------------------------------------------
 */

void
ConvertColor(unsigned char l,
	     unsigned char cr, 
	     unsigned char cb,
	     unsigned char *r,
	     unsigned char *g,
	     unsigned char *b)
{
  double fl, fcr, fcb, fr, fg, fb;

  fl = (double) l;
  fcr =  ((double) cr) - 128.0;
  fcb =  ((double) cb) - 128.0;


  fr = fl + (1.40200 * fcb);
  fg = fl - (0.71414 * fcb) - (0.34414 * fcr);
  fb = fl + (1.77200 * fcr);

  if (fr < 0.0) fr = 0.0;
  else if (fr > 255.0) fr = 255.0;

  if (fg < 0.0) fg = 0.0;
  else if (fg > 255.0) fg = 255.0;

  if (fb < 0.0) fb = 0.0;
  else if (fb > 255.0) fb = 255.0;

  *r = (unsigned char) fr;
  *g = (unsigned char) fg;
  *b = (unsigned char) fb;

}




