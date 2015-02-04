#ifdef __STDC__
# define	P(s) s
#else
# define P(s) ()
#endif


/* util.c */
void correct_underflow P((void ));
int next_bits P((int num , unsigned int mask ));
char *get_ext_data P((void ));
int next_start_code P((void ));
char *get_extra_bit_info P((void ));

/* video.c */
void init_stats P((void ));
void PrintAllStats P((void ));
/*double ReadSysClock P((void ));*/
void PrintTimeInfo P((void ));
VidStream *NewVidStream P((int bufLength ));
void DestroyVidStream P((VidStream *astream ));
PictImage *NewPictImage P((unsigned int width , unsigned int height ));
void DestroyPictImage P((PictImage *apictimage ));
int mpegVidRsrc P((TimeStamp time_stamp , VidStream *vid_stream ));
void ToggleBFlag P((void ));
void TogglePFlag P((void ));

/* parseblock.c */
void ParseReconBlock P((int n ));
void ParseAwayBlock P((int n ));

/* motionvector.c */
void ComputeForwVector P((int *recon_right_for_ptr , int *recon_down_for_ptr ));
void ComputeBackVector P((int *recon_right_back_ptr , int *recon_down_back_ptr ));

/* decoders.c */
void init_tables P((void ));
void decodeDCTDCSizeLum P((unsigned int *value ));
void decodeDCTDCSizeChrom P((unsigned int *value ));
void decodeDCTCoeffFirst P((unsigned int *run , int *level ));
void decodeDCTCoeffNext P((unsigned int *run , int *level ));

/* main.c */
int get_more_data P((unsigned int *buf_start , int max_length , int *length_ptr , unsigned int **buf_ptr ));
void int_handler P((void ));
void main P((int argc , char **argv ));
void usage P((char *s ));
void DoDitherImage P((unsigned char *l , unsigned char *Cr , unsigned char *Cb , unsigned char *disp , int h , int w ));

/* gdith.c */
void InitColor P((void ));
/* int HandleXError P((Display *dpy , XErrorEvent *event ));*/
void InstallXErrorHandler P((void ));
void DeInstallXErrorHandler P((void ));
void ResizeDisplay P((int w , int h ));
void InitDisplay P((char *name ));
void InitGrayDisplay P((char *name ));
void InitMonoDisplay P((char *name ));
void InitColorDisplay P((char *name ));
void ExecuteDisplay P((VidStream *vid_stream ));


/* ordered2.c */
void InitOrdered2Dither P((void ));
void Ordered2DitherImage P((unsigned char *lum , unsigned char *cr , unsigned char *cb , unsigned char *out , int h , int w ));

#undef P
