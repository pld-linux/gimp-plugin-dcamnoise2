/*
 * Copyright (C) 2005 Peter Heckert
 *                    <peter /dot/ heckert /at/ arcor /dot/ de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

 /*
  * =========Version History=============
  *
  *
  * 13. July, Version 0.62
  *
  * All radius,lookahead, and Phase Jitter damping params are internally processed with pixel granularity.
  * (This is unvisible at user interface level, because the step_increment param for gimp_scale_entry_new()
  *  apparently doesnt work as advertised)
  *
  * This prepares for an important change: I will go from double float to single float and I will use
  * an FIR quad boxfilter instead of an IIR gauss. This should give considerable speedup.
  *
  * Introduced new param "Erosion". The new filter gives better sharpness and this also gives problems
  * with spike noise. The "Erosion" param erodes singular spikes and it has a smooth effect to edges, and sharpens
  * edges by erosion, so noise at edges is eroded.  
  * The effect is dependant from sharpness,phase-jitter damping and lookahead. Use it only as last ressort.
  *
  * Using this params at the above mentioned ISO 800 Image I got better results than noise ninja:
  * Radius 4, Thresh. 0.15 , Texture -0.42, Sharpness 0.5, Erosion 3, gamma 1.9, other params default.
  *
  * Updated usage instructions.
  *
  * 5. Aug, Version 0.63
  * Better Filter effect. Internaly a noise clipping threshold is used.
  * If you used previous versions then jugde yourself.
  * Possibly I will make this param adjudstable, it depend on future experiments and on user-feedback 
  *
  * Speeded up the gaussian filter. Now it should be the fastest and the best gauss filter in the west ;-)
  * Therefore I dropped any plans to use boxfilters.
  * Found out that the main speed-bottleneck is gimps plugin-interface.
  * Need for 16 Bit. The filter can not only remove sensor noise, it can remove 8-bit quantization noise.
  * So it should be possible to convert an 8-Bit image into a pseudo-16 Bit image which would allow
  * color and gradadation manipulations which are impossible now.
  * Thinking about porting the whole thing to cinepaint, to overcome the speed-bottleneck and the 8-Bit
  * bottleneck 
  *
  * 6. Aug Version 0.63
  * Corrected minor error in code, uploaded again, without version change
  */


  
  
  
/* ========= SHORT DESCRIPTION AND REQUIREMENTS ==================================================================
  * This is a gimp plugin to remove typical digital camera noise.  
 *
 * Gimp >= 2.2 is required. 
 *
 * If you use windows, then possibly you cannot compile it.
 *
 * By courtesy of Michael Schumacher, there is a windows binary for dcamnoise2-0.54.
 *
 * For windows there also are commercial alternatives:
 * Look at http://www.imagenomic.com, they have a very nice denoising tool
 * for windows, in a commercial and in a light free version.
 * (This is my personal opinion, there are other tools (neatimage,noiseninja....)
 * Of course they all have disadvantages: They are not gimp-plugins.
 *
 * And possibly my tool does a better job in preserving sharpness without artifacts and segmentation, because I
 * do not use something artificial like edge recognition. 
 *
 * There is another tool "Picture cooler" at <http://beam.to/picturecooler> that looks very 
 * interesting. Presently it is beta and free.
 *
 * Code for dcamnoise2.c is derived from common unsharp-plugin for gimp-2.2.
 * I used this as a sceleton, but idea and core code are basing on my own
 * ideas and work.
 * It is still very experimental code and a little bit chaotic.
 * It has undergone a lot of changes in the last weeks.
 *
 * ======== INSTALLATION =========================================================================================
 * Installation hints for  Linux systems:
 * INSTALL with: $ gimptool  --install dcamnoise2-#.##.c (#.## is the version number)
 * (gimp-devel must be installed before)
 * (This is for a rpm-based system, if you use Debian or else, I dont know exactly)
 *
 * It will install in Gimp under "Filter->Verbessern"  
 * (This should be "Filters->Enhance" in the english version)
 * To uninstall, simply delete it from  ~/.gimp-2.2/plugins
 * or type gimptool --uninstall-bin dcamnoise2-#.## 
 * Or type gimptool --help and learn about it ;-)
 * On some sytems you must use "gimptool-2.2" instead of "gimptool" to compile and install it.
 *
 * =========USAGE ================================================================================================
 * Let me explain, how the filter works, some understanding is necessary to use it:
 *
 * Hint for the novice user:
 * In most cases only Filter Max Radius, Filter treshold and Texture Detail are needed and the other
 * params can be left at their default setting.
 *
 *...... Main Filter (Preprocessing) ........................................................................
 * First, a filtered template is generated, using an adaptive filter.
 * To see this template, we must set _Luminance tolerance,  _Color tolerance to 1.0.
 *------------------------------------------------------------------------------------------------------------
 *
 * "Filter max. Radius" is preset to 5.0 
 * This is good for most noise situations. 
 * In any case it must be about the same size as noise granularity ore somewhat more.
 * If it is set higher than necessary, then it can cause unwanted blur.
 *------------------------------------------------------------------------------------------------------------
 * "Filter Threshold" should be set so that edges are clearly visible and noise is  smoothed out. 
 * This threshold value is not bound to any intensity value, it is bound to the second derivative of
 * intensity values.
 * Simply adjust it and watch the preview. Adjustment must be made carefully, because the gap
 * between "noisy", "smooth", and "blur" is very small. Adjust it as carefully as you would adjust
 * the focus of a camera.
 *------------------------------------------------------------------------------------------------------------
 *
 * "Lookahead" defines the pixel distance in which the filter looks ahead for luminance variations
 * Normally the default value should do.
 * When _Lookahead is increased, then spikenoise is erased. 
 * Eventually readjust Filter treshold, when you changed lookahead.
 * When the value is to high, then the adaptive filter cannot longer accurately track image details, and
 * noise can reappear or blur can occur.
 *
 * Minimum value is 1.0,this gives best accuracy when blurring very weak noise.
 *
 * I never had good success with other values than 2.0.
 * However, for images with extemely high or low resolution another value possibly is better.
 * Use it only as a last ressort.
 * It can destroy your computer and can kill cats :-)
 *------------------------------------------------------------------------------------------------------------
 *
 * "Phase Jitter Damping" defines how fast the adaptive filter-radius reacts to luminance variations.
 * I have preset a value, that should do in most cases.
 * If increased, then edges appear smoother, if too high, then blur may occur.
 * If at minimum then noise and phase jitter at edges can occur.
 * It can supress Spike noise when increased and this is the preferred method to remove spike noise.
 *------------------------------------------------------------------------------------------------------------
 *
 * "Sharpness" does just what it says, it improves sharpness. It improves the frequency response for the filter.
 * When it is too strong then not all noise can be removed, or spike noise may appear.
 * Set it near to maximum, if you want to remove weak noise or JPEG-artifacts, without loosing detail.
 *
 *------------------------------------------------------------------------------------------------------------
 *
 * Introduced new param "Erosion". The new filter gives better sharpness and this also gives problems
 * with spike noise. The Erosion param erodes singular spikes and it has a smooth effect to edges, and sharpens
 * edges by erosion, so noise at edges is eroded.  
 * The effect is dependant from sharpness,phase-jitter damping and lookahead.
 * Set it to minimum (zero), if you want to remove weak noise or JPEG-artifacts.
 * When "Erosion" is increased, then also increasing "Phase Jitter Damping" is often useful 
 *
 * It works nicely. Apart from removing spike noise it has a sharpening and antialiasing effect to edges 
 * (Sharpening occurs by erosion, not by deconvolution) 
 *------------------------------------------------------------------------------------------------------------
 * 
 * "Texture Detail" can be used, to get more or less texture accuracy.
 * When decreased, then noise and texture are blurred out, when increased then texture is
 * amplified, but also noise will increase.
 * It has almost no effect to image edges, opposed to Filter theshold, which would blur edges, when increased.  
 *
 * E.g. if Threshold is adjusted in away so that edges are sharp, and there is still too much area noise, then
 * Texture detail could be used to reduce noise without blurring edges.
 * (Another way would be to decrease radius and to increase threshold)
 *
 * It can make beautyful skin :-)
 *------------------------------------------------------------------------------------------------------------
 *------------------------------------------------------------------------------------------------------------
 * The filtered image that is now seen in the preview, is used as template for the following processing steps,
 * therefore it is  important to do this adjustment in first place and to do it as good as possible.
 *------------------------------------------------------------------------------------------------------------
 *------------------------------------------------------------------------------------------------------------
 * 
 *..... Combining original image and filtered image, using tolerance thresholds (Postprocessing).............. 
 * This can give a final touch of sharpness to your image. 
 * It is not necessary to do this, if you want to reduce JPEG-artifacts or weak noise.
 * It's purpose is to master strong noise without loosing too much sharpness.
 *
 * Note, that this all is done in one filter invocation. Preprocessing and postprocessing is done in one run,
 * but logically and in the algorithm they are different and ordered processes.
 *
 *
 * Adjust _Color tolerance or/and Luminance tolerance, (if necessary) so that you get the final image.
 * I recommend to use only one, either _Color  or _Luminance. 
 * These settings dont influence the main smoothing process. What they really do is this:
 *
 * The tolerance values are used as error-thresholds to compare the filtered template with the original
 * image. The plugin algorithm uses them to combine the filtered template with the original image
 * so that  noise and filter errors (blur) are thrown out.
 * A filtered pixel, that is too far away from the original pixel will be overriden by original image content. 
 *
 * Hint:
 * If you cange other sliders, like lookahead or Texture Detail, then you should set color tolerance and 
 * luminance tolerance to 1.0 (right end), because otherwise the filtered template is partially hidden 
 * and e.g. the effects for the damping filter cant be seen clearly and cant be optimized. 
 *------------------------------------------------------------------------------------------------------------
 *
 * _Gamma can be used to increase the tolerance values for darker areas (which commonly are more noisy)
 * This results in more blur for shadow areas.
 *
 * Hint for users of previous versions:
 * Gamma also influences the main-filter process. While the previous version did not have this feature,
 * I have reimplemented it, however, the algorithm used is totally new.
 * 
 *
 * Keep in mind, how the filter works, then usage should be easy!  
 *
 *
 * ============== THANKS ======================================================================================
 * 
 * Thanks go to the gimp-developers for their very fine work!
 * I got a lot of positive feedback for this plugin, thanks!
 * 
 */ 

/* 
 * ================ THEORY AND TECHNIC =======================================================================
 *
 * Some interesting things (theoretic and technic)
 * This plugin bases on the assumption, that noise has no 2-dimensional correlation and therefore
 * can be removed in a 1-dimensional process.
 * To remove noise, I use a four-times  boxfilter with variable radius.
 *
 * The radius is calculated from 2nd derivative of pixeldata.
 * A gauss filter is used to calculte 2nd derivative.
 * The filter has some inbuilt features to clip low amplitude noise to clip very high values that would
 * slow down response time.
 * The 2nd derivative is lowpassfiltered and then radius is calculated as (Filter Treshold)/2nd_derivative. 
 * The radius modulation data is precalulated and buffered an is used to steer filter radius when
 * the actual filtering occurs.
 *
 * Noise and texture can be further supressed by nonlinear distortion before adaptive filtering.
 * To make this possible I subtract low frequency from image data before denoising, so that I get a 
 * bipolar, zerosymmetric image signal.
 *
 * The filter works in a /one-dimensional/ way. It is applied to x and then to y axis.
 *
 * After filtering a zerodimensional point operator  (pixel by pixel comparison)  is used, where 
 * filter-errors are thrown out.
 * This is meant to limit and control filter errors,it can give "final touch" to the image, but it has 
 * nothing to do with the main filter process. 
 *
 *...............................................................................................
 *
 * I dont know if something like this filter already exists.
 * It is all based on my own ideas and experiments.
 * Possibly a separable adaptive gauss-filter is a new thing.
 * Also it is an impossible thing, from a mathemathical point of view ;-)
 * It is possible only for bandwidth limited images.
 * Happyly most photographic images are bandwidth limited, or when they are noisy then we want
 * to limit banwith locally. And this is, what the filter does: It limits bandwidth locally, dependent
 * from (approximately) 2nd derivative of intensity.
 * 
 * Because gauss filtering is essentially linear diffusion, and because this filter uses a variable
 * nonlinear modulated gaussfilter (four box passes are almost gauss) we could say, that this filter
 * implements a special subclass of nonlinear adaptive diffusion, which is separable, and indeed, 
 * results are very similar to nonlinear diffusion filters.
 * However, because the filter is separable, it is much faster and needs less memory.
 *
 */

#define STANDALONE


#ifndef STANDALONE
#include "config.h"
#endif


#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <gtk/gtk.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>


#ifdef STANDALONE

#define INIT_I18N voidproc
void voidproc(void){};
#define _(x) (x)
#define N_(x)  (x)

#else

#include "libgimp/stdplugins-intl.h"

#endif


#define PLUG_IN_VERSION "0.63"

#define SCALE_WIDTH   150
#define ENTRY_WIDTH     4

// This is for testing only!
//#define float double
//#define gfloat gdouble

/* Uncomment this line to get a rough estimate of how long the plug-in
 * takes to run.
 */

/*  #define TIMER  */


typedef struct
{
  gdouble  radius;
  gdouble  lsmooth;
  gdouble  csmooth;
  gdouble  effect;
  gdouble  lookahead;
  gdouble  gamma;
  gdouble  damping;
  gdouble phase;
  gdouble texture;
  gdouble sharp;
  gboolean update_preview;

} UnsharpMaskParams;

typedef struct
{
  gboolean  run;
} UnsharpMaskInterface;

/* local function prototypes */
static void      query (void);
static void      run   (const gchar      *name,
                        gint              nparams,
                        const GimpParam  *param,
                        gint             *nreturn_vals,
                        GimpParam       **return_vals);

static void      blur_line           (gfloat  * const data,
                                      gfloat  * const data2,
                                      gfloat  * const buffer,
                                      gfloat *rbuf,
                                      gfloat *tbuf,
                                      const guchar   *src,
                                      guchar         *dest,
                                      gint            len,
                                      glong           bytes);

static void      unsharp_region      (GimpPixelRgn   *srcPTR,
                                      GimpPixelRgn   *dstPTR,
                                      gint            bytes,
                                      gdouble         radius,
                                      gdouble         lsmooth,
                                      gint            x1,
                                      gint            x2,
                                      gint            y1,
                                      gint            y2,
                                      gboolean        show_progress);

static void      unsharp_mask        (GimpDrawable   *drawable,
                                      gdouble         radius,
                                      gdouble         lsmooth);

static gboolean  unsharp_mask_dialog (GimpDrawable   *drawable);
static void      preview_update      (GimpPreview    *preview);


/* create a few globals, set default values */
static UnsharpMaskParams unsharp_params =
  {
    5.0, /* default radius = 5 */
    1.0, /* Luminance Tolerance */
    1.0, /* RGB Tolerance  */
    0.08, /* Adaptive filter-effect threshold */
    2.0, /* Lookahead */ 
    1.0, /* Filter gamma */
    5.0,  /* Phase jitter Damping  */
    1.0,  /* Area Noise Clip */
    0.0,  /* Texture Detail */
    0.25,  /* Sharpness factor */  
    TRUE /* default is to update the preview */
  /* These values are values that I used for a test image with some success. 
     They are not optimal for every image*/
  
  };

/* Setting PLUG_IN_INFO */
GimpPlugInInfo PLUG_IN_INFO =
  {
    NULL,  /* init_proc  */
    NULL,  /* quit_proc  */
    query, /* query_proc */
    run,   /* run_proc   */
  };


gfloat lut_gamma, inv_gamma;  
  
gfloat lut[256];  
  
static void 
lut_init(gfloat g)
{
 int i;
 
 for (i = 1; i< 255; i++) lut[i] = pow((gfloat) i/255.0, g);
 lut[0]= 0.0;
 lut[255] = 1.0;

 lut_gamma = g;
 inv_gamma = 1.0/g;

}




MAIN ()

static void
query (void)
{
  static GimpParamDef args[] =
    {
      { GIMP_PDB_INT32,    "run_mode",  "Interactive, non-interactive" },
      { GIMP_PDB_IMAGE,    "image",     "(unused)" },
      { GIMP_PDB_DRAWABLE, "drawable",  "Drawable to draw on" },
      { GIMP_PDB_FLOAT,    "radius",    "Radius of gaussian blur (in pixels)" },
      { GIMP_PDB_FLOAT,    "lsmooth",    "Luminance Tolerance" },
      { GIMP_PDB_FLOAT,    "csmooth", "Color Tolerance" },
      { GIMP_PDB_FLOAT,    "effect", "Threshold for 2nd derivative of luminance" },      
      { GIMP_PDB_FLOAT,    "lookahead", "Sharpness" } ,      
      
      { GIMP_PDB_FLOAT,    "gamma", "Gamma" } ,      
      
      { GIMP_PDB_FLOAT,    "damping", "Phase jitter damping" },      
      { GIMP_PDB_FLOAT,    "phase", "Phase shift for edges" },      
      { GIMP_PDB_FLOAT,    "texture", "Texture accuracy" },      
      { GIMP_PDB_FLOAT,    "sharp", "Edge accuracy" }      
    
    
    };

  gimp_install_procedure ("plug_in_dcamnoise2-"PLUG_IN_VERSION,
                          "A Digital Camera Noise filter",
                          "It is commonly "
                          "used on photographic images, and is provides a much "
                          "more pleasing result than the standard denoising "
                          "filters.",

                          "This is an experimental version. ",
			  "",
			  "",
			    

			  N_("_Dcam Noise 2 "PLUG_IN_VERSION" ..."),
                          "GRAY*, RGB*",
                          GIMP_PLUGIN,
                          G_N_ELEMENTS (args), 0,
                          args, NULL);

  gimp_plugin_menu_register ("plug_in_dcamnoise2-"PLUG_IN_VERSION,
                             "<Image>/Filters/Enhance");
}

static void
run (const gchar      *name,
     gint              nparams,
     const GimpParam  *param,
     gint             *nreturn_vals,
     GimpParam       **return_vals)
{
  static GimpParam   values[1];
  GimpPDBStatusType  status = GIMP_PDB_SUCCESS;
  GimpDrawable      *drawable;
  GimpRunMode        run_mode;
#ifdef TIMER
  GTimer            *timer = g_timer_new ();
#endif

  run_mode = param[0].data.d_int32;

  *return_vals  = values;
  *nreturn_vals = 1;

  values[0].type          = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  INIT_I18N ();

  /*
   * Get drawable information...
   */
  drawable = gimp_drawable_get (param[2].data.d_drawable);
  gimp_tile_cache_ntiles (2 * (drawable->width / gimp_tile_width () + 1));

  switch (run_mode)
    {
    case GIMP_RUN_INTERACTIVE:
      gimp_get_data ("plug_in_dcamnoise2-"PLUG_IN_VERSION, &unsharp_params);
      /* Reset default values show preview unmodified */

      /* initialize pixel regions and buffer */
      if (! unsharp_mask_dialog (drawable))
        return;

      break;

    case GIMP_RUN_NONINTERACTIVE:
      if (nparams != 13)
        {
          status = GIMP_PDB_CALLING_ERROR;
        }
      else
        {
          unsharp_params.radius = param[3].data.d_float;
          unsharp_params.lsmooth = param[4].data.d_float;
          unsharp_params.csmooth = param[5].data.d_float;
          unsharp_params.effect = param[6].data.d_float;
          unsharp_params.lookahead = param[7].data.d_float;
          unsharp_params.gamma = param[8].data.d_float;
          unsharp_params.damping = param[9].data.d_float; 
          unsharp_params.phase = param[10].data.d_float; 
          unsharp_params.texture = param[11].data.d_float; 
          unsharp_params.sharp = param[12].data.d_float; 
          
          /* make sure there are legal values */
          if ((unsharp_params.radius < 0.0) ||
              (unsharp_params.lsmooth < 0.0))
            status = GIMP_PDB_CALLING_ERROR;
        }
      break;

    case GIMP_RUN_WITH_LAST_VALS:
      gimp_get_data ("plug_in_dcamnoise2-"PLUG_IN_VERSION, &unsharp_params);
      break;

    default:
      break;
    }

  if (status == GIMP_PDB_SUCCESS)
    {
      drawable = gimp_drawable_get (param[2].data.d_drawable);

      /* here we go */
      unsharp_mask (drawable, unsharp_params.radius, unsharp_params.lsmooth);

      gimp_displays_flush ();

      /* set data for next use of filter */
      gimp_set_data ("plug_in_dcamnoise2-"PLUG_IN_VERSION, &unsharp_params,
                     sizeof (UnsharpMaskParams));

      gimp_drawable_detach(drawable);
      values[0].data.d_status = status;
    }

#ifdef TIMER
  g_printerr ("%f seconds\n", g_timer_elapsed (timer, NULL));
  g_timer_destroy (timer);
#endif
}


static struct iir_param
{
 gdouble B,b1,b2,b3,b0,r,q;
 gdouble *p;
} iir; 
#define GAUSS
#define PI 3.141592653  
static void iir_init(double r)
{
  if (iir.r == r) return;
  iir.r = r; // = unsharp_params.damping;
  gdouble q,l;
  
  if ( r >= 2.5) q = 0.98711 * r - 0.96330;
  else q = 3.97156 - 4.14554 * sqrt(1.0-0.26891 * r);
  
  iir.q = q;
  iir.b0 = 1.57825 + ((0.422205 * q  + 1.4281) * q + 2.44413) *  q;
  iir.b1 = ((1.26661 * q +2.85619) * q + 2.44413) * q / iir.b0;
  iir.b2 = - ((1.26661*q +1.4281) * q * q ) / iir.b0;
  iir.b3 = 0.422205 * q * q * q / iir.b0;
    
  iir.B = 1.0 - (iir.b1 + iir.b2 + iir.b3);


}

static inline gdouble bsqrt(gdouble val)
{
 if (val >= 0.0) return sqrt(val);
 else return -sqrt(-val);
}

static inline gdouble sq(gdouble val)
{
  return val*val;
}

static inline gdouble bsq(gdouble val)
{
  return fabs(val) * val;
}     

static inline gdouble mypow( gdouble val,gdouble ex)
{
 //return val;
 if (fabs(val) < 1e-16) return 0;
 if (val > 0.0) return exp(log(val)*ex);
 return -exp(log(-val)*ex);
} 




static void box_filter(gdouble *src, gdouble *end, gdouble * dest, gdouble radius)
/* src and dest must be different */
{
    gfloat fbw = 2.0 * radius;
    if (fbw < 1.0) fbw = 1.0;
    gfloat box = (*src); 
    
    gint boxwidth=1;

    while(boxwidth+2 <= (int) fbw) boxwidth+=2, box += (src[boxwidth/2]) + (src[-boxwidth/2]);  
    gdouble frac = (fbw - (gdouble) boxwidth)/2.0;
    gint  bh = boxwidth/2, bh1 =boxwidth/2+1;
    
     
    for ( ; src <= end; src++, dest++){
    
      *dest = (box + frac * ((src[bh1])+(src[-bh1])))/fbw;
      box = box - (src[-bh]) + (src[bh1]);
    
    }

}

static void iir_filter(gfloat * const start, gfloat * const end, gfloat * dstart,  gdouble radius, const gint type)
/* Bidirectional IIR-filter, speed optimized */
{
  if (!dstart) dstart = start;
  gfloat *src = start;
  gfloat *dest = dstart;
  gfloat *dend = dstart + (end - start);
  gint width;
  radius = floor((radius+0.1)/0.5)*0.5;
//  gfloat boxwidth = radius * 2.0;
//  gint bw = (gint) boxwidth;
  
  gint ofs = radius;
  if (ofs < 1) ofs=1;
  
  gdouble d1,d2,d3;
  
  width = end-start+1;
         

  if (radius < 0.25){ 
    if ( start != dest ){ 
       memcpy(dest,start,width*sizeof(*dest));
       return;       
    }
  }
  iir_init(radius);
  const gdouble b1 = iir.b1;
  const gdouble b2 = iir.b2/iir.b1;
  const gdouble b3 = iir.b3/iir.b2;
  const gdouble b  = iir.B/iir.b3;
  
  switch(type){
  

#define IIR1(dest,src)  (dest) = (d3 = ((((src) * b + d3) * b3 + d2) * b2 + d1) * b1)
#define IIR2(dest,src)  (dest) = (d2 = ((((src) * b + d2) * b3 + d1) * b2 + d3) * b1)
#define IIR3(dest,src)  (dest) = (d1 = ((((src) * b + d1) * b3 + d3) * b2 + d2) * b1)


  case 0: /*Gauss*/
       d1=d2=d3 = *dest; 

       dend -= 6; 
       src--,dest--;
       while (dest < dend)  {
         IIR1(*(++dest),*(++src));
         IIR2(*(++dest),*(++src));
         IIR3(*(++dest),*(++src));
         IIR1(*(++dest),*(++src));
         IIR2(*(++dest),*(++src));
         IIR3(*(++dest),*(++src));
       }
       dend += 6;
       while (1){
        if (++dest > dend) break; 
        IIR1(*dest,*(++src));
        if (++dest > dend) break; 
        IIR2(*dest,*(++src));
        if (++dest > dend) break; 
        IIR3(*dest,*(++src));
       }

       d1=d2=d3 = dest[-1];
       dstart += 6;
       while (dest > dstart)  {
        --dest, IIR1(*dest,*dest);
        --dest, IIR2(*dest,*dest);
        --dest, IIR3(*dest,*dest);
        --dest, IIR1(*dest,*dest);
        --dest, IIR2(*dest,*dest);
        --dest, IIR3(*dest,*dest);
       }
       dstart -= 6;
       while (1){
        if (--dest < dstart) break; 
        IIR1(*dest,*dest);
        if (--dest < dstart) break; 
        IIR2(*dest,*dest);
        if (--dest < dstart) break; 
        IIR3(*dest,*dest);
       }

       break;
       
   case 2: /* rectified and filtered second derivative, source and dest may be equal */
       
       
       d1=d2=d3 =0.0; 
       dest[0]=dest[ofs]=0.0;
       dend -= 6; 
       dest--; src--;
       while (dest < dend)  {
         ++src, IIR1(*(++dest), src[ofs]-src[0]);
         ++src, IIR2(*(++dest), src[ofs]-src[0]);
         ++src, IIR3(*(++dest), src[ofs]-src[0]);
         ++src, IIR1(*(++dest), src[ofs]-src[0]);
         ++src, IIR2(*(++dest), src[ofs]-src[0]);
         ++src, IIR3(*(++dest), src[ofs]-src[0]);
       }
       dend += 6; 
       while (1)  {
         if (++dest > dend) break; 
         ++src, IIR1(*dest, src[ofs]-src[0]);
         if (++dest > dend) break; 
         ++src, IIR2(*dest, src[ofs]-src[0]);
         if (++dest > dend) break; 
         ++src, IIR3(*dest, src[ofs]-src[0]);
       }
       
       d1=d2=d3 = 0.0;
       dest[-1]=dest[-ofs-1]=0.0;
       dstart += 6;
        
#define IIR1A(dest,src)  (dest) = fabs(d3 = ((((src) * b + d3) * b3 + d2) * b2 + d1) * b1)
#define IIR2A(dest,src)  (dest) = fabs(d2 = ((((src) * b + d2) * b3 + d1) * b2 + d3) * b1)
#define IIR3A(dest,src)  (dest) = fabs(d1 = ((((src) * b + d1) * b3 + d3) * b2 + d2) * b1)

       while (dest > dstart){
        --dest, IIR1A(*dest, dest[0]-dest[-ofs]);
        --dest, IIR2A(*dest, dest[0]-dest[-ofs]);
        --dest, IIR3A(*dest, dest[0]-dest[-ofs]);
        --dest, IIR1A(*dest, dest[0]-dest[-ofs]);
        --dest, IIR2A(*dest, dest[0]-dest[-ofs]);
        --dest, IIR3A(*dest, dest[0]-dest[-ofs]);
        }
        dstart -= 6;
        while (1){
         if (--dest < dstart) break; 
         IIR1A(*dest, dest[0]-dest[-ofs]);
         if (--dest < dstart) break; 
         IIR2A(*dest, dest[0]-dest[-ofs]);
         if (--dest < dstart) break; 
         IIR3A(*dest, dest[0]-dest[-ofs]);
        }
      break;
      }
   
}
   
   
   



//define FR 0.3
//define FG 0.59
//#define FB 0.11
          
#define FR 0.212671
#define FG 0.715160
#define FB 0.072169



static void filter(gfloat *buffer, gfloat *data, gfloat *data2, gfloat *rbuf, gfloat *tbuf, gint width, gint color)
/* 
 * A  forward-backward box filter is used here and the radius is adapted to luminance jump.
 * Radius is calculated fron 1st and 2nd derivative of intensity values.
 * (Its not exactly 2nd derivative, but something similar, optimized by experiment) 
 * The radius variations are filtered. This reduces spatial phase jitter. 
 *
 */

{
  gfloat *lp = data, *rp = data + width-1;
  gfloat *lp2 = data2, *rp2 = data2 +width-1;
  gfloat *blp = buffer, *brp = buffer + width-1;
  gfloat *rbuflp = rbuf, *rbufrp = rbuf + width-1;
  gfloat *p1,*p2;
  gfloat fboxwidth = unsharp_params.radius*2.0;
  gfloat fradius = unsharp_params.radius;
  if (fboxwidth < 1.0) fboxwidth = 1.0 ;
  if (fradius < 0.5) fradius = 0.5;
  gint i;
  gdouble box,lbox,rbox,frac;
  gint boxwidth,ofs,ofs2;
  gfloat maxrad;
  gfloat fbw;
  gfloat val,val2,lval,rval;
  gdouble rfact = sq(unsharp_params.effect);
  gdouble sharp=unsharp_params.sharp;  
  ofs2 = floor(unsharp_params.damping*2.0+0.1);
  
    
  ofs = floor(unsharp_params.lookahead*2.0+0.1);
  gint pass,w=fboxwidth+unsharp_params.damping+unsharp_params.lookahead+unsharp_params.phase + 2.0;
  

    
  for (i=1; i <= w; i++) blp[-i]=blp[i]; /* Mirror image edges */
  for (i=1; i <= w; i++) brp[i] = brp[-i];     
  
  if (color < 0){  
  
    /* Calc 2nd derivative */  
    
    for (p1 = blp,p2=rbuflp;p1<= brp;p1++,p2++){ /* boost high frequency in rbuf */
      *p2 = (sharp+1.0) * p1[0] - sharp * 0.5 * (p1[-ofs]+p1[ofs]);
    }  

    iir_filter(rbuflp-w,rbufrp+w,blp-w,unsharp_params.lookahead,2);
    for (i=1; i <= w; i++) blp[-i]=blp[i]; /* Mirror image edges */
    for (i=1; i <= w; i++) brp[i] = brp[-i];     
    
    for (p1 = blp,p2=rbuflp;p1<= brp;p1++,p2++){ /* boost high frequency in rbuf */
      *p2 = ((sharp+1.0) * (p1[0]) - sharp * 0.5 * ((p1[-ofs2])+(p1[ofs2])));
    }  
  
    for (i=1; i <= w; i++) rbuflp[-i]=rbuflp[i]; /* Mirror rbuf edges */
    for (i=1; i <= w; i++) rbufrp[i] = rbufrp[-i];     
    
    /* Lowpass (gauss) filter rbuf, remove phase jitter */
    iir_filter(rbuflp-w+5,rbufrp+w-5,rbuflp-w+5,unsharp_params.damping,0);

    for (i=-w+5; i< width-1+w-5;i++){
//      val = rbuflp[i];
      val = rbuflp[i]-rfact;
      if (val < rfact/fradius) val=rfact/fradius; /* Avoid division by zero, clip negative filter overshoot */
      val = rfact/val;
     // val = pow(val/fradius,unsharp_params.phase)*fradius;
      if (val < 0.5) val = 0.5;
      
      rbuflp[i] = val*2.0;  
    } 
 
    for (i=1; i <= w; i++) rbuflp[-i]=rbuflp[i]; /* Mirror rbuf edges */
    for (i=1; i <= w; i++) rbufrp[i] = rbufrp[-i];     
    return;
  } /* if color < 0 */
  
    
  
  /* Calc lowpass filtered input signal */
  iir_filter(blp-w+1,brp+w-1,lp2-w+1,unsharp_params.radius,0); 
  
  
  /* Subtract low frequency from input signal (aka original image data) 
   * and predistort this signal
   */
  val = unsharp_params.texture+1.0;  
  for (i = -w+1;i <= width-1+w-1;i++){
     blp[i] = mypow(blp[i]- lp2[i],val);
  }

  
  gfloat *src, *dest;
  val = unsharp_params.texture+1.0;
  
  pass = 2;
  
  while (pass--) {
    
    gint ibw;
    src = blp;
    dest =lp;
    gfloat sum;
    
    maxrad = 0.0;
    
    for (i=1; i <= w; i++) src[-i]=src[i]; /* Mirror left edge */

    sum =  (src[-1] += src[-2]);
  
    /* forward pass */
    for (rbuf = rbuflp-(int) unsharp_params.phase ; rbuf <= rbufrp; src++, dest++, rbuf++){
    
      //fbw = fabs( rbuf[-ofs2]*ll2+rbuf[-ofs2-1]*rl2);
      fbw = *rbuf;
      
      if (fbw > (maxrad += 1.0)) fbw = maxrad;
      else if (fbw < maxrad) maxrad = fbw;
      ibw = (gint) fbw;

      *src = sum += *src;
      *dest = (sum-src[-ibw]+(src[-ibw]-src[-ibw-1])*(fbw-ibw))/fbw;
    }
    
    
    
    src = rp; 
    dest = brp; 
    
    maxrad = 0.0;    
    
    for (i=1; i <= w; i++) src[i] = src[-i];  /* Mirror right edge */   

    sum =   (src[1] += src[2]);
    /* backward pass */
    for ( rbuf = rbufrp +(int) unsharp_params.phase ; rbuf >= rbuflp; src--, dest--, rbuf--){
    
      //fbw = fabs( rbuf[ofs2]*ll2+rbuf[ofs2+1]*rl2);
      fbw = *rbuf;
      
      if (fbw > (maxrad +=1.0)) fbw = maxrad;
      else if (fbw < maxrad) maxrad = fbw;
      
      ibw = (gint) fbw;

      *src = sum += *src;
      *dest = (sum-src[ibw]+(src[ibw]-src[ibw+1])*(fbw-ibw))/fbw;

    }
    
    
  } /* Next pass */
    
    val = 1.0/(unsharp_params.texture+1.0);
    for (i = -w+1;i <= width-1+w-1;i++){
      /* Undo  predistortion */
      
      blp[i]= mypow(blp[i],val);
      blp[i] += lp2[i]; /* Add in low frequency */

//      if (blp[i] >= 0.0) blp[i] = pow(blp[i],val);
//      else blp[i] = 0.0;      
      
    }
 
}

/* This function is written as if it is blurring a column at a time,
 * even though it can operate on rows, too.  There is no difference
 * in the processing of the lines, at least to the blur_line function.
 */
static void
blur_line (gfloat * const data,
           gfloat * const data2,
           gfloat * const buffer,
           gfloat * rbuf,
           gfloat * tbuf,
           const guchar  *src,
           guchar        *dest,
           gint           len,   /* length of src and dest */
           glong          bytes) /* Bits per plane */  
{
  gfloat scale;
  gfloat sum;
  gint    b,i = 0;
  gint    j = 0;
  gint    row;
  gint colors = 3; 
  if (bytes < 3) colors = 1;

     gint idx;
     /* Calculate radius factors */
     if (colors < 3){
         for (row =0, idx=0 ; idx < len; row +=bytes, idx++) data[idx] = lut[dest[row]]; 
         filter(data, data2, buffer,rbuf,tbuf, len, -1);
     }
     else
     {
         for (row =0, idx=0 ; idx < len; row +=bytes, idx++){
            /* Color weigths are choosen proportional to Bayer Sensor pixel count */             
            data[idx] = (gfloat) dest[row] / 255.0 * 0.25; /* Red color */ 
            data[idx] +=  (gfloat) dest[row+1] / 255.0 * 0.5; /* Green color */
            data[idx] +=  (gfloat) dest[row+2] / 255.0 * 0.25;  /* Blue color */
            data[idx] = mypow(data[idx],lut_gamma);
         }
         filter(data, data2, buffer,rbuf,tbuf, len, -1);
     }
     
     /* Do actual filtering */
     for (b = 0; b<colors; b++){
       for (row =b, idx=0 ; idx < len; row +=bytes, idx++) data[idx] = (gfloat) src[row]/255.0; 
        
       filter(data, data2, buffer,rbuf,tbuf, len, b);
     
       for (row =b, idx=0; idx < len; row +=bytes, idx++){
          gint value = data[idx]*255.0+0.5;
          dest[row] = CLAMP( value, 0, 255); 
       }
    }
  
  
}

static void
unsharp_mask (GimpDrawable *drawable,
              gdouble       radius,
              gdouble       lsmooth)
{
  GimpPixelRgn srcPR, destPR;
  gint         x1, y1, x2, y2;

  /* initialize pixel regions */
  gimp_pixel_rgn_init (&srcPR, drawable,
                       0, 0, drawable->width, drawable->height, FALSE, FALSE);
  gimp_pixel_rgn_init (&destPR, drawable,
                       0, 0, drawable->width, drawable->height, TRUE, TRUE);

  /* Get the input */
  gimp_drawable_mask_bounds (drawable->drawable_id, &x1, &y1, &x2, &y2);

  unsharp_region (&srcPR, &destPR, drawable->bpp,
                  radius, lsmooth,
                  x1, x2, y1, y2,
                  TRUE);

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id, x1, y1, x2 - x1, y2 - y1);
}


/* Remove noise on the region, given a source region, dest.
 * region, width and height of the regions, and corner coordinates of
 * a subregion to act upon.  Everything outside the subregion is unaffected.
 */
static void
unsharp_region (GimpPixelRgn *srcPR,
                GimpPixelRgn *destPR,
                gint          bytes, /* Bytes per pixel */
                gdouble       radius,
                gdouble       lsmooth,
                gint          x1,    /* Corners of subregion */
                gint          x2,
                gint          y1,
                gint          y2,
                gboolean      show_progress)
{
  guchar  *src;
  guchar  *dest;
  gint     width   = x2 - x1;
  gint     height  = y2 - y1;
  gfloat *data, *data2;
  gfloat *buffer=NULL;
  gint     row, col,i;
  gfloat prob = 0.0;
  gint w = (radius+unsharp_params.lookahead+unsharp_params.damping+unsharp_params.phase) * 4.0 + 40.0;
//  if (radius < unsharp_params.lookahead) w = unsharp_params.lookahead * 4.0 + 40.0;
  
    
  lut_init(unsharp_params.gamma);
  
  gfloat     csmooth =  unsharp_params.csmooth;
  if (csmooth >= 0.99) csmooth = 1.0; /* Raw Filter preview */  
     
  if (show_progress)
    gimp_progress_init (_("Blurring..."));

  
  /* allocate buffers */
  src  = g_new (guchar, MAX (width, height) * bytes);
  dest = g_new (guchar, MAX (width, height) * bytes);
  data = g_new(gfloat, MAX (width, height) + 2*w);
  data2 = g_new(gfloat, MAX (width, height) + 2*w);
  buffer = g_new(gfloat, MAX (width, height) + 2*w);
  gfloat *rbuf = g_new(gfloat,MAX (width, height) + 2*w);
  gfloat *tbuf = g_new(gfloat,MAX (width, height) + 2*w);
  
  for (i=0;i<MAX(width,height)+2*w-1;i++) data[i]=data2[i]=buffer[i]=rbuf[i]=tbuf[i]=0.0;

  /* Initialize the damping filter coefficients */  
  iir_init(unsharp_params.radius);
  
  /* blur the rows */
  for (row = 0; row < height; row++)
    {
      gimp_pixel_rgn_get_row (srcPR, src, x1, y1 + row, width);
      memcpy(dest,src,width*bytes);
//      gimp_pixel_rgn_get_row (srcPR, dest, x1, y1 + row, width);
      blur_line (data+w, data2+w, buffer+w, rbuf+w, tbuf+w, src, dest, width, bytes);
      gimp_pixel_rgn_set_row (destPR, dest, x1, y1 + row, width);

      if (show_progress && row % 8 == 0)
        gimp_progress_update ((gdouble) row / (3 * height));
    }

  /* blur the cols */

  for (col = 0; col < width; col++)
    {
      gimp_pixel_rgn_get_col (destPR, src, x1 + col, y1, height);
      gimp_pixel_rgn_get_col (srcPR, dest, x1 + col, y1, height);
      blur_line (data+w, data2+w, buffer+w, rbuf+w, tbuf+w, src, dest, height, bytes);
      gimp_pixel_rgn_set_col (destPR, dest, x1 + col, y1, height);

      if (show_progress && col % 8 == 0)
        gimp_progress_update ((gdouble) col / (3 * width) + 0.33);
    }

  if (show_progress)
    gimp_progress_init (_("Merging..."));

#define MERGE
#ifdef MERGE
  /* merge the source and destination (which currently contains
     the blurred version) images */
  for (row = 0; row < height; row++)
    {
      guchar *s = src;
      guchar       *d = dest;
      gfloat        value;
      gint          u, v;

      /* get source row */
      gimp_pixel_rgn_get_row (srcPR, src, x1, y1 + row, width);

      gimp_pixel_rgn_get_row (destPR, dest, x1, y1 + row, width);

      /* get dest row */
      /* combine the two */
      
      gfloat t=csmooth;
      gfloat t2 = lsmooth;
      
      t*=t;  /* Easier adjustment for small values */
      t2*=t2;
              
      for (u = 0; u < width; u++)
          {
          gfloat dpix[3],spix[3];      

          gfloat lum,red,green,blue;
          gfloat lum2,red2,green2,blue2; 
          
          red = (gfloat) s[0]/255.0;
          if (bytes > 2) 
            {
            green = (gfloat) s[1]/255.0;
            blue =  (gfloat) s[2]/255.0;
            }
          else  green = blue = red;

          spix[0] = red;
          spix[1] = green;
          spix[2] = blue;


                            
          lum = (FR*red + FG*green + FB*blue);
          
          red2 = (gfloat) d[0]/255.0;
          if (bytes > 2) 
            {
            green2 = (gfloat) d[1]/255.0;
            blue2 =  (gfloat) d[2]/255.0;
            }
          else  green2 = blue2 = red2;
          
          
                  
          lum2 = (FR*red2 + FG*green2 + FB*blue2);
          /*
           * Calculate luminance error (contrast error) for filtered template.
           * This error is biggest, where edges are. Edges anyway cannot be filtered.
           * Therefore we can correct luminance error in edges without increasing noise.
           * Should be adjusted carefully, or not so carefully if you intentionally want to add noise.
           * Noise, if not colorized, /can/ look good, so this makes sense.
           */ 
          gfloat dl = lum - lum2;            
            
          /* Multiply dl with first derivative of gamma curve divided by derivative value for midtone 0.5 
           * So bright tones will be corrected more (get more luminance noise and -information) than darker values 
           * Because bright parts of image generally are less noisy, this is what we want.
           */  
          dl *= pow(lum2/0.5,lut_gamma-1.0); 
  
          if (t2 >= 0.0) dl *= (1.0 - exp(-dl*dl/(2.0*t2*t2)));        
          
   //       if (dl > p) dl = p;
   //       if (dl < -p) dl = -p;
          
          dpix[0] =   red2 + dl;
          dpix[1] = green2 + dl;
          dpix[2] =  blue2 + dl;
          
          for (v = 0; v < bytes; v++)
            {
       
            gfloat value = spix[v];
            gfloat fvalue = dpix[v];
            gfloat mvalue = (value + fvalue)/2.0;
            
            gfloat diff = (value) - (fvalue);
             
            /* Multiply diff with first derivative of gamma curve divided by derivative value for midtone 0.5 
             * So midtones will stay unchanged, darker values get more blur and brighter values get less blur
             * when we increase gamma.
             */  
            diff *= pow(mvalue/0.5,lut_gamma-1.0); 

            /* Calculate noise probability for pixel                 
             * Ok, probably it's not probability but an
             * arbitrary curve ;-)
             * Probably we should provide a GUI-interface for this
             */
            
             
            if (t > 0.0) prob = exp(-diff*diff/(2.0*t*t));
            else prob = 0.0;
            
            
            if (t >= 0.99) prob = 1.0; /* Allow viewing of raw filter output */
            dpix[v] = value = fvalue * prob + value * (1.0 - prob);    
            
            }

          
          value = dpix[0]*255.0+0.5;
          d[0] = CLAMP(value,0,255);
          if (bytes > 2)
             {
             value = dpix[1]*255.0+0.5;
             d[1] = CLAMP(value,0,255);
             value = dpix[2]*255.0+0.5;
             d[2] = CLAMP(value,0,255);
             } 
          d += bytes;
          s +=bytes;
          }

      if (show_progress && row % 8 == 0)
        gimp_progress_update ((gdouble) row / (3 * height) + 0.67);

      gimp_pixel_rgn_set_row (destPR, dest, x1, y1 + row, width);
    }
#endif
  if (show_progress)
    gimp_progress_update (0.0);
  g_free (data);
  g_free (data2);
  g_free (buffer);
  g_free (rbuf);
  g_free (tbuf);
  g_free (dest);
  g_free (src);
}

static gboolean
unsharp_mask_dialog (GimpDrawable *drawable)
{
  GtkWidget *dialog;
  GtkWidget *main_vbox;
  GtkWidget *preview;
  GtkWidget *table;
  GtkObject *adj;
  gboolean   run;

  #define MAX_RAD 0
  #define F_THRESH 1
  #define TEXTURE 2
  #define SHARP 3
  #define CLIP 4
  #define L_AHEAD 5
  #define PHASE 6
  #define GAMMA 7
  #define LUM_TOL 8
  #define COL_TOL 9
  
  
  gimp_ui_init ("unsharp", TRUE);

  dialog = gimp_dialog_new (_("Dcam Noise 2 V "PLUG_IN_VERSION), "dcamnoise2-"PLUG_IN_VERSION,
                            NULL, 0,
                            gimp_standard_help_func, "plug-in-dcamnoise2-"PLUG_IN_VERSION,

                            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                            GTK_STOCK_OK,     GTK_RESPONSE_OK,

                            NULL);

  main_vbox = gtk_vbox_new (FALSE, 12);
  gtk_container_set_border_width (GTK_CONTAINER (main_vbox), 12);
  gtk_container_add (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), main_vbox);
  gtk_widget_show (main_vbox);

  preview = gimp_drawable_preview_new (drawable,
                                       &unsharp_params.update_preview);
  gtk_box_pack_start (GTK_BOX (main_vbox), preview, TRUE, TRUE, 0);
  gtk_widget_show (preview);

  g_signal_connect (preview, "invalidated",
                    G_CALLBACK (preview_update),
                    NULL);

  table = gtk_table_new (3, 3, FALSE);
  gtk_table_set_col_spacings (GTK_TABLE (table), 6);
  gtk_table_set_row_spacings (GTK_TABLE (table), 6);
  gtk_box_pack_start (GTK_BOX (main_vbox), table, FALSE, FALSE, 0);
  gtk_widget_show (table);

  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, MAX_RAD,
                              _("_Filter max. Radius:"), SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.radius, 0.5, 50.0, 0.5, 0.5, 1,
                              TRUE, 0, 0,
                              NULL, NULL);
 // gimp_scale_entry_set_logarithmic(adj,TRUE);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.radius);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);

  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, LUM_TOL,
                              _("_Luminance tolerance:"), SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.lsmooth, 0.0, 1.0, 0.1, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);
  
  gimp_scale_entry_set_logarithmic(adj,TRUE);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.lsmooth);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);

  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, COL_TOL,
                              _("_Color tolerance:"), SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.csmooth,
                              0.0, 1.0, 0.1, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);

  gimp_scale_entry_set_logarithmic(adj,TRUE);
  
  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.csmooth);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);


  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, F_THRESH,
                              "_Filter Threshold:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.effect,
                              0.0, 1.0, 0.1, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);
  gimp_scale_entry_set_logarithmic(adj,TRUE);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.effect);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);


  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, L_AHEAD,
                              "_Lookahead:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.lookahead,
                              0.5, 20.0, 0.5, 0.5, 1,
                              TRUE, 0, 0,
                              NULL, NULL);
//gimp_scale_entry_set_logarithmic(adj,TRUE);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.lookahead);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);


  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, GAMMA,
                              "_Gamma:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.gamma,
                              0.02, 5.0, 0.01, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.gamma);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);


  
  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, PHASE,
                              "_Phase Jitter damping:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.damping,
                              0.5, 20.0, 0.5, 0.5, 1,
                              TRUE, 0, 0,
                              NULL, NULL);
  
  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.damping);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);
#if 1
  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, CLIP,
                              "_Erosion:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.phase,
                              0.0, 10.0, 1.0, 1.0, 0,
                              TRUE, 0, 0,
                              NULL, NULL);
  
  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.phase);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);
//  gimp_scale_entry_set_logarithmic(adj,TRUE);
#endif
  
  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, TEXTURE,
                              "_Texture Detail:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.texture,
                              -1.0, 1.0, 0.01, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);
  
  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.texture);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);

  
  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, SHARP,
                              "_Sharpness:", SCALE_WIDTH, ENTRY_WIDTH,
                              unsharp_params.sharp,
                              0.0, 1.0, 0.01, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);
  
  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &unsharp_params.sharp);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);
 // gimp_scale_entry_set_logarithmic(adj,TRUE);



  gtk_widget_show (dialog);

  run = (gimp_dialog_run (GIMP_DIALOG (dialog)) == GTK_RESPONSE_OK);

  gtk_widget_destroy (dialog);

  return run;
}

static void
preview_update (GimpPreview *preview)
{
  GimpDrawable *drawable;
  gint          x1, x2;
  gint          y1, y2;
  gint          x, y;
  gint          width, height;
  gint          border;
  GimpPixelRgn  srcPR;
  GimpPixelRgn  destPR;

  
  drawable =
    gimp_drawable_preview_get_drawable (GIMP_DRAWABLE_PREVIEW (preview));

  gimp_pixel_rgn_init (&srcPR, drawable,
                       0, 0, drawable->width, drawable->height, FALSE, FALSE);
  gimp_pixel_rgn_init (&destPR, drawable,
                       0, 0, drawable->width, drawable->height, TRUE, TRUE);

  gimp_preview_get_position (preview, &x, &y);
  gimp_preview_get_size (preview, &width, &height);

  /* enlarge the region to avoid artefacts at the edges of the preview */
  border = 2.0 * unsharp_params.radius + 0.5;
  
  x1 = MAX (0, x - border);
  y1 = MAX (0, y - border);
  x2 = MIN (x + width  + border, drawable->width);
  y2 = MIN (y + height + border, drawable->height);

  unsharp_region (&srcPR, &destPR, drawable->bpp,
                  unsharp_params.radius, unsharp_params.lsmooth,
                  x1, x2, y1, y2,
                  FALSE);

  gimp_pixel_rgn_init (&destPR, drawable, x, y, width, height, FALSE, TRUE);
  gimp_drawable_preview_draw_region (GIMP_DRAWABLE_PREVIEW (preview), &destPR);
 
}

