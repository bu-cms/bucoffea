# Utility functions for dealing with ROOT plots
import ROOT as r
import string
import random

def create_tdr_style(want_title=False):
   tdrStyle =r.TStyle("tdrStyle","Style for P-TDR");

   tdrStyle.SetCanvasBorderMode(0);
   tdrStyle.SetCanvasColor(r.kWhite);
   tdrStyle.SetCanvasDefH(600); #Height of canvas
   tdrStyle.SetCanvasDefW(600); #Width of canvas
   tdrStyle.SetCanvasDefX(0);   #POsition on screen
   tdrStyle.SetCanvasDefY(0);

   # For the Pad:
   tdrStyle.SetPadBorderMode(0);
   tdrStyle.SetPadColor(r.kWhite);
   tdrStyle.SetPadGridX(False);
   tdrStyle.SetPadGridY(False);
   tdrStyle.SetGridColor(0);
   tdrStyle.SetGridStyle(3);
   tdrStyle.SetGridWidth(1);

   # For the frame:
   tdrStyle.SetFrameBorderMode(0);
   tdrStyle.SetFrameBorderSize(1);
   tdrStyle.SetFrameFillColor(0);
   tdrStyle.SetFrameFillStyle(0);
   tdrStyle.SetFrameLineColor(1);
   tdrStyle.SetFrameLineStyle(1);
   tdrStyle.SetFrameLineWidth(1);

   # For the histo:
   tdrStyle.SetHistLineColor(1);
   tdrStyle.SetHistLineStyle(0);
   tdrStyle.SetHistLineWidth(2);
   tdrStyle.SetErrorX(0.5);


   #For the fit/function:
   tdrStyle.SetOptFit(0);
   tdrStyle.SetFitFormat("5.4g");
   tdrStyle.SetFuncColor(2);
   tdrStyle.SetFuncStyle(1);
   tdrStyle.SetFuncWidth(1);

   #For the date:
   tdrStyle.SetOptDate(0);


   # For the statistics box:
   tdrStyle.SetOptFile(0);
   tdrStyle.SetOptStat(""); # To display the mean and RMS:   SetOptStat("mr");
   tdrStyle.SetStatColor(r.kWhite);
   tdrStyle.SetStatFont(42);
   tdrStyle.SetStatFontSize(0.04);
   tdrStyle.SetStatFormat("6.4g");
   tdrStyle.SetStatBorderSize(0);


   # Margins:
   tdrStyle.SetPadTopMargin(0.07);
   tdrStyle.SetPadBottomMargin(0.13);
   tdrStyle.SetPadLeftMargin(0.13);
   tdrStyle.SetPadRightMargin(0.06);

   # For the Global title:
   if( want_title ):
      tdrStyle.SetOptTitle(1);
      tdrStyle.SetTitleFont(42);
      tdrStyle.SetTitleColor(1);
      tdrStyle.SetTitleTextColor(r.kBlack);
      tdrStyle.SetTitleFillColor(r.kWhite);
      tdrStyle.SetTitleFontSize(0.14);
      tdrStyle.SetTitleBorderSize(0);
      tdrStyle.SetTitleAlign(23);
      tdrStyle.SetTitleH(0.1); # Set the height of the title box
      tdrStyle.SetTitleW(0.8); # Set the height of the title box
      tdrStyle.SetTitleX(0.5); # Set the position of the title box
      tdrStyle.SetTitleY(1.); # Set the position of the title box
      tdrStyle.SetPadTopMargin(0.15);

   else:
      tdrStyle.SetOptTitle(0);
      tdrStyle.SetTitleFont(42);
      tdrStyle.SetTitleColor(1);
      tdrStyle.SetTitleTextColor(1);
      tdrStyle.SetTitleFillColor(10);
      tdrStyle.SetTitleFontSize(0.05);
      tdrStyle.SetTitleH(0.04); # Set the height of the title box
      tdrStyle.SetTitleX(0.1); # Set the position of the title box
      tdrStyle.SetTitleY(1.); # Set the position of the title box

   # For the axis titles:

   tdrStyle.SetTitleColor(1, "XYZ");
   tdrStyle.SetTitleFont(42, "XYZ");
   tdrStyle.SetTitleSize(0.05, "XYZ");

   # For the axis labels:

   tdrStyle.SetLabelColor(1, "XYZ");
   tdrStyle.SetLabelFont(42, "XYZ");
   tdrStyle.SetLabelOffset(0.007, "XYZ");
   tdrStyle.SetLabelSize(0.04, "XYZ");

   # For the axis:

   tdrStyle.SetAxisColor(1, "XYZ");
   tdrStyle.SetStripDecimals(r.kTRUE);
   tdrStyle.SetTickLength(0.03, "XYZ");
   tdrStyle.SetNdivisions(510, "XYZ");
   tdrStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame
   tdrStyle.SetPadTickY(1);

   # Change for log plots:
   tdrStyle.SetOptLogx(0);
   tdrStyle.SetOptLogy(0);
   tdrStyle.SetOptLogz(0);

   # Postscript options:
   tdrStyle.SetPaperSize(20.,20.);

   tdrStyle.SetLegendBorderSize(0)
   tdrStyle.cd();
   return tdrStyle;


def setup_canvas( want_ratio = False, width_x=None, width_y=None ):
   ID="".join(random.sample(string.ascii_uppercase+string.digits,8))

   if(want_ratio):
      c1 = r.TCanvas( 'c1'+ID,'c1'+ID, 200, 10, width_x if width_x is not None else 500, width_y if width_y is not None else 500 )

      t1 = r.TPad("t1"+ID,"t1"+ID, 0.0, 0.3, 1.0, 1.0)
      t2 = r.TPad("t2"+ID,"t2"+ID, 0.0, 0.0, 1.0, 0.3)
      t1.SetLeftMargin(0.2)
      t1.SetBottomMargin(0.065)
      t2.SetBottomMargin(0.35)
      t2.SetLeftMargin(0.2)
      t2.SetTopMargin(0.05)
      c1.cd()
      t1.Draw()
      t2.Draw()
      return c1,t1,t2
   else:
      c1 = r.TCanvas( 'c1'+ID,'c1'+ID, 200, 10, width_x if width_x is not None else 700, width_y if width_y is not None else 500 )
      return c1


def apply_style_to_axis( histogram, is_ratio, ymin=None, ymax=None, xtitle="", ytitle="" ):
   x_axis = histogram.GetXaxis()
   y_axis = histogram.GetYaxis()

   title_size = 0.07
   label_size = 0.06
   title_offset = 1.0
   tick_length = 0.03
   if( is_ratio ):
      yscale = (1.0-0.3)/(0.3-0);
      title_size *= yscale
      label_size *= yscale
      tick_length *= yscale

      if (ymin and ymax):
            y_axis.SetRangeUser(ymin, ymax)
            histogram.SetMinimum(ymin)
            histogram.SetMaximum(ymax)
      else:
            y_axis.SetRangeUser(0.3,1.8)
      if(len(xtitle)): x_axis.SetTitle(xtitle)
      if(len(ytitle)): y_axis.SetTitle(ytitle) # Whitespace is important so that centering does not cause the title to be otuside the pad
      y_axis.CenterTitle(True)
      y_axis.SetNdivisions(5)

      y_axis.SetTitleOffset(0.4)
   else:
      if(ymin and ymax):
         histogram.SetMinimum(ymin)
         histogram.SetMaximum(ymax)
      if(len(ytitle)): y_axis.SetTitle(ytitle)
      if(len(ytitle)): y_axis.SetTitleOffset(title_offset)

   x_axis.SetTitleOffset(title_offset)
   x_axis.SetTickLength(tick_length)
   x_axis.SetNdivisions(505,True)
   for axis in [x_axis, y_axis]:
      axis.SetTitleSize(title_size)
      axis.SetLabelSize(label_size)
