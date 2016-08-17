#!MC 1400
# Created by Tecplot 360 build 14.0.1.26249
$!VarSet |FILEIND|    = 'FIGIFILEINDEX'
$!VarSet |PATH|       = 'FIGIPATH'
$!VarSet |MYRE|       = 'FIGIMYRE'
# Polydeg
$!VarSet |MYN1|       = 'FIGIMYN1'
$!VarSet |MYN2|       = 'FIGIMYN2'
$!VarSet |MYN3|       = 'FIGIMYN3'
# Path of datasets
$!VarSet |INPATH1|    = 'FIGIINPATH1'
$!VarSet |INPATH2|    = 'FIGIINPATH2'
$!VarSet |INPATH3|    = 'FIGIINPATH3'
# Path of references
$!VarSet |REFPATH|    = 'FIGIREFPATH'
$!VarSet |MYTYPE|     = 'FIGITYPE'
$!VarSet |MYRANGEMIN| = FIGIRANGEMIN
$!VarSet |MYRANGEMAX| = FIGIRANGEMAX
$!VarSet |MYLEGENDX|  = FIGILEGENDX
$!VarSet |MYLEGENDY|  = FIGILEGENDY
$!VarSet |MYXPOS|     = 'FIGIXPOS'
$!VarSet |MYVARNAME|  = 'FIGIVARNAME'
$!VarSet |MYTPVARNAME|= 'FIGITPVARNAME'
$!DRAWGRAPHICS FALSE
$!VarSet |LFDSFN1| = '"|REFPATH|/phill_ref_lesocc_|MYRE|_|FILEIND|.dat" "|REFPATH|/phill_ref_mglet_|MYRE|_|FILEIND|.dat" "|INPATH1|/phill_|MYRE|_|MYN1|_|MYTYPE|_|FILEIND|.dat" "|INPATH2|/phill_|MYRE|_|MYN2|_|MYTYPE|_|FILEIND|.dat" "|INPATH3|/phill_|MYRE|_|MYN3|_|MYTYPE|_|FILEIND|.dat"'
#$!VarSet |LFDSFN1| = '"|REFPATH|/phill_ref_lesocc_|MYRE|_|FILEIND|.dat" "|REFPATH|/phill_ref_mglet_|MYRE|_|FILEIND|.dat" "|INPATH1|/phill_|MYRE|_|MYN1|_|MYTYPE|_|FILEIND|.dat"'
$!VarSet |LFDSVL1| = '"CoordinateY" "|MYTPVARNAME|"'
$!SETSTYLEBASE FACTORY
$!GLOBALLINKING 
  LINKCOLORMAPS = YES
$!GLOBALCOLORMAP  1
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL 1 RESETTOFACTORY
$!GLOBALCOLORMAP  1
  SMRAINBOW
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.25
      LEADRGB
        {
        R = 0
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 0
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.75
      LEADRGB
        {
        R = 255
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    }
$!GLOBALCOLORMAP  2
  CONTOURCOLORMAP = LGRAINBOW
$!COLORMAPCONTROL 2 RESETTOFACTORY
$!GLOBALCOLORMAP  2
  LGRAINBOW
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.1667
      LEADRGB
        {
        R = 0
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.3333
      LEADRGB
        {
        R = 0
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 255
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 0.6667
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    CONTROLPOINT 6
      {
      COLORMAPFRACTION = 0.8333
      LEADRGB
        {
        R = 255
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 7
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 255
        }
      }
    }
$!GLOBALCOLORMAP  3
  CONTOURCOLORMAP = MODERN
$!COLORMAPCONTROL 3 RESETTOFACTORY
$!GLOBALCOLORMAP  3
  MODERN
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 100
        G = 35
        B = 100
        }
      TRAILRGB
        {
        R = 100
        G = 35
        B = 100
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.1429
      LEADRGB
        {
        R = 35
        G = 35
        B = 100
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.2857
      LEADRGB
        {
        R = 35
        G = 100
        B = 100
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.4286
      LEADRGB
        {
        R = 35
        G = 100
        B = 35
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 0.5714
      LEADRGB
        {
        R = 100
        G = 100
        B = 35
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 6
      {
      COLORMAPFRACTION = 0.7143
      LEADRGB
        {
        R = 85
        G = 46
        B = 10
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 7
      {
      COLORMAPFRACTION = 0.8571
      LEADRGB
        {
        R = 100
        G = 35
        B = 35
        }
      TRAILRGB
        {
        R = 217
        G = 117
        B = 26
        }
      }
    CONTROLPOINT 8
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    }
$!GLOBALCOLORMAP  4
  CONTOURCOLORMAP = GRAYSCALE
$!COLORMAPCONTROL 4 RESETTOFACTORY
$!GLOBALCOLORMAP  4
  GRAYSCALE
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 0
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.125
      LEADRGB
        {
        R = 32
        G = 32
        B = 32
        }
      TRAILRGB
        {
        R = 32
        G = 32
        B = 32
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.25
      LEADRGB
        {
        R = 64
        G = 64
        B = 64
        }
      TRAILRGB
        {
        R = 64
        G = 64
        B = 64
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.375
      LEADRGB
        {
        R = 96
        G = 96
        B = 96
        }
      TRAILRGB
        {
        R = 96
        G = 96
        B = 96
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 128
        G = 128
        B = 128
        }
      TRAILRGB
        {
        R = 128
        G = 128
        B = 128
        }
      }
    CONTROLPOINT 6
      {
      COLORMAPFRACTION = 0.625
      LEADRGB
        {
        R = 160
        G = 160
        B = 160
        }
      TRAILRGB
        {
        R = 160
        G = 160
        B = 160
        }
      }
    CONTROLPOINT 7
      {
      COLORMAPFRACTION = 0.75
      LEADRGB
        {
        R = 192
        G = 192
        B = 192
        }
      TRAILRGB
        {
        R = 192
        G = 192
        B = 192
        }
      }
    CONTROLPOINT 8
      {
      COLORMAPFRACTION = 0.875
      LEADRGB
        {
        R = 224
        G = 224
        B = 224
        }
      TRAILRGB
        {
        R = 224
        G = 224
        B = 224
        }
      }
    CONTROLPOINT 9
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 255
        }
      }
    }
$!GLOBALCOLORMAP  5
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL 5 RESETTOFACTORY
$!GLOBALCOLORMAP  5
  SMRAINBOW
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.25
      LEADRGB
        {
        R = 0
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 0
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.75
      LEADRGB
        {
        R = 255
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    }
$!GLOBALCOLORMAP  6
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL 6 RESETTOFACTORY
$!GLOBALCOLORMAP  6
  SMRAINBOW
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.25
      LEADRGB
        {
        R = 0
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 0
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.75
      LEADRGB
        {
        R = 255
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    }
$!GLOBALCOLORMAP  7
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL 7 RESETTOFACTORY
$!GLOBALCOLORMAP  7
  SMRAINBOW
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.25
      LEADRGB
        {
        R = 0
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 0
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.75
      LEADRGB
        {
        R = 255
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    }
$!GLOBALCOLORMAP  8
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL 8 RESETTOFACTORY
$!GLOBALCOLORMAP  8
  SMRAINBOW
    {
    CONTROLPOINT 1
      {
      COLORMAPFRACTION = 0
      LEADRGB
        {
        R = 0
        G = 0
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 0
        B = 255
        }
      }
    CONTROLPOINT 2
      {
      COLORMAPFRACTION = 0.25
      LEADRGB
        {
        R = 0
        G = 255
        B = 255
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 255
        }
      }
    CONTROLPOINT 3
      {
      COLORMAPFRACTION = 0.5
      LEADRGB
        {
        R = 0
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 0
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 4
      {
      COLORMAPFRACTION = 0.75
      LEADRGB
        {
        R = 255
        G = 255
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 255
        B = 0
        }
      }
    CONTROLPOINT 5
      {
      COLORMAPFRACTION = 1
      LEADRGB
        {
        R = 255
        G = 0
        B = 0
        }
      TRAILRGB
        {
        R = 255
        G = 0
        B = 0
        }
      }
    }
$!PLOTOPTIONS 
  SUBDIVIDEALLCELLS = NO
$!GLOBALPAPER 
  PAPERSIZEINFO
    {
    LETTER
      {
      WIDTH = 8.5
      HEIGHT = 11
      LEFTHARDCLIPOFFSET = 0.125
      RIGHTHARDCLIPOFFSET = 0.125
      TOPHARDCLIPOFFSET = 0.125
      BOTTOMHARDCLIPOFFSET = 0.125
      }
    }
$!PAGE 
  NAME = ''
  PAPERATTRIBUTES
    {
    BACKGROUNDCOLOR = WHITE
    ISTRANSPARENT = YES
    ORIENTPORTRAIT = NO
    SHOWGRID = YES
    SHOWRULER = NO
    SHOWPAPER = NO
    PAPERSIZE = LETTER
    RULERSPACING = ONEINCH
    PAPERGRIDSPACING = HALFINCH
    REGIONINWORKAREA
      {
      X1 = 1
      Y1 = 0.25
      X2 = 10
      Y2 = 8.25
      }
    }
### Frame Number 1 ###
$!READDATASET  '|LFDSFN1|'
  INITIALPLOTTYPE = XYLINE
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  ASSIGNSTRANDIDS = YES
  VARLOADMODE = BYNAME
  VARNAMELIST = '|LFDSVL1|'
$!REMOVEVAR |LFDSVL1|
$!REMOVEVAR |LFDSFN1|

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Festlegen der Zonen
$!RENAMEDATASETZONE 
  ZONE = 1
  NAME = 'LESOCC'
$!RENAMEDATASETZONE
  ZONE = 2
  NAME = 'MGLET'
$!RENAMEDATASETZONE
  ZONE = 3
  NAME = 'DG-|MYN1|'
$!RENAMEDATASETZONE
  ZONE = 4
  NAME = 'DG-|MYN2|'
$!RENAMEDATASETZONE
  ZONE = 5
  NAME = 'DG-|MYN3|'
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

$!FRAMELAYOUT 
  SHOWHEADER = NO
  HEADERCOLOR = RED
  XYPOS
    {
    X = 1
    Y = 0.25
    }
  WIDTH = 9
  HEIGHT = 8
$!THREEDAXIS 
  ASPECTRATIOLIMIT = 25
  BOXASPECTRATIOLIMIT = 25
$!PLOTTYPE  = XYLINE
$!FRAMENAME  = 'Frame 001'
$!GLOBALTIME 
  SOLUTIONTIME = 0
$!DELETELINEMAPS 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
$!ACTIVELINEMAPS  =  [1-5]   #! TODO: increase number of zones for additional plots
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
$!GLOBALLINEPLOT 
  DATALABELS
    {
    DISTANCESKIP = 5
    }
  LEGEND
    {
    SHOW = YES
    TEXTSHAPE
      {
      HEIGHT = 2.5
      }
    BOX
      {
      BOXTYPE = NONE
      }
    XYPOS
      {
      X = |MYLEGENDX|
      Y = |MYLEGENDY|
      }
    }
# Plot LESOCC and MGLET as dots and quads
# LESOCC
$!LINEMAP  [1]
  NAME = '&ZN&'
  ASSIGN
    {
    ZONE = 1
    XAXISVAR = 2
    YAXISVAR = 1
    }
  LINES
    {
    SHOW = NO
    COLOR = BLACK
    LINETHICKNESS = 0.3
    }
  SYMBOLS
    {
    SHOW = YES
    SYMBOLSHAPE
      {
      GEOMSHAPE = CIRCLE
      }
    COLOR = BLACK
    FILLCOLOR = BLACK
    SIZE = 1.5
    LINETHICKNESS = 0.2
    SKIPMODE = BYFRAMEUNITS
    SKIPPING = 6
    }
  BARCHARTS
    {
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  ERRORBARS
    {
    COLOR = BLACK
    }
# MGLET
$!LINEMAP  [2]
  NAME = '&ZN&'
  ASSIGN
    {
    ZONE = 2
    XAXISVAR = 2
    YAXISVAR = 1
    }
  LINES
    {
    SHOW = NO
    COLOR = BLACK
    LINETHICKNESS = 0.3
    }
  SYMBOLS
    {
    SHOW = YES
    SYMBOLSHAPE
      {
      GEOMSHAPE = DEL
      }
    COLOR = BLACK
    FILLCOLOR = BLACK
    SIZE = 1.5
    LINETHICKNESS = 0.2
    SKIPMODE = BYFRAMEUNITS
    SKIPPING = 6
    }
  BARCHARTS
    {
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  ERRORBARS
    {
    COLOR = BLACK
    }
# Plot Flexi results as lines
#FLEXI N4
$!LINEMAP  [3]
  NAME = '&ZN&'
  ASSIGN
    {
    ZONE = 3
    XAXISVAR = 2
    YAXISVAR = 1
    }
  LINES
    {
    COLOR = BLACK
    LINEPATTERN = DASHDOT
    PATTERNLENGTH = 1.5
    LINETHICKNESS = 0.2
    }
  SYMBOLS
    {
    SHOW = NO
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  BARCHARTS
    {
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  ERRORBARS
    {
    COLOR = BLACK
    }
##FLEXI N6
$!LINEMAP  [4]
  NAME = '&ZN&'
  ASSIGN
    {
    ZONE = 4
    XAXISVAR = 2
    YAXISVAR = 1
    }
  LINES
    {
    SHOW = YES
    COLOR = BLACK
    LINEPATTERN = DASHED
    PATTERNLENGTH = 1.0
    LINETHICKNESS = 0.2
    }
  SYMBOLS
    {
    SHOW = NO
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  BARCHARTS
    {
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  ERRORBARS
    {
    COLOR = BLACK
    }
#FLEXI N9
$!LINEMAP  [5]
  NAME = '&ZN&'
  ASSIGN
    {
    ZONE = 5
    XAXISVAR = 2
    YAXISVAR = 1
    }
  LINES
    {
    SHOW = YES
    COLOR = BLACK
    LINEPATTERN = SOLID
    PATTERNLENGTH = 0.8
    LINETHICKNESS = 0.2
    }
  SYMBOLS
    {
    SHOW = NO
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  BARCHARTS
    {
    COLOR = BLACK
    FILLCOLOR = BLACK
    }
  ERRORBARS
    {
    COLOR = BLACK
    }
$!XYLINEAXIS 
  DEPXTOYRATIO = 1
$!XYLINEAXIS 
  GRIDAREA{DRAWBORDER = YES}
$!XYLINEAXIS 
  YDETAIL 1
    {
    RANGEMIN = 0.00
    RANGEMAX = 3.035
    GRSPACING = 0.5
    TICKS
      {
      TICKDIRECTION = CENTERED
      }
    TITLE
      {
      TITLEMODE = USETEXT
      OFFSET = 6
      TEXT = 'y/h'
      }
    }
$!XYLINEAXIS 
  XDETAIL 1
    {
    RANGEMIN = |MYRANGEMIN|
    RANGEMAX = |MYRANGEMAX|
    TICKS
      {
      TICKDIRECTION = CENTERED
      }
    TITLE
      {
      TITLEMODE = USETEXT
      OFFSET = 6
      TEXT = 'FIGIAXIS'
      }
    }
$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 12.66618738404446
    Y = 2.715909090909063
    }
  TEXTSHAPE
    {
    HEIGHT = 20
    }
  TEXT = 'Position: x=|MYXPOS|'
$!LINEPLOTLAYERS 
  SHOWSYMBOLS = YES
$!FRAMECONTROL ACTIVATEBYNUMBER
  FRAME = 1
$!SETSTYLEBASE CONFIG

# Export as EPS, LPK and PNG
$!EXPORTSETUP EXPORTFORMAT = EPS
$!EXPORTSETUP IMAGEWIDTH = 1275
$!EXPORTSETUP EPSPREVIEWIMAGE{IMAGETYPE = NONE}
#$!EXPORTSETUP EXPORTFNAME = "|PATH|/phill_|MYRE|_|MYN|_|MYVARNAME|_|FILEIND|.eps"
$!EXPORTSETUP EXPORTFNAME = "|PATH|/phill_|MYRE|_|MYVARNAME|_|FILEIND|.eps"
$!EXPORT 
  EXPORTREGION = CURRENTFRAME

#$!SAVELAYOUT = "|PATH|/phill_|MYRE|_|MYN|_|MYVARNAME|_|FILEIND|.lpk"
$!SAVELAYOUT = "|PATH|/phill_|MYRE|_|MYVARNAME|_|FILEIND|.lpk"
  INCLUDEDATA = YES
  INCLUDEPREVIEW = NO

$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH = 1024
$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
#$!EXPORTSETUP EXPORTFNAME = "|PATH|/phill_|MYRE|_|MYN|_|MYVARNAME|_|FILEIND|.png"
$!EXPORTSETUP EXPORTFNAME = "|PATH|/phill_|MYRE|_|MYVARNAME|_|FILEIND|.png"
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
