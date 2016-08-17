#!MC 1100
$!VarSet |MFBD| = '/univ_2/ws3/ws/iagtbole-phill_10595-0/postproc'
$!OPENLAYOUT  '/univ_2/ws3/ws/iagtbole-phill_10595-0/postproc/video_velslices_x.lay'

$!REDRAWALL
$!EXPORTSETUP EXPORTFORMAT = JPEG
$!EXPORTSETUP IMAGEWIDTH = 2000
$!EXPORTSETUP QUALITY = 100
$!EXPORTSETUP EXPORTFNAME = '/univ_2/ws3/ws/iagtbole-phill_10595-0/postproc/picture/video_velslices_x/phill_video_velslices_x_0001060.0000000.jpg'
$!EXPORT
  EXPORTREGION = WORKAREA
$!QUIT

$!RemoveVar |MFBD|
