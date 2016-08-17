#!MC 1100
$!VarSet |MFBD| = '/univ_2/ws3/ws/iagtbole-phill_10595-0/postproc/qcrit'
$!OPENLAYOUT  '/univ_2/ws3/ws/iagtbole-phill_10595-0/postproc/qcrit/video_q.lay'

$!REDRAWALL
$!EXPORTSETUP EXPORTFORMAT = JPEG
$!EXPORTSETUP IMAGEWIDTH = 2000
$!EXPORTSETUP QUALITY = 100
$!EXPORTSETUP EXPORTFNAME = '/univ_2/ws3/ws/iagtbole-phill_10595-0/postproc/qcrit/picture/video_q/../phill_video_q_0000950.0000000.jpg'
$!EXPORT
  EXPORTREGION = CURRENTFRAME
$!QUIT

$!RemoveVar |MFBD|
