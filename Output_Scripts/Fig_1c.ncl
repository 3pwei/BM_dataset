;*********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;**********************************************************************
begin
;Mtype
;1 闊葉樹林型  2 竹林            3 針葉樹林型
;4 竹闊混淆林  5 針闊葉樹混淆林  6 竹針混淆林	
;--------------------------------------------------------------------------
;all Stype
; 1 肖楠人工林          2 樟樹人工林      3 柳杉人工林      4 檜木人工林
; 5 松類人工林          6 櫸木人工林      7 杉木人工林      8 大葉桃花心木
; 9 松類天然林         10 叢生狀竹       11 單桿狀竹       12 雲杉林
;13 鐵杉林             14 相思樹人工林   15 肖楠人工林     16 台灣櫸人工林
;17 杉木人工林         18 桃花心木人工林 19 光臘樹人工林   20 台灣杉人工林
;21 赤楊人工林         22 其他闊葉人工林 23 其他人工針葉林 24 冷杉林
;25 天然檜木林         26 天然松類       27 其他針葉林     28 闊葉樹混淆林
;29 天然闊葉混淆林     30 竹闊混淆林     31 竹針混淆林     32 針闊葉樹混淆林
;33 人工針闊葉樹混淆林 34 其他人工闊葉樹混淆林
;--------------------------------------------------------------------------
  msg_val = -999.
;---------------------------------------------------------------------
  dir = "/work/flyeh/DBH/Data/"
  filename = "DBHin2.nc"
  infile = addfile(dir+filename,"r")
  in_SampleID = infile->SampleID
  in_Mtype = infile->Mtype
  in_Stype = infile->Stype
  in_allStype = infile->allStype
  in_DBH = infile->DBH
  in_VOL = infile->VOL

  filename2 = "DBHout2.nc"
  infile2 = addfile(dir+filename2,"r")
  out_SampleID = infile2->SampleID
  out_Mtype = infile2->Mtype
  out_Stype = infile2->Stype
  out_allStype = infile2->allStype
  out_DBH = infile2->DBH
  out_VOL = infile2->VOL

;printVarSummary(in_Mtype)
;printVarSummary(out_Mtype)
;選取要的年份資料-------------------------------------------
;過濾事業區內外樹種--------------------------------------------------------------------------------
  ;事業區內-------------------------------------------------------------------------
  ;needle leaf
  Mtype_in_needle = where(in_Mtype.eq.3,in_Mtype,in_Mtype@_FillValue)
  DBH_in_needle = where(in_Mtype.eq.3,in_DBH,in_DBH@_FillValue)

  ;事業區外-------------------------------------------------------------------------
  ;needle leaf
  Mtype_out_needle = where(out_Mtype.eq.3,out_Mtype,out_Mtype@_FillValue)
  DBH_out_needle = where(out_Mtype.eq.3,out_DBH,out_DBH@_FillValue)

;計算DBH區間的個數---------------------------------------------------------------------------------
; 0<  <=5   ,   5<  <=10  ,  10<  <=15  ,  15<  <=20  ,  20<  <=25
;25<  <=30  ,  30<  <=35  ,  35<  <=40  ,  40<  <=45  ,  45<  <=50
;50<  <=55  ,  55<  <=60  ,  60<  <=65  ,  65<  <=70  ,  >70

  n_in  = dimsizes(in_Mtype)
  n_out = dimsizes(out_Mtype)
  n = n_in+n_out

  broad   = new((/(70/5)+1/),"integer")
  needle  = new((/(70/5)+1/),"integer")
  do i = 0,(70/5)+1-1,1
    if(i.lt.(70/5)+1-1)then
      print(""+(0+i*5)+"<  <="+(5+i*5)+"")
      needle(i)     = num(DBH_in_needle.gt.(0+i*5).and.DBH_in_needle.le.(5+i*5))+\
                      num(DBH_out_needle.gt.(0+i*5).and.DBH_out_needle.le.(5+i*5))
    end if
    if(i.eq.(70/5)+1-1)then
      print(">"+(i*5)+"")
      needle(i)     = num(DBH_in_needle.gt.i*5)+num(DBH_out_needle.gt.i*5)
    end if
  end do

  needle_sum = dim_sum_n(needle(1:14),0)
  needle_per = (tofloat(needle(1:14))/tofloat(needle_sum))*100
print(needle_per)
print(sum(needle_per))
;modle----------------------------------------------------------------------------
  year = ispan(1974,1978,1)
  ny = dimsizes(year)
  DIAMETER = new((/ny,28,86,41/),"float")
  dir2 = "/lfs/home/whuang/RE_EXPT/EXPT1/stomate/"
  do yyyy = year(0),year(ny-1),1
    filename2 = "re.exp.t1_"+yyyy+"0101_"+yyyy+"1231_1M_stomate_history.nc"
    infile2 = addfile(dir2+filename2,"r")
    DIAMETER(yyyy-year(0),:,:,:) = infile2->DIAMETER(0,:,:,:)*100
  end do
  DIAMETER!0 = "year"
  printVarSummary(DIAMETER)

  DIAMETER_needle_od = ndtooned(DIAMETER(:,3,:,:)) ;PFT4常綠針葉
  qsort(DIAMETER_needle_od)
  nm = num(.not.ismissing(DIAMETER_needle_od))
  DIAMETER_needle_oned = DIAMETER_needle_od(0:nm-1)
  delete(nm)

  nm = new((/ny/),"integer")
  yDIAMETER_needle_od = new((/ny,dimsizes(DIAMETER(0,0,:,0))*dimsizes(DIAMETER(0,0,0,:))/),"float")
  do i = 0,ny-1,1
    yDIAMETER_needle_od(i,:) = ndtooned(DIAMETER(i,3,:,:)) ;PFT4常綠針葉
    qsort(yDIAMETER_needle_od(i,:))
    nm(i) = num(.not.ismissing(yDIAMETER_needle_od(i,:)))
  end do
  yDIAMETER_needle_oned = yDIAMETER_needle_od(:,:nm(0)-1)

  DIAMETER_needle  = new((/(70/5)+1/),"integer")
  yDIAMETER_needle = new((/ny,(70/5)+1/),"integer")
  do yy = 0,ny-1,1
  do i = 0,(70/5)+1-1,1
    if(i.lt.(70/5)+1-1)then
      print(""+(0+i*5)+"<  <="+(5+i*5)+"")
      DIAMETER_needle(i)  = num(DIAMETER_needle_oned.gt.(0+i*5).and.DIAMETER_needle_oned.le.(5+i*5))
      yDIAMETER_needle(yy,i)  = num(yDIAMETER_needle_oned(yy,:).gt.(0+i*5).and.yDIAMETER_needle_oned(yy,:).le.(5+i*5))
    end if
    if(i.eq.(70/5)+1-1)then
      print(">"+(i*5)+"")
      DIAMETER_needle(i)  = num(DIAMETER_needle_oned.gt.i*5)
      yDIAMETER_needle(yy,i)  = num(yDIAMETER_needle_oned(yy,:).gt.i*5)
    end if
  end do
  end do

  DIAMETER_needle_sum = dim_sum_n(DIAMETER_needle(1:14),0)
  DIAMETER_needle_per = (tofloat(DIAMETER_needle(1:14))/tofloat(DIAMETER_needle_sum))*100
  print(DIAMETER_needle_per)

  line = new((/2,dimsizes(needle)/),"float",msg_val)
  line(0,1:dimsizes(needle)-1) = (/needle_per/)
  line(1,1:dimsizes(needle)-1) = (/DIAMETER_needle_per/)

  yDIAMETER_needle_sum = dim_sum_n(yDIAMETER_needle(:,1:14),1)
  yDIAMETER_needle_per = new((/ny,14/),"float")
  do yy = 0,ny-1,1
    yDIAMETER_needle_per(yy,:) = (tofloat(yDIAMETER_needle(yy,1:14))/tofloat(yDIAMETER_needle_sum(yy)))*100
  end do
  print(sum(DIAMETER_needle_per))
  do i = 0,ny-1,1
    print(sum(yDIAMETER_needle_per(i,:)))
  end do
;print(max(yDIAMETER_needle_per(:,1)))
;calculate 四分位距-------------------------------------------
  nn = dimsizes(needle)
  qsyDIAMETER_needle_per = yDIAMETER_needle_per
  do i = 0,13,1
    qsort(qsyDIAMETER_needle_per(:,i))
  end do

  mmin = new((/nn/),"float")
  do i = 1,nn-1,1
    mmin(i) = min(qsyDIAMETER_needle_per(:,i-1))
  end do
  delete(i)

  mmax = new((/nn/),"float")
  do i = 1,nn-1,1
    mmax(i) = max(qsyDIAMETER_needle_per(:,i-1))
  end do
  delete(i)

  med = new((/nn/),"float")
  do i = 1,nn-1,1
    med(i) = dim_median_n(qsyDIAMETER_needle_per(:,i-1),0)
  end do
  delete(i)

  Q1 = new((/nn/),"float")
  do i = 1,nn-1,1
    Q1(i) = qsyDIAMETER_needle_per(round(0.25*ny,3)-1,i-1) ;NCL數字要-1
  end do
  delete(i)

  Q3 = new((/nn/),"float")
  do i = 1,nn-1,1
    Q3(i) = qsyDIAMETER_needle_per(round(0.75*ny,3)-1,i-1) ;NCL數字要-1
  end do
  delete(i)

  iarr = new((/dimsizes(needle),5/),"float")
  iarr(:,0) = (/mmin/)
  iarr(:,1) = (/Q1/)
  iarr(:,2) = (/med/)
  iarr(:,3) = (/Q3/)
  iarr(:,4) = (/mmax/)

;exit
;plot**********************************************************************
  type = "ps"
;  type@wkHeight = 2048
;  type@wkWidth  = 2048
  wks = gsn_open_wks(type,"DBH_modle_needle_box_EXPT1_903")

  res                             = True
  res@gsnDraw                     = False
  res@gsnFrame                    = False
  res@gsnMaximize                 = True
  res@gsnLeftString               = "Diameter at breast height (%)"
  res@gsnRightString              = "PFT4 Needle leaf";"~F34~0~F~C year~S~-1~N~"+"~S~o~N~C /year"
  res@gsnStringFontHeightF        = 0.023
  res@gsnYRefLine                 = 0.0
  res@gsnYRefLineThicknessF       = 3

  ;---Make legend smaller and move into plot
  res@pmLegendDisplayMode         = "Always"            ; turn on legend
  res@pmLegendSide                = "Top"               ; Change location of
  res@pmLegendParallelPosF        = .78                 ; move units right
  res@pmLegendOrthogonalPosF      = -0.68              ; move units down
  res@pmLegendWidthF              = 0.085               ; Change width and
  res@pmLegendHeightF             = 0.09               ; height of legend.

  res@lgPerimOn                   = False                ; turn off/on box around
  res@lgLabelFontHeightF          = .018                ; label font height

  res@xyExplicitLegendLabels      = (/"  Observation","  ORCHIDEE DIAMETER"/)

  res@vpOn                        = True
  res@vpXF                        = 0.1
  res@vpYF                        = 0.7
  res@vpWidthF                    = 0.9;0.6
  res@vpHeightF                   = 0.3;0.3

  res@xyLineThicknessF            = 5.0
  res@xyLineColors                = (/"darkgreen","orange3"/) ;chartreuse3
;  res@xyLineOpacities             = (/1,1/)
  res@xyDashPatterns              = 0;(/0,1,1/)

  res@tmXTOn                      = False
  res@tmBorderThicknessF          = 3.0
  res@tmXBMajorThicknessF         = 3.0
  res@tmXBMinorThicknessF         = 3.0
  res@tmYLMajorThicknessF         = 3.0
;  res@tmYLMajorLineColor          = "red"
  res@tmYLMinorThicknessF         = 3.0
;  res@tmYLMinorLineColor          = "red"
;  res@tmYRMajorThicknessF         = 5.0
;  res@tmYRMinorThicknessF         = 5.0

  res@tmXBMode                    = "Explicit"
  res@tmXBValues                  = ispan(0,70,5)
  res@tmXBLabels                  = (/"~F18~0<HK5   ","~F18~5<HK10   ","~F18~10<HK15   ",\
                                      "~F18~15<HK20   ","~F18~20<HK25   ","~F18~25<HK30   ",\
                                      "~F18~30<HK35   ","~F18~35<HK40   ","~F18~40<HK45   ",\
                                      "~F18~45<HK50   ","~F18~50<HK55   ","~F18~55<HK60   ",\
                                      "~F18~60<HK65   ","~F18~65<HK70   ","~F18~70<H   "/)
  res@tmXBLabelAngleF             = 25
  res@tmXBLabelFontHeightF        = 0.015
  res@tmXBLabelFontThicknessF     = 2.5

;  res@tmXBMinorValues             = ispan(1,101,1)

;  res@tmYLMode                    = "Explicit"                ; explicit labels
;  res@tmYLValues                  = fspan(-5,5,11)
;  res@tmYLLabels                  = ""+res@tmYLValues         ; make strings
;  res@tmYLMinorValues             = fspan(-5,5,51)
  res@tmYLLabelFontHeightF        = 0.018

  res@tmXMajorGrid                = True          ; implement x grid
  res@tmXMajorGridLineColor       = (/"gray60"/)
  res@tmXMajorGridThicknessF      = 1;5.0           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines

  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridLineColor       = (/"gray60"/)
  res@tmYMajorGridThicknessF      = 1;5.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
;  res@tmYLFormat                  = "#^se"
;  res@tmGridDrawOrder             = "PreDraw"


  res@trXMinF                     = 0
  res@trXMaxF                     = 70
  res@trYMaxF                     = 30
  res@trYMinF                     = 0

  res@tiYAxisOn                   = False

  x = ispan(0,70,5)
  plot = gsn_csm_xy(wks,x,line,res)

;************************************************
; add error bars
;************************************************
  polyres                      = True                       ; marker resources
  polyres@gsMarkerIndex        = 1                          ; index for circle
  polyres@gsMarkerSizeF        = .02                        ; size
  polyres@gsLineThicknessF     = 2                        ; size

  centers   = new(14,graphic)
  ycenters   = new(14,graphic)

  c = needle_per
  cy = DIAMETER_needle_per
  do t = 0,14-1,1
; add a polymarker as the center of the error bar located at each data
; point
    polyres@gsMarkerColor     = "skyblue4"                      ; index for circle
    centers(t)    = gsn_add_polymarker(wks,plot,x(t+1),c(t),polyres)

    polyres@gsMarkerColor     = "black"                      ; index for circle
    ycenters(t)   = gsn_add_polymarker(wks,plot,x(t+1),cy(t),polyres)
  end do
  bres                             = True

  bres@vpOn                        = True
  bres@vpXF                        = 0.1
  bres@vpYF                        = 0.7
  bres@vpWidthF                    = 0.9;0.6
  bres@vpHeightF                   = 0.3;0.3

  bres@tmXTOn                      = False
  bres@tmXBOn                      = False
  bres@tmYLOn                      = False
  bres@tmYROn                      = False

  bplot = boxplot(wks,x,iarr,False,bres,False)

  overlay(plot,bplot)

draw(plot)
frame(wks)
end
