load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;=================================
begin
;=================================
; Settings.
;---------------------------------
;nz = 16
;nx = 128

;nz = 32
;nx = 256

nz = 64
nx = 512

;nz = 128
;nx = 1024


;nt = 4
;nt = 3
;nt = 2
;nt = 1
;nt = 11
;nt = 17
;nt = 26
;nt = 28
;nt = 32
;nt = 41
;nt = 47
;nt = 49
;nt = 51
;nt = 56
;nt = 64
;nt = 69
;nt = 85
;nt = 86
nt = 91
;nt = 101
;nt = 1001
;nt = 201
;nt = 501
;nt = 1001
;nt = 5001
;nt = 90
;nt = 601
;nt = 501
;nt = 1
;nt = 1501
var = "theta_M_1"
;var = "theta_M"
;var = "theta"
;var = "pi_1"
;var = "u"
;var = "w"
halo = 4
;=================================
; Read data.
;---------------------------------
setfileoption ("bin", "ReadByteOrder", "BigEndian")
;data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz,nx/),"float")
data = cbinread("../output/modelvar_"+var+".bin",(/nt,nz+2*halo,nx+2*halo/),"double")
;print(data)
;=================================
; Plot.
;---------------------------------
wks = gsn_open_wks ("pdf","dc")                  ; open workstation
;gsn_define_colormap(wks,"temp_diff_18lev")
;gsn_define_colormap(wks,"WhiteYellowOrangeRed")
gsn_define_colormap(wks,"WhiteBlue")

res = True

res@gsnMaximize = True
res@vpWidthF  = 1.0                    ; change the aspect ratio
res@vpHeightF = 0.23
;res@vpWidthF  = 0.5                    ; change the aspect ratio
;res@vpHeightF = 1.0
;res@vpXF      = .1                     ; location of where plot starts
;res@vpYF      = .7
res@tmXBLabelFontHeightF = 0.015       ; label font height
res@tmYLLabelFontHeightF = 0.015

res@tmXTOn = False
res@tmYROn = False


; LabelBar
res@cnFillOn = True
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
;res@pmLabelBarDisplayMode   = "Always" ; turn on color bar
res@pmLabelBarSide          = "Bottom" ; put it on the bottom of plot
res@pmLabelBarWidthF        = 1.0      ; control size of colorbar
res@pmLabelBarHeightF       = 0.1      ;
;res@pmLabelBarOrthogonalPosF = -0.02   ; position wrt plot
res@lbLabelFontHeightF      = 0.015     ; label bar font
res@lbOrientation           = "Horizontal" ; label bar orientation
res@lbPerimOn               = False    ; no box around label bar
res@lbAutoManage            = False	 ; we control how drawn not plot 
;res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
;res@cnLevels            = (//)

res@gsnSpreadColors = True
res@gsnSpreadColorStart = 250
res@gsnSpreadColorEnd = 2

res@cnLevelSelectionMode = "ManualLevels"   ; set explicit contour levels
res@cnMinLevelValF          = -14.
res@cnMaxLevelValF          = -1.
res@cnLevelSpacingF         = 1.
;res@cnMinLevelValF          = 0.10
;res@cnMaxLevelValF          = 0.20
;res@cnLevelSpacingF         = 0.001

;res@cnFillColors = ispan(250,10,16)

; Line
;res@cnFillOn = False
;res@cnLinesOn = True
;res@cnLineLabelsOn = True
;res@cnMonoLevelFlag = False
;res@cnLevelFlags = "LineAndLabel"

res@tiMainString = "Density Current"
;res@tiXAxisString      = "X"  
;res@tiYAxisString      = "Z (height)"

res@tiMainFontHeightF  = 0.015                    ; change font heights
;res@tiXAxisFontHeightF  = 0.015                    ; change font heights
;res@tiYAxisFontHeightF  = 0.015                    ; change font heights

;res@gsnLeftString   = "dx = 100m; dz = 100m; u = 50m/s; w = 0m/s; dt = 1s"               ; add the gsn titles
;res@gsnLeftString   = "dx = 100m; dz = 100m; u = 50m/s; w = 0m/s; dt = 0.1s"               ; add the gsn titles
;res@gsnLeftString   = "dx = 100m; dz = 100m; dt = 1s"               ; add the gsn titles
res@gsnLeftString   = "dx = 100m; dz = 100m; dt = 0.1s"               ; add the gsn titles
;res@gsnLeftString   = "dx = 50m; dz = 50m; dt = 0.1s"               ; add the gsn titles
;res@gsnLeftString   = "dx = 100m; dz = 100m; dt = 0.05s"               ; add the gsn titles
;res@gsnLeftString   = "dx = 100m; dz = 100m; dt = 0.01s"               ; add the gsn titles
res@gsnLeftStringFontHeightF   = 0.015

;plot = gsn_contour(wks,pi_data,res)
plots = new(nt,graphic)
;do i = 0, nt-1, 100
do i = 0, nt-1
	print("Step: "+i+" plotting...")
	;res@gsnCenterString = i
	;res@gsnCenterStringFontHeightF   = 0.015
	k = i*100
	;k = i*1000
	res@gsnRightString  = "Step: " + k
	;res@gsnRightString  = "Step: " + i
	res@gsnRightStringFontHeightF   = 0.015
	;plots(0) = gsn_csm_contour(wks,data(i,:,:),res)
	plots(0) = gsn_csm_contour(wks,data(i,halo:halo+nz-1,halo:halo+nx-1),res)
end do

;=================================
end
;=================================
